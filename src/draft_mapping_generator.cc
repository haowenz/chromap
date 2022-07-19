#include "draft_mapping_generator.h"

#include <vector>

#include "alignment.h"

namespace chromap {

void DraftMappingGenerator::GenerateDraftMappings(
    const SequenceBatch &read_batch, uint32_t read_index,
    const SequenceBatch &reference, MappingMetadata &mapping_metadata) {
  mapping_metadata.SetMinNumErrors(error_threshold_ + 1);
  mapping_metadata.SetNumBestMappings(0);
  mapping_metadata.SetSecondMinNumErrors(error_threshold_ + 1);
  mapping_metadata.SetNumSecondBestMappings(0);

  // Directly obtain the non-split mapping in ideal case and return without
  // running actual verification.
  const bool is_mapping_generated =
      GenerateNonSplitDraftMappingSupportedByAllMinimizers(
          read_batch, read_index, reference, mapping_metadata);
  if (is_mapping_generated) {
    return;
  }

  // Use more sophicated approach to obtain the mapping.
  // Sort the candidates by their count in descending order.
  // TODO: check if this sorting is necessary.
  mapping_metadata.SortCandidates();

  // For split alignments, SIMD cannot be used.
  if (split_alignment_) {
    GenerateDraftMappingsOnOneDirection(kPositive, read_index, read_batch,
                                        reference, mapping_metadata);

    GenerateDraftMappingsOnOneDirection(kNegative, read_index, read_batch,
                                        reference, mapping_metadata);
    return;
  }

  // For non-split alignments, use SIMD when possible.
  if (mapping_metadata.GetNumPositiveCandidates() < (size_t)num_vpu_lanes_) {
    GenerateDraftMappingsOnOneDirection(kPositive, read_index, read_batch,
                                        reference, mapping_metadata);
  } else {
    GenerateDraftMappingsOnOneDirectionUsingSIMD(
        kPositive, read_index, read_batch, reference, mapping_metadata);
  }

  if (mapping_metadata.GetNumNegativeCandidates() < (size_t)num_vpu_lanes_) {
    GenerateDraftMappingsOnOneDirection(kNegative, read_index, read_batch,
                                        reference, mapping_metadata);
  } else {
    GenerateDraftMappingsOnOneDirectionUsingSIMD(
        kNegative, read_index, read_batch, reference, mapping_metadata);
  }
}

bool DraftMappingGenerator::IsValidCandidate(uint32_t rid, uint32_t position,
                                             uint32_t read_length,
                                             const SequenceBatch &reference) {
  const uint32_t reference_length = reference.GetSequenceLengthAt(rid);

  if (position < (uint32_t)error_threshold_ || position >= reference_length ||
      position + read_length + (uint32_t)error_threshold_ >= reference_length) {
    return false;
  }

  return true;
}

bool DraftMappingGenerator::
    GenerateNonSplitDraftMappingSupportedByAllMinimizers(
        const SequenceBatch &read_batch, uint32_t read_index,
        const SequenceBatch &reference, MappingMetadata &mapping_metadata) {
  if (split_alignment_) {
    return false;
  }

  const bool has_one_candidate = (mapping_metadata.GetNumCandidates() == 1);

  if (!has_one_candidate) {
    return false;
  }

  const std::vector<Candidate> &positive_candidates =
      mapping_metadata.positive_candidates_;
  const std::vector<Candidate> &negative_candidates =
      mapping_metadata.negative_candidates_;

  const std::vector<std::pair<uint64_t, uint64_t>> &minimizers =
      mapping_metadata.minimizers_;

  std::vector<DraftMapping> &positive_mappings =
      mapping_metadata.positive_mappings_;
  std::vector<DraftMapping> &negative_mappings =
      mapping_metadata.negative_mappings_;

  uint32_t num_all_minimizer_candidates = 0;
  uint32_t all_minimizer_candidate_index = 0;
  Direction all_minimizer_candidate_direction = kPositive;

  for (uint32_t i = 0; i < positive_candidates.size(); ++i) {
#ifdef LI_DEBUG
    printf("%s + %u %u %d:%d\n", __func__, i, positive_candidates[i].count,
           (int)(positive_candidates[i].position >> 32),
           (int)positive_candidates[i].position);
#endif
    if (positive_candidates[i].count == minimizers.size()) {
      all_minimizer_candidate_index = i;
      ++num_all_minimizer_candidates;
    }
  }

  for (uint32_t i = 0; i < negative_candidates.size(); ++i) {
#ifdef LI_DEBUG
    printf("%s - %u %u %d:%d\n", __func__, i, negative_candidates[i].count,
           (int)(negative_candidates[i].position >> 32),
           (int)negative_candidates[i].position);
#endif
    if (negative_candidates[i].count == minimizers.size()) {
      all_minimizer_candidate_index = i;
      all_minimizer_candidate_direction = kNegative;
      ++num_all_minimizer_candidates;
    }
  }

  if (num_all_minimizer_candidates != 1) {
    return false;
  }

  mapping_metadata.SetMinNumErrors(0);
  mapping_metadata.SetNumBestMappings(1);
  mapping_metadata.SetNumSecondBestMappings(0);

  uint32_t rid = 0;
  uint32_t position = 0;
  const uint32_t read_length = read_batch.GetSequenceLengthAt(read_index);

  if (all_minimizer_candidate_direction == kPositive) {
    rid = positive_candidates[all_minimizer_candidate_index]
              .GetReferenceSequenceIndex();
    position = positive_candidates[all_minimizer_candidate_index]
                   .GetReferenceSequencePosition();
  } else {
    rid = negative_candidates[all_minimizer_candidate_index]
              .GetReferenceSequenceIndex();
    position = negative_candidates[all_minimizer_candidate_index]
                   .GetReferenceSequencePosition() -
               read_length + 1;
  }

  const bool is_valid_candidate =
      IsValidCandidate(rid, position, read_length, reference);
  if (is_valid_candidate) {
    if (all_minimizer_candidate_direction == kPositive) {
      positive_mappings.emplace_back(
          0, positive_candidates[all_minimizer_candidate_index].position +
                 read_length - 1);
    } else {
      negative_mappings.emplace_back(
          0, negative_candidates[all_minimizer_candidate_index].position);
    }
    return true;
  }

  return false;
}

void DraftMappingGenerator::GenerateDraftMappingsOnOneDirectionUsingSIMD(
    Direction candidate_direction, uint32_t read_index,
    const SequenceBatch &read_batch, const SequenceBatch &reference,
    MappingMetadata &mapping_metadata) {
  const char *read = read_batch.GetSequenceAt(read_index);
  const uint32_t read_length = read_batch.GetSequenceLengthAt(read_index);
  const std::string &negative_read =
      read_batch.GetNegativeSequenceAt(read_index);

  const std::vector<Candidate> &candidates =
      candidate_direction == kPositive ? mapping_metadata.positive_candidates_
                                       : mapping_metadata.negative_candidates_;
  std::vector<DraftMapping> &mappings =
      candidate_direction == kPositive ? mapping_metadata.positive_mappings_
                                       : mapping_metadata.negative_mappings_;
  int &min_num_errors = mapping_metadata.min_num_errors_;
  int &num_best_mappings = mapping_metadata.num_best_mappings_;
  int &second_min_num_errors = mapping_metadata.second_min_num_errors_;
  int &num_second_best_mappings = mapping_metadata.num_second_best_mappings_;

  Candidate valid_candidates[num_vpu_lanes_];
  const char *valid_candidate_starts[num_vpu_lanes_];
  uint32_t valid_candidate_index = 0;
  size_t candidate_index = 0;
  uint32_t candidate_count_threshold = 0;

  while (candidate_index < candidates.size()) {
    if (candidates[candidate_index].count < candidate_count_threshold) {
      break;
    }

    uint32_t rid = candidates[candidate_index].GetReferenceSequenceIndex();
    uint32_t position =
        candidates[candidate_index].GetReferenceSequencePosition();

    if (candidate_direction == kNegative) {
      position = position - read_length + 1;
    }

    if (!IsValidCandidate(rid, position, read_length, reference)) {
      ++candidate_index;
      continue;
    }

    valid_candidates[valid_candidate_index] = candidates[candidate_index];
    valid_candidate_starts[valid_candidate_index] =
        reference.GetSequenceAt(rid) + position - error_threshold_;
    ++valid_candidate_index;
    ++candidate_index;

    if (valid_candidate_index < (uint32_t)num_vpu_lanes_) {
      continue;
    }

    if (num_vpu_lanes_ == 8) {
      int16_t mapping_edit_distances[num_vpu_lanes_];
      int16_t mapping_end_positions[num_vpu_lanes_];
      for (int li = 0; li < num_vpu_lanes_; ++li) {
        mapping_end_positions[li] = read_length - 1;
      }
      if (candidate_direction == kPositive) {
        BandedAlign8PatternsToText(error_threshold_, valid_candidate_starts,
                                   read, read_length, mapping_edit_distances,
                                   mapping_end_positions);
      } else {
        BandedAlign8PatternsToText(
            error_threshold_, valid_candidate_starts, negative_read.data(),
            read_length, mapping_edit_distances, mapping_end_positions);
      }
      for (int mi = 0; mi < num_vpu_lanes_; ++mi) {
        if (mapping_edit_distances[mi] <= error_threshold_) {
          if (mapping_edit_distances[mi] < min_num_errors) {
            second_min_num_errors = min_num_errors;
            num_second_best_mappings = num_best_mappings;
            min_num_errors = mapping_edit_distances[mi];
            num_best_mappings = 1;
          } else if (mapping_edit_distances[mi] == min_num_errors) {
            num_best_mappings++;
          } else if (mapping_edit_distances[mi] == second_min_num_errors) {
            num_second_best_mappings++;
          } else if (mapping_edit_distances[mi] < second_min_num_errors) {
            num_second_best_mappings = 1;
            second_min_num_errors = mapping_edit_distances[mi];
          }
          if (candidate_direction == kPositive) {
            mappings.emplace_back(mapping_edit_distances[mi],
                                  valid_candidates[mi].position -
                                      error_threshold_ +
                                      mapping_end_positions[mi]);
          } else {
            mappings.emplace_back(mapping_edit_distances[mi],
                                  valid_candidates[mi].position - read_length +
                                      1 - error_threshold_ +
                                      mapping_end_positions[mi]);
          }
        } else {
          candidate_count_threshold = valid_candidates[mi].count;
        }
      }
    } else if (num_vpu_lanes_ == 4) {
      int32_t mapping_edit_distances[num_vpu_lanes_];
      int32_t mapping_end_positions[num_vpu_lanes_];
      for (int li = 0; li < num_vpu_lanes_; ++li) {
        mapping_end_positions[li] = read_length - 1;
      }
      if (candidate_direction == kPositive) {
        BandedAlign4PatternsToText(error_threshold_, valid_candidate_starts,
                                   read, read_length, mapping_edit_distances,
                                   mapping_end_positions);
      } else {
        BandedAlign4PatternsToText(
            error_threshold_, valid_candidate_starts, negative_read.data(),
            read_length, mapping_edit_distances, mapping_end_positions);
      }
      for (int mi = 0; mi < num_vpu_lanes_; ++mi) {
        if (mapping_edit_distances[mi] <= error_threshold_) {
          if (mapping_edit_distances[mi] < min_num_errors) {
            second_min_num_errors = min_num_errors;
            num_second_best_mappings = num_best_mappings;
            min_num_errors = mapping_edit_distances[mi];
            num_best_mappings = 1;
          } else if (mapping_edit_distances[mi] == min_num_errors) {
            num_best_mappings++;
          } else if (mapping_edit_distances[mi] == second_min_num_errors) {
            num_second_best_mappings++;
          } else if (mapping_edit_distances[mi] < second_min_num_errors) {
            num_second_best_mappings = 1;
            second_min_num_errors = mapping_edit_distances[mi];
          }
          if (candidate_direction == kPositive) {
            mappings.emplace_back(mapping_edit_distances[mi],
                                  valid_candidates[mi].position -
                                      error_threshold_ +
                                      mapping_end_positions[mi]);
          } else {
            mappings.emplace_back(mapping_edit_distances[mi],
                                  valid_candidates[mi].position - read_length +
                                      1 - error_threshold_ +
                                      mapping_end_positions[mi]);
          }
        } else {
          candidate_count_threshold = valid_candidates[mi].count;
        }
      }
    }

    valid_candidate_index = 0;
  }

  for (uint32_t ci = 0; ci < valid_candidate_index; ++ci) {
    uint32_t rid = valid_candidates[ci].GetReferenceSequenceIndex();
    uint32_t position = valid_candidates[ci].GetReferenceSequencePosition();
    if (candidate_direction == kNegative) {
      position = position - read_length + 1;
    }

    if (!IsValidCandidate(rid, position, read_length, reference)) {
      continue;
    }

    int mapping_end_position;
    int num_errors;
    if (candidate_direction == kPositive) {
      num_errors = BandedAlignPatternToText(
          error_threshold_,
          reference.GetSequenceAt(rid) + position - error_threshold_, read,
          read_length, &mapping_end_position);
    } else {
      num_errors = BandedAlignPatternToText(
          error_threshold_,
          reference.GetSequenceAt(rid) + position - error_threshold_,
          negative_read.data(), read_length, &mapping_end_position);
    }
    if (num_errors <= error_threshold_) {
      if (num_errors < min_num_errors) {
        second_min_num_errors = min_num_errors;
        num_second_best_mappings = num_best_mappings;
        min_num_errors = num_errors;
        num_best_mappings = 1;
      } else if (num_errors == min_num_errors) {
        num_best_mappings++;
      } else if (num_errors == second_min_num_errors) {
        num_second_best_mappings++;
      } else if (num_errors < second_min_num_errors) {
        num_second_best_mappings = 1;
        second_min_num_errors = num_errors;
      }
      if (candidate_direction == kPositive) {
        mappings.emplace_back(num_errors, valid_candidates[ci].position -
                                              error_threshold_ +
                                              mapping_end_position);
      } else {
        mappings.emplace_back(num_errors,
                              valid_candidates[ci].position - read_length + 1 -
                                  error_threshold_ + mapping_end_position);
      }
    }
  }
}

void DraftMappingGenerator::GenerateDraftMappingsOnOneDirection(
    Direction candidate_direction, uint32_t read_index,
    const SequenceBatch &read_batch, const SequenceBatch &reference,
    MappingMetadata &mapping_metadata) {
  const char *read = read_batch.GetSequenceAt(read_index);
  const uint32_t read_length = read_batch.GetSequenceLengthAt(read_index);
  const std::string &negative_read =
      read_batch.GetNegativeSequenceAt(read_index);

  const std::vector<Candidate> &candidates =
      candidate_direction == kPositive ? mapping_metadata.positive_candidates_
                                       : mapping_metadata.negative_candidates_;
  std::vector<DraftMapping> &mappings =
      candidate_direction == kPositive ? mapping_metadata.positive_mappings_
                                       : mapping_metadata.negative_mappings_;
  std::vector<int> &split_sites = candidate_direction == kPositive
                                      ? mapping_metadata.positive_split_sites_
                                      : mapping_metadata.negative_split_sites_;
  int &min_num_errors = mapping_metadata.min_num_errors_;
  int &num_best_mappings = mapping_metadata.num_best_mappings_;
  int &second_min_num_errors = mapping_metadata.second_min_num_errors_;
  int &num_second_best_mappings = mapping_metadata.num_second_best_mappings_;

  uint32_t candidate_count_threshold = 0;

  for (uint32_t ci = 0; ci < candidates.size(); ++ci) {
    if (candidates[ci].count < candidate_count_threshold) {
      break;
    }

    uint32_t rid = candidates[ci].GetReferenceSequenceIndex();
    uint32_t position = candidates[ci].GetReferenceSequencePosition();
    if (candidate_direction == kNegative) {
      position = position - read_length + 1;
    }

    if (!IsValidCandidate(rid, position, read_length, reference)) {
      continue;
    }

    int mapping_end_position = read_length;
    int gap_beginning = 0;
    int num_errors = 0;
    const int allow_gap_beginning_ = 20;
    const int mapping_length_threshold = 30;
    int allow_gap_beginning = allow_gap_beginning_ - error_threshold_;
    int actual_num_errors = 0;
    int read_mapping_length = 0;
    int best_mapping_longest_match = 0;
    int longest_match = 0;

    if (split_alignment_) {
      if (candidate_direction == kPositive) {
        num_errors = BandedAlignPatternToTextWithDropOff(
            error_threshold_,
            reference.GetSequenceAt(rid) + position - error_threshold_, read,
            read_length, &mapping_end_position, &read_mapping_length);
        if (mapping_end_position < 0 && allow_gap_beginning > 0) {
          int backup_num_errors = num_errors;
          int backup_mapping_end_position = -mapping_end_position;
          int backup_read_mapping_length = read_mapping_length;
          num_errors = BandedAlignPatternToTextWithDropOff(
              error_threshold_,
              reference.GetSequenceAt(rid) + position - error_threshold_ +
                  allow_gap_beginning,
              read + allow_gap_beginning, read_length - allow_gap_beginning,
              &mapping_end_position, &read_mapping_length);
          if (num_errors > error_threshold_ || mapping_end_position < 0) {
            num_errors = backup_num_errors;
            mapping_end_position = backup_mapping_end_position;
            read_mapping_length = backup_read_mapping_length;
          } else {
            gap_beginning = allow_gap_beginning;
            // Realign the mapping end position as it is the alignment from the
            // whole read.
            mapping_end_position += gap_beginning;
            // I use this adjustment since "position" is based on the whole
            // read, and it will be more consistent with no gap beginning case.
            read_mapping_length += gap_beginning;
          }
        }
      } else {
        num_errors = BandedAlignPatternToTextWithDropOffFrom3End(
            error_threshold_,
            reference.GetSequenceAt(rid) + position - error_threshold_,
            negative_read.data(), read_length, &mapping_end_position,
            &read_mapping_length);
        if (mapping_end_position < 0 && allow_gap_beginning > 0) {
          int backup_num_errors = num_errors;
          int backup_mapping_end_position = -mapping_end_position;
          int backup_read_mapping_length = read_mapping_length;
          num_errors = BandedAlignPatternToTextWithDropOffFrom3End(
              error_threshold_,
              reference.GetSequenceAt(rid) + position - error_threshold_,
              negative_read.data(), read_length - allow_gap_beginning,
              &mapping_end_position, &read_mapping_length);
          if (num_errors > error_threshold_ || mapping_end_position < 0) {
            num_errors = backup_num_errors;
            mapping_end_position = backup_mapping_end_position;
            read_mapping_length = backup_read_mapping_length;
          } else {
            gap_beginning = allow_gap_beginning;
            mapping_end_position += gap_beginning;
            read_mapping_length += gap_beginning;
          }
        }
      }

      if (mapping_end_position + 1 - error_threshold_ - num_errors -
              gap_beginning >=
          mapping_length_threshold) {
        actual_num_errors = num_errors;
        num_errors = -(mapping_end_position - error_threshold_ - num_errors -
                       gap_beginning);

        if (candidates.size() > 200) {
          if (candidate_direction == kPositive) {
            longest_match = GetLongestMatchLength(
                reference.GetSequenceAt(rid) + position, read, read_length);
          } else {
            longest_match =
                GetLongestMatchLength(reference.GetSequenceAt(rid) + position,
                                      negative_read.data(), read_length);
          }
        }
      } else {
        num_errors = error_threshold_ + 1;
        actual_num_errors = error_threshold_ + 1;
      }
    } else {
      if (candidate_direction == kPositive) {
        num_errors = BandedAlignPatternToText(
            error_threshold_,
            reference.GetSequenceAt(rid) + position - error_threshold_, read,
            read_length, &mapping_end_position);
      } else {
        num_errors = BandedAlignPatternToText(
            error_threshold_,
            reference.GetSequenceAt(rid) + position - error_threshold_,
            negative_read.data(), read_length, &mapping_end_position);
      }
    }

    if (num_errors <= error_threshold_) {
      if (num_errors < min_num_errors) {
        second_min_num_errors = min_num_errors;
        num_second_best_mappings = num_best_mappings;
        min_num_errors = num_errors;
        num_best_mappings = 1;
        if (split_alignment_) {
          if (candidates.size() > 50) {
            candidate_count_threshold = candidates[ci].count;
          } else {
            candidate_count_threshold = candidates[ci].count / 2;
          }
          if (second_min_num_errors < min_num_errors + error_threshold_ / 2 &&
              best_mapping_longest_match > longest_match &&
              candidates.size() > 200) {
            second_min_num_errors = min_num_errors;
          }
        }
        best_mapping_longest_match = longest_match;
      } else if (num_errors == min_num_errors) {
        num_best_mappings++;
      } else if (num_errors == second_min_num_errors) {
        num_second_best_mappings++;
      } else if (num_errors < second_min_num_errors) {
        num_second_best_mappings = 1;
        second_min_num_errors = num_errors;
      }

      if (candidate_direction == kPositive) {
        mappings.emplace_back(
            num_errors,
            candidates[ci].position - error_threshold_ + mapping_end_position);
      } else {
        if (split_alignment_ && mapping_output_format_ != MAPPINGFORMAT_SAM) {
          // TODO: this if condition is suspicious. Check this later.
          mappings.emplace_back(num_errors,
                                candidates[ci].position - gap_beginning);
        } else {
          // Need to minus gap_beginning because mapping_end_position is
          // adjusted by it, but read_length is not.
          // printf("%d %d %d\n", candidates[ci].position, mapping_end_position,
          // gap_beginning);
          mappings.emplace_back(num_errors,
                                candidates[ci].position - read_length + 1 -
                                    error_threshold_ + mapping_end_position);
        }
      }

      if (split_alignment_) {
        split_sites.emplace_back(((actual_num_errors & 0xff) << 24) |
                                 ((gap_beginning & 0xff) << 16) |
                                 (read_mapping_length & 0xffff));
      }
    }
  }
}

}  // namespace chromap
