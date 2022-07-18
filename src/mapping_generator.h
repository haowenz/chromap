#ifndef MAPPING_GENERATOR_H_
#define MAPPING_GENERATOR_H_

#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>
#include <tuple>
#include <vector>

#include "alignment.h"
#include "bed_mapping.h"
#include "ksw.h"
#include "mapping.h"
#include "mapping_in_memory.h"
#include "mapping_metadata.h"
#include "mapping_parameters.h"
#include "paf_mapping.h"
#include "paired_end_mapping_metadata.h"
#include "pairs_mapping.h"
#include "sam_mapping.h"
#include "sequence_batch.h"
#include "utils.h"

namespace chromap {

// Class to generate mappings from candidates. It supports multi-threadidng as
// only the parameters are owned by the class.
template <typename MappingRecord>
class MappingGenerator {
 public:
  MappingGenerator() = delete;

  MappingGenerator(const MappingParameters &mapping_parameters,
                   const std::vector<int> &pairs_custom_rid_rank)
      : mapping_parameters_(mapping_parameters),
        NUM_VPU_LANES_(mapping_parameters.GetNumVPULanes()),
        pairs_custom_rid_rank_(pairs_custom_rid_rank) {}

  ~MappingGenerator() = default;

  void VerifyCandidates(const SequenceBatch &read_batch, uint32_t read_index,
                        const SequenceBatch &reference,
                        MappingMetadata &mapping_metadata);

  void GenerateBestMappingsForSingleEndRead(
      const SequenceBatch &read_batch, uint32_t read_index,
      const SequenceBatch &reference, const SequenceBatch &barcode_batch,
      MappingMetadata &mapping_metadata,
      std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs);

  // When the number of supplemented candidates is greater than 0, the
  // force_mapq will be 0, and thereby setting the mapq to 0 (not mapq1 or
  // mapq2). Split alignment won't run candidate supplement.
  void GenerateBestMappingsForPairedEndRead(
      uint32_t pair_index, const SequenceBatch &read_batch1,
      const SequenceBatch &read_batch2, const SequenceBatch &barcode_batch,
      const SequenceBatch &reference, std::vector<int> &best_mapping_indices,
      std::mt19937 &generator, int force_mapq,
      PairedEndMappingMetadata &paired_end_mapping_metadata,
      std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs);

 private:
  bool IsValidCandidate(uint32_t rid, uint32_t position, uint32_t read_length,
                        const SequenceBatch &reference);

  bool DirectlyGenerateOneNonSplitMappingSupportedByAllMinimizers(
      const SequenceBatch &read_batch, uint32_t read_index,
      const SequenceBatch &reference, MappingMetadata &mapping_metadata);

  void VerifyCandidatesOnOneDirectionUsingSIMD(
      Direction candidate_direction, uint32_t read_index,
      const SequenceBatch &read_batch, const SequenceBatch &reference,
      MappingMetadata &mapping_metadata);

  void VerifyCandidatesOnOneDirection(Direction candidate_direction,
                                      uint32_t read_index,
                                      const SequenceBatch &read_batch,
                                      const SequenceBatch &reference,
                                      MappingMetadata &mapping_metadata);

  void ProcessBestMappingsForSingleEndRead(
      Direction mapping_direction, uint32_t read_index,
      const SequenceBatch &read_batch, const SequenceBatch &barcode_batch,
      const SequenceBatch &reference, const MappingMetadata &mapping_metadata,
      const std::vector<int> &best_mapping_indices, int &best_mapping_index,
      int &num_best_mappings_reported,
      std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs);

  void GenerateBestMappingsForPairedEndReadOnOneDirection(
      Direction first_read_direction, Direction second_read_direction,
      uint32_t pair_index, const SequenceBatch &read_batch1,
      const SequenceBatch &read_batch2, const SequenceBatch &reference,
      PairedEndMappingMetadata &paired_end_mapping_metadata);

  void ProcessBestMappingsForPairedEndReadOnOneDirection(
      Direction first_read_direction, Direction second_read_direction,
      uint32_t pair_index, const SequenceBatch &read_batch1,
      const SequenceBatch &read_batch2, const SequenceBatch &barcode_batch,
      const SequenceBatch &reference,
      const std::vector<int> &best_mapping_indices, int &best_mapping_index,
      int &num_best_mappings_reported, int force_mapq,
      const PairedEndMappingMetadata &paired_end_mapping_metadata,
      std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs);

  void GetRefStartEndPositionForReadFromMapping(
      const std::pair<int, uint64_t> &mapping, const SequenceBatch &reference,
      MappingInMemory &mapping_in_memory);

  // For single-end. It should be fully specialized.
  void EmplaceBackSingleEndMappingRecord(
      MappingInMemory &mapping_in_memory,
      std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs);

  // For paired-end. It should be fully specialized.
  void EmplaceBackPairedEndMappingRecord(
      PairedEndMappingInMemory &paired_mapping_in_memory,
      std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs);

  uint8_t GetMAPQForSingleEndRead(Direction direction, int num_errors,
                                  uint16_t alignment_length, int read_length,
                                  int max_num_error_difference,
                                  const MappingMetadata &mapping_metadata);

  uint8_t GetMAPQForPairedEndRead(
      Direction first_read_direction, Direction second_read_direction,
      int read1_num_errors, int read2_num_errors,
      uint16_t read1_alignment_length, uint16_t read2_alignment_length,
      int read1_length, int read2_length, int force_mapq,
      const PairedEndMappingMetadata &paired_end_mapping_metadata,
      uint8_t &mapq1, uint8_t &mapq2);

  const MappingParameters mapping_parameters_;
  const int NUM_VPU_LANES_;
  const std::vector<int> pairs_custom_rid_rank_;
};

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::VerifyCandidates(
    const SequenceBatch &read_batch, uint32_t read_index,
    const SequenceBatch &reference, MappingMetadata &mapping_metadata) {
  int &min_num_errors = mapping_metadata.min_num_errors_;
  int &num_best_mappings = mapping_metadata.num_best_mappings_;
  int &second_min_num_errors = mapping_metadata.second_min_num_errors_;
  int &num_second_best_mappings = mapping_metadata.num_second_best_mappings_;

  min_num_errors = mapping_parameters_.error_threshold + 1;
  num_best_mappings = 0;
  second_min_num_errors = mapping_parameters_.error_threshold + 1;
  num_second_best_mappings = 0;

  // Directly obtain the non-split mapping in ideal case and return without
  // running actual verification.
  const bool is_mapping_generated =
      DirectlyGenerateOneNonSplitMappingSupportedByAllMinimizers(
          read_batch, read_index, reference, mapping_metadata);
  if (is_mapping_generated) {
    return;
  }

  // Use more sophicated approach to obtain the mapping.
  // Sort the candidates by their count in descending order.
  // TODO: check if this sorting is necessary.
  mapping_metadata.SortCandidates();

  // For split alignments, SIMD cannot be used.
  if (mapping_parameters_.split_alignment) {
    VerifyCandidatesOnOneDirection(kPositive, read_index, read_batch, reference,
                                   mapping_metadata);

    VerifyCandidatesOnOneDirection(kNegative, read_index, read_batch, reference,
                                   mapping_metadata);
    return;
  }

  const size_t num_positive_candidates =
      mapping_metadata.positive_candidates_.size();
  const size_t num_negative_candidates =
      mapping_metadata.negative_candidates_.size();

  // For non-split alignments, use SIMD when possible.
  if (num_positive_candidates < (size_t)NUM_VPU_LANES_) {
    VerifyCandidatesOnOneDirection(kPositive, read_index, read_batch, reference,
                                   mapping_metadata);
  } else {
    VerifyCandidatesOnOneDirectionUsingSIMD(kPositive, read_index, read_batch,
                                            reference, mapping_metadata);
  }

  if (num_negative_candidates < (size_t)NUM_VPU_LANES_) {
    VerifyCandidatesOnOneDirection(kNegative, read_index, read_batch, reference,
                                   mapping_metadata);
  } else {
    VerifyCandidatesOnOneDirectionUsingSIMD(kNegative, read_index, read_batch,
                                            reference, mapping_metadata);
  }
}

// Return true if the candidate position is valid on the reference with rid.
// Note only the position is checked and the input rid is not checked in this
// function. So the input rid must be valid.
template <typename MappingRecord>
bool MappingGenerator<MappingRecord>::IsValidCandidate(
    uint32_t rid, uint32_t position, uint32_t read_length,
    const SequenceBatch &reference) {
  const uint32_t reference_length = reference.GetSequenceLengthAt(rid);

  if (position < (uint32_t)mapping_parameters_.error_threshold ||
      position >= reference_length ||
      position + read_length + mapping_parameters_.error_threshold >=
          reference_length) {
    return false;
  }

  return true;
}

// Return true when there is one non-split mapping generated and the mapping is
// supported by all the minimizers.
template <typename MappingRecord>
bool MappingGenerator<MappingRecord>::
    DirectlyGenerateOneNonSplitMappingSupportedByAllMinimizers(
        const SequenceBatch &read_batch, uint32_t read_index,
        const SequenceBatch &reference, MappingMetadata &mapping_metadata) {
  if (mapping_parameters_.split_alignment) {
    return false;
  }

  const std::vector<Candidate> &positive_candidates =
      mapping_metadata.positive_candidates_;
  const std::vector<Candidate> &negative_candidates =
      mapping_metadata.negative_candidates_;

  const bool has_one_candidate =
      (positive_candidates.size() + negative_candidates.size() == 1);

  if (!has_one_candidate) {
    return false;
  }

  const std::vector<std::pair<uint64_t, uint64_t>> &minimizers =
      mapping_metadata.minimizers_;

  std::vector<std::pair<int, uint64_t>> &positive_mappings =
      mapping_metadata.positive_mappings_;
  std::vector<std::pair<int, uint64_t>> &negative_mappings =
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

  mapping_metadata.min_num_errors_ = 0;
  mapping_metadata.num_best_mappings_ = 1;
  mapping_metadata.num_second_best_mappings_ = 0;

  uint32_t rid = 0;
  uint32_t position = 0;
  uint32_t read_length = read_batch.GetSequenceLengthAt(read_index);

  if (all_minimizer_candidate_direction == kPositive) {
    rid = positive_candidates[all_minimizer_candidate_index].position >> 32;
    position = positive_candidates[all_minimizer_candidate_index].position;
  } else {
    rid = negative_candidates[all_minimizer_candidate_index].position >> 32;
    position =
        (uint32_t)negative_candidates[all_minimizer_candidate_index].position -
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

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::VerifyCandidatesOnOneDirectionUsingSIMD(
    Direction candidate_direction, uint32_t read_index,
    const SequenceBatch &read_batch, const SequenceBatch &reference,
    MappingMetadata &mapping_metadata) {
  const char *read = read_batch.GetSequenceAt(read_index);
  uint32_t read_length = read_batch.GetSequenceLengthAt(read_index);
  const std::string &negative_read =
      read_batch.GetNegativeSequenceAt(read_index);

  const std::vector<Candidate> &candidates =
      candidate_direction == kPositive ? mapping_metadata.positive_candidates_
                                       : mapping_metadata.negative_candidates_;
  std::vector<std::pair<int, uint64_t>> &mappings =
      candidate_direction == kPositive ? mapping_metadata.positive_mappings_
                                       : mapping_metadata.negative_mappings_;
  int &min_num_errors = mapping_metadata.min_num_errors_;
  int &num_best_mappings = mapping_metadata.num_best_mappings_;
  int &second_min_num_errors = mapping_metadata.second_min_num_errors_;
  int &num_second_best_mappings = mapping_metadata.num_second_best_mappings_;

  size_t num_candidates = candidates.size();
  Candidate valid_candidates[NUM_VPU_LANES_];
  const char *valid_candidate_starts[NUM_VPU_LANES_];
  uint32_t valid_candidate_index = 0;
  size_t candidate_index = 0;
  uint32_t candidate_count_threshold = 0;
  while (candidate_index < num_candidates) {
    if (candidates[candidate_index].count < candidate_count_threshold) break;
    uint32_t rid = candidates[candidate_index].position >> 32;
    uint32_t position = candidates[candidate_index].position;
    if (candidate_direction == kNegative) {
      position = position - read_length + 1;
    }

    if (!IsValidCandidate(rid, position, read_length, reference)) {
      ++candidate_index;
      continue;
    } else {
      valid_candidates[valid_candidate_index] =
          candidates[candidate_index];  // reference.GetSequenceAt(rid) +
                                        // position -
                                        // mapping_parameters_.error_threshold;
      valid_candidate_starts[valid_candidate_index] =
          reference.GetSequenceAt(rid) + position -
          mapping_parameters_.error_threshold;
      ++valid_candidate_index;
    }
    if (valid_candidate_index == (uint32_t)NUM_VPU_LANES_) {
      if (NUM_VPU_LANES_ == 8) {
        int16_t mapping_edit_distances[NUM_VPU_LANES_];
        int16_t mapping_end_positions[NUM_VPU_LANES_];
        for (int li = 0; li < NUM_VPU_LANES_; ++li) {
          mapping_end_positions[li] = read_length - 1;
        }
        if (candidate_direction == kPositive) {
          BandedAlign8PatternsToText(
              mapping_parameters_.error_threshold, valid_candidate_starts, read,
              read_length, mapping_edit_distances, mapping_end_positions);
        } else {
          BandedAlign8PatternsToText(
              mapping_parameters_.error_threshold, valid_candidate_starts,
              negative_read.data(), read_length, mapping_edit_distances,
              mapping_end_positions);
        }
        for (int mi = 0; mi < NUM_VPU_LANES_; ++mi) {
          if (mapping_edit_distances[mi] <=
              mapping_parameters_.error_threshold) {
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
              mappings.emplace_back((uint8_t)mapping_edit_distances[mi],
                                    valid_candidates[mi].position -
                                        mapping_parameters_.error_threshold +
                                        mapping_end_positions[mi]);
            } else {
              mappings.emplace_back((uint8_t)mapping_edit_distances[mi],
                                    valid_candidates[mi].position -
                                        read_length + 1 -
                                        mapping_parameters_.error_threshold +
                                        mapping_end_positions[mi]);
            }
          } else {
            candidate_count_threshold = valid_candidates[mi].count;
          }
        }
      } else if (NUM_VPU_LANES_ == 4) {
        int32_t mapping_edit_distances[NUM_VPU_LANES_];
        int32_t mapping_end_positions[NUM_VPU_LANES_];
        for (int li = 0; li < NUM_VPU_LANES_; ++li) {
          mapping_end_positions[li] = read_length - 1;
        }
        if (candidate_direction == kPositive) {
          BandedAlign4PatternsToText(
              mapping_parameters_.error_threshold, valid_candidate_starts, read,
              read_length, mapping_edit_distances, mapping_end_positions);
        } else {
          BandedAlign4PatternsToText(
              mapping_parameters_.error_threshold, valid_candidate_starts,
              negative_read.data(), read_length, mapping_edit_distances,
              mapping_end_positions);
        }
        for (int mi = 0; mi < NUM_VPU_LANES_; ++mi) {
          if (mapping_edit_distances[mi] <=
              mapping_parameters_.error_threshold) {
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
              mappings.emplace_back((uint8_t)mapping_edit_distances[mi],
                                    valid_candidates[mi].position -
                                        mapping_parameters_.error_threshold +
                                        mapping_end_positions[mi]);
            } else {
              mappings.emplace_back((uint8_t)mapping_edit_distances[mi],
                                    valid_candidates[mi].position -
                                        read_length + 1 -
                                        mapping_parameters_.error_threshold +
                                        mapping_end_positions[mi]);
            }
          } else {
            candidate_count_threshold = valid_candidates[mi].count;
          }
        }
      }
      valid_candidate_index = 0;
    }
    ++candidate_index;
  }

  for (uint32_t ci = 0; ci < valid_candidate_index; ++ci) {
    uint32_t rid = valid_candidates[ci].position >> 32;
    uint32_t position = valid_candidates[ci].position;
    if (candidate_direction == kNegative) {
      position = position - read_length + 1;
    }

    if (!IsValidCandidate(rid, position, read_length, reference)) {
      continue;
    }

    int mapping_end_position;
    int num_errors;
    if (candidate_direction == kPositive) {
      num_errors =
          BandedAlignPatternToText(mapping_parameters_.error_threshold,
                                   reference.GetSequenceAt(rid) + position -
                                       mapping_parameters_.error_threshold,
                                   read, read_length, &mapping_end_position);
    } else {
      num_errors = BandedAlignPatternToText(
          mapping_parameters_.error_threshold,
          reference.GetSequenceAt(rid) + position -
              mapping_parameters_.error_threshold,
          negative_read.data(), read_length, &mapping_end_position);
    }
    if (num_errors <= mapping_parameters_.error_threshold) {
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
        mappings.emplace_back(num_errors,
                              valid_candidates[ci].position -
                                  mapping_parameters_.error_threshold +
                                  mapping_end_position);
      } else {
        mappings.emplace_back(num_errors,
                              valid_candidates[ci].position - read_length + 1 -
                                  mapping_parameters_.error_threshold +
                                  mapping_end_position);
      }
    }
  }
}

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::VerifyCandidatesOnOneDirection(
    Direction candidate_direction, uint32_t read_index,
    const SequenceBatch &read_batch, const SequenceBatch &reference,
    MappingMetadata &mapping_metadata) {
  const char *read = read_batch.GetSequenceAt(read_index);
  uint32_t read_length = read_batch.GetSequenceLengthAt(read_index);
  const std::string &negative_read =
      read_batch.GetNegativeSequenceAt(read_index);

  const std::vector<Candidate> &candidates =
      candidate_direction == kPositive ? mapping_metadata.positive_candidates_
                                       : mapping_metadata.negative_candidates_;
  std::vector<std::pair<int, uint64_t>> &mappings =
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
    if (candidates[ci].count < candidate_count_threshold) break;
    uint32_t rid = candidates[ci].position >> 32;
    uint32_t position = candidates[ci].position;
    if (candidate_direction == kNegative) {
      position = position - read_length + 1;
    }

    if (!IsValidCandidate(rid, position, read_length, reference)) {
      continue;
    }

    int mapping_end_position = read_length;
    int gap_beginning = 0;
    int num_errors;
    int allow_gap_beginning_ = 20;
    int mapping_length_threshold = 30;
    int allow_gap_beginning =
        allow_gap_beginning_ - mapping_parameters_.error_threshold;
    int actual_num_errors = 0;
    int read_mapping_length = 0;
    int best_mapping_longest_match = 0;
    int longest_match = 0;

    if (mapping_parameters_.split_alignment) {
      if (candidate_direction == kPositive) {
        num_errors = BandedAlignPatternToTextWithDropOff(
            mapping_parameters_.error_threshold,
            reference.GetSequenceAt(rid) + position -
                mapping_parameters_.error_threshold,
            read, read_length, &mapping_end_position, &read_mapping_length);
        if (mapping_end_position < 0 && allow_gap_beginning > 0) {
          int backup_num_errors = num_errors;
          int backup_mapping_end_position = -mapping_end_position;
          int backup_read_mapping_length = read_mapping_length;
          num_errors = BandedAlignPatternToTextWithDropOff(
              mapping_parameters_.error_threshold,
              reference.GetSequenceAt(rid) + position -
                  mapping_parameters_.error_threshold + allow_gap_beginning,
              read + allow_gap_beginning, read_length - allow_gap_beginning,
              &mapping_end_position, &read_mapping_length);
          if (num_errors > mapping_parameters_.error_threshold ||
              mapping_end_position < 0) {
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
            mapping_parameters_.error_threshold,
            reference.GetSequenceAt(rid) + position -
                mapping_parameters_.error_threshold,
            negative_read.data(), read_length, &mapping_end_position,
            &read_mapping_length);
        if (mapping_end_position < 0 && allow_gap_beginning > 0) {
          int backup_num_errors = num_errors;
          int backup_mapping_end_position = -mapping_end_position;
          int backup_read_mapping_length = read_mapping_length;
          num_errors = BandedAlignPatternToTextWithDropOffFrom3End(
              mapping_parameters_.error_threshold,
              reference.GetSequenceAt(rid) + position -
                  mapping_parameters_.error_threshold,
              negative_read.data(), read_length - allow_gap_beginning,
              &mapping_end_position, &read_mapping_length);
          if (num_errors > mapping_parameters_.error_threshold ||
              mapping_end_position < 0) {
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

      if (mapping_end_position + 1 - mapping_parameters_.error_threshold -
              num_errors - gap_beginning >=
          mapping_length_threshold) {
        actual_num_errors = num_errors;
        num_errors =
            -(mapping_end_position - mapping_parameters_.error_threshold -
              num_errors - gap_beginning);

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
        num_errors = mapping_parameters_.error_threshold + 1;
        actual_num_errors = mapping_parameters_.error_threshold + 1;
      }
    } else {
      if (candidate_direction == kPositive) {
        num_errors =
            BandedAlignPatternToText(mapping_parameters_.error_threshold,
                                     reference.GetSequenceAt(rid) + position -
                                         mapping_parameters_.error_threshold,
                                     read, read_length, &mapping_end_position);
      } else {
        num_errors = BandedAlignPatternToText(
            mapping_parameters_.error_threshold,
            reference.GetSequenceAt(rid) + position -
                mapping_parameters_.error_threshold,
            negative_read.data(), read_length, &mapping_end_position);
      }
    }

    if (num_errors <= mapping_parameters_.error_threshold) {
      if (num_errors < min_num_errors) {
        second_min_num_errors = min_num_errors;
        num_second_best_mappings = num_best_mappings;
        min_num_errors = num_errors;
        num_best_mappings = 1;
        if (mapping_parameters_.split_alignment) {
          if (candidates.size() > 50) {
            candidate_count_threshold = candidates[ci].count;
          } else {
            candidate_count_threshold = candidates[ci].count / 2;
          }
          if (second_min_num_errors <
                  min_num_errors + mapping_parameters_.error_threshold / 2 &&
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
        mappings.emplace_back(num_errors,
                              candidates[ci].position -
                                  mapping_parameters_.error_threshold +
                                  mapping_end_position);
      } else {
        if (mapping_parameters_.split_alignment &&
            mapping_parameters_.mapping_output_format != MAPPINGFORMAT_SAM) {
          mappings.emplace_back(num_errors,
                                candidates[ci].position - gap_beginning);
        } else {
          // Need to minus gap_beginning because mapping_end_position is
          // adjusted by it, but read_length is not.
          // printf("%d %d %d\n", candidates[ci].position, mapping_end_position,
          // gap_beginning);
          mappings.emplace_back(num_errors,
                                candidates[ci].position - read_length + 1 -
                                    mapping_parameters_.error_threshold +
                                    mapping_end_position);
        }
      }

      if (mapping_parameters_.split_alignment) {
        split_sites.emplace_back(((actual_num_errors & 0xff) << 24) |
                                 ((gap_beginning & 0xff) << 16) |
                                 (read_mapping_length & 0xffff));
      }
    }
  }
}

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::GenerateBestMappingsForSingleEndRead(
    const SequenceBatch &read_batch, uint32_t read_index,
    const SequenceBatch &reference, const SequenceBatch &barcode_batch,
    MappingMetadata &mapping_metadata,
    std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs) {
  int num_best_mappings = mapping_metadata.num_best_mappings_;

  // We will use reservoir sampling.
  std::vector<int> best_mapping_indices(
      mapping_parameters_.max_num_best_mappings);
  std::iota(best_mapping_indices.begin(), best_mapping_indices.end(), 0);
  if (num_best_mappings > mapping_parameters_.max_num_best_mappings) {
    std::mt19937 generator(11);
    for (int i = mapping_parameters_.max_num_best_mappings;
         i < num_best_mappings; ++i) {
      // important: inclusive range
      std::uniform_int_distribution<int> distribution(0, i);
      int j = distribution(generator);
      if (j < mapping_parameters_.max_num_best_mappings) {
        best_mapping_indices[j] = i;
      }
    }
    std::sort(best_mapping_indices.begin(), best_mapping_indices.end());
  }

  int best_mapping_index = 0;
  int num_best_mappings_reported = 0;

  ProcessBestMappingsForSingleEndRead(
      kPositive, read_index, read_batch, barcode_batch, reference,
      mapping_metadata, best_mapping_indices, best_mapping_index,
      num_best_mappings_reported, mappings_on_diff_ref_seqs);

  if (num_best_mappings_reported !=
      std::min(num_best_mappings, mapping_parameters_.max_num_best_mappings)) {
    ProcessBestMappingsForSingleEndRead(
        kNegative, read_index, read_batch, barcode_batch, reference,
        mapping_metadata, best_mapping_indices, best_mapping_index,
        num_best_mappings_reported, mappings_on_diff_ref_seqs);
  }
}

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::GenerateBestMappingsForPairedEndRead(
    uint32_t pair_index, const SequenceBatch &read_batch1,
    const SequenceBatch &read_batch2, const SequenceBatch &barcode_batch,
    const SequenceBatch &reference, std::vector<int> &best_mapping_indices,
    std::mt19937 &generator, int force_mapq,
    PairedEndMappingMetadata &paired_end_mapping_metadata,
    std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs) {
  int &min_sum_errors = paired_end_mapping_metadata.min_sum_errors_;
  int &num_best_mappings = paired_end_mapping_metadata.num_best_mappings_;
  int &second_min_sum_errors =
      paired_end_mapping_metadata.second_min_sum_errors_;
  int &num_second_best_mappings =
      paired_end_mapping_metadata.num_second_best_mappings_;

  min_sum_errors = 2 * mapping_parameters_.error_threshold + 1;
  num_best_mappings = 0;
  second_min_sum_errors = min_sum_errors;
  num_second_best_mappings = 0;

  GenerateBestMappingsForPairedEndReadOnOneDirection(
      kPositive, kNegative, pair_index, read_batch1, read_batch2, reference,
      paired_end_mapping_metadata);
  GenerateBestMappingsForPairedEndReadOnOneDirection(
      kNegative, kPositive, pair_index, read_batch1, read_batch2, reference,
      paired_end_mapping_metadata);

  if (mapping_parameters_.split_alignment) {
    GenerateBestMappingsForPairedEndReadOnOneDirection(
        kPositive, kPositive, pair_index, read_batch1, read_batch2, reference,
        paired_end_mapping_metadata);
    GenerateBestMappingsForPairedEndReadOnOneDirection(
        kNegative, kNegative, pair_index, read_batch1, read_batch2, reference,
        paired_end_mapping_metadata);
  }

  if (num_best_mappings <= mapping_parameters_.drop_repetitive_reads) {
    // We will use reservoir sampling.
    // std::vector<int>
    // best_mapping_indices(mapping_parameters_.max_num_best_mappings);
    std::iota(best_mapping_indices.begin(), best_mapping_indices.end(), 0);
    if (num_best_mappings > mapping_parameters_.max_num_best_mappings) {
      // std::mt19937 generator(11);
      for (int i = mapping_parameters_.max_num_best_mappings;
           i < num_best_mappings; ++i) {
        // Important: inclusive range.
        std::uniform_int_distribution<int> distribution(0, i);
        int j = distribution(generator);
        // int j = distribution(tmp_generator);
        if (j < mapping_parameters_.max_num_best_mappings) {
          best_mapping_indices[j] = i;
        }
      }
      std::sort(best_mapping_indices.begin(), best_mapping_indices.end());
    }

    int best_mapping_index = 0;
    int num_best_mappings_reported = 0;
    ProcessBestMappingsForPairedEndReadOnOneDirection(
        kPositive, kNegative, pair_index, read_batch1, read_batch2,
        barcode_batch, reference, best_mapping_indices, best_mapping_index,
        num_best_mappings_reported, force_mapq, paired_end_mapping_metadata,
        mappings_on_diff_ref_seqs);

    if (num_best_mappings_reported !=
        std::min(mapping_parameters_.max_num_best_mappings,
                 num_best_mappings)) {
      ProcessBestMappingsForPairedEndReadOnOneDirection(
          kNegative, kPositive, pair_index, read_batch1, read_batch2,
          barcode_batch, reference, best_mapping_indices, best_mapping_index,
          num_best_mappings_reported, force_mapq, paired_end_mapping_metadata,
          mappings_on_diff_ref_seqs);
    }

    if (mapping_parameters_.split_alignment &&
        num_best_mappings_reported !=
            std::min(mapping_parameters_.max_num_best_mappings,
                     num_best_mappings)) {
      ProcessBestMappingsForPairedEndReadOnOneDirection(
          kPositive, kPositive, pair_index, read_batch1, read_batch2,
          barcode_batch, reference, best_mapping_indices, best_mapping_index,
          num_best_mappings_reported, force_mapq, paired_end_mapping_metadata,
          mappings_on_diff_ref_seqs);
    }

    if (mapping_parameters_.split_alignment &&
        num_best_mappings_reported !=
            std::min(mapping_parameters_.max_num_best_mappings,
                     num_best_mappings)) {
      ProcessBestMappingsForPairedEndReadOnOneDirection(
          kNegative, kNegative, pair_index, read_batch1, read_batch2,
          barcode_batch, reference, best_mapping_indices, best_mapping_index,
          num_best_mappings_reported, force_mapq, paired_end_mapping_metadata,
          mappings_on_diff_ref_seqs);
    }
  }
}

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::ProcessBestMappingsForSingleEndRead(
    Direction mapping_direction, uint32_t read_index,
    const SequenceBatch &read_batch, const SequenceBatch &barcode_batch,
    const SequenceBatch &reference, const MappingMetadata &mapping_metadata,
    const std::vector<int> &best_mapping_indices, int &best_mapping_index,
    int &num_best_mappings_reported,
    std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs) {
  const std::vector<std::pair<int, uint64_t>> &mappings =
      mapping_direction == kPositive ? mapping_metadata.positive_mappings_
                                     : mapping_metadata.negative_mappings_;
  const std::vector<int> &split_sites =
      mapping_direction == kPositive ? mapping_metadata.positive_split_sites_
                                     : mapping_metadata.negative_split_sites_;

  const char *read = read_batch.GetSequenceAt(read_index);
  const uint32_t read_id = read_batch.GetSequenceIdAt(read_index);
  const char *read_name = read_batch.GetSequenceNameAt(read_index);
  const uint32_t read_length = read_batch.GetSequenceLengthAt(read_index);
  const std::string &negative_read =
      read_batch.GetNegativeSequenceAt(read_index);

  MappingInMemory mapping_in_memory;
  mapping_in_memory.read_id = read_id;
  mapping_in_memory.read_name = read_name;
  // const uint8_t is_unique = mapping_metadata.num_best_mappings_ == 1 ? 1 : 0;
  mapping_in_memory.is_unique = (mapping_metadata.num_best_mappings_ == 1);

  uint64_t barcode_key = 0;
  if (!mapping_parameters_.is_bulk_data) {
    barcode_key = barcode_batch.GenerateSeedFromSequenceAt(
        read_index, /*start_position=*/0,
        barcode_batch.GetSequenceLengthAt(read_index));
  }
  mapping_in_memory.barcode_key = barcode_key;

  mapping_in_memory.direction = mapping_direction;
  mapping_in_memory.read_sequence =
      mapping_direction == kPositive ? read : negative_read.data();
  mapping_in_memory.read_length = read_length;

  for (uint32_t mi = 0; mi < mappings.size(); ++mi) {
    if (mappings[mi].first > mapping_metadata.min_num_errors_) {
      continue;
    }

    if (best_mapping_index ==
        best_mapping_indices[num_best_mappings_reported]) {
      mapping_in_memory.rid = mappings[mi].second >> 32;

      if (mapping_parameters_.split_alignment) {
        mapping_in_memory.read_split_site = split_sites[mi];
      }

      GetRefStartEndPositionForReadFromMapping(mappings[mi], reference,
                                               mapping_in_memory);

      const uint16_t alignment_length = mapping_in_memory.ref_end_position -
                                        mapping_in_memory.ref_start_position +
                                        1;
      const uint8_t mapq = GetMAPQForSingleEndRead(
          mapping_direction, /*num_errors=*/mappings[mi].first,
          alignment_length, read_length,
          /*max_num_error_difference=*/mapping_parameters_.error_threshold,
          mapping_metadata);
      mapping_in_memory.mapq = mapq;

      if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_SAM) {
        uint16_t flag = mapping_direction == kPositive ? 0 : BAM_FREVERSE;
        if (num_best_mappings_reported >= 1) {
          flag |= BAM_FSECONDARY;
        }
        mapping_in_memory.SAM_flag = flag;
        mapping_in_memory.qual_sequence =
            read_batch.GetSequenceQualAt(read_index);
      }

      EmplaceBackSingleEndMappingRecord(mapping_in_memory,
                                        mappings_on_diff_ref_seqs);

      num_best_mappings_reported++;
      if (num_best_mappings_reported ==
          std::min(mapping_parameters_.max_num_best_mappings,
                   mapping_metadata.num_best_mappings_)) {
        break;
      }
    }

    best_mapping_index++;
  }
}

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::
    GenerateBestMappingsForPairedEndReadOnOneDirection(
        Direction first_read_direction, Direction second_read_direction,
        uint32_t pair_index, const SequenceBatch &read_batch1,
        const SequenceBatch &read_batch2, const SequenceBatch &reference,
        PairedEndMappingMetadata &paired_end_mapping_metadata) {
  uint32_t i1 = 0;
  uint32_t i2 = 0;
  uint32_t min_overlap_length = mapping_parameters_.min_read_length;
  uint32_t read1_length = read_batch1.GetSequenceLengthAt(pair_index);
  uint32_t read2_length = read_batch2.GetSequenceLengthAt(pair_index);

  const std::vector<std::pair<int, uint64_t>> &mappings1 =
      first_read_direction == kPositive
          ? paired_end_mapping_metadata.mapping_metadata1_.positive_mappings_
          : paired_end_mapping_metadata.mapping_metadata1_.negative_mappings_;
  const std::vector<std::pair<int, uint64_t>> &mappings2 =
      second_read_direction == kPositive
          ? paired_end_mapping_metadata.mapping_metadata2_.positive_mappings_
          : paired_end_mapping_metadata.mapping_metadata2_.negative_mappings_;

  std::vector<std::pair<uint32_t, uint32_t>> &best_mappings =
      paired_end_mapping_metadata.GetBestMappings(first_read_direction,
                                                  second_read_direction);
  int &min_sum_errors = paired_end_mapping_metadata.min_sum_errors_;
  int &num_best_mappings = paired_end_mapping_metadata.num_best_mappings_;
  int &second_min_sum_errors =
      paired_end_mapping_metadata.second_min_sum_errors_;
  int &num_second_best_mappings =
      paired_end_mapping_metadata.num_second_best_mappings_;

#ifdef LI_DEBUG
  for (int i = 0; i < mappings1.size(); ++i)
    printf("mappings1 %d %d:%d\n", i, (int)(mappings1[i].second >> 32),
           (int)mappings1[i].second);
  for (int i = 0; i < mappings1.size(); ++i)
    printf("mappings2 %d %d:%d\n", i, (int)(mappings2[i].second >> 32),
           (int)mappings2[i].second);
#endif

  if (mapping_parameters_.split_alignment) {
    if (mappings1.size() == 0 || mappings2.size() == 0) {
      return;
    }
    // For split alignment, selecting the pairs whose both single-end are the
    // best.
    for (i1 = 0; i1 < mappings1.size(); ++i1) {
      if (mappings1[i1].first !=
          paired_end_mapping_metadata.mapping_metadata1_.min_num_errors_) {
        continue;
      }
      for (i2 = 0; i2 < mappings2.size(); ++i2) {
        if (mappings2[i2].first !=
            paired_end_mapping_metadata.mapping_metadata2_.min_num_errors_) {
          continue;
        }
        best_mappings.emplace_back(i1, i2);
        min_sum_errors =
            paired_end_mapping_metadata.mapping_metadata1_.min_num_errors_ +
            paired_end_mapping_metadata.mapping_metadata2_.min_num_errors_;
        //*second_min_sum_errors = min_num_errors1 + min_num_errors2 + 1;
        num_best_mappings++;
      }
    }

    return;
  }

  while (i1 < mappings1.size() && i2 < mappings2.size()) {
    if ((first_read_direction == kNegative &&
         mappings1[i1].second > mappings2[i2].second +
                                    mapping_parameters_.max_insert_size -
                                    read1_length) ||
        (first_read_direction == kPositive &&
         mappings1[i1].second >
             mappings2[i2].second + read2_length - min_overlap_length)) {
      ++i2;
    } else if ((first_read_direction == kPositive &&
                mappings2[i2].second > mappings1[i1].second +
                                           mapping_parameters_.max_insert_size -
                                           read2_length) ||
               (first_read_direction == kNegative &&
                mappings2[i2].second >
                    mappings1[i1].second + read1_length - min_overlap_length)) {
      ++i1;
    } else {
      // ok, find a pair, we store current ni2 somewhere and keep looking until
      // we go out of the range, then we go back and then move to next pi1 and
      // keep doing the similar thing.
      uint32_t current_i2 = i2;
      while (current_i2 < mappings2.size() &&
             ((first_read_direction == kPositive &&
               mappings2[current_i2].second <=
                   mappings1[i1].second + mapping_parameters_.max_insert_size -
                       read2_length) ||
              (first_read_direction == kNegative &&
               mappings2[current_i2].second <=
                   mappings1[i1].second + read1_length - min_overlap_length))) {
#ifdef LI_DEBUG
        printf("%s passed: %llu %d %llu %d: %d %d %d\n", __func__,
               mappings1[i1].second >> 32, int(mappings1[i1].second),
               mappings2[current_i2].second >> 32,
               int(mappings2[current_i2].second),
               mappings1[i1].first + mappings2[current_i2].first,
               mappings1[i1].first, mappings2[current_i2].first);
#endif

        int current_sum_errors =
            mappings1[i1].first + mappings2[current_i2].first;
        if (current_sum_errors < min_sum_errors) {
          second_min_sum_errors = min_sum_errors;
          num_second_best_mappings = num_best_mappings;
          min_sum_errors = current_sum_errors;
          num_best_mappings = 1;
          best_mappings.clear();
          best_mappings.emplace_back(i1, current_i2);
        } else if (current_sum_errors == min_sum_errors) {
          num_best_mappings++;
          best_mappings.emplace_back(i1, current_i2);
        } else if (current_sum_errors == second_min_sum_errors) {
          num_second_best_mappings++;
        } else if (current_sum_errors < second_min_sum_errors) {
          second_min_sum_errors = current_sum_errors;
          num_second_best_mappings = 1;
        }
        ++current_i2;
      }
      ++i1;
    }
  }
}

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::
    ProcessBestMappingsForPairedEndReadOnOneDirection(
        Direction first_read_direction, Direction second_read_direction,
        uint32_t pair_index, const SequenceBatch &read_batch1,
        const SequenceBatch &read_batch2, const SequenceBatch &barcode_batch,
        const SequenceBatch &reference,
        const std::vector<int> &best_mapping_indices, int &best_mapping_index,
        int &num_best_mappings_reported, int force_mapq,
        const PairedEndMappingMetadata &paired_end_mapping_metadata,
        std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs) {
  PairedEndMappingInMemory paired_end_mapping_in_memory;

  paired_end_mapping_in_memory.mapping_in_memory1.direction =
      first_read_direction;
  paired_end_mapping_in_memory.mapping_in_memory2.direction =
      second_read_direction;

  const char *read1 = read_batch1.GetSequenceAt(pair_index);
  const char *read2 = read_batch2.GetSequenceAt(pair_index);
  const uint32_t read1_length = read_batch1.GetSequenceLengthAt(pair_index);
  const uint32_t read2_length = read_batch2.GetSequenceLengthAt(pair_index);
  const char *read1_name = read_batch1.GetSequenceNameAt(pair_index);
  const char *read2_name = read_batch2.GetSequenceNameAt(pair_index);
  const std::string &negative_read1 =
      read_batch1.GetNegativeSequenceAt(pair_index);
  const std::string &negative_read2 =
      read_batch2.GetNegativeSequenceAt(pair_index);
  const uint32_t read_id = read_batch1.GetSequenceIdAt(pair_index);

  paired_end_mapping_in_memory.mapping_in_memory1.read_id = read_id;
  paired_end_mapping_in_memory.mapping_in_memory2.read_id = read_id;

  paired_end_mapping_in_memory.mapping_in_memory1.read_name = read1_name;
  paired_end_mapping_in_memory.mapping_in_memory2.read_name = read2_name;

  paired_end_mapping_in_memory.mapping_in_memory1.read_length = read1_length;
  paired_end_mapping_in_memory.mapping_in_memory2.read_length = read2_length;

  const MappingMetadata &mapping_metadata1 =
      paired_end_mapping_metadata.mapping_metadata1_;
  const MappingMetadata &mapping_metadata2 =
      paired_end_mapping_metadata.mapping_metadata2_;

  const std::vector<std::pair<int, uint64_t>> &mappings1 =
      first_read_direction == kPositive ? mapping_metadata1.positive_mappings_
                                        : mapping_metadata1.negative_mappings_;
  const std::vector<std::pair<int, uint64_t>> &mappings2 =
      second_read_direction == kPositive ? mapping_metadata2.positive_mappings_
                                         : mapping_metadata2.negative_mappings_;

  const std::vector<int> &split_sites1 =
      first_read_direction == kPositive
          ? mapping_metadata1.positive_split_sites_
          : mapping_metadata1.negative_split_sites_;
  const std::vector<int> &split_sites2 =
      second_read_direction == kPositive
          ? mapping_metadata2.positive_split_sites_
          : mapping_metadata2.negative_split_sites_;

  const std::vector<std::pair<uint32_t, uint32_t>> &best_mappings =
      paired_end_mapping_metadata.GetBestMappings(first_read_direction,
                                                  second_read_direction);

  const uint8_t is_unique =
      (paired_end_mapping_metadata.num_best_mappings_ == 1 ||
       mapping_metadata1.num_best_mappings_ == 1 ||
       mapping_metadata2.num_best_mappings_ == 1)
          ? 1
          : 0;
  paired_end_mapping_in_memory.is_unique = is_unique;

  uint64_t barcode_key = 0;
  if (!mapping_parameters_.is_bulk_data) {
    barcode_key = barcode_batch.GenerateSeedFromSequenceAt(
        pair_index, /*start_position=*/0,
        barcode_batch.GetSequenceLengthAt(pair_index));
  }
  paired_end_mapping_in_memory.mapping_in_memory1.barcode_key = barcode_key;
  paired_end_mapping_in_memory.mapping_in_memory2.barcode_key = barcode_key;

  for (uint32_t mi = 0; mi < best_mappings.size(); ++mi) {
    const uint32_t i1 = best_mappings[mi].first;
    const uint32_t i2 = best_mappings[mi].second;
    const int current_sum_errors = mappings1[i1].first + mappings2[i2].first;

    if (current_sum_errors > paired_end_mapping_metadata.min_sum_errors_) {
      continue;
    }

    if (best_mapping_index ==
        best_mapping_indices[num_best_mappings_reported]) {
      const uint32_t rid1 = mappings1[i1].second >> 32;
      const uint32_t rid2 = mappings2[i2].second >> 32;

      paired_end_mapping_in_memory.mapping_in_memory1.rid = rid1;
      paired_end_mapping_in_memory.mapping_in_memory2.rid = rid2;

      paired_end_mapping_in_memory.mapping_in_memory1.read_sequence =
          first_read_direction == kPositive ? read1 : negative_read1.data();
      paired_end_mapping_in_memory.mapping_in_memory2.read_sequence =
          second_read_direction == kPositive ? read2 : negative_read2.data();

      if (mapping_parameters_.split_alignment) {
        paired_end_mapping_in_memory.mapping_in_memory1.read_split_site =
            split_sites1[i1];
        paired_end_mapping_in_memory.mapping_in_memory2.read_split_site =
            split_sites2[i2];
      }

      GetRefStartEndPositionForReadFromMapping(
          mappings1[i1], reference,
          paired_end_mapping_in_memory.mapping_in_memory1);
      GetRefStartEndPositionForReadFromMapping(
          mappings2[i2], reference,
          paired_end_mapping_in_memory.mapping_in_memory2);

      uint8_t mapq1 = 0;
      uint8_t mapq2 = 0;
      const uint8_t mapq = GetMAPQForPairedEndRead(
          first_read_direction, second_read_direction,
          /*read1_num_errors=*/mappings1[i1].first,
          /*read2_num_errors=*/mappings2[i2].first,
          paired_end_mapping_in_memory.mapping_in_memory1.GetFragmentLength(),
          paired_end_mapping_in_memory.mapping_in_memory2.GetFragmentLength(),
          read1_length, read2_length, force_mapq, paired_end_mapping_metadata,
          mapq1, mapq2);
      paired_end_mapping_in_memory.mapq = mapq;
      paired_end_mapping_in_memory.mapping_in_memory1.mapq = mapq;
      paired_end_mapping_in_memory.mapping_in_memory2.mapq = mapq;

      if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_SAM) {
        uint16_t flag1 = 3;
        uint16_t flag2 = 3;
        if (first_read_direction == kNegative) {
          flag1 |= BAM_FREVERSE;
          flag2 |= BAM_FMREVERSE;
        }
        if (second_read_direction == kNegative) {
          flag1 |= BAM_FMREVERSE;
          flag2 |= BAM_FREVERSE;
        }
        flag1 |= BAM_FREAD1;
        flag2 |= BAM_FREAD2;
        if (num_best_mappings_reported >= 1) {
          flag1 |= BAM_FSECONDARY;
          flag2 |= BAM_FSECONDARY;
        }
        paired_end_mapping_in_memory.mapping_in_memory1.SAM_flag = flag1;
        paired_end_mapping_in_memory.mapping_in_memory2.SAM_flag = flag2;
        paired_end_mapping_in_memory.mapping_in_memory1.qual_sequence =
            read_batch1.GetSequenceQualAt(pair_index);
        paired_end_mapping_in_memory.mapping_in_memory2.qual_sequence =
            read_batch2.GetSequenceQualAt(pair_index);
        paired_end_mapping_in_memory.mapping_in_memory1.is_unique = is_unique;
        paired_end_mapping_in_memory.mapping_in_memory2.is_unique = is_unique;
      }

      EmplaceBackPairedEndMappingRecord(paired_end_mapping_in_memory,
                                        mappings_on_diff_ref_seqs);

      num_best_mappings_reported++;
      if (num_best_mappings_reported ==
          std::min(mapping_parameters_.max_num_best_mappings,
                   paired_end_mapping_metadata.num_best_mappings_)) {
        break;
      }
    }

    best_mapping_index++;
  }
}

// The computed ref start and end coordinates are left closed and right closed.
template <typename MappingRecord>
void MappingGenerator<MappingRecord>::GetRefStartEndPositionForReadFromMapping(
    const std::pair<int, uint64_t> &mapping, const SequenceBatch &reference,
    MappingInMemory &mapping_in_memory) {
  // For now this mat is only used by ksw to generate mappings in SAM format.
  int8_t mat[25];
  // if (output_mapping_in_SAM_) {
  int i, j, k;
  for (i = k = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j)
      mat[k++] = i == j ? mapping_parameters_.match_score
                        : -mapping_parameters_.mismatch_penalty;
    mat[k++] = 0;  // ambiguous base
  }
  for (j = 0; j < 5; ++j) mat[k++] = 0;
  //}

  const uint32_t rid = mapping.second >> 32;
  const uint32_t ref_position = mapping.second;

  const int full_read_length = mapping_in_memory.read_length;
  int read_length = mapping_in_memory.read_length;

  const int min_num_errors = mapping.first;

  int split_site = mapping_in_memory.direction == kPositive
                       ? 0
                       : mapping_in_memory.read_length;

  int gap_beginning = 0;
  int actual_num_errors = 0;

  if (mapping_parameters_.split_alignment) {
    split_site = mapping_in_memory.read_split_site & 0xffff;
    // Beginning means the 5' end of the read.
    gap_beginning = (mapping_in_memory.read_split_site >> 16) & 0xff;
    // In split alignment, -num_errors is the number of matches.
    actual_num_errors = (mapping_in_memory.read_split_site >> 24) & 0xff;
    read_length = split_site - gap_beginning;
  }

  uint32_t verification_window_start_position =
      ref_position + 1 >
              (uint32_t)(read_length + mapping_parameters_.error_threshold)
          ? ref_position + 1 - read_length - mapping_parameters_.error_threshold
          : 0;

  if (ref_position + mapping_parameters_.error_threshold >=
      reference.GetSequenceLengthAt(rid)) {
    verification_window_start_position = reference.GetSequenceLengthAt(rid) -
                                         mapping_parameters_.error_threshold -
                                         read_length;
  }

  if (verification_window_start_position < 0) {
    verification_window_start_position = 0;
  }

  if (mapping_parameters_.split_alignment) {
    if (split_site < full_read_length &&
        mapping_parameters_.mapping_output_format == MAPPINGFORMAT_SAM &&
        split_site > 3 * mapping_parameters_.error_threshold) {
      split_site -= 3 * mapping_parameters_.error_threshold;
    }
    read_length = split_site - gap_beginning;
  }

  if (mapping_in_memory.direction == kPositive) {
    if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_SAM) {
      mapping_in_memory.n_cigar = 0;

      int mapping_start_position = 0;
      int mapping_end_position = 0;

      ksw_semi_global3(
          read_length + 2 * mapping_parameters_.error_threshold,
          reference.GetSequenceAt(rid) + verification_window_start_position,
          read_length, mapping_in_memory.read_sequence + gap_beginning, 5, mat,
          mapping_parameters_.gap_open_penalties[0],
          mapping_parameters_.gap_extension_penalties[0],
          mapping_parameters_.gap_open_penalties[1],
          mapping_parameters_.gap_extension_penalties[1],
          mapping_parameters_.error_threshold * 2 + 1,
          &(mapping_in_memory.n_cigar), &(mapping_in_memory.cigar),
          &mapping_start_position, &mapping_end_position);

      if (gap_beginning > 0) {
        int new_ref_start_position = AdjustGapBeginning(
            mapping_in_memory.direction, reference.GetSequenceAt(rid),
            mapping_in_memory.read_sequence, &gap_beginning, read_length - 1,
            verification_window_start_position + mapping_start_position,
            verification_window_start_position + mapping_end_position - 1,
            &(mapping_in_memory.n_cigar), &(mapping_in_memory.cigar));
        mapping_start_position =
            new_ref_start_position - verification_window_start_position;
      }

      GenerateNMAndMDTag(
          reference.GetSequenceAt(rid),
          mapping_in_memory.read_sequence + gap_beginning,
          verification_window_start_position + mapping_start_position,
          mapping_in_memory);

      mapping_in_memory.ref_start_position =
          verification_window_start_position + mapping_start_position;
      mapping_in_memory.ref_end_position =
          verification_window_start_position + mapping_end_position - 1;
    } else {
      int mapping_start_position = 0;
      if (!mapping_parameters_.split_alignment) {
        BandedTraceback(
            mapping_parameters_.error_threshold, min_num_errors,
            reference.GetSequenceAt(rid) + verification_window_start_position,
            mapping_in_memory.read_sequence, read_length,
            &mapping_start_position);
      } else {
        BandedTraceback(
            mapping_parameters_.error_threshold, actual_num_errors,
            reference.GetSequenceAt(rid) + verification_window_start_position,
            mapping_in_memory.read_sequence + gap_beginning, read_length,
            &mapping_start_position);
      }

      if (gap_beginning > 0) {
        int new_ref_start_position = AdjustGapBeginning(
            mapping_in_memory.direction, reference.GetSequenceAt(rid),
            mapping_in_memory.read_sequence, &gap_beginning, read_length - 1,
            verification_window_start_position + mapping_start_position,
            ref_position, &(mapping_in_memory.n_cigar),
            &(mapping_in_memory.cigar));

        mapping_start_position =
            new_ref_start_position - verification_window_start_position;
      }

      mapping_in_memory.ref_start_position =
          verification_window_start_position + mapping_start_position;
      mapping_in_memory.ref_end_position = ref_position;
    }

    return;
  }

  //  reversed read looks like:
  //
  //      veri_start_pos       ref_position
  //  ref   --|-------------------|------------------->
  //  read     <-|--read_length---|--gap_beginning--
  //          split_site
  //

  const int read_start_site = full_read_length - split_site;
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_SAM) {
    mapping_in_memory.n_cigar = 0;

    int mapping_start_position = 0;
    int mapping_end_position = 0;

    ksw_semi_global3(read_length + 2 * mapping_parameters_.error_threshold,
                     reference.GetSequenceAt(rid) +
                         verification_window_start_position + read_start_site,
                     read_length,
                     mapping_in_memory.read_sequence + read_start_site, 5, mat,
                     mapping_parameters_.gap_open_penalties[0],
                     mapping_parameters_.gap_extension_penalties[0],
                     mapping_parameters_.gap_open_penalties[1],
                     mapping_parameters_.gap_extension_penalties[1],
                     mapping_parameters_.error_threshold * 2 + 1,
                     &(mapping_in_memory.n_cigar), &(mapping_in_memory.cigar),
                     &mapping_start_position, &mapping_end_position);

    if (gap_beginning > 0) {
      int new_ref_end_position = AdjustGapBeginning(
          mapping_in_memory.direction, reference.GetSequenceAt(rid),
          mapping_in_memory.read_sequence + read_start_site, &gap_beginning,
          read_length - 1,
          verification_window_start_position + mapping_start_position,
          verification_window_start_position + mapping_end_position - 1,
          &(mapping_in_memory.n_cigar), &(mapping_in_memory.cigar));

      // The returned position is right-closed, so need to plus one to match
      // bed convention
      mapping_end_position = new_ref_end_position + 1 -
                             verification_window_start_position -
                             read_start_site;
      read_length = split_site - gap_beginning;
    }

    GenerateNMAndMDTag(reference.GetSequenceAt(rid),
                       mapping_in_memory.read_sequence + read_start_site,
                       verification_window_start_position + read_start_site +
                           mapping_start_position,
                       mapping_in_memory);

    mapping_in_memory.ref_start_position = verification_window_start_position +
                                           read_start_site +
                                           mapping_start_position;
    mapping_in_memory.ref_end_position = verification_window_start_position +
                                         read_start_site +
                                         mapping_end_position - 1;
  } else {
    int mapping_start_position = mapping_parameters_.error_threshold;
    int mapping_end_position =
        ref_position - verification_window_start_position + 1;
    // int n_cigar = 0;
    // uint32_t *cigar;
    // ksw_semi_global3(read_length + 2 * mapping_parameters_.error_threshold,
    // reference.GetSequenceAt(rid) + verification_window_start_position,
    // read_length, negative_read.data() + split_sites[mi], 5, mat,
    // mapping_parameters_.gap_open_penalties[0],
    // mapping_parameters_.gap_extension_penalties[0],
    // mapping_parameters_.gap_open_penalties[1],
    // mapping_parameters_.gap_extension_penalties[1],
    // mapping_parameters_.error_threshold * 2 + 1, &n_cigar, &cigar,
    // &mapping_start_position, &mapping_end_position); mapq =
    // GetMAPQForSingleEndRead(mapping_parameters_.error_threshold, 0, 0,
    // mapping_end_position - mapping_start_position + 1, min_num_errors,
    // num_best_mappings, second_min_num_errors, num_second_best_mappings);
    // uint32_t fragment_start_position = verification_window_start_position +
    // mapping_start_position; uint16_t fragment_length = mapping_end_position
    // - mapping_start_position + 1;
    if (!mapping_parameters_.split_alignment) {
      BandedTraceback(
          mapping_parameters_.error_threshold, min_num_errors,
          reference.GetSequenceAt(rid) + verification_window_start_position,
          mapping_in_memory.read_sequence + read_start_site, read_length,
          &mapping_start_position);
    } else {
      // BandedTracebackToEnd(mapping_parameters_.error_threshold,actual_num_errors,
      // reference.GetSequenceAt(rid)
      // + verification_window_start_position, read + read_start_site,
      // read_length, &mapping_end_position);
      BandedAlignPatternToText(
          mapping_parameters_.error_threshold,
          reference.GetSequenceAt(rid) + verification_window_start_position,
          mapping_in_memory.read_sequence + read_start_site, read_length,
          &mapping_end_position);
      // seems banded align's mapping end position is included?
      mapping_end_position += 1;
    }

    if (gap_beginning > 0) {
      int new_ref_end_position = AdjustGapBeginning(
          mapping_in_memory.direction, reference.GetSequenceAt(rid),
          mapping_in_memory.read_sequence + read_start_site, &gap_beginning,
          read_length - 1,
          verification_window_start_position + mapping_start_position,
          verification_window_start_position + mapping_end_position,
          &(mapping_in_memory.n_cigar), &(mapping_in_memory.cigar));

      // The returned position is right-closed, so need to plus one to match
      // bed convention.
      mapping_end_position =
          new_ref_end_position - verification_window_start_position + 1;
      read_length = split_site - gap_beginning;
    }

    mapping_in_memory.ref_start_position =
        verification_window_start_position + mapping_start_position;
    mapping_in_memory.ref_end_position =
        verification_window_start_position + mapping_end_position - 1;
  }
}

template <typename MappingRecord>
uint8_t MappingGenerator<MappingRecord>::GetMAPQForSingleEndRead(
    Direction direction, int num_errors, uint16_t alignment_length,
    int read_length, int max_num_error_difference,
    const MappingMetadata &mapping_metadata) {
  int mapq_coef_length = 50;
  int mapq_coef_fraction = log(mapq_coef_length);

  if (!mapping_parameters_.split_alignment) {
    alignment_length =
        alignment_length > read_length ? alignment_length : read_length;
  }

  double alignment_identity = 1 - (double)num_errors / alignment_length;

  if (mapping_parameters_.split_alignment) {
    alignment_identity = (double)(-num_errors) / alignment_length;
    if (alignment_identity > 1) alignment_identity = 1;
  }

  int mapq = 0;
  int second_min_num_errors = mapping_metadata.second_min_num_errors_;

  if (mapping_metadata.num_best_mappings_ > 1) {
    // mapq = -4.343 * log(1 - 1.0 / num_best_mappings);
    // if (num_best_mappings == 2) {
    //  mapq = 3;
    //} else if (num_best_mappings == 3) {
    //  mapq = 2;
    //} else if (num_best_mappings < 10) {
    //  mapq = 1;
    //} else {
    //  mapq = 0;
    //}
  } else {
    if (second_min_num_errors > num_errors + max_num_error_difference) {
      second_min_num_errors = num_errors + max_num_error_difference;
    }

    double tmp = alignment_length < mapq_coef_length
                     ? 1.0
                     : mapq_coef_fraction / log(alignment_length);
    tmp *= alignment_identity * alignment_identity;
    mapq = 5 * 6.02 * (second_min_num_errors - num_errors) * tmp * tmp + 0.499;
  }

  if (mapping_metadata.num_second_best_mappings_ > 0) {
    mapq -= (int)(4.343 * log(mapping_metadata.num_second_best_mappings_ + 1) +
                  0.499);
  }

  if (mapq > 60) {
    mapq = 60;
  }
  if (mapq < 0) {
    mapq = 0;
  }

  if (mapping_metadata.repetitive_seed_length_ > 0) {
    double frac_rep =
        (mapping_metadata.repetitive_seed_length_) / (double)read_length;
    if (mapping_metadata.repetitive_seed_length_ >= (uint32_t)read_length) {
      frac_rep = 0.999;
    }
    if (alignment_identity <= 0.95) {
      mapq = mapq * (1 - sqrt(frac_rep)) + 0.499;
    } else if (alignment_identity <= 0.97) {
      mapq = mapq * (1 - frac_rep) + 0.499;
    } else if (alignment_identity >= 0.999) {
      mapq = mapq * (1 - frac_rep * frac_rep * frac_rep * frac_rep) + 0.499;
    } else {
      mapq = mapq * (1 - frac_rep * frac_rep) + 0.499;
    }
  }

  if (mapping_parameters_.split_alignment &&
      alignment_length < read_length - mapping_parameters_.error_threshold &&
      second_min_num_errors != num_errors) {
    if (mapping_metadata.repetitive_seed_length_ >= alignment_length &&
        mapping_metadata.repetitive_seed_length_ < (uint32_t)read_length &&
        alignment_length < read_length / 3) {
      mapq = 0;
    }
    const int diff = second_min_num_errors - num_errors;
    const uint32_t num_candidates =
        direction == kPositive ? mapping_metadata.positive_candidates_.size()
                               : mapping_metadata.negative_candidates_.size();
    if (second_min_num_errors - num_errors <=
            mapping_parameters_.error_threshold * 3 / 4 &&
        num_candidates >= 5) {
      mapq -= (num_candidates / 5 / diff);
    }
    if (mapq < 0) {
      mapq = 0;
    }
    if (mapping_metadata.num_second_best_mappings_ > 0 &&
        second_min_num_errors - num_errors <=
            mapping_parameters_.error_threshold * 3 / 4) {
      mapq /= (mapping_metadata.num_second_best_mappings_ / diff + 1);
    }
  }

  return (uint8_t)mapq;
}

#define raw_mapq(diff, a) ((int)(5 * 6.02 * (diff) / (a) + .499))

template <typename MappingRecord>
uint8_t MappingGenerator<MappingRecord>::GetMAPQForPairedEndRead(
    Direction first_read_direction, Direction second_read_direction,
    int read1_num_errors, int read2_num_errors, uint16_t read1_alignment_length,
    uint16_t read2_alignment_length, int read1_length, int read2_length,
    int force_mapq, const PairedEndMappingMetadata &paired_end_mapping_metadata,
    uint8_t &mapq1, uint8_t &mapq2) {
  const MappingMetadata &mapping_metadata1 =
      paired_end_mapping_metadata.mapping_metadata1_;
  const MappingMetadata &mapping_metadata2 =
      paired_end_mapping_metadata.mapping_metadata2_;

#ifdef CHROMAP_DEBUG
  std::cerr
      << " rl1:"
      << (int)paired_end_mapping_metadata.mapping_metadata1_
             .repetitive_seed_length_
      << " rl2:"
      << (int)paired_end_mapping_metadata.mapping_metadata2_
             .repetitive_seed_length_
      << " pal:" << (int)read1_alignment_length
      << " nal:" << (int)read2_alignment_length
      << " me:" << paired_end_mapping_metadata.min_sum_errors_
      << " #bm:" << paired_end_mapping_metadata.num_best_mappings_
      << " sme:" << paired_end_mapping_metadata.second_min_sum_errors_
      << " #sbm:" << paired_end_mapping_metadata.num_second_best_mappings_
      << " ne1:" << read1_num_errors << " ne2:" << read2_num_errors << " me1:"
      << paired_end_mapping_metadata.mapping_metadata1_.min_num_errors_
      << " me2:"
      << paired_end_mapping_metadata.mapping_metadata2_.min_num_errors_
      << " #bm1:"
      << paired_end_mapping_metadata.mapping_metadata1_.num_best_mappings_
      << " #bm2:"
      << paired_end_mapping_metadata.mapping_metadata2_.num_best_mappings_
      << " sme1:"
      << paired_end_mapping_metadata.mapping_metadata1_.second_min_num_errors_
      << " sme2:"
      << paired_end_mapping_metadata.mapping_metadata2_.second_min_num_errors_
      << " #sbm1:"
      << paired_end_mapping_metadata.mapping_metadata1_
             .num_second_best_mappings_
      << " #sbm2:"
      << paired_end_mapping_metadata.mapping_metadata2_
             .num_second_best_mappings_
      << "\n";
#endif

  uint8_t mapq_pe = 0;
  int min_num_unpaired_sum_errors =
      mapping_metadata1.min_num_errors_ + mapping_metadata2.min_num_errors_ + 3;

  if (paired_end_mapping_metadata.num_best_mappings_ <= 1) {
    int adjusted_second_min_sum_errors =
        paired_end_mapping_metadata.second_min_sum_errors_ <
                min_num_unpaired_sum_errors
            ? paired_end_mapping_metadata.second_min_sum_errors_
            : min_num_unpaired_sum_errors;

    mapq_pe = raw_mapq(adjusted_second_min_sum_errors -
                           paired_end_mapping_metadata.min_sum_errors_,
                       1);

#ifdef CHROMAP_DEBUG
    std::cerr << "mapqpe: " << (int)mapq_pe << "\n";
#endif

    if (paired_end_mapping_metadata.num_second_best_mappings_ > 0) {
      mapq_pe -=
          (int)(4.343 *
                    log(paired_end_mapping_metadata.num_second_best_mappings_ +
                        1) +
                0.499);
    }

    if (mapq_pe > 60) {
      mapq_pe = 60;
    }
    if (mapq_pe < 0) {
      mapq_pe = 0;
    }

#ifdef CHROMAP_DEBUG
    std::cerr << "mapqpe: " << (int)mapq_pe << "\n";
#endif

    int repetitive_seed_length = mapping_metadata1.repetitive_seed_length_ +
                                 mapping_metadata2.repetitive_seed_length_;

    if (repetitive_seed_length > 0) {
      double total_read_length = read1_length + read2_length;
      double frac_rep = (double)repetitive_seed_length / total_read_length;
      if (repetitive_seed_length >= total_read_length) {
        frac_rep = 0.999;
      }

      double alignment_identity1 =
          1 - (double)read1_num_errors / (read1_length > read1_alignment_length
                                              ? read1_length
                                              : read1_alignment_length);

      double alignment_identity2 =
          1 - (double)read2_num_errors / (read2_length > read2_alignment_length
                                              ? read2_length
                                              : read2_alignment_length);

      double alignment_identity = alignment_identity1 < alignment_identity2
                                      ? alignment_identity1
                                      : alignment_identity2;

      if (alignment_identity <= 0.95) {
        mapq_pe = mapq_pe * (1 - sqrt(frac_rep)) + 0.499;
      } else if (alignment_identity <= 0.97) {
        mapq_pe = mapq_pe * (1 - frac_rep) + 0.499;
      } else if (alignment_identity >= 0.999) {
        mapq_pe =
            mapq_pe * (1 - frac_rep * frac_rep * frac_rep * frac_rep) + 0.499;
      } else {
        mapq_pe = mapq_pe * (1 - frac_rep * frac_rep) + 0.499;
      }
    }
  }

  mapq1 = GetMAPQForSingleEndRead(first_read_direction, read1_num_errors,
                                  read1_alignment_length, read1_length,
                                  /*max_num_error_difference=*/2,
                                  mapping_metadata1);

  mapq2 = GetMAPQForSingleEndRead(second_read_direction, read2_num_errors,
                                  read2_alignment_length, read2_length,
                                  /*max_num_error_difference=*/2,
                                  mapping_metadata2);

#ifdef CHROMAP_DEBUG
  std::cerr << " 1:" << (int)mapq1 << " 2:" << (int)mapq2
            << " mapq_pe:" << (int)mapq_pe << "\n";
#endif

  if (!mapping_parameters_.split_alignment) {
    mapq1 = mapq1 > mapq_pe                    ? mapq1
            : mapq_pe < mapq1 + mapq_pe * 0.65 ? mapq_pe
                                               : mapq1 + mapq_pe * 0.65;
    mapq2 = mapq2 > mapq_pe                    ? mapq2
            : mapq_pe < mapq2 + mapq_pe * 0.65 ? mapq_pe
                                               : mapq2 + mapq_pe * 0.65;
  }

  mapq1 *= 1.2;
  if (mapq1 > 60) {
    mapq1 = 60;
  }

  mapq2 *= 1.2;
  if (mapq2 > 60) {
    mapq2 = 60;
  }

#ifdef CHROMAP_DEBUG
  std::cerr << " 1:" << (int)mapq1 << " 2:" << (int)mapq2 << "\n\n";
#endif

  uint8_t mapq = mapq1 < mapq2 ? mapq1 : mapq2;

  if (mapq < 60 && force_mapq >= 0 && force_mapq < mapq) {
    mapq = force_mapq;
  }

  return (uint8_t)mapq;
}

}  // namespace chromap

#endif  // MAPPING_GENERATOR_H_
