#include "candidate_processor.h"

#include <cinttypes>
#include <cstring>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

namespace chromap {

void CandidateProcessor::GenerateCandidates(
    int error_threshold, const Index &index,
    MappingMetadata &mapping_metadata) const {
  const std::vector<Minimizer> &minimizers = mapping_metadata.minimizers_;
  std::vector<uint64_t> &positive_hits = mapping_metadata.positive_hits_;
  std::vector<uint64_t> &negative_hits = mapping_metadata.negative_hits_;
  std::vector<Candidate> &positive_candidates =
      mapping_metadata.positive_candidates_;
  std::vector<Candidate> &negative_candidates =
      mapping_metadata.negative_candidates_;
  uint32_t &repetitive_seed_length = mapping_metadata.repetitive_seed_length_;

  const CandidatePositionGeneratingConfig first_round_generating_config(
      /*max_seed_frequency=*/max_seed_frequencies_[0],
      /*repetitive_seed_frequency=*/max_seed_frequencies_[0],
      /*use_heap_merge=*/false);

  repetitive_seed_length = 0;
  int repetitive_seed_count = index.GenerateCandidatePositions(
      first_round_generating_config, mapping_metadata);

  bool use_high_frequency_minimizers = false;
  if (positive_hits.size() + negative_hits.size() == 0) {
    positive_hits.clear();
    negative_hits.clear();
    repetitive_seed_length = 0;

    const CandidatePositionGeneratingConfig second_round_generating_config(
        /*max_seed_frequency=*/max_seed_frequencies_[1],
        /*repetitive_seed_frequency=*/max_seed_frequencies_[0],
        /*use_heap_merge=*/true);

    repetitive_seed_count = index.GenerateCandidatePositions(
        second_round_generating_config, mapping_metadata);
    use_high_frequency_minimizers = true;
    if (positive_hits.size() == 0 || negative_hits.size() == 0) {
      use_high_frequency_minimizers = false;
    }
  }

  int num_required_seeds = minimizers.size() - repetitive_seed_count;
  num_required_seeds = num_required_seeds > 1 ? num_required_seeds : 1;
  num_required_seeds = num_required_seeds > min_num_seeds_required_for_mapping_
                           ? min_num_seeds_required_for_mapping_
                           : num_required_seeds;
  if (use_high_frequency_minimizers) {
    num_required_seeds = min_num_seeds_required_for_mapping_;
  }

  // std::cerr << "Normal positive gen on one dir\n";
  GenerateCandidatesOnOneStrand(error_threshold, num_required_seeds,
                                minimizers.size(), positive_hits,
                                positive_candidates);
  // std::cerr << "Normal negative gen on one dir\n";
  GenerateCandidatesOnOneStrand(error_threshold, num_required_seeds,
                                minimizers.size(), negative_hits,
                                negative_candidates);
  // fprintf(stderr, "p+n: %d\n", positive_candidates->size() +
  // negative_candidates->size()) ;
}

// Return 0 if it supplements normally. Return 1 if the supplement could be too
// aggressive, and MAPQ needs setting to 0.
int CandidateProcessor::SupplementCandidates(
    int error_threshold, uint32_t search_range, const Index &index,
    PairedEndMappingMetadata &paired_end_mapping_metadata) const {
  std::vector<Candidate> augment_positive_candidates1;
  std::vector<Candidate> augment_positive_candidates2;
  std::vector<Candidate> augment_negative_candidates1;
  std::vector<Candidate> augment_negative_candidates2;

  int ret = 0;

  for (int mate = 0; mate <= 1; ++mate) {
    std::vector<Minimizer> *minimizers;
    std::vector<uint64_t> *positive_hits;
    std::vector<uint64_t> *negative_hits;
    std::vector<Candidate> *positive_candidates;
    std::vector<Candidate> *negative_candidates;
    std::vector<Candidate> *mate_positive_candidates;
    std::vector<Candidate> *mate_negative_candidates;
    std::vector<Candidate> *augment_positive_candidates;
    std::vector<Candidate> *augment_negative_candidates;
    uint32_t *repetitive_seed_length;

    if (mate == 0) {
      minimizers = &paired_end_mapping_metadata.mapping_metadata1_.minimizers_;
      positive_hits =
          &paired_end_mapping_metadata.mapping_metadata1_.positive_hits_;
      negative_hits =
          &paired_end_mapping_metadata.mapping_metadata1_.negative_hits_;
      positive_candidates =
          &paired_end_mapping_metadata.mapping_metadata1_.positive_candidates_;
      negative_candidates =
          &paired_end_mapping_metadata.mapping_metadata1_.negative_candidates_;
      mate_positive_candidates =
          &paired_end_mapping_metadata.mapping_metadata2_.positive_candidates_;
      mate_negative_candidates =
          &paired_end_mapping_metadata.mapping_metadata2_.negative_candidates_;
      augment_positive_candidates = &augment_positive_candidates1;
      augment_negative_candidates = &augment_negative_candidates1;
      repetitive_seed_length = &paired_end_mapping_metadata.mapping_metadata1_
                                    .repetitive_seed_length_;
    } else {
      minimizers = &paired_end_mapping_metadata.mapping_metadata2_.minimizers_;
      positive_hits =
          &paired_end_mapping_metadata.mapping_metadata2_.positive_hits_;
      negative_hits =
          &paired_end_mapping_metadata.mapping_metadata2_.negative_hits_;
      positive_candidates =
          &paired_end_mapping_metadata.mapping_metadata2_.positive_candidates_;
      negative_candidates =
          &paired_end_mapping_metadata.mapping_metadata2_.negative_candidates_;
      mate_positive_candidates =
          &paired_end_mapping_metadata.mapping_metadata1_.positive_candidates_;
      mate_negative_candidates =
          &paired_end_mapping_metadata.mapping_metadata1_.negative_candidates_;
      augment_positive_candidates = &augment_positive_candidates2;
      augment_negative_candidates = &augment_negative_candidates2;
      repetitive_seed_length = &paired_end_mapping_metadata.mapping_metadata2_
                                    .repetitive_seed_length_;
    }

    uint32_t mm_count = minimizers->size();
    bool augment_flag = true;
    uint32_t candidate_num = positive_candidates->size();

    for (uint32_t i = 0; i < candidate_num; ++i) {
      if ((*positive_candidates)[i].count >= mm_count / 2) {
        augment_flag = false;
        break;
      }
    }

    candidate_num = negative_candidates->size();
    if (augment_flag) {
      for (uint32_t i = 0; i < candidate_num; ++i) {
        if ((*negative_candidates)[i].count >= mm_count / 2) {
          augment_flag = false;
          break;
        }
      }
    }

    if (augment_flag) {
      positive_hits->clear();
      negative_hits->clear();
      positive_hits->reserve(max_seed_frequencies_[0]);
      negative_hits->reserve(max_seed_frequencies_[0]);
      int positive_rescue_result = 0;
      int negative_rescue_result = 0;
      if (mate_positive_candidates->size() > 0) {
        positive_rescue_result =
            GenerateCandidatesFromRepetitiveReadWithMateInfoOnOneStrand(
                kNegative, search_range, error_threshold, index, *minimizers,
                *mate_positive_candidates, *repetitive_seed_length,
                *negative_hits, *augment_negative_candidates);
      }

      if (mate_negative_candidates->size() > 0) {
        negative_rescue_result =
            GenerateCandidatesFromRepetitiveReadWithMateInfoOnOneStrand(
                kPositive, search_range, error_threshold, index, *minimizers,
                *mate_negative_candidates, *repetitive_seed_length,
                *positive_hits, *augment_positive_candidates);
      }

      // If one of the strand did not supplement due to too many best candidate,
      // and the filtered strand have better best candidates,
      // and there is no candidate directly from minimizers,
      // then we remove the supplement
      if (((positive_rescue_result < 0 && negative_rescue_result > 0 &&
            -positive_rescue_result >= negative_rescue_result) ||
           (positive_rescue_result > 0 && negative_rescue_result < 0 &&
            positive_rescue_result <= -negative_rescue_result)) &&
          positive_candidates->size() + negative_candidates->size() == 0) {
        // augment_positive_candidates->clear();
        // augment_negative_candidates->clear();
        ret = 1;
      }
    }
  }

  if (augment_positive_candidates1.size() > 0) {
    MergeCandidates(
        error_threshold,
        paired_end_mapping_metadata.mapping_metadata1_.positive_candidates_,
        augment_positive_candidates1,
        paired_end_mapping_metadata.mapping_metadata1_
            .positive_candidates_buffer_);
  }

  if (augment_negative_candidates1.size() > 0) {
    MergeCandidates(
        error_threshold,
        paired_end_mapping_metadata.mapping_metadata1_.negative_candidates_,
        augment_negative_candidates1,
        paired_end_mapping_metadata.mapping_metadata1_
            .negative_candidates_buffer_);
  }

  if (augment_positive_candidates2.size() > 0) {
    MergeCandidates(
        error_threshold,
        paired_end_mapping_metadata.mapping_metadata2_.positive_candidates_,
        augment_positive_candidates2,
        paired_end_mapping_metadata.mapping_metadata2_
            .positive_candidates_buffer_);
  }

  if (augment_negative_candidates2.size() > 0) {
    MergeCandidates(
        error_threshold,
        paired_end_mapping_metadata.mapping_metadata2_.negative_candidates_,
        augment_negative_candidates2,
        paired_end_mapping_metadata.mapping_metadata2_
            .negative_candidates_buffer_);
  }
  return ret;
}

void CandidateProcessor::ReduceCandidatesForPairedEndRead(
    uint32_t mapping_positions_distance,
    PairedEndMappingMetadata &paired_end_mapping_metadata) const {
  const std::vector<Candidate> &positive_candidates1 =
      paired_end_mapping_metadata.mapping_metadata1_
          .positive_candidates_buffer_;
  const std::vector<Candidate> &negative_candidates1 =
      paired_end_mapping_metadata.mapping_metadata1_
          .negative_candidates_buffer_;
  const std::vector<Candidate> &positive_candidates2 =
      paired_end_mapping_metadata.mapping_metadata2_
          .positive_candidates_buffer_;
  const std::vector<Candidate> &negative_candidates2 =
      paired_end_mapping_metadata.mapping_metadata2_
          .negative_candidates_buffer_;
  std::vector<Candidate> &filtered_positive_candidates1 =
      paired_end_mapping_metadata.mapping_metadata1_.positive_candidates_;
  std::vector<Candidate> &filtered_negative_candidates1 =
      paired_end_mapping_metadata.mapping_metadata1_.negative_candidates_;
  std::vector<Candidate> &filtered_positive_candidates2 =
      paired_end_mapping_metadata.mapping_metadata2_.positive_candidates_;
  std::vector<Candidate> &filtered_negative_candidates2 =
      paired_end_mapping_metadata.mapping_metadata2_.negative_candidates_;

  ReduceCandidatesForPairedEndReadOnOneDirection(
      mapping_positions_distance, positive_candidates1, negative_candidates2,
      filtered_positive_candidates1, filtered_negative_candidates2);
  ReduceCandidatesForPairedEndReadOnOneDirection(
      mapping_positions_distance, negative_candidates1, positive_candidates2,
      filtered_negative_candidates1, filtered_positive_candidates2);
}

int CandidateProcessor::
    GenerateCandidatesFromRepetitiveReadWithMateInfoOnOneStrand(
        const Strand strand, uint32_t search_range, int error_threshold,
        const Index &index, const std::vector<Minimizer> &minimizers,
        const std::vector<Candidate> &mate_candidates,
        uint32_t &repetitive_seed_length, std::vector<uint64_t> &hits,
        std::vector<Candidate> &candidates) const {
  int max_seed_count =
      index.GenerateCandidatePositionsFromRepetitiveReadWithMateInfoOnOneStrand(
          strand, search_range, min_num_seeds_required_for_mapping_,
          max_seed_frequencies_[0], error_threshold, minimizers,
          mate_candidates, repetitive_seed_length, hits);

  GenerateCandidatesOnOneStrand(error_threshold, /*num_seeds_required=*/1,
                                minimizers.size(), hits, candidates);
  return max_seed_count;
}

void CandidateProcessor::GenerateCandidatesOnOneStrand(
    int error_threshold, int num_seeds_required, uint32_t num_minimizers,
    std::vector<uint64_t> &hits, std::vector<Candidate> &candidates) const {
  hits.emplace_back(UINT64_MAX);
  if (hits.size() > 0) {
    int minimizer_count = 1;
    // The number of seeds with the exact same reference position.
    int equal_count = 1;
    int best_equal_count = 1;
    uint64_t previous_hit = hits[0];
    uint32_t previous_reference_id = previous_hit >> 32;
    uint32_t previous_reference_position = previous_hit;
    uint64_t best_local_hit = hits[0];
    for (uint32_t pi = 1; pi < hits.size(); ++pi) {
      uint32_t current_reference_id = hits[pi] >> 32;
      uint32_t current_reference_position = hits[pi];
#ifdef LI_DEBUG
      printf("%s: %d %d\n", __func__, current_reference_id,
             current_reference_position);
#endif
      if (current_reference_id != previous_reference_id ||
          current_reference_position >
              previous_reference_position + error_threshold ||
          ((uint32_t)minimizer_count >= num_minimizers &&
           current_reference_position >
               (uint32_t)best_local_hit + error_threshold)) {
        if (minimizer_count >= num_seeds_required) {
          Candidate candidate;
          candidate.position = best_local_hit;
          candidate.count = best_equal_count;
          candidates.push_back(candidate);
        }

        minimizer_count = 1;
        equal_count = 1;
        best_equal_count = 1;
        best_local_hit = hits[pi];
      } else {
        if (hits[pi] == best_local_hit) {
          ++equal_count;
          ++best_equal_count;
        } else if (hits[pi] == previous_hit) {
          ++equal_count;
          if (equal_count > best_equal_count) {
            best_local_hit = previous_hit;
            best_equal_count = equal_count;
          }
        } else {
          equal_count = 1;
        }

        ++minimizer_count;
      }

      previous_hit = hits[pi];
      previous_reference_id = current_reference_id;
      previous_reference_position = current_reference_position;
    }
  }
}

// Merge c1 and c2 into buffer and then swap the results into c1.
void CandidateProcessor::MergeCandidates(int error_threshold,
                                         std::vector<Candidate> &c1,
                                         std::vector<Candidate> &c2,
                                         std::vector<Candidate> &buffer) const {
  if (c1.size() == 0) {
    c1.swap(c2);
    return;
  }

  uint32_t i, j;
  uint32_t size1, size2;
  size1 = c1.size();
  size2 = c2.size();
  buffer.clear();

#ifdef LI_DEBUG
  for (i = 0; i < size1; ++i)
    printf("c1: %d %d %d\n", (int)(c1[i].position >> 32), (int)c1[i].position,
           c1[i].count);
  for (i = 0; i < size2; ++i)
    printf("c2: %d %d %d\n", (int)(c2[i].position >> 32), (int)c2[i].position,
           c2[i].count);
#endif

  i = 0;
  j = 0;
  while (i < size1 && j < size2) {
    if (c1[i].position == c2[j].position) {
      if (buffer.empty() ||
          c1[i].position > buffer.back().position + error_threshold) {
        if (c1[i].count > c2[j].count) {
          buffer.push_back(c1[i]);
        } else {
          buffer.push_back(c2[j]);
        }
      }
      ++i, ++j;
    } else if (c1[i].position < c2[j].position) {
      if (buffer.empty() ||
          c1[i].position > buffer.back().position + error_threshold) {
        buffer.push_back(c1[i]);
      }
      ++i;
    } else {
      if (buffer.empty() ||
          c2[j].position > buffer.back().position + error_threshold) {
        buffer.push_back(c2[j]);
      }
      ++j;
    }
  }

  while (i < size1) {
    if (buffer.empty() ||
        c1[i].position > buffer.back().position + error_threshold) {
      buffer.push_back(c1[i]);
    }
    ++i;
  }

  while (j < size2) {
    if (buffer.empty() ||
        c2[j].position > buffer.back().position + error_threshold) {
      buffer.push_back(c2[j]);
    }
    ++j;
  }

  c1.swap(buffer);
}

void CandidateProcessor::ReduceCandidatesForPairedEndReadOnOneDirection(
    uint32_t mapping_positions_distance,
    const std::vector<Candidate> &candidates1,
    const std::vector<Candidate> &candidates2,
    std::vector<Candidate> &filtered_candidates1,
    std::vector<Candidate> &filtered_candidates2) const {
  uint32_t i1 = 0;
  uint32_t i2 = 0;
  int num_unpaired_candidate1 = 0;
  int num_unpaired_candidate2 = 0;
  int num_unpaired_candidate_threshold = 5;
  int max_candidate_count1 = 6;
  int max_candidate_count2 = 6;
  uint32_t previous_end_i2 = i2;
#ifdef LI_DEBUG
  for (uint32_t i = 0; i < candidates1.size(); ++i)
    printf("%s 0: %d %d:%d\n", __func__, i,
           (int)(candidates1[i].position >> 32), (int)candidates1[i].position);
  for (uint32_t i = 0; i < candidates2.size(); ++i)
    printf("%s 1: %d %d:%d\n", __func__, i,
           (int)(candidates2[i].position >> 32), (int)candidates2[i].position);
#endif
  while (i1 < candidates1.size() && i2 < candidates2.size()) {
    if (candidates1[i1].position >
        candidates2[i2].position + mapping_positions_distance) {
      if (i2 >= previous_end_i2 &&
          num_unpaired_candidate2 < num_unpaired_candidate_threshold &&
          (candidates1[i1].position >> 32) ==
              (candidates2[i2].position >> 32) &&
          candidates2[i2].count >= max_candidate_count2) {
        filtered_candidates2.emplace_back(candidates2[i2]);
        ++num_unpaired_candidate2;
      }
      ++i2;
    } else if (candidates2[i2].position >
               candidates1[i1].position + mapping_positions_distance) {
      if (num_unpaired_candidate1 < num_unpaired_candidate_threshold &&
          (candidates1[i1].position >> 32) ==
              (candidates2[i2].position >> 32) &&
          candidates1[i1].count >= max_candidate_count1) {
        filtered_candidates1.emplace_back(candidates1[i1]);
        ++num_unpaired_candidate1;
      }
      ++i1;
    } else {
      // ok, find a pair, we store current ni2 somewhere and keep looking until
      // we go out of the range, then we go back and then move to next pi1 and
      // keep doing the similar thing.
      filtered_candidates1.emplace_back(candidates1[i1]);
      if (candidates1[i1].count > max_candidate_count1) {
        max_candidate_count1 = candidates1[i1].count;
      }
      uint32_t current_i2 = i2;
      while (current_i2 < candidates2.size() &&
             candidates2[current_i2].position <=
                 candidates1[i1].position + mapping_positions_distance) {
        if (current_i2 >= previous_end_i2) {
          filtered_candidates2.emplace_back(candidates2[current_i2]);
          if (candidates2[current_i2].count > max_candidate_count2) {
            max_candidate_count2 = candidates2[current_i2].count;
          }
        }
        ++current_i2;
      }
      previous_end_i2 = current_i2;
      ++i1;
    }
  }
}

}  // namespace chromap
