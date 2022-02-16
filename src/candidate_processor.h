#ifndef CANDIDATE_PROCESSOR_H_
#define CANDIDATE_PROCESSOR_H_

#include <cinttypes>
#include <cstring>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "candidate.h"
#include "index.h"
#include "mapping_metadata.h"
#include "paired_end_mapping_metadata.h"
#include "sequence_batch.h"

namespace chromap {

class CandidateProcessor {
 public:
  CandidateProcessor() = delete;

  CandidateProcessor(int min_num_seeds_required_for_mapping,
                     const std::vector<int> max_seed_frequencies)
      : min_num_seeds_required_for_mapping_(min_num_seeds_required_for_mapping),
        max_seed_frequencies_(max_seed_frequencies) {}

  ~CandidateProcessor() = default;

  void GenerateCandidates(int error_threshold, const Index &index,
                          MappingMetadata &mapping_metadata) const;

  int SupplementCandidates(
      int error_threshold, uint32_t search_range, const Index &index,
      PairedEndMappingMetadata &paired_end_mapping_metadata) const;

  void ReduceCandidatesForPairedEndRead(
      uint32_t mapping_positions_distance,
      PairedEndMappingMetadata &paired_end_mapping_metadata) const;

 private:
  void GenerateCandidatesOnOneDirection(
      int error_threshold, int num_seeds_required, uint32_t num_minimizers,
      std::vector<uint64_t> &hits, std::vector<Candidate> &candidates) const;

  int GenerateCandidatesFromRepetitiveReadWithMateInfo(
      int error_threshold, const Index &index,
      const std::vector<std::pair<uint64_t, uint64_t> > &minimizers,
      uint32_t &repetitive_seed_length, std::vector<uint64_t> &hits,
      std::vector<Candidate> &candidates,
      const std::vector<Candidate> &mate_candidates, const Direction direction,
      uint32_t search_range) const;

  void MergeCandidates(int error_threshold, std::vector<Candidate> &c1,
                       std::vector<Candidate> &c2,
                       std::vector<Candidate> &buffer) const;

  void ReduceCandidatesForPairedEndReadOnOneDirection(
      uint32_t mapping_positions_distance,
      const std::vector<Candidate> &candidates1,
      const std::vector<Candidate> &candidates2,
      std::vector<Candidate> &filtered_candidates1,
      std::vector<Candidate> &filtered_candidates2) const;

  const int min_num_seeds_required_for_mapping_;
  // Vector of size 2. The first element is the frequency threshold, and the
  // second element is the frequency threshold to run rescue. The second element
  // should always larger than the first one.
  // TODO(Haowen): add an error check.
  const std::vector<int> max_seed_frequencies_;
};

}  // namespace chromap

#endif  // CANDIDATE_PROCESSOR_H_
