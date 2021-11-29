#ifndef MAPPINGBUFFER_H_
#define MAPPINGBUFFER_H_

#include <utility>
#include <vector>

namespace chromap {
struct MappingBuffer {
  std::vector<std::pair<uint64_t, uint64_t>> minimizers1;
  std::vector<std::pair<uint64_t, uint64_t>> minimizers2;
  std::vector<uint64_t> positive_hits1;
  std::vector<uint64_t> positive_hits2;
  std::vector<uint64_t> negative_hits1;
  std::vector<uint64_t> negative_hits2;
  std::vector<Candidate> positive_candidates1;
  std::vector<Candidate> positive_candidates2;
  std::vector<Candidate> negative_candidates1;
  std::vector<Candidate> negative_candidates2;
  std::vector<Candidate> positive_candidates1_buffer;
  std::vector<Candidate> positive_candidates2_buffer;
  std::vector<Candidate> negative_candidates1_buffer;
  std::vector<Candidate> negative_candidates2_buffer;
  std::vector<std::pair<int, uint64_t>> positive_mappings1;
  std::vector<std::pair<int, uint64_t>> positive_mappings2;
  std::vector<std::pair<int, uint64_t>> negative_mappings1;
  std::vector<std::pair<int, uint64_t>> negative_mappings2;
  std::vector<int> positive_split_sites1;
  std::vector<int> negative_split_sites1;
  std::vector<int> positive_split_sites2;
  std::vector<int> negative_split_sites2;
  std::vector<std::pair<uint32_t, uint32_t>> F1R2_best_mappings;
  std::vector<std::pair<uint32_t, uint32_t>> F2R1_best_mappings;
  std::vector<std::pair<uint32_t, uint32_t>> F1F2_best_mappings;
  std::vector<std::pair<uint32_t, uint32_t>> R1R2_best_mappings;
};

}  // namespace chromap

#endif  // MAPPINGBUFFER_H_
