#ifndef MAPPING_METADATA_H_
#define MAPPING_METADATA_H_

#include <algorithm>
#include <utility>
#include <vector>

namespace chromap {

struct Candidate {
  uint64_t position;
  uint8_t count;
  bool operator<(const Candidate &c) const {
    if (count != c.count)
      return count > c.count;
    else
      return position < c.position;
  }
};

class mm_cache;
class Index;
class PairedEndMappingMetadata;
template <typename MappingRecord>
class Chromap;

class MappingMetadata {
 public:
  inline void PrepareForMappingNextRead(int reserve_size) {
    Clear();
    ReserveSpace(reserve_size);
    repetitive_seed_length_ = 0;
  }

  inline size_t GetNumCandidates() const {
    return positive_candidates_.size() + negative_candidates_.size();
  }

  inline size_t GetNumMappings() const {
    return positive_mappings_.size() + negative_mappings_.size();
  }

  inline void MoveCandidiatesToBuffer() {
    positive_candidates_.swap(positive_candidates_buffer_);
    positive_candidates_.clear();
    negative_candidates_.swap(negative_candidates_buffer_);
    negative_candidates_.clear();
  }

  // Callback function to update all candidates.
  inline void UpdateCandidates(void (*Update)(std::vector<Candidate> &)) {
    Update(positive_candidates_);
    Update(negative_candidates_);
  }

  inline void SortMappingsByPositions() {
    auto compare_function = [](const std::pair<int, uint64_t> &a,
                               const std::pair<int, uint64_t> &b) {
      return a.second < b.second;
    };
    std::sort(positive_mappings_.begin(), positive_mappings_.end(),
              compare_function);
    std::sort(negative_mappings_.begin(), negative_mappings_.end(),
              compare_function);
  }

  inline int GetMinNumErrors() const { return min_num_errors_; }
  inline int GetSecondMinNumErrors() const { return second_min_num_errors_; }
  inline int GetNumBestMappings() const { return num_best_mappings_; }
  inline int GetNumSecondBestMappings() const {
    return num_second_best_mappings_;
  }

  inline void SetMinNumErrors(int min_num_errors) {
    min_num_errors_ = min_num_errors;
  }
  inline void SetSecondMinNumErrors(int second_min_num_errors) {
    second_min_num_errors_ = second_min_num_errors;
  }
  inline void SetNumBestMappings(int num_best_mappings) {
    num_best_mappings_ = num_best_mappings;
  }
  inline void SetNumSecondBestMappings(int num_second_best_mappings) {
    num_second_best_mappings_ = num_second_best_mappings;
  }

 protected:
  inline void ReserveSpace(int reserve_size) {
    minimizers_.reserve(reserve_size);
    positive_hits_.reserve(reserve_size);
    negative_hits_.reserve(reserve_size);
    positive_candidates_.reserve(reserve_size);
    negative_candidates_.reserve(reserve_size);
    positive_candidates_buffer_.reserve(reserve_size);
    negative_candidates_buffer_.reserve(reserve_size);
    positive_mappings_.reserve(reserve_size);
    negative_mappings_.reserve(reserve_size);
    positive_split_sites_.reserve(reserve_size);
    negative_split_sites_.reserve(reserve_size);
  }

  inline void Clear() {
    minimizers_.clear();
    positive_hits_.clear();
    negative_hits_.clear();
    positive_candidates_.clear();
    negative_candidates_.clear();
    positive_candidates_buffer_.clear();
    negative_candidates_buffer_.clear();
    positive_mappings_.clear();
    negative_mappings_.clear();
    positive_split_sites_.clear();
    negative_split_sites_.clear();
  }

  int min_num_errors_, second_min_num_errors_;
  int num_best_mappings_, num_second_best_mappings_;

  uint32_t repetitive_seed_length_;

  std::vector<std::pair<uint64_t, uint64_t>> minimizers_;

  std::vector<uint64_t> positive_hits_;
  std::vector<uint64_t> negative_hits_;

  std::vector<Candidate> positive_candidates_;
  std::vector<Candidate> negative_candidates_;

  std::vector<Candidate> positive_candidates_buffer_;
  std::vector<Candidate> negative_candidates_buffer_;

  std::vector<std::pair<int, uint64_t>> positive_mappings_;
  std::vector<std::pair<int, uint64_t>> negative_mappings_;

  std::vector<int> positive_split_sites_;
  std::vector<int> negative_split_sites_;

  friend class mm_cache;
  friend class Index;
  friend class PairedEndMappingMetadata;
  template <typename MappingRecord>
  friend class Chromap;
};

}  // namespace chromap

#endif  // MAPPING_METADATA_H_
