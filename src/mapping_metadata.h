#ifndef MAPPING_METADATA_H_
#define MAPPING_METADATA_H_

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <utility>
#include <vector>

#include "minimizer.h"
#include "candidate.h"
#include "draft_mapping.h"

namespace chromap {

class mm_cache;
class Index;
class CandidateProcessor;
class PairedEndMappingMetadata;
class DraftMappingGenerator;
template <typename MappingRecord>
class MappingGenerator;
class Chromap;

class MappingMetadata {
 public:
  inline void PrepareForMappingNextRead(int reserve_size) {
    Clear();
    ReserveSpace(reserve_size);
    repetitive_seed_length_ = 0;
  }

  inline size_t GetNumMinimizers() const {
    return minimizers_.size();
  }

  inline size_t GetNumPositiveCandidates() const {
    return positive_candidates_.size();
  }

  inline size_t GetNumNegativeCandidates() const {
    return negative_candidates_.size();
  }

  inline size_t GetNumCandidates() const {
    return positive_candidates_.size() + negative_candidates_.size();
  }

  inline size_t GetNumDraftMappings() const {
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

  inline void SortCandidates() {
    SortCandidatesByCount(positive_candidates_);
    SortCandidatesByCount(negative_candidates_);
  }

  inline void SortMappingsByPositions() {
    auto compare_function = [](const DraftMapping &a, const DraftMapping &b) {
      return a.position < b.position;
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

  // For debug only.
  inline void PrintCandidates(FILE *fp) {
    uint32_t i;
    for (i = 0; i < positive_candidates_.size(); ++i)
      fprintf(fp, "+ %d %d %d %d\n", i,
              (int)(positive_candidates_[i].position >> 32),
              (int)(positive_candidates_[i].position),
              positive_candidates_[i].count);
    for (i = 0; i < negative_candidates_.size(); ++i)
      fprintf(fp, "- %d %d %d %d\n", i,
              (int)(negative_candidates_[i].position >> 32),
              (int)(negative_candidates_[i].position),
              negative_candidates_[i].count);
  }

 protected:
  // Sorts the candidates into the order of Candidate::operator< (count
  // descending, then position ascending). GenerateCandidatesOnOneStrand and
  // MergeCandidates emit candidates in ascending position order, so after one
  // pass that checks this, a stable counting sort on count gives that order;
  // candidates equal in both count and position are identical, so the result
  // equals std::sort's. Any other input falls back to std::sort.
  static inline void SortCandidatesByCount(std::vector<Candidate> &candidates) {
    const size_t num_candidates = candidates.size();
    if (num_candidates < 2) {
      return;
    }

    Candidate *c = candidates.data();
    bool is_position_sorted = true;
    bool is_count_sorted = true;
    uint8_t min_count = c[0].count;
    uint8_t max_count = c[0].count;
    for (size_t i = 1; i < num_candidates; ++i) {
      is_position_sorted &= c[i - 1].position <= c[i].position;
      is_count_sorted &= c[i - 1].count >= c[i].count;
      min_count = std::min(min_count, c[i].count);
      max_count = std::max(max_count, c[i].count);
    }

    if (!is_position_sorted) {
      std::sort(candidates.begin(), candidates.end());
      return;
    }

    if (is_count_sorted) {
      return;
    }

    // Bucket b holds the candidates with count max_count - b.
    const uint32_t num_buckets = (uint32_t)max_count - min_count + 1;
    uint32_t bucket_starts[256];
    std::memset(bucket_starts, 0, num_buckets * sizeof(uint32_t));
    for (size_t i = 0; i < num_candidates; ++i) {
      ++bucket_starts[max_count - c[i].count];
    }

    uint32_t num_before = 0;
    for (uint32_t b = 0; b < num_buckets; ++b) {
      const uint32_t bucket_size = bucket_starts[b];
      bucket_starts[b] = num_before;
      num_before += bucket_size;
    }

    static thread_local std::vector<Candidate> sorted_candidates;
    if (sorted_candidates.size() < num_candidates) {
      sorted_candidates.resize(num_candidates);
    }

    Candidate *sorted = sorted_candidates.data();
    for (size_t i = 0; i < num_candidates; ++i) {
      sorted[bucket_starts[max_count - c[i].count]++] = c[i];
    }
    std::copy(sorted, sorted + num_candidates, c);
  }

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

  std::vector<Minimizer> minimizers_;

  std::vector<uint64_t> positive_hits_;
  std::vector<uint64_t> negative_hits_;

  std::vector<Candidate> positive_candidates_;
  std::vector<Candidate> negative_candidates_;

  std::vector<Candidate> positive_candidates_buffer_;
  std::vector<Candidate> negative_candidates_buffer_;

  // The first element is ed, and the second element is position.
  std::vector<DraftMapping> positive_mappings_;
  std::vector<DraftMapping> negative_mappings_;

  std::vector<int> positive_split_sites_;
  std::vector<int> negative_split_sites_;

  friend class mm_cache;
  friend class Index;
  friend class CandidateProcessor;
  friend class PairedEndMappingMetadata;
  friend class DraftMappingGenerator;
  template <typename MappingRecord>
  friend class MappingGenerator;
  friend class Chromap;
};

}  // namespace chromap

#endif  // MAPPING_METADATA_H_
