#ifndef PAIRED_END_MAPPING_METADATA_H_
#define PAIRED_END_MAPPING_METADATA_H_

#include <algorithm>
#include <utility>
#include <vector>

#include "mapping_metadata.h"

namespace chromap {

class PairedEndMappingMetadata {
 public:
  inline void PreparedForMappingNextReadPair(int reserve_size) {
    mapping_metadata1_.PrepareForMappingNextRead(reserve_size);
    mapping_metadata2_.PrepareForMappingNextRead(reserve_size);

    F1R2_best_mappings_.clear();
    F2R1_best_mappings_.clear();
    F1F2_best_mappings_.clear();
    R1R2_best_mappings_.clear();

    F1R2_best_mappings_.reserve(reserve_size);
    F2R1_best_mappings_.reserve(reserve_size);
    F1F2_best_mappings_.reserve(reserve_size);
    R1R2_best_mappings_.reserve(reserve_size);
  }

  inline void MoveCandidiatesToBuffer() {
    mapping_metadata1_.MoveCandidiatesToBuffer();
    mapping_metadata2_.MoveCandidiatesToBuffer();
  }

  // Callback function to update all candidates.
  inline void UpdateCandidates(void (*Update)(std::vector<Candidate> &)) {
    mapping_metadata1_.UpdateCandidates(Update);
    mapping_metadata2_.UpdateCandidates(Update);
  }

  inline void SortMappingsByPositions() {
    mapping_metadata1_.SortMappingsByPositions();
    mapping_metadata2_.SortMappingsByPositions();
  }
  // inline void ClearAndReserveMinimizers(int reserve_size) {
  //  mapping_metadata1_.minimizers_.clear();
  //  mapping_metadata2_.minimizers_.clear();
  //  mapping_metadata1_.minimizers_.reserve(reserve_size);
  //  mapping_metadata2_.minimizers_.reserve(reserve_size);
  //}

  inline bool BothEndsHaveMinimizers() const {
    return !mapping_metadata1_.minimizers_.empty() &&
           !mapping_metadata2_.minimizers_.empty();
  }

  inline int GetMinSumErrors() const { return min_sum_errors_; }
  inline int GetSecondMinSumErrors() const { return second_min_sum_errors_; }
  inline int GetNumBestMappings() const { return num_best_mappings_; }
  inline int GetNumSecondBestMappings() const {
    return num_second_best_mappings_;
  }

  inline void SetMinSumErrors(int min_sum_errors) {
    min_sum_errors_ = min_sum_errors;
  }
  inline void SetSecondMinSumErrors(int second_min_sum_errors) {
    second_min_sum_errors_ = second_min_sum_errors;
  }
  inline void SetNumBestMappings(int num_best_mappings) {
    num_best_mappings_ = num_best_mappings;
  }
  inline void SetNumSecondBestMappings(int num_second_best_mappings) {
    num_second_best_mappings_ = num_second_best_mappings;
  }

 protected:
  MappingMetadata mapping_metadata1_;
  MappingMetadata mapping_metadata2_;

  int min_sum_errors_, second_min_sum_errors_;
  int num_best_mappings_, num_second_best_mappings_;

  std::vector<std::pair<uint32_t, uint32_t>> F1R2_best_mappings_;
  std::vector<std::pair<uint32_t, uint32_t>> F2R1_best_mappings_;
  std::vector<std::pair<uint32_t, uint32_t>> F1F2_best_mappings_;
  std::vector<std::pair<uint32_t, uint32_t>> R1R2_best_mappings_;

  friend class CandidateProcessor;
  template <typename MappingRecord>
  friend class MappingGenerator;
  template <typename MappingRecord>
  friend class Chromap;
};

}  // namespace chromap

#endif  // PAIRED_END_MAPPING_METADATA_H_
