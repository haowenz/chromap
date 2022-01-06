#ifndef FEATURE_BARCODE_MATRIX_H_
#define FEATURE_BARCODE_MATRIX_H_

#include <assert.h>

#include <algorithm>
#include <cinttypes>
#include <cstring>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "bed_mapping.h"
#include "feature_barcode_matrix_writer.h"
#include "khash.h"
#include "utils.h"

namespace chromap {

struct Peak {
  uint32_t start_position;
  uint16_t length;
  uint32_t index;
};

class FeatureBarcodeMatrix {
 public:
  FeatureBarcodeMatrix(bool cell_by_bin, int bin_size, int overlap_distance,
                       uint16_t depth_cutoff_to_call_peak)
      : cell_by_bin_(cell_by_bin),
        bin_size_(bin_size),
        overlap_distance_(overlap_distance),
        depth_cutoff_to_call_peak_(depth_cutoff_to_call_peak) {
    barcode_index_table_ = kh_init(k64_seq);
  }

  ~FeatureBarcodeMatrix() {
    if (barcode_index_table_ != NULL) {
      kh_destroy(k64_seq, barcode_index_table_);
    }
  }

  void OutputFeatureMatrix(
      uint32_t num_sequences, const SequenceBatch &reference,
      const std::vector<std::vector<PairedEndMappingWithBarcode>> &mappings,
      const std::string &matrix_output_prefix);

 private:
  void BuildAugmentedTreeForPeaks(uint32_t ref_id);

  uint32_t GetNumOverlappedPeaks(
      uint32_t ref_id, const PairedEndMappingWithBarcode &mapping,
      std::vector<uint32_t> &overlapped_peak_indices);

  void GetNumOverlappedBins(uint32_t rid, uint32_t start_position,
                            uint16_t mapping_length, uint32_t num_sequences,
                            const SequenceBatch &reference,
                            std::vector<uint32_t> &overlapped_peak_indices);

  uint32_t CallPeaks(
      uint16_t coverage_threshold, uint32_t num_reference_sequences,
      const SequenceBatch &reference,
      const std::vector<std::vector<PairedEndMappingWithBarcode>> &mappings);

  const bool cell_by_bin_;
  const int bin_size_;
  const int overlap_distance_;
  const uint16_t depth_cutoff_to_call_peak_;

  khash_t(k64_seq) * barcode_index_table_;
  // (max_level, # nodes)
  std::vector<std::pair<int, uint32_t>> tree_info_on_diff_ref_seqs_;

  // max
  std::vector<std::vector<uint32_t>> tree_extras_on_diff_ref_seqs_;

  // For peak calling.
  std::vector<std::vector<uint16_t>> pileup_on_diff_ref_seqs_;
  std::vector<std::vector<Peak>> peaks_on_diff_ref_seqs_;

  FeatureBarcodeMatrixWriter feature_barcode_matrix_writer_;
};

}  // namespace chromap

#endif  // FEATURE_BARCODE_MATRIX_H_
