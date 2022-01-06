#include "feature_barcode_matrix.h"

#include <cinttypes>
#include <cstring>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

namespace chromap {

void FeatureBarcodeMatrix::BuildAugmentedTreeForPeaks(uint32_t ref_id) {
  // std::sort(mappings.begin(), mappings.end(), IntervalLess());
  int max_level = 0;
  size_t i, last_i = 0;  // last_i points to the rightmost node in the tree
  uint32_t last = 0;     // last is the max value at node last_i
  int k;
  std::vector<Peak> &peaks = peaks_on_diff_ref_seqs_[ref_id];
  std::vector<uint32_t> &extras = tree_extras_on_diff_ref_seqs_[ref_id];
  if (peaks.size() == 0) {
    max_level = -1;
  }

  for (i = 0; i < peaks.size(); i += 2) {
    last_i = i;
    // last = mappings[i].max = mappings[i].en; // leaves (i.e. at level 0)
    last = extras[i] =
        peaks[i].start_position + peaks[i].length;  // leaves (i.e. at level 0)
  }

  for (k = 1; 1LL << k <= (int64_t)peaks.size();
       ++k) {  // process internal nodes in the bottom-up order
    size_t x = 1LL << (k - 1);
    size_t i0 = (x << 1) - 1;
    size_t step = x << 2;  // i0 is the first node
    for (i = i0; i < peaks.size();
         i += step) {               // traverse all nodes at level k
      uint32_t el = extras[i - x];  // max value of the left child
      uint32_t er =
          i + x < peaks.size() ? extras[i + x] : last;  // of the right child
      uint32_t e = peaks[i].start_position + peaks[i].length;
      e = e > el ? e : el;
      e = e > er ? e : er;
      extras[i] = e;  // set the max value for node i
    }
    last_i =
        last_i >> k & 1
            ? last_i - x
            : last_i +
                  x;  // last_i now points to the parent of the original last_i
    if (last_i < peaks.size() &&
        extras[last_i] > last)  // update last accordingly
      last = extras[last_i];
  }

  max_level = k - 1;
  tree_info_on_diff_ref_seqs_.emplace_back(max_level, peaks.size());
}

uint32_t FeatureBarcodeMatrix::CallPeaks(
    uint16_t coverage_threshold, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const std::vector<std::vector<PairedEndMappingWithBarcode>> &mappings) {
  double real_start_time = GetRealTime();
  // std::vector<std::vector<PairedEndMappingWithBarcode>> &mappings =
  //    allocate_multi_mappings_
  //        ? allocated_mappings_on_diff_ref_seqs_
  //        : (remove_pcr_duplicates_ ? deduped_mappings_on_diff_ref_seqs_
  //                                  : mappings_on_diff_ref_seqs_);
  // Build pileup.
  for (uint32_t ri = 0; ri < num_reference_sequences; ++ri) {
    pileup_on_diff_ref_seqs_.emplace_back(std::vector<uint16_t>());
    pileup_on_diff_ref_seqs_[ri].assign(reference.GetSequenceLengthAt(ri), 0);
    for (size_t mi = 0; mi < mappings[ri].size(); ++mi) {
      for (uint16_t pi = 0; pi < mappings[ri][mi].fragment_length_; ++pi) {
        ++pileup_on_diff_ref_seqs_[ri]
                                  [mappings[ri][mi].GetStartPosition() + pi];
      }
    }
  }
  std::cerr << "Built pileup in " << GetRealTime() - real_start_time << "s.\n";

  real_start_time = GetRealTime();
  // Call and save peaks.
  tree_extras_on_diff_ref_seqs_.clear();
  tree_info_on_diff_ref_seqs_.clear();
  tree_extras_on_diff_ref_seqs_.reserve(num_reference_sequences);
  tree_info_on_diff_ref_seqs_.reserve(num_reference_sequences);
  uint32_t peak_count = 0;
  for (uint32_t ri = 0; ri < num_reference_sequences; ++ri) {
    tree_extras_on_diff_ref_seqs_.emplace_back(std::vector<uint32_t>());
    tree_extras_on_diff_ref_seqs_[ri].reserve(
        reference.GetSequenceLengthAt(ri) / 100);
    peaks_on_diff_ref_seqs_.emplace_back(std::vector<Peak>());
    uint32_t peak_start_position = 0;
    uint16_t peak_length = 0;
    for (size_t pi = 0; pi < reference.GetSequenceLengthAt(ri); ++pi) {
      if (pileup_on_diff_ref_seqs_[ri][pi] >= coverage_threshold) {
        if (peak_length == 0) {  // start a new peak
          peak_start_position = pi;
        }
        ++peak_length;               // extend the peak
      } else if (peak_length > 0) {  // save the previous peak
        // TODO(Haowen): improve peak calling
        peaks_on_diff_ref_seqs_[ri].emplace_back(
            Peak{peak_start_position, peak_length, peak_count});
        tree_extras_on_diff_ref_seqs_[ri].emplace_back(0);
        feature_barcode_matrix_writer_.OutputPeaks(peak_start_position,
                                                   peak_length, ri, reference);
        ++peak_count;
        peak_length = 0;
      }
    }
    BuildAugmentedTreeForPeaks(ri);
  }
  std::cerr << "Call peaks and built peak augmented tree in "
            << GetRealTime() - real_start_time << "s.\n";
  // Output feature matrix
  return peak_count;
}

void FeatureBarcodeMatrix::OutputFeatureMatrix(
    uint32_t num_sequences, const SequenceBatch &reference,
    const std::vector<std::vector<PairedEndMappingWithBarcode>> &mappings,
    const std::string &matrix_output_prefix) {
  feature_barcode_matrix_writer_.InitializeMatrixOutput(matrix_output_prefix);

  uint32_t num_peaks = 0;
  if (cell_by_bin_) {
    feature_barcode_matrix_writer_.OutputPeaks(bin_size_, num_sequences,
                                               reference);
    for (uint32_t i = 0; i < num_sequences; ++i) {
      uint32_t ref_seq_length = reference.GetSequenceLengthAt(i);
      num_peaks += ref_seq_length / bin_size_;
      if (ref_seq_length % bin_size_ != 0) {
        ++num_peaks;
      }
    }
  } else {
    num_peaks = CallPeaks(depth_cutoff_to_call_peak_, num_sequences, reference,
                          mappings);
  }

  // std::vector<std::vector<PairedEndMappingWithBarcode>> &mappings =
  //    allocate_multi_mappings_
  //        ? allocated_mappings_on_diff_ref_seqs_
  //        : (remove_pcr_duplicates_ ? deduped_mappings_on_diff_ref_seqs_
  //                                  : mappings_on_diff_ref_seqs_);
  double real_start_time = GetRealTime();
  // First pass to index barcodes
  uint32_t barcode_index = 0;
  for (uint32_t rid = 0; rid < num_sequences; ++rid) {
    for (uint32_t mi = 0; mi < mappings[rid].size(); ++mi) {
      uint64_t barcode_key = mappings[rid][mi].cell_barcode_;
      khiter_t barcode_index_table_iterator =
          kh_get(k64_seq, barcode_index_table_, barcode_key);
      if (barcode_index_table_iterator == kh_end(barcode_index_table_)) {
        int khash_return_code;
        barcode_index_table_iterator = kh_put(k64_seq, barcode_index_table_,
                                              barcode_key, &khash_return_code);
        assert(khash_return_code != -1 && khash_return_code != 0);
        kh_value(barcode_index_table_, barcode_index_table_iterator) =
            barcode_index;
        ++barcode_index;
        feature_barcode_matrix_writer_.AppendBarcodeOutput(barcode_key);
      }
    }
  }
  std::cerr << "Index and output barcodes in "
            << GetRealTime() - real_start_time << "s.\n";

  real_start_time = GetRealTime();
  // Second pass to generate matrix
  khash_t(kmatrix) *matrix = kh_init(kmatrix);
  std::vector<uint32_t> overlapped_peak_indices;
  for (uint32_t rid = 0; rid < num_sequences; ++rid) {
    for (uint32_t mi = 0; mi < mappings[rid].size(); ++mi) {
      uint64_t barcode_key = mappings[rid][mi].cell_barcode_;
      khiter_t barcode_index_table_iterator =
          kh_get(k64_seq, barcode_index_table_, barcode_key);
      uint64_t barcode_index =
          kh_value(barcode_index_table_, barcode_index_table_iterator);
      overlapped_peak_indices.clear();
      if (cell_by_bin_) {
        GetNumOverlappedBins(rid, mappings[rid][mi].GetStartPosition(),
                             mappings[rid][mi].GetEndPosition() -
                                 mappings[rid][mi].GetStartPosition(),
                             num_sequences, reference, overlapped_peak_indices);
      } else {
        GetNumOverlappedPeaks(rid, mappings[rid][mi], overlapped_peak_indices);
      }
      size_t num_overlapped_peaks = overlapped_peak_indices.size();
      for (size_t pi = 0; pi < num_overlapped_peaks; ++pi) {
        uint32_t peak_index = overlapped_peak_indices[pi];
        uint64_t entry_index = (barcode_index << 32) | peak_index;
        khiter_t matrix_iterator = kh_get(kmatrix, matrix, entry_index);
        if (matrix_iterator == kh_end(matrix)) {
          int khash_return_code;
          matrix_iterator =
              kh_put(kmatrix, matrix, entry_index, &khash_return_code);
          assert(khash_return_code != -1 && khash_return_code != 0);
          kh_value(matrix, matrix_iterator) = 1;
        } else {
          kh_value(matrix, matrix_iterator) += 1;
        }
      }
    }
  }
  std::cerr << "Generate feature matrix in " << GetRealTime() - real_start_time
            << "s.\n";
  // Output matrix
  real_start_time = GetRealTime();
  feature_barcode_matrix_writer_.WriteMatrixOutputHead(
      num_peaks, kh_size(barcode_index_table_), kh_size(matrix));
  uint64_t key;
  uint32_t value;
  std::vector<std::pair<uint64_t, uint32_t>> feature_matrix;
  feature_matrix.reserve(kh_size(matrix));
  // kh_foreach(matrix, key, value,
  // output_tools_->AppendMatrixOutput((uint32_t)key, (uint32_t)(key >> 32),
  // value));
  kh_foreach(matrix, key, value, feature_matrix.emplace_back(key, value));
  kh_destroy(kmatrix, matrix);
  std::sort(feature_matrix.begin(), feature_matrix.end());
  for (size_t i = 0; i < feature_matrix.size(); ++i) {
    feature_barcode_matrix_writer_.AppendMatrixOutput(
        (uint32_t)feature_matrix[i].first,
        (uint32_t)(feature_matrix[i].first >> 32), feature_matrix[i].second);
  }

  feature_barcode_matrix_writer_.FinalizeMatrixOutput();
  std::cerr << "Output feature matrix in " << GetRealTime() - real_start_time
            << "s.\n";
}

void FeatureBarcodeMatrix::GetNumOverlappedBins(
    uint32_t rid, uint32_t start_position, uint16_t mapping_length,
    uint32_t num_sequences, const SequenceBatch &reference,
    std::vector<uint32_t> &overlapped_peak_indices) {
  uint32_t bin_index = 0;
  for (uint32_t i = 0; i < rid; ++i) {
    uint32_t ref_seq_length = reference.GetSequenceLengthAt(i);
    bin_index += ref_seq_length / bin_size_;
    if (ref_seq_length % bin_size_ != 0) {
      ++bin_index;
    }
  }
  bin_index += start_position / bin_size_;
  overlapped_peak_indices.emplace_back(bin_index);
  uint32_t max_num_overlapped_bins = mapping_length / bin_size_ + 2;
  for (uint32_t i = 0; i < max_num_overlapped_bins; ++i) {
    if (start_position + mapping_length - 1 >=
        (bin_index + 1 + i) * bin_size_) {
      overlapped_peak_indices.emplace_back(bin_index + 1 + i);
    }
  }
}

uint32_t FeatureBarcodeMatrix::GetNumOverlappedPeaks(
    uint32_t ref_id, const PairedEndMappingWithBarcode &mapping,
    std::vector<uint32_t> &overlapped_peak_indices) {
  int t = 0;
  StackCell stack[64];
  // out.clear();
  overlapped_peak_indices.clear();
  int num_overlapped_peaks = 0;
  int max_level = tree_info_on_diff_ref_seqs_[ref_id].first;
  uint32_t num_tree_nodes = tree_info_on_diff_ref_seqs_[ref_id].second;
  std::vector<Peak> &peaks = peaks_on_diff_ref_seqs_[ref_id];
  std::vector<uint32_t> &extras = tree_extras_on_diff_ref_seqs_[ref_id];
  // uint32_t interval_start = mapping.fragment_start_position;
  uint32_t interval_start =
      mapping.GetStartPosition() > (uint32_t)overlap_distance_
          ? mapping.GetStartPosition() - overlap_distance_
          : 0;
  uint32_t interval_end =
      mapping.GetEndPosition() + (uint32_t)overlap_distance_;
  stack[t++] = StackCell(max_level, (1LL << max_level) - 1,
                         0);  // push the root; this is a top down traversal
  while (
      t) {  // the following guarantees that numbers in out[] are always sorted
    StackCell z = stack[--t];
    if (z.k <=
        3) {  // we are in a small subtree; traverse every node in this subtree
      size_t i, i0 = z.x >> z.k << z.k, i1 = i0 + (1LL << (z.k + 1)) - 1;
      if (i1 >= num_tree_nodes) {
        i1 = num_tree_nodes;
      }
      for (i = i0; i < i1 && peaks[i].start_position < interval_end; ++i) {
        if (interval_start <
            peaks[i].start_position +
                peaks[i].length) {  // if overlap, append to out[]
          // out.push_back(i);
          overlapped_peak_indices.emplace_back(peaks[i].index);
          ++num_overlapped_peaks;
        }
      }
    } else if (z.w == 0) {  // if left child not processed
      size_t y =
          z.x - (1LL << (z.k - 1));  // the left child of z.x; NB: y may be out
                                     // of range (i.e. y>=a.size())
      stack[t++] = StackCell(
          z.k, z.x,
          1);  // re-add node z.x, but mark the left child having been processed
      if (y >= num_tree_nodes ||
          extras[y] > interval_start)  // push the left child if y is out of
                                       // range or may overlap with the query
        stack[t++] = StackCell(z.k - 1, y, 0);
    } else if (z.x < num_tree_nodes &&
               peaks[z.x].start_position <
                   interval_end) {  // need to push the right child
      if (interval_start < peaks[z.x].start_position + peaks[z.x].length) {
        // out.push_back(z.x); // test if z.x overlaps the query; if yes, append
        // to out[]
        overlapped_peak_indices.emplace_back(peaks[z.x].index);
        ++num_overlapped_peaks;
      }
      stack[t++] = StackCell(z.k - 1, z.x + (1LL << (z.k - 1)),
                             0);  // push the right child
    }
  }
  return num_overlapped_peaks;
}

}  // namespace chromap
