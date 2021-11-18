#include "chromap.h"

#include <assert.h>
#include <math.h>
#include <omp.h>
#include <smmintrin.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>

#include "cxxopts.hpp"
#include "ksw.h"
#include "mmcache.hpp"

namespace chromap {
template <typename MappingRecord>
void Chromap<MappingRecord>::TrimAdapterForPairedEndRead(
    uint32_t pair_index, SequenceBatch *read_batch1,
    SequenceBatch *read_batch2) {
  const char *read1 = read_batch1->GetSequenceAt(pair_index);
  uint32_t read2_length = read_batch2->GetSequenceLengthAt(pair_index);
  const std::string &negative_read2 =
      read_batch2->GetNegativeSequenceAt(pair_index);
  int min_overlap_length = min_read_length_;
  int seed_length = min_overlap_length / 2;
  int error_threshold_for_merging = 1;
  bool is_merged = false;
  for (int si = 0; si < error_threshold_for_merging + 1; ++si) {
    int seed_start_position =
        negative_read2.find(read1 + si * seed_length, 0, seed_length);
    while ((uint32_t)seed_start_position != std::string::npos &&
           read2_length - seed_start_position + seed_length * si >=
               (uint32_t)min_overlap_length &&
           seed_start_position >= si * seed_length) {
      bool can_merge = true;
      int num_errors = 0;
      for (int i = 0; i < seed_length * si; ++i) {
        if (negative_read2[seed_start_position - si * seed_length + i] !=
            read1[i]) {
          ++num_errors;
        }
        if (num_errors > error_threshold_for_merging) {
          can_merge = false;
          break;
        }
      }
      for (uint32_t i = seed_length; i + seed_start_position < read2_length;
           ++i) {
        if (negative_read2[seed_start_position + i] !=
            read1[si * seed_length + i]) {
          ++num_errors;
        }
        if (num_errors > error_threshold_for_merging) {
          can_merge = false;
          break;
        }
      }
      if (can_merge) {
        // Trim adapters and TODO: fix sequencing errors
        int overlap_length =
            read2_length - seed_start_position + si * seed_length;
        read_batch1->TrimSequenceAt(pair_index, overlap_length);
        read_batch2->TrimSequenceAt(pair_index, overlap_length);
        is_merged = true;
        // std::cerr << "Trimed! overlap length: " << overlap_length << ", " <<
        // read1.GetLength() << " " << read2.GetLength() << "\n";
        break;
      }
      seed_start_position = negative_read2.find(
          read1 + si * seed_length, seed_start_position + 1, seed_length);
    }
    if (is_merged) {
      break;
    }
  }
}

template <typename MappingRecord>
bool Chromap<MappingRecord>::PairedEndReadWithBarcodeIsDuplicate(
    uint32_t pair_index, const SequenceBatch &barcode_batch,
    const SequenceBatch &read_batch1, const SequenceBatch &read_batch2) {
  int dedupe_seed_length = 16;
  uint32_t barcode_length = barcode_batch.GetSequenceLengthAt(pair_index);
  uint64_t barcode_key =
      barcode_batch.GenerateSeedFromSequenceAt(pair_index, 0, barcode_length);
  uint64_t read1_seed1 =
      read_batch1.GenerateSeedFromSequenceAt(pair_index, 0, dedupe_seed_length);
  uint64_t read2_seed1 =
      read_batch2.GenerateSeedFromSequenceAt(pair_index, 0, dedupe_seed_length);
  uint64_t read_seed_key =
      (read1_seed1 << (dedupe_seed_length * 2)) | read2_seed1;
  uint64_t read1_seed2 = read_batch1.GenerateSeedFromSequenceAt(
      pair_index, dedupe_seed_length, dedupe_seed_length * 2);
  uint64_t read2_seed2 = read_batch2.GenerateSeedFromSequenceAt(
      pair_index, dedupe_seed_length, dedupe_seed_length * 2);
  khiter_t barcode_table_iterator =
      kh_get(k64_seq, barcode_lookup_table_, barcode_key);
  if (barcode_table_iterator != kh_end(barcode_lookup_table_)) {
    uint32_t read_lookup_table_index =
        kh_value(barcode_lookup_table_, barcode_table_iterator);
    // std::cerr << "Have barcode, try to check read. " <<
    // read_lookup_table_index << "\n";
    khash_t(k128) *read_lookup_table =
        read_lookup_tables_[read_lookup_table_index];
    khiter_t read_lookup_table_iterator =
        kh_get(k128, read_lookup_table, read_seed_key);
    if (read_lookup_table_iterator != kh_end(read_lookup_table)) {
      // std::cerr << "Have barcode, have read, try whether match.\n";
      uint128_t read_seeds =
          kh_value(read_lookup_table, read_lookup_table_iterator);
      if (read_seeds.first == read1_seed2 && read_seeds.second == read2_seed2) {
        // std::cerr << "Have barcode, have read, and match.\n";
        return true;
      } else {
        // std::cerr << "Have barcode, have read, but don't match.\n";
        return false;
      }
    } else {
      // std::cerr << "Have barcode, no read.\n";
      uint128_t read_seeds = {.first = read1_seed2, .second = read2_seed2};
      int khash_return_code;
      khiter_t read_lookup_table_insert_iterator =
          kh_put(k128, read_lookup_table, read_seed_key, &khash_return_code);
      assert(khash_return_code != -1 && khash_return_code != 0);
      kh_value(read_lookup_table, read_lookup_table_insert_iterator) =
          read_seeds;
      // std::cerr << "Have barcode, no read.\n";
      return false;
    }
  } else {
    // insert the barcode and append a new read hash table to tables and then
    // insert the reads
    // std::cerr << "No barcode, no read.\n";
    int khash_return_code;
    khiter_t barcode_table_insert_iterator =
        kh_put(k64_seq, barcode_lookup_table_, barcode_key, &khash_return_code);
    assert(khash_return_code != -1 && khash_return_code != 0);
    kh_value(barcode_lookup_table_, barcode_table_insert_iterator) =
        read_lookup_tables_.size();
    khash_t(k128) *read_lookup_table = kh_init(k128);
    khiter_t read_lookup_table_iterator =
        kh_put(k128, read_lookup_table, read_seed_key, &khash_return_code);
    assert(khash_return_code != -1 && khash_return_code != 0);
    uint128_t read_seeds = {.first = read1_seed2, .second = read2_seed2};
    kh_value(read_lookup_table, read_lookup_table_iterator) = read_seeds;
    read_lookup_tables_.push_back(read_lookup_table);
    if (kh_size(barcode_lookup_table_) >=
        (uint32_t)allocated_barcode_lookup_table_size_) {
      allocated_barcode_lookup_table_size_ <<= 1;
      kh_resize(k64_seq, barcode_lookup_table_,
                allocated_barcode_lookup_table_size_);
    }
    // std::cerr << "No barcode, no read.\n";
    return false;
  }
}

template <typename MappingRecord>
bool Chromap<MappingRecord>::CorrectBarcodeAt(
    uint32_t barcode_index, SequenceBatch *barcode_batch,
    uint64_t *num_barcode_in_whitelist, uint64_t *num_corrected_barcode) {
  uint32_t barcode_length = barcode_batch->GetSequenceLengthAt(barcode_index);
  uint64_t barcode_key = barcode_batch->GenerateSeedFromSequenceAt(
      barcode_index, 0, barcode_length);
  khiter_t barcode_whitelist_lookup_table_iterator =
      kh_get(k64_seq, barcode_whitelist_lookup_table_, barcode_key);
  if (barcode_whitelist_lookup_table_iterator !=
      kh_end(barcode_whitelist_lookup_table_)) {
    // Correct barcode
    ++(*num_barcode_in_whitelist);
    return true;
  } else if (barcode_correction_error_threshold_ > 0) {
    // Need to correct this barcode
    // const char *barcode = barcode_batch->GetSequenceAt(barcode_index);
    // std::cerr << barcode_index << " barcode " << barcode << " needs
    // correction\n";
    const char *barcode_qual = barcode_batch->GetSequenceQualAt(barcode_index);
    std::vector<BarcodeWithQual> corrected_barcodes_with_quals;
    uint64_t mask = (uint64_t)3;
    for (uint32_t i = 0; i < barcode_length; ++i) {
      uint64_t barcode_key_to_change = mask << (2 * i);
      barcode_key_to_change = ~barcode_key_to_change;
      barcode_key_to_change &= barcode_key;
      uint64_t base_to_change1 = (barcode_key >> (2 * i)) & mask;
      for (uint32_t ti = 0; ti < 3; ++ti) {
        // change the base
        base_to_change1 += 1;
        base_to_change1 &= mask;
        // generate the corrected key
        uint64_t corrected_barcode_key =
            barcode_key_to_change | (base_to_change1 << (2 * i));
        barcode_whitelist_lookup_table_iterator = kh_get(
            k64_seq, barcode_whitelist_lookup_table_, corrected_barcode_key);
        if (barcode_whitelist_lookup_table_iterator !=
            kh_end(barcode_whitelist_lookup_table_)) {
          // find one possible corrected barcode
          double barcode_abundance =
              kh_value(barcode_whitelist_lookup_table_,
                       barcode_whitelist_lookup_table_iterator) /
              (double)num_sample_barcodes_;
          int qual_offset = 33;
          int adjusted_qual =
              barcode_qual[barcode_length - 1 - i] - qual_offset;
          adjusted_qual = adjusted_qual > 40 ? 40 : adjusted_qual;
          adjusted_qual = adjusted_qual < 3 ? 3 : adjusted_qual;
          double score =
              pow(10.0, ((-adjusted_qual) / 10.0)) * barcode_abundance;
          corrected_barcodes_with_quals.emplace_back(BarcodeWithQual{
              barcode_length - 1 - i,
              SequenceBatch::Uint8ToChar(base_to_change1), 0, 0, score});
          // std::cerr << "1score: " << score << " pos1: " << barcode_length - 1
          // - i << " b1: " << base_to_change1 << " pos2: " << 0 << " b2: " <<
          // (char)0 << "\n";
        }
        if (barcode_correction_error_threshold_ == 2) {
          for (uint32_t j = i + 1; j < barcode_length; ++j) {
            uint64_t barcode_key_to_change2 = mask << (2 * i);
            barcode_key_to_change2 = mask << (2 * j);
            barcode_key_to_change2 = ~barcode_key_to_change2;
            barcode_key_to_change2 &= corrected_barcode_key;
            uint64_t base_to_change2 =
                (corrected_barcode_key >> (2 * j)) & mask;
            for (uint32_t ti2 = 0; ti2 < 3; ++ti2) {
              // change the base
              base_to_change2 += 1;
              base_to_change2 &= mask;
              // generate the corrected key
              uint64_t corrected_barcode_key2 =
                  barcode_key_to_change2 | (base_to_change2 << (2 * j));
              barcode_whitelist_lookup_table_iterator =
                  kh_get(k64_seq, barcode_whitelist_lookup_table_,
                         corrected_barcode_key2);
              if (barcode_whitelist_lookup_table_iterator !=
                  kh_end(barcode_whitelist_lookup_table_)) {
                // find one possible corrected barcode
                double barcode_abundance =
                    kh_value(barcode_whitelist_lookup_table_,
                             barcode_whitelist_lookup_table_iterator) /
                    (double)num_sample_barcodes_;
                int qual_offset = 33;
                int adjusted_qual =
                    barcode_qual[barcode_length - 1 - j] - qual_offset;
                adjusted_qual = adjusted_qual > 40 ? 40 : adjusted_qual;
                adjusted_qual = adjusted_qual < 3 ? 3 : adjusted_qual;
                int adjusted_qual1 =
                    barcode_qual[barcode_length - 1 - i] - qual_offset;
                adjusted_qual1 = adjusted_qual1 > 40 ? 40 : adjusted_qual1;
                adjusted_qual1 = adjusted_qual1 < 3 ? 3 : adjusted_qual1;
                adjusted_qual += adjusted_qual1;
                double score =
                    pow(10.0, ((-adjusted_qual) / 10.0)) * barcode_abundance;
                corrected_barcodes_with_quals.emplace_back(BarcodeWithQual{
                    barcode_length - 1 - i,
                    SequenceBatch::Uint8ToChar(base_to_change1),
                    barcode_length - 1 - j,
                    SequenceBatch::Uint8ToChar(base_to_change2), score});
                // std::cerr << "2score: " << score << " pos1: " <<
                // barcode_length - 1 - i << " b1: " << base_to_change1 << "
                // pos2: " << barcode_length - 1 -j << " b2: " <<
                // base_to_change2
                // << "\n";
              }
            }
          }
        }
      }
    }
    size_t num_possible_corrected_barcodes =
        corrected_barcodes_with_quals.size();
    if (num_possible_corrected_barcodes == 0) {
      // Barcode cannot be corrected, leave it for downstream
      return false;
    } else if (num_possible_corrected_barcodes == 1) {
      // Just correct it
      // std::cerr << "Corrected the barcode from " << barcode << " to ";
      barcode_batch->CorrectBaseAt(
          barcode_index, corrected_barcodes_with_quals[0].corrected_base_index1,
          corrected_barcodes_with_quals[0].correct_base1);
      if (corrected_barcodes_with_quals[0].correct_base2 != 0) {
        barcode_batch->CorrectBaseAt(
            barcode_index,
            corrected_barcodes_with_quals[0].corrected_base_index2,
            corrected_barcodes_with_quals[0].correct_base2);
      }
      // std::cerr << barcode << "\n";
      // std::cerr << "score: " << corrected_barcodes_with_quals[0].score <<
      // "\n"; std::cerr << "score: " << corrected_barcodes_with_quals[0].score
      // << " pos1: " << corrected_barcodes_with_quals[0].corrected_base_index1
      // << " b1: " << corrected_barcodes_with_quals[0].correct_base1 << " pos2:
      // " << corrected_barcodes_with_quals[0].corrected_base_index2 << " b2: "
      // << corrected_barcodes_with_quals[0].correct_base2 << "\n";
      ++(*num_corrected_barcode);
      return true;
    } else {
      // Select the best correction
      std::sort(corrected_barcodes_with_quals.begin(),
                corrected_barcodes_with_quals.end(),
                std::greater<BarcodeWithQual>());
      // int num_ties = 0;
      double sum_score = 0;
      for (size_t ci = 0; ci < num_possible_corrected_barcodes; ++ci) {
        sum_score += corrected_barcodes_with_quals[ci].score;
        // std::cerr << ci << " score: " <<
        // corrected_barcodes_with_quals[ci].score << " pos1: " <<
        // corrected_barcodes_with_quals[ci].corrected_base_index1 << " b1: " <<
        // corrected_barcodes_with_quals[ci].correct_base1 << " pos2: " <<
        // corrected_barcodes_with_quals[ci].corrected_base_index2 << " b2: " <<
        // corrected_barcodes_with_quals[ci].correct_base2 << "\n"; if
        // (corrected_barcodes_with_quals[ci].qual ==
        // corrected_barcodes_with_quals[0].qual) {
        //  ++num_ties;
        //}
      }
      int best_corrected_barcode_index = 0;
      // if (num_ties > 0) {
      //  std::mt19937 tmp_generator(11);
      //  std::uniform_int_distribution<int> distribution(0, num_ties); //
      //  important: inclusive range best_corrected_barcode_index =
      //  distribution(tmp_generator);
      //}
      // std::cerr << "Corrected the barcode from " << barcode << " to ";
      double confidence_threshold = barcode_correction_probability_threshold_;
      if (corrected_barcodes_with_quals[best_corrected_barcode_index].score /
              sum_score >
          confidence_threshold) {
        barcode_batch->CorrectBaseAt(
            barcode_index,
            corrected_barcodes_with_quals[best_corrected_barcode_index]
                .corrected_base_index1,
            corrected_barcodes_with_quals[best_corrected_barcode_index]
                .correct_base1);
        if (corrected_barcodes_with_quals[best_corrected_barcode_index]
                .correct_base2 != 0) {
          barcode_batch->CorrectBaseAt(
              barcode_index,
              corrected_barcodes_with_quals[best_corrected_barcode_index]
                  .corrected_base_index2,
              corrected_barcodes_with_quals[best_corrected_barcode_index]
                  .correct_base2);
        }
        // std::cerr << barcode << "\n";
        // std::cerr << "score: " <<
        // corrected_barcodes_with_quals[best_corrected_barcode_index].score <<
        // "\n"; std::cerr << "best score: " <<
        // corrected_barcodes_with_quals[best_corrected_barcode_index].score <<
        // " sum score: " << sum_score << "\n";
        ++(*num_corrected_barcode);
        return true;
      } else {
        // std::cerr << "Didnt pass filter: " <<
        // corrected_barcodes_with_quals[best_corrected_barcode_index].score /
        // sum_score << "\n"; std::cerr << "best score: " <<
        // corrected_barcodes_with_quals[best_corrected_barcode_index].score <<
        // " sum score: " << sum_score << "\n";
        return false;
      }
    }
  } else {
    return false;
  }
}

template <typename MappingRecord>
uint32_t Chromap<MappingRecord>::CallPeaks(uint16_t coverage_threshold,
                                           uint32_t num_reference_sequences,
                                           const SequenceBatch &reference) {
  return 0;
}

template <>
uint32_t Chromap<PairedEndMappingWithBarcode>::CallPeaks(
    uint16_t coverage_threshold, uint32_t num_reference_sequences,
    const SequenceBatch &reference) {
  double real_start_time = GetRealTime();
  std::vector<std::vector<PairedEndMappingWithBarcode>> &mappings =
      allocate_multi_mappings_
          ? allocated_mappings_on_diff_ref_seqs_
          : (remove_pcr_duplicates_ ? deduped_mappings_on_diff_ref_seqs_
                                    : mappings_on_diff_ref_seqs_);
  // Build pileup
  for (uint32_t ri = 0; ri < num_reference_sequences; ++ri) {
    pileup_on_diff_ref_seqs_.emplace_back(std::vector<uint16_t>());
    pileup_on_diff_ref_seqs_[ri].assign(reference.GetSequenceLengthAt(ri), 0);
    for (size_t mi = 0; mi < mappings[ri].size(); ++mi) {
      for (uint16_t pi = 0; pi < mappings[ri][mi].fragment_length; ++pi) {
        ++pileup_on_diff_ref_seqs_[ri]
                                  [mappings[ri][mi].GetStartPosition() + pi];
      }
    }
  }
  std::cerr << "Built pileup in " << Chromap<>::GetRealTime() - real_start_time
            << "s.\n";
  real_start_time = GetRealTime();
  // Call and save peaks
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
        output_tools_.OutputPeaks(peak_start_position, peak_length, ri,
                                  reference);
        ++peak_count;
        peak_length = 0;
      }
    }
    BuildAugmentedTreeForPeaks(ri);
  }
  std::cerr << "Call peaks and built peak augmented tree in "
            << Chromap<>::GetRealTime() - real_start_time << "s.\n";
  // Output feature matrix
  return peak_count;
}

template <typename MappingRecord>
void Chromap<MappingRecord>::OutputFeatureMatrix(
    uint32_t num_sequences, const SequenceBatch &reference) {}

template <>
void Chromap<PairedEndMappingWithBarcode>::OutputFeatureMatrix(
    uint32_t num_sequences, const SequenceBatch &reference) {
  uint32_t num_peaks = 0;
  if (cell_by_bin_) {
    output_tools_.OutputPeaks(bin_size_, num_sequences, reference);
    for (uint32_t i = 0; i < num_sequences; ++i) {
      uint32_t ref_seq_length = reference.GetSequenceLengthAt(i);
      num_peaks += ref_seq_length / bin_size_;
      if (ref_seq_length % bin_size_ != 0) {
        ++num_peaks;
      }
    }
  } else {
    num_peaks = CallPeaks(depth_cutoff_to_call_peak_, num_sequences, reference);
  }
  std::vector<std::vector<PairedEndMappingWithBarcode>> &mappings =
      allocate_multi_mappings_
          ? allocated_mappings_on_diff_ref_seqs_
          : (remove_pcr_duplicates_ ? deduped_mappings_on_diff_ref_seqs_
                                    : mappings_on_diff_ref_seqs_);
  double real_start_time = GetRealTime();
  // First pass to index barcodes
  uint32_t barcode_index = 0;
  for (uint32_t rid = 0; rid < num_sequences; ++rid) {
    for (uint32_t mi = 0; mi < mappings[rid].size(); ++mi) {
      uint64_t barcode_key = mappings[rid][mi].cell_barcode;
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
        output_tools_.AppendBarcodeOutput(barcode_key);
      }
    }
  }
  std::cerr << "Index and output barcodes in "
            << Chromap<>::GetRealTime() - real_start_time << "s.\n";
  real_start_time = GetRealTime();
  // Second pass to generate matrix
  khash_t(kmatrix) *matrix = kh_init(kmatrix);
  std::vector<uint32_t> overlapped_peak_indices;
  for (uint32_t rid = 0; rid < num_sequences; ++rid) {
    for (uint32_t mi = 0; mi < mappings[rid].size(); ++mi) {
      uint64_t barcode_key = mappings[rid][mi].cell_barcode;
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
  std::cerr << "Generate feature matrix in "
            << Chromap<>::GetRealTime() - real_start_time << "s.\n";
  // Output matrix
  real_start_time = GetRealTime();
  output_tools_.WriteMatrixOutputHead(num_peaks, kh_size(barcode_index_table_),
                                      kh_size(matrix));
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
    output_tools_.AppendMatrixOutput((uint32_t)feature_matrix[i].first,
                                     (uint32_t)(feature_matrix[i].first >> 32),
                                     feature_matrix[i].second);
  }
  std::cerr << "Output feature matrix in "
            << Chromap<>::GetRealTime() - real_start_time << "s.\n";
}

template <typename MappingRecord>
void Chromap<MappingRecord>::GetNumOverlappedBins(
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

template <typename MappingRecord>
void Chromap<MappingRecord>::BuildAugmentedTreeForPeaks(uint32_t ref_id) {
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

template <typename MappingRecord>
uint32_t Chromap<MappingRecord>::GetNumOverlappedPeaks(
    uint32_t ref_id, const MappingRecord &mapping,
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
      mapping.GetStartPosition() > (uint32_t)multi_mapping_allocation_distance_
          ? mapping.GetStartPosition() - multi_mapping_allocation_distance_
          : 0;
  uint32_t interval_end =
      mapping.GetEndPosition() + (uint32_t)multi_mapping_allocation_distance_;
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

template <typename MappingRecord>
uint32_t Chromap<MappingRecord>::LoadPairedEndReadsWithBarcodes(
    SequenceBatch *read_batch1, SequenceBatch *read_batch2,
    SequenceBatch *barcode_batch) {
  // double real_start_time = Chromap<>::GetRealTime();
  uint32_t num_loaded_pairs = 0;
  while (num_loaded_pairs < read_batch_size_) {
    bool no_more_read1 =
        read_batch1->LoadOneSequenceAndSaveAt(num_loaded_pairs);
    bool no_more_read2 =
        read_batch2->LoadOneSequenceAndSaveAt(num_loaded_pairs);
    bool no_more_barcode = no_more_read2;
    if (!is_bulk_data_) {
      no_more_barcode =
          barcode_batch->LoadOneSequenceAndSaveAt(num_loaded_pairs);
    }
    if ((!no_more_read1) && (!no_more_read2) && (!no_more_barcode)) {
      if (read_batch1->GetSequenceLengthAt(num_loaded_pairs) <
              (uint32_t)min_read_length_ ||
          read_batch2->GetSequenceLengthAt(num_loaded_pairs) <
              (uint32_t)min_read_length_) {
        continue;  // reads are too short, just drop.
      }
      // if (PairedEndReadWithBarcodeIsDuplicate(num_loaded_pairs,
      // (*barcode_batch), (*read_batch1), (*read_batch2))) {
      //  num_duplicated_reads_ += 2;
      //  continue;
      //}
    } else if (no_more_read1 && no_more_read2 && no_more_barcode) {
      break;
    } else {
      Chromap<>::ExitWithMessage("Numbers of reads and barcodes don't match!");
    }
    ++num_loaded_pairs;
  }
  // if (num_loaded_pairs > 0) {
  //  std::cerr << "Loaded " << num_loaded_pairs << " pairs in "<<
  //  Chromap<>::GetRealTime() - real_start_time << "s. ";
  //} else {
  //  std::cerr << "No more reads.\n";
  //}
  return num_loaded_pairs;
}

template <typename MappingRecord>
void Chromap<MappingRecord>::ComputeBarcodeAbundance(
    uint64_t max_num_sample_barcodes) {
  double real_start_time = Chromap<>::GetRealTime();
  SequenceBatch barcode_batch(read_batch_size_);
  for (size_t read_file_index = 0; read_file_index < read_file1_paths_.size();
       ++read_file_index) {
    barcode_batch.InitializeLoading(barcode_file_paths_[read_file_index]);
    uint32_t num_loaded_barcodes = barcode_batch.LoadBatch();
    while (num_loaded_barcodes > 0) {
      for (uint32_t barcode_index = 0; barcode_index < num_loaded_barcodes;
           ++barcode_index) {
        uint32_t barcode_length =
            barcode_batch.GetSequenceLengthAt(barcode_index);
        uint64_t barcode_key = barcode_batch.GenerateSeedFromSequenceAt(
            barcode_index, 0, barcode_length);
        khiter_t barcode_whitelist_lookup_table_iterator =
            kh_get(k64_seq, barcode_whitelist_lookup_table_, barcode_key);
        if (barcode_whitelist_lookup_table_iterator !=
            kh_end(barcode_whitelist_lookup_table_)) {
          // Correct barcode
          kh_value(barcode_whitelist_lookup_table_,
                   barcode_whitelist_lookup_table_iterator) += 1;
          ++num_sample_barcodes_;
        }
      }

      if (!skip_barcode_check_  
        && num_sample_barcodes_ * 20 < num_loaded_barcodes) {
        // Since num_loaded_pairs is a constant, this if is actuaclly only effective in the first iteration
        Chromap<>::ExitWithMessage("Less than 5\% barcodes can be found or corrected based on the barcode whitelist.\nPlease check whether the barcode whitelist matches the data, e.g. length, reverse-complement. If this is a false positive warning, please run Chromap with the option --skip-barcode-check.");
      }

      if (num_sample_barcodes_ >= max_num_sample_barcodes) {
        break;
      }
      num_loaded_barcodes = barcode_batch.LoadBatch();
    }
    barcode_batch.FinalizeLoading();
    if (num_sample_barcodes_ >= max_num_sample_barcodes) {
      break;
    }
  }

  std::cerr << "Compute barcode abundance using " << num_sample_barcodes_
            << " in " << Chromap<>::GetRealTime() - real_start_time << "s.\n";
}

template <typename MappingRecord>
void Chromap<MappingRecord>::UpdateBarcodeAbundance(
    uint32_t num_loaded_barcodes, const SequenceBatch &barcode_batch) {
  double real_start_time = Chromap<>::GetRealTime();
  for (uint32_t barcode_index = 0; barcode_index < num_loaded_barcodes;
       ++barcode_index) {
    uint32_t barcode_length = barcode_batch.GetSequenceLengthAt(barcode_index);
    uint64_t barcode_key = barcode_batch.GenerateSeedFromSequenceAt(
        barcode_index, 0, barcode_length);
    khiter_t barcode_whitelist_lookup_table_iterator =
        kh_get(k64_seq, barcode_whitelist_lookup_table_, barcode_key);
    if (barcode_whitelist_lookup_table_iterator !=
        kh_end(barcode_whitelist_lookup_table_)) {
      // Correct barcode
      kh_value(barcode_whitelist_lookup_table_,
               barcode_whitelist_lookup_table_iterator) += 1;
      ++num_sample_barcodes_;
    }
  }
  std::cerr << "Update barcode abundance using " << num_sample_barcodes_
            << " in " << Chromap<>::GetRealTime() - real_start_time << "s.\n";
}

// 0-supplement normally
// 1-the supplement could be too aggressive, need to set MAPQ to 0.
template <typename MappingRecord>
int Chromap<MappingRecord>::SupplementCandidates(
    const Index &index, uint32_t repetitive_seed_length1,
    uint32_t repetitive_seed_length2,
    std::vector<std::pair<uint64_t, uint64_t>> &minimizers1,
    std::vector<std::pair<uint64_t, uint64_t>> &minimizers2,
    std::vector<uint64_t> &positive_hits1,
    std::vector<uint64_t> &positive_hits2,
    std::vector<Candidate> &positive_candidates1,
    std::vector<Candidate> &positive_candidates2,
    std::vector<Candidate> &positive_candidates1_buffer,
    std::vector<Candidate> &positive_candidates2_buffer,
    std::vector<uint64_t> &negative_hits1,
    std::vector<uint64_t> &negative_hits2,
    std::vector<Candidate> &negative_candidates1,
    std::vector<Candidate> &negative_candidates2,
    std::vector<Candidate> &negative_candidates1_buffer,
    std::vector<Candidate> &negative_candidates2_buffer) {
  std::vector<Candidate> augment_positive_candidates1;
  std::vector<Candidate> augment_positive_candidates2;
  std::vector<Candidate> augment_negative_candidates1;
  std::vector<Candidate> augment_negative_candidates2;
  int ret = 0;
  for (int mate = 0; mate <= 1; ++mate) {
    std::vector<std::pair<uint64_t, uint64_t>> *minimizers;
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
      minimizers = &minimizers1;
      positive_hits = &positive_hits1;
      negative_hits = &negative_hits1;
      positive_candidates = &positive_candidates1;
      negative_candidates = &negative_candidates1;
      mate_positive_candidates = &positive_candidates2;
      mate_negative_candidates = &negative_candidates2;
      augment_positive_candidates = &augment_positive_candidates1;
      augment_negative_candidates = &augment_negative_candidates1;
      repetitive_seed_length = &repetitive_seed_length1;
    } else {
      minimizers = &minimizers2;
      positive_hits = &positive_hits2;
      negative_hits = &negative_hits2;
      positive_candidates = &positive_candidates2;
      negative_candidates = &negative_candidates2;
      mate_positive_candidates = &positive_candidates1;
      mate_negative_candidates = &negative_candidates1;
      augment_positive_candidates = &augment_positive_candidates2;
      augment_negative_candidates = &augment_negative_candidates2;
      repetitive_seed_length = &repetitive_seed_length2;
    }
    uint32_t mm_count = minimizers->size();
    bool augment_flag = true;
    uint32_t candidate_num = positive_candidates->size();
    // if (positive_candidates->size() >= (uint32_t)max_seed_frequencies_[0]) {
    //  augment_flag = false;
    //} else {
    for (uint32_t i = 0; i < candidate_num; ++i) {
      if (positive_candidates->at(i).count >= mm_count / 2) {
        augment_flag = false;
        break;
      }
    }
    //}
    candidate_num = negative_candidates->size();
    if (augment_flag) {
      // if (negative_candidates->size() >= (uint32_t)max_seed_frequencies_[0])
      // {
      //  augment_flag = false;
      //} else {
      for (uint32_t i = 0; i < candidate_num; ++i) {
        if (negative_candidates->at(i).count >= mm_count / 2) {
          augment_flag = false;
          break;
        }
      }
      //}
    }

    if (augment_flag) {
      positive_hits->clear();
      negative_hits->clear();
      int positive_rescue_result = 0;
      int negative_rescue_result = 0;
      if (mate_positive_candidates->size() >
          0) {  // && mate_positive_candidates->size() <=
                // (uint32_t)max_seed_frequencies_[0]) {
        // std::cerr << "Supplement positive" << "\n";
        positive_rescue_result =
            index.GenerateCandidatesFromRepetitiveReadWithMateInfo(
                error_threshold_, *minimizers, repetitive_seed_length,
                negative_hits, augment_negative_candidates,
                mate_positive_candidates, kNegative, 2 * max_insert_size_);
      }
      if (mate_negative_candidates->size() >
          0) {  // && mate_negative_candidates->size() <=
                // (uint32_t)max_seed_frequencies_[0]) {
        // std::cerr << "Supplement negative" << "\n";
        negative_rescue_result =
            index.GenerateCandidatesFromRepetitiveReadWithMateInfo(
                error_threshold_, *minimizers, repetitive_seed_length,
                positive_hits, augment_positive_candidates,
                mate_negative_candidates, kPositive, 2 * max_insert_size_);
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
    MergeCandidates(positive_candidates1, augment_positive_candidates1,
                    positive_candidates1_buffer);
  }
  if (augment_negative_candidates1.size() > 0) {
    MergeCandidates(negative_candidates1, augment_negative_candidates1,
                    negative_candidates1_buffer);
  }
  if (augment_positive_candidates2.size() > 0) {
    MergeCandidates(positive_candidates2, augment_positive_candidates2,
                    positive_candidates2_buffer);
  }
  if (augment_negative_candidates2.size() > 0) {
    MergeCandidates(negative_candidates2, augment_negative_candidates2,
                    negative_candidates2_buffer);
  }
  return ret;
}

template <typename MappingRecord>
void Chromap<MappingRecord>::PostProcessingInLowMemory(
    uint32_t num_mappings_in_mem, uint32_t num_reference_sequences,
    const SequenceBatch &reference) {
  if (num_mappings_in_mem > 0) {
    TempMappingFileHandle<MappingRecord> temp_mapping_file_handle;
    temp_mapping_file_handle.file_path =
        mapping_output_file_path_ + ".temp" +
        std::to_string(temp_mapping_file_handles_.size());
    temp_mapping_file_handles_.emplace_back(temp_mapping_file_handle);
    SortOutputMappings(num_reference_sequences, &mappings_on_diff_ref_seqs_);
    // double output_temp_mapping_start_time = Chromap<>::GetRealTime();
    output_tools_.OutputTempMapping(temp_mapping_file_handle.file_path,
                                    num_reference_sequences,
                                    mappings_on_diff_ref_seqs_);
    // std::cerr << "Output temp mappings in " << Chromap<>::GetRealTime() -
    // output_temp_mapping_start_time << "s.\n";
    num_mappings_in_mem = 0;
    for (uint32_t i = 0; i < num_reference_sequences; ++i) {
      mappings_on_diff_ref_seqs_[i].clear();
    }
  }
  if (temp_mapping_file_handles_.size() == 0) {
    return;
  }
  double sort_and_dedupe_start_time = Chromap<>::GetRealTime();
  // Calculate block size and initialize
  uint64_t max_mem_size = 10 * ((uint64_t)1 << 30);
  if (mapping_output_format_ == MAPPINGFORMAT_SAM ||
      mapping_output_format_ == MAPPINGFORMAT_PAIRS ||
      mapping_output_format_ == MAPPINGFORMAT_PAF) {
    max_mem_size = (uint64_t)1 << 30;
  }
  for (size_t hi = 0; hi < temp_mapping_file_handles_.size(); ++hi) {
    temp_mapping_file_handles_[hi].block_size =
        max_mem_size / temp_mapping_file_handles_.size() /
        sizeof(MappingRecord);
    temp_mapping_file_handles_[hi].InitializeTempMappingLoading(
        num_reference_sequences);
    temp_mapping_file_handles_[hi].LoadTempMappingBlock(
        num_reference_sequences);
  }
  // Merge and dedupe
  bool all_merged = false;
  uint32_t last_rid = std::numeric_limits<uint32_t>::max();
  MappingRecord last_mapping;
  uint32_t dup_count = 0;
  uint64_t num_uni_mappings = 0;
  uint64_t num_multi_mappings = 0;
  uint64_t num_mappings_passing_filters = 0;
  std::vector<MappingRecord> temp_dups_for_bulk_level_dedup;
  temp_dups_for_bulk_level_dedup.reserve(255);
  while (!all_merged) {
    // Merge, dedupe and output
    // Find min first (sorted by rid and then barcode and then positions)
    size_t min_handle_index = temp_mapping_file_handles_.size();
    uint32_t min_rid = std::numeric_limits<uint32_t>::max();
    for (size_t hi = 0; hi < temp_mapping_file_handles_.size(); ++hi) {
      TempMappingFileHandle<MappingRecord> &current_handle =
          temp_mapping_file_handles_[hi];
      if (!current_handle.all_loaded) {
        if (min_handle_index == temp_mapping_file_handles_.size() ||
            current_handle.current_rid < min_rid ||
            (current_handle.current_rid == min_rid &&
             current_handle.mappings[current_handle.current_mapping_index] <
                 temp_mapping_file_handles_[min_handle_index]
                     .mappings[temp_mapping_file_handles_[min_handle_index]
                                   .current_mapping_index])) {
          min_handle_index = hi;
          min_rid = current_handle.current_rid;
        }
      }
    }
    // Append current min to mappings if not a duplicate
    if (!remove_pcr_duplicates_ ||
        min_handle_index != temp_mapping_file_handles_.size()) {
      MappingRecord &current_min_mapping =
          temp_mapping_file_handles_[min_handle_index]
              .mappings[temp_mapping_file_handles_[min_handle_index]
                            .current_mapping_index];
      if (remove_pcr_duplicates_ && last_rid == min_rid &&
          (current_min_mapping == last_mapping ||
           (remove_pcr_duplicates_at_bulk_level_ &&
            current_min_mapping.HasSamePosition(last_mapping)))) {
        ++dup_count;
        if (!is_bulk_data_ && remove_pcr_duplicates_at_bulk_level_) {
          if (!temp_dups_for_bulk_level_dedup.empty() &&
              current_min_mapping == temp_dups_for_bulk_level_dedup.back()) {
            current_min_mapping.num_dups =
                temp_dups_for_bulk_level_dedup.back().num_dups + 1;
            temp_dups_for_bulk_level_dedup.back() = current_min_mapping;
          } else {
            temp_dups_for_bulk_level_dedup.push_back(current_min_mapping);
            temp_dups_for_bulk_level_dedup.back().num_dups = 1;
          }
        }
      } else {
        if (dup_count > 0) {
          if (!is_bulk_data_ && remove_pcr_duplicates_at_bulk_level_ &&
              temp_dups_for_bulk_level_dedup.size() > 0) {
            // Find the best barcode, break ties first by the number of the
            // barcodes in the dups, then by the barcode abundance
            last_mapping = temp_dups_for_bulk_level_dedup[0];
            khiter_t barcode_whitelist_lookup_table_iterator =
                kh_get(k64_seq, barcode_whitelist_lookup_table_,
                       last_mapping.GetBarcode());
            double last_mapping_barcode_abundance =
                kh_value(barcode_whitelist_lookup_table_,
                         barcode_whitelist_lookup_table_iterator) /
                (double)num_sample_barcodes_;
            for (uint32_t bulk_dup_i = 1;
                 bulk_dup_i < temp_dups_for_bulk_level_dedup.size();
                 ++bulk_dup_i) {
              barcode_whitelist_lookup_table_iterator = kh_get(
                  k64_seq, barcode_whitelist_lookup_table_,
                  temp_dups_for_bulk_level_dedup[bulk_dup_i].GetBarcode());
              double current_mapping_barcode_abundance =
                  kh_value(barcode_whitelist_lookup_table_,
                           barcode_whitelist_lookup_table_iterator) /
                  (double)num_sample_barcodes_;
              if (temp_dups_for_bulk_level_dedup[bulk_dup_i].num_dups >
                      last_mapping.num_dups ||
                  (temp_dups_for_bulk_level_dedup[bulk_dup_i].num_dups ==
                       last_mapping.num_dups &&
                   current_mapping_barcode_abundance >
                       last_mapping_barcode_abundance)) {
                last_mapping = temp_dups_for_bulk_level_dedup[bulk_dup_i];
                last_mapping_barcode_abundance =
                    current_mapping_barcode_abundance;
              }
            }
            temp_dups_for_bulk_level_dedup.clear();
          }
          if (last_mapping.mapq >= mapq_threshold_) {
            // if (allocate_multi_mappings_ || (only_output_unique_mappings_ &&
            // last_mapping.is_unique == 1)) {
            last_mapping.num_dups = std::min(
                (uint32_t)std::numeric_limits<uint8_t>::max(), dup_count);
            if (Tn5_shift_) {
              // last_mapping.fragment_start_position += 4;
              // last_mapping.positive_alignment_length -= 4;
              // last_mapping.fragment_length -= 9;
              // last_mapping.negative_alignment_length -= 5;
              last_mapping.Tn5Shift();
            }
            output_tools_.AppendMapping(last_rid, reference, last_mapping);
            ++num_mappings_passing_filters;
            //}
          }
          if (last_mapping.is_unique == 1) {
            ++num_uni_mappings;
          } else {
            ++num_multi_mappings;
          }
        }
        last_mapping = current_min_mapping;
        last_rid = min_rid;
        dup_count = 1;
        if (remove_pcr_duplicates_at_bulk_level_) {
          temp_dups_for_bulk_level_dedup.push_back(current_min_mapping);
          temp_dups_for_bulk_level_dedup.back().num_dups = 1;
        }
      }
      temp_mapping_file_handles_[min_handle_index].Next(
          num_reference_sequences);
    }
    // Check if all are merged.
    all_merged = true;
    for (size_t hi = 0; hi < temp_mapping_file_handles_.size(); ++hi) {
      if (!temp_mapping_file_handles_[hi].all_loaded) {
        all_merged = false;
      }
    }
  }
  if (last_mapping.mapq >= mapq_threshold_) {
    if (!is_bulk_data_ && remove_pcr_duplicates_at_bulk_level_ &&
        temp_dups_for_bulk_level_dedup.size() > 0) {
      // Find the best barcode, break ties first by the number of the barcodes
      // in the dups, then by the barcode abundance
      last_mapping = temp_dups_for_bulk_level_dedup[0];
      khiter_t barcode_whitelist_lookup_table_iterator = kh_get(
          k64_seq, barcode_whitelist_lookup_table_, last_mapping.GetBarcode());
      double last_mapping_barcode_abundance =
          kh_value(barcode_whitelist_lookup_table_,
                   barcode_whitelist_lookup_table_iterator) /
          (double)num_sample_barcodes_;
      for (uint32_t bulk_dup_i = 1;
           bulk_dup_i < temp_dups_for_bulk_level_dedup.size(); ++bulk_dup_i) {
        barcode_whitelist_lookup_table_iterator =
            kh_get(k64_seq, barcode_whitelist_lookup_table_,
                   temp_dups_for_bulk_level_dedup[bulk_dup_i].GetBarcode());
        double current_mapping_barcode_abundance =
            kh_value(barcode_whitelist_lookup_table_,
                     barcode_whitelist_lookup_table_iterator) /
            (double)num_sample_barcodes_;
        if (temp_dups_for_bulk_level_dedup[bulk_dup_i].num_dups >
                last_mapping.num_dups ||
            (temp_dups_for_bulk_level_dedup[bulk_dup_i].num_dups ==
                 last_mapping.num_dups &&
             current_mapping_barcode_abundance >
                 last_mapping_barcode_abundance)) {
          last_mapping = temp_dups_for_bulk_level_dedup[bulk_dup_i];
          last_mapping_barcode_abundance = current_mapping_barcode_abundance;
        }
      }
      temp_dups_for_bulk_level_dedup.clear();
    }
    // if (allocate_multi_mappings_ || (only_output_unique_mappings_ &&
    // last_mapping.is_unique == 1)) {
    last_mapping.num_dups =
        std::min((uint32_t)std::numeric_limits<uint8_t>::max(), dup_count);
    if (Tn5_shift_) {
      // last_mapping.fragment_start_position += 4;
      // last_mapping.positive_alignment_length -= 4;
      // last_mapping.fragment_length -= 9;
      // last_mapping.negative_alignment_length -= 5;
      last_mapping.Tn5Shift();
    }
    output_tools_.AppendMapping(last_rid, reference, last_mapping);
    ++num_mappings_passing_filters;
    //}
  }
  if (last_mapping.is_unique == 1) {
    ++num_uni_mappings;
  } else {
    ++num_multi_mappings;
  }
  // Delete temp files
  for (size_t hi = 0; hi < temp_mapping_file_handles_.size(); ++hi) {
    temp_mapping_file_handles_[hi].FinalizeTempMappingLoading();
    remove(temp_mapping_file_handles_[hi].file_path.c_str());
  }
  if (remove_pcr_duplicates_) {
    std::cerr << "Sorted, deduped and outputed mappings in "
              << Chromap<>::GetRealTime() - sort_and_dedupe_start_time
              << "s.\n";
  } else {
    std::cerr << "Sorted and outputed mappings in "
              << Chromap<>::GetRealTime() - sort_and_dedupe_start_time
              << "s.\n";
  }
  std::cerr << "# uni-mappings: " << num_uni_mappings
            << ", # multi-mappings: " << num_multi_mappings
            << ", total: " << num_uni_mappings + num_multi_mappings << ".\n";
  std::cerr << "Number of output mappings (passed filters): "
            << num_mappings_passing_filters << "\n";
}

template <typename MappingRecord>
void Chromap<MappingRecord>::GetRidRank(const std::string rid_order_path,
                                        const SequenceBatch &reference,
                                        std::vector<int> &rid_rank) {
  int ref_size = reference.GetSequenceBatchSize();
  int i = 0;
  rid_rank.resize(ref_size);
  for (i = 0; i < ref_size; ++i) {
    rid_rank[i] = i;
  }

  if (rid_order_path.length() == 0) {
    return;
  }

  std::map<std::string, int> rname_to_rank;
  std::ifstream file_stream(rid_order_path);
  std::string line;
  i = 0;
  while (getline(file_stream, line)) {
    rname_to_rank[line] = i;
    i += 1;
  }
  file_stream.close();

  // First put the chrosomes in the list provided by user
  for (i = 0; i < ref_size; ++i) {
    std::string rname(reference.GetSequenceNameAt(i));
    if (rname_to_rank.find(rname) != rname_to_rank.end()) {
      rid_rank[i] = rname_to_rank[rname];
    } else {
      rid_rank[i] = -1;
    }
  }

  // we may have some rank without any rid associated with, this helps if
  // cutstom list contains rid not in the reference`
  int k = rname_to_rank.size();
  // Put the remaining chrosomes
  for (i = 0; i < ref_size; ++i) {
    if (rid_rank[i] == -1) {
      rid_rank[i] = k;
      ++k;
    }
  }

  if (k > ref_size) {
    chromap::Chromap<>::ExitWithMessage(
        "Unknown chromsome names found in chromosome order file");
  }

  /*for (i = 0 ; i < ref_size; ++i) {
          std::cerr<<rid_rank_[i]<<"\n";
  }*/
}

template <typename MappingRecord>
void Chromap<MappingRecord>::RerankCandidatesRid(
    std::vector<Candidate> &candidates) {
  int i;
  int size = candidates.size();
  for (i = 0; i < size; ++i) {
    uint64_t rid = (uint32_t)(candidates[i].position >> 32);
    rid = custom_rid_rank_[rid];
    candidates[i].position =
        (candidates[i].position & (uint64_t)0xffffffff) | (rid << 32);
  }
}

template <typename MappingRecord>
void Chromap<MappingRecord>::MapPairedEndReads() {
  double real_start_time = Chromap<>::GetRealTime();
  // Load reference
  SequenceBatch reference;
  reference.InitializeLoading(reference_file_path_);
  uint32_t num_reference_sequences = reference.LoadAllSequences();
  if (custom_rid_order_path_.length() > 0) {
    GetRidRank(custom_rid_order_path_, reference, custom_rid_rank_);
    reference.ReorderSequences(custom_rid_rank_);
  }
  if (mapping_output_format_ == MAPPINGFORMAT_PAIRS) {
    GetRidRank(pairs_custom_rid_order_path_, reference, pairs_custom_rid_rank_);
  }

  // Load index
  Index index(min_num_seeds_required_for_mapping_, max_seed_frequencies_,
              index_file_path_);
  index.Load();
  kmer_size_ = index.GetKmerSize();
  window_size_ = index.GetWindowSize();
  // index.Statistics(num_sequences, reference);
  // Initialize read batches
  SequenceBatch read_batch1(read_batch_size_);
  SequenceBatch read_batch2(read_batch_size_);
  SequenceBatch barcode_batch(read_batch_size_);
  SequenceBatch read_batch1_for_loading(read_batch_size_);
  SequenceBatch read_batch2_for_loading(read_batch_size_);
  SequenceBatch barcode_batch_for_loading(read_batch_size_);
  read_batch1.SetSeqEffectiveRange(read1_format_[0], read1_format_[1]);
  read_batch2.SetSeqEffectiveRange(read2_format_[0], read2_format_[1]);
  barcode_batch.SetSeqEffectiveRange(barcode_format_[0], barcode_format_[1]);
  read_batch1_for_loading.SetSeqEffectiveRange(read1_format_[0],
                                               read1_format_[1]);
  read_batch2_for_loading.SetSeqEffectiveRange(read2_format_[0],
                                               read2_format_[1]);
  barcode_batch_for_loading.SetSeqEffectiveRange(barcode_format_[0],
                                                 barcode_format_[1]);
  // Initialize cache
  mm_cache mm_to_candidates_cache(2000003);
  mm_to_candidates_cache.SetKmerLength(kmer_size_);
  struct _mm_history *mm_history1 = new struct _mm_history[read_batch_size_];
  struct _mm_history *mm_history2 = new struct _mm_history[read_batch_size_];
  // Initialize mapping container
  mappings_on_diff_ref_seqs_.reserve(num_reference_sequences);
  deduped_mappings_on_diff_ref_seqs_.reserve(num_reference_sequences);
  for (uint32_t i = 0; i < num_reference_sequences; ++i) {
    mappings_on_diff_ref_seqs_.emplace_back(std::vector<MappingRecord>());
    deduped_mappings_on_diff_ref_seqs_.emplace_back(
        std::vector<MappingRecord>());
  }
  // Preprocess barcodes for single cell data
  if (!is_bulk_data_) {
    if (!barcode_whitelist_file_path_.empty()) {
      LoadBarcodeWhitelist();
      ComputeBarcodeAbundance(initial_num_sample_barcodes_);
    }
  }

  // Initialize output tools
  if (mapping_output_format_ == MAPPINGFORMAT_PAIRS) {
    output_tools_.SetPairsCustomRidRank(pairs_custom_rid_rank_);
  }

  output_tools_.InitializeMappingOutput(
      barcode_length_, mapping_output_file_path_, mapping_output_format_);
  output_tools_.OutputHeader(num_reference_sequences, reference);

  uint32_t num_mappings_in_mem = 0;
  uint64_t max_num_mappings_in_mem =
      1 * ((uint64_t)1 << 30) / sizeof(MappingRecord);
  if (mapping_output_format_ == MAPPINGFORMAT_SAM ||
      mapping_output_format_ == MAPPINGFORMAT_PAF ||
      mapping_output_format_ == MAPPINGFORMAT_PAIRS) {
    max_num_mappings_in_mem = 1 * ((uint64_t)1 << 29) / sizeof(MappingRecord);
  }
  static uint64_t thread_num_candidates = 0;
  static uint64_t thread_num_mappings = 0;
  static uint64_t thread_num_mapped_reads = 0;
  static uint64_t thread_num_uniquely_mapped_reads = 0;
  static uint64_t thread_num_barcode_in_whitelist = 0;
  static uint64_t thread_num_corrected_barcode = 0;
#pragma omp threadprivate(                                               \
    thread_num_candidates, thread_num_mappings, thread_num_mapped_reads, \
    thread_num_uniquely_mapped_reads, thread_num_barcode_in_whitelist,   \
    thread_num_corrected_barcode)
  double real_start_mapping_time = Chromap<>::GetRealTime();
  for (size_t read_file_index = 0; read_file_index < read_file1_paths_.size();
       ++read_file_index) {
    read_batch1_for_loading.InitializeLoading(
        read_file1_paths_[read_file_index]);
    read_batch2_for_loading.InitializeLoading(
        read_file2_paths_[read_file_index]);
    if (!is_bulk_data_) {
      barcode_batch_for_loading.InitializeLoading(
          barcode_file_paths_[read_file_index]);
    }
    uint32_t num_loaded_pairs_for_loading = 0;
    uint32_t num_loaded_pairs = LoadPairedEndReadsWithBarcodes(
        &read_batch1_for_loading, &read_batch2_for_loading,
        &barcode_batch_for_loading);
    read_batch1_for_loading.SwapSequenceBatch(read_batch1);
    read_batch2_for_loading.SwapSequenceBatch(read_batch2);
    barcode_batch_for_loading.SwapSequenceBatch(barcode_batch);
    std::vector<std::vector<std::vector<MappingRecord>>>
        mappings_on_diff_ref_seqs_for_diff_threads;
    std::vector<std::vector<std::vector<MappingRecord>>>
        mappings_on_diff_ref_seqs_for_diff_threads_for_saving;
    mappings_on_diff_ref_seqs_for_diff_threads.reserve(num_threads_);
    mappings_on_diff_ref_seqs_for_diff_threads_for_saving.reserve(num_threads_);
    for (int ti = 0; ti < num_threads_; ++ti) {
      mappings_on_diff_ref_seqs_for_diff_threads.emplace_back(
          std::vector<std::vector<MappingRecord>>(num_reference_sequences));
      mappings_on_diff_ref_seqs_for_diff_threads_for_saving.emplace_back(
          std::vector<std::vector<MappingRecord>>(num_reference_sequences));
      for (uint32_t i = 0; i < num_reference_sequences; ++i) {
        mappings_on_diff_ref_seqs_for_diff_threads[ti][i].reserve(
            (num_loaded_pairs +
             num_loaded_pairs / 1000 * max_num_best_mappings_) /
            num_threads_ / num_reference_sequences);
        mappings_on_diff_ref_seqs_for_diff_threads_for_saving[ti][i].reserve(
            (num_loaded_pairs +
             num_loaded_pairs / 1000 * max_num_best_mappings_) /
            num_threads_ / num_reference_sequences);
      }
    }
#pragma omp parallel default(none) shared(output_mappings_not_in_whitelist_, barcode_whitelist_file_path_, trim_adapters_, window_size_, custom_rid_order_path_, split_alignment_, error_threshold_, max_num_best_mappings_, num_threads_, num_reads_, low_memory_mode_, reference, index, read_batch1, read_batch2, barcode_batch, read_batch1_for_loading, read_batch2_for_loading, barcode_batch_for_loading, output_tools_, std::cerr, num_loaded_pairs_for_loading, num_loaded_pairs, num_reference_sequences, mappings_on_diff_ref_seqs_for_diff_threads, mappings_on_diff_ref_seqs_for_diff_threads_for_saving, mappings_on_diff_ref_seqs_, mapping_output_file_path_, num_mappings_in_mem, max_num_mappings_in_mem, temp_mapping_file_handles_, mm_to_candidates_cache, mm_history1, mm_history2) num_threads(num_threads_) reduction(+:num_candidates_, num_mappings_, num_mapped_reads_, num_uniquely_mapped_reads_, num_barcode_in_whitelist_, num_corrected_barcode_)
    {
      thread_num_candidates = 0;
      thread_num_mappings = 0;
      thread_num_mapped_reads = 0;
      thread_num_uniquely_mapped_reads = 0;
      thread_num_barcode_in_whitelist = 0;
      thread_num_corrected_barcode = 0;
      std::vector<std::pair<uint64_t, uint64_t>> minimizers1;
      std::vector<std::pair<uint64_t, uint64_t>> minimizers2;
      std::vector<uint64_t> positive_hits1;
      std::vector<uint64_t> positive_hits2;
      std::vector<uint64_t> negative_hits1;
      std::vector<uint64_t> negative_hits2;
      positive_hits1.reserve(max_seed_frequencies_[0]);
      positive_hits2.reserve(max_seed_frequencies_[0]);
      negative_hits1.reserve(max_seed_frequencies_[0]);
      negative_hits2.reserve(max_seed_frequencies_[0]);
      std::vector<Candidate> positive_candidates1;
      std::vector<Candidate> positive_candidates2;
      std::vector<Candidate> negative_candidates1;
      std::vector<Candidate> negative_candidates2;
      positive_candidates1.reserve(max_seed_frequencies_[0]);
      positive_candidates2.reserve(max_seed_frequencies_[0]);
      negative_candidates1.reserve(max_seed_frequencies_[0]);
      negative_candidates2.reserve(max_seed_frequencies_[0]);
      std::vector<Candidate> positive_candidates1_buffer;
      std::vector<Candidate> positive_candidates2_buffer;
      std::vector<Candidate> negative_candidates1_buffer;
      std::vector<Candidate> negative_candidates2_buffer;
      positive_candidates1_buffer.reserve(max_seed_frequencies_[0]);
      positive_candidates2_buffer.reserve(max_seed_frequencies_[0]);
      negative_candidates1_buffer.reserve(max_seed_frequencies_[0]);
      negative_candidates2_buffer.reserve(max_seed_frequencies_[0]);
      std::vector<std::pair<int, uint64_t>> positive_mappings1;
      std::vector<std::pair<int, uint64_t>> positive_mappings2;
      std::vector<std::pair<int, uint64_t>> negative_mappings1;
      std::vector<std::pair<int, uint64_t>> negative_mappings2;
      positive_mappings1.reserve(max_seed_frequencies_[0]);
      positive_mappings2.reserve(max_seed_frequencies_[0]);
      negative_mappings1.reserve(max_seed_frequencies_[0]);
      negative_mappings2.reserve(max_seed_frequencies_[0]);
      std::vector<int> positive_split_sites1;
      std::vector<int> negative_split_sites1;
      std::vector<int> positive_split_sites2;
      std::vector<int> negative_split_sites2;
      positive_split_sites1.reserve(max_seed_frequencies_[0]);
      negative_split_sites1.reserve(max_seed_frequencies_[0]);
      positive_split_sites2.reserve(max_seed_frequencies_[0]);
      negative_split_sites2.reserve(max_seed_frequencies_[0]);
      std::vector<std::pair<uint32_t, uint32_t>> F1R2_best_mappings;
      std::vector<std::pair<uint32_t, uint32_t>> F2R1_best_mappings;
      std::vector<std::pair<uint32_t, uint32_t>> F1F2_best_mappings;
      std::vector<std::pair<uint32_t, uint32_t>> R1R2_best_mappings;
      F1R2_best_mappings.reserve(max_seed_frequencies_[0]);
      F2R1_best_mappings.reserve(max_seed_frequencies_[0]);
      if (split_alignment_) {
        F1F2_best_mappings.reserve(max_seed_frequencies_[0]);
        R1R2_best_mappings.reserve(max_seed_frequencies_[0]);
      }
      // we will use reservoir sampling
      std::vector<int> best_mapping_indices(max_num_best_mappings_);
      std::mt19937 generator(11);
#pragma omp single
      {
        double real_batch_start_time = Chromap<>::GetRealTime();
        while (num_loaded_pairs > 0) {
          num_reads_ += num_loaded_pairs;
          num_reads_ += num_loaded_pairs;
#pragma omp task
          {
            num_loaded_pairs_for_loading = LoadPairedEndReadsWithBarcodes(
                &read_batch1_for_loading, &read_batch2_for_loading,
                &barcode_batch_for_loading);
          }  // end of openmp loading task
          int grain_size = 5000;
#pragma omp taskloop grainsize(grain_size)  // num_tasks(num_threads_* 50)
          //#pragma omp taskloop num_tasks(num_threads_* num_threads_)
          for (uint32_t pair_index = 0; pair_index < num_loaded_pairs;
               ++pair_index) {
            bool current_barcode_is_whitelisted = true;
            if (!barcode_whitelist_file_path_.empty()) {
              current_barcode_is_whitelisted = CorrectBarcodeAt(
                  pair_index, &barcode_batch, &thread_num_barcode_in_whitelist,
                  &thread_num_corrected_barcode);
            }
            if (current_barcode_is_whitelisted ||
                output_mappings_not_in_whitelist_) {
              read_batch1.PrepareNegativeSequenceAt(pair_index);
              read_batch2.PrepareNegativeSequenceAt(pair_index);
              // std::cerr << pair_index<<"
              // "<<read_batch1.GetSequenceNameAt(pair_index) << "\n";
              if (trim_adapters_) {
                TrimAdapterForPairedEndRead(pair_index, &read_batch1,
                                            &read_batch2);
              }
              minimizers1.clear();
              minimizers2.clear();
              minimizers1.reserve(read_batch1.GetSequenceLengthAt(pair_index) /
                                  window_size_ * 2);
              minimizers2.reserve(read_batch2.GetSequenceLengthAt(pair_index) /
                                  window_size_ * 2);
              index.GenerateMinimizerSketch(read_batch1, pair_index,
                                            &minimizers1);
              index.GenerateMinimizerSketch(read_batch2, pair_index,
                                            &minimizers2);
              // std::cerr << "m1" << " " << minimizers1.size() << "\n";
              // for (auto &mi : minimizers1) {
              //  std::cerr << (mi.second >> 33) << " " << (uint32_t) (mi.second
              //  >> 1) << "\n";
              //}
              // std::cerr << "m2" << " " << minimizers2.size() << "\n";
              // for (auto &mi : minimizers2) {
              //  std::cerr << (mi.second >> 33) << " " << (uint32_t) (mi.second
              //  >> 1) << "\n";
              //}
              if (minimizers1.size() != 0 && minimizers2.size() != 0) {
                positive_hits1.clear();
                positive_hits2.clear();
                negative_hits1.clear();
                negative_hits2.clear();
                positive_candidates1.clear();
                positive_candidates2.clear();
                negative_candidates1.clear();
                negative_candidates2.clear();
                positive_candidates1_buffer.clear();
                positive_candidates2_buffer.clear();
                negative_candidates1_buffer.clear();
                negative_candidates2_buffer.clear();
                uint32_t repetitive_seed_length1 = 0;
                uint32_t repetitive_seed_length2 = 0;
                // Generate candidates
                if (mm_to_candidates_cache.Query(
                        minimizers1, positive_candidates1, negative_candidates1,
                        repetitive_seed_length1,
                        read_batch1.GetSequenceLengthAt(pair_index)) == -1) {
                  index.GenerateCandidates(
                      error_threshold_, minimizers1, &repetitive_seed_length1,
                      &positive_hits1, &negative_hits1, &positive_candidates1,
                      &negative_candidates1);
                }
                uint32_t current_num_candidates1 =
                    positive_candidates1.size() + negative_candidates1.size();
                if (mm_to_candidates_cache.Query(
                        minimizers2, positive_candidates2, negative_candidates2,
                        repetitive_seed_length2,
                        read_batch2.GetSequenceLengthAt(pair_index)) == -1) {
                  index.GenerateCandidates(
                      error_threshold_, minimizers2, &repetitive_seed_length2,
                      &positive_hits2, &negative_hits2, &positive_candidates2,
                      &negative_candidates2);
                }
                // std::cerr << "p1" << "\n";
                // for (auto &ci : positive_candidates1) {
                //  std::cerr << (ci.position >> 32) << " " << (uint32_t)
                //  ci.position << " " << (int)ci.count << "\n";
                //}
                // std::cerr << "n1" << "\n";
                // for (auto &ci : negative_candidates1) {
                //  std::cerr << (ci.position >> 32) << " " << (uint32_t)
                //  ci.position << " " << (int)ci.count << "\n";
                //}
                // std::cerr << "p2" << "\n";
                // for (auto &ci : positive_candidates2) {
                //  std::cerr << (ci.position >> 32) << " " << (uint32_t)
                //  ci.position << " " << (int)ci.count << "\n";
                //}
                // std::cerr << "n2" << "\n";
                // for (auto &ci : negative_candidates2) {
                //  std::cerr << (ci.position >> 32) << " " << (uint32_t)
                //  ci.position << " " << (int)ci.count << "\n";
                //}
                // std::cerr << "After generation, #pc1: " <<
                // positive_candidates1.size() << ", #nc1: " <<
                // negative_candidates1.size() << ", #pc2: " <<
                // positive_candidates2.size() << ", #nc2: " <<
                // negative_candidates2.size() << "\n";
                uint32_t current_num_candidates2 =
                    positive_candidates2.size() + negative_candidates2.size();
                if (pair_index < num_loaded_pairs &&
                    (pair_index < num_loaded_pairs / num_threads_ ||
                     num_reads_ <= 5000000)) {
                  mm_history1[pair_index].timestamp =
                      mm_history2[pair_index].timestamp = num_reads_;
                  mm_history1[pair_index].minimizers = minimizers1;
                  mm_history1[pair_index].positive_candidates =
                      positive_candidates1;
                  mm_history1[pair_index].negative_candidates =
                      negative_candidates1;
                  mm_history1[pair_index].repetitive_seed_length =
                      repetitive_seed_length1;
                  mm_history2[pair_index].minimizers = minimizers2;
                  mm_history2[pair_index].positive_candidates =
                      positive_candidates2;
                  mm_history2[pair_index].negative_candidates =
                      negative_candidates2;
                  mm_history2[pair_index].repetitive_seed_length =
                      repetitive_seed_length2;
                }
                // Test whether we need to augment the candidate list with mate
                // information.
                // std::cerr << "before supplement" << "\n";
                // std::cerr << "p1" << "\n";
                // for (auto &ci : positive_candidates1) {
                //  std::cerr << (ci.position >> 32) << " " << (uint32_t)
                //  ci.position << " " << (int)ci.count << "\n";
                //}
                // std::cerr << "n1" << "\n";
                // for (auto &ci : negative_candidates1) {
                //  std::cerr << (ci.position >> 32) << " " << (uint32_t)
                //  ci.position << " " << (int)ci.count << "\n";
                //}
                // std::cerr << "p2" << "\n";
                // for (auto &ci : positive_candidates2) {
                //  std::cerr << (ci.position >> 32) << " " << (uint32_t)
                //  ci.position << " " << (int)ci.count << "\n";
                //}
                // std::cerr << "n2" << "\n";
                // for (auto &ci : negative_candidates2) {
                //  std::cerr << (ci.position >> 32) << " " << (uint32_t)
                //  ci.position << " " << (int)ci.count << "\n";
                //}
                int supplementCandidateResult = 0;
                if (!split_alignment_) {
                  supplementCandidateResult = SupplementCandidates(
                      index, repetitive_seed_length1, repetitive_seed_length2,
                      minimizers1, minimizers2, positive_hits1, positive_hits2,
                      positive_candidates1, positive_candidates2,
                      positive_candidates1_buffer, positive_candidates2_buffer,
                      negative_hits1, negative_hits2, negative_candidates1,
                      negative_candidates2, negative_candidates1_buffer,
                      negative_candidates2_buffer);
                  current_num_candidates1 =
                      positive_candidates1.size() + negative_candidates1.size();
                  current_num_candidates2 =
                      positive_candidates2.size() + negative_candidates2.size();
                }
                if (current_num_candidates1 > 0 &&
                    current_num_candidates2 > 0 && !split_alignment_) {
                  positive_candidates1.swap(positive_candidates1_buffer);
                  negative_candidates1.swap(negative_candidates1_buffer);
                  positive_candidates2.swap(positive_candidates2_buffer);
                  negative_candidates2.swap(negative_candidates2_buffer);
                  positive_candidates1.clear();
                  positive_candidates2.clear();
                  negative_candidates1.clear();
                  negative_candidates2.clear();
                  // Paired-end filter
                  // std::cerr << "p1" << "\n";
                  // for (auto &ci : positive_candidates1_buffer) {
                  //  std::cerr << (ci.position >> 32) << " " << (uint32_t)
                  //  ci.position << " " << (int)ci.count << "\n";
                  //}
                  // std::cerr << "n1" << "\n";
                  // for (auto &ci : negative_candidates1_buffer) {
                  //  std::cerr << (ci.position >> 32) << " " << (uint32_t)
                  //  ci.position << " " << (int)ci.count << "\n";
                  //}
                  // std::cerr << "p2" << "\n";
                  // for (auto &ci : positive_candidates2_buffer) {
                  //  std::cerr << (ci.position >> 32) << " " << (uint32_t)
                  //  ci.position << " " << (int)ci.count << "\n";
                  //}
                  // std::cerr << "n2" << "\n";
                  // for (auto &ci : negative_candidates2_buffer) {
                  //  std::cerr << (ci.position >> 32) << " " << (uint32_t)
                  //  ci.position << " " << (int)ci.count << "\n";
                  //}
                  // std::cerr << "#pc1: " << positive_candidates1_buffer.size()
                  // << ", #nc1: " << negative_candidates1_buffer.size() << ",
                  // #pc2: " << positive_candidates2_buffer.size() << ", #nc2: "
                  // << negative_candidates2_buffer.size() << "\n";
                  ReduceCandidatesForPairedEndRead(
                      positive_candidates1_buffer, negative_candidates1_buffer,
                      positive_candidates2_buffer, negative_candidates2_buffer,
                      &positive_candidates1, &negative_candidates1,
                      &positive_candidates2, &negative_candidates2);
                  // std::cerr << "p1" << "\n";
                  // for (auto &ci : positive_candidates1) {
                  //  std::cerr << (ci.position >> 32) << " " << (uint32_t)
                  //  ci.position << " " << (int)ci.count << "\n";
                  //}
                  // std::cerr << "n1" << "\n";
                  // for (auto &ci : negative_candidates1) {
                  //  std::cerr << (ci.position >> 32) << " " << (uint32_t)
                  //  ci.position << " " << (int)ci.count << "\n";
                  //}
                  // std::cerr << "p2" << "\n";
                  // for (auto &ci : positive_candidates2) {
                  //  std::cerr << (ci.position >> 32) << " " << (uint32_t)
                  //  ci.position << " " << (int)ci.count << "\n";
                  //}
                  // std::cerr << "n2" << "\n";
                  // for (auto &ci : negative_candidates2) {
                  //  std::cerr << (ci.position >> 32) << " " << (uint32_t)
                  //  ci.position << " " << (int)ci.count << "\n";
                  //}
                  // std::cerr << "After pe filter, #pc1: " <<
                  // positive_candidates1.size() << ", #nc1: " <<
                  // negative_candidates1.size() << ", #pc2: " <<
                  // positive_candidates2.size() << ", #nc2: " <<
                  // negative_candidates2.size() << "\n";
                  current_num_candidates1 =
                      positive_candidates1.size() + negative_candidates1.size();
                  current_num_candidates2 =
                      positive_candidates2.size() + negative_candidates2.size();
                }
                // Verify candidates
                if (current_num_candidates1 > 0 &&
                    current_num_candidates2 > 0) {
                  thread_num_candidates += positive_candidates1.size() +
                                           positive_candidates2.size() +
                                           negative_candidates1.size() +
                                           negative_candidates2.size();
                  if (custom_rid_order_path_.length() > 0) {
                    RerankCandidatesRid(positive_candidates1);
                    RerankCandidatesRid(negative_candidates1);
                    RerankCandidatesRid(positive_candidates2);
                    RerankCandidatesRid(negative_candidates2);
                  }
                  positive_mappings1.clear();
                  positive_mappings2.clear();
                  negative_mappings1.clear();
                  negative_mappings2.clear();
                  positive_split_sites1.clear();
                  negative_split_sites1.clear();
                  positive_split_sites2.clear();
                  negative_split_sites2.clear();
                  int min_num_errors1, second_min_num_errors1;
                  int num_best_mappings1, num_second_best_mappings1;
                  int min_num_errors2, second_min_num_errors2;
                  int num_best_mappings2, num_second_best_mappings2;
                  VerifyCandidates(read_batch1, pair_index, reference,
                                   minimizers1, positive_candidates1,
                                   negative_candidates1, &positive_mappings1,
                                   &positive_split_sites1, &negative_mappings1,
                                   &negative_split_sites1, &min_num_errors1,
                                   &num_best_mappings1, &second_min_num_errors1,
                                   &num_second_best_mappings1);
                  uint32_t current_num_mappings1 =
                      positive_mappings1.size() + negative_mappings1.size();
                  VerifyCandidates(read_batch2, pair_index, reference,
                                   minimizers2, positive_candidates2,
                                   negative_candidates2, &positive_mappings2,
                                   &positive_split_sites2, &negative_mappings2,
                                   &negative_split_sites2, &min_num_errors2,
                                   &num_best_mappings2, &second_min_num_errors2,
                                   &num_second_best_mappings2);
                  // std::cerr << "me1: " << min_num_errors1 << " me2: " <<
                  // min_num_errors2 << "\n";
                  uint32_t current_num_mappings2 =
                      positive_mappings2.size() + negative_mappings2.size();
                  if (current_num_mappings1 > 0 && current_num_mappings2 > 0) {
                    int min_sum_errors, second_min_sum_errors;
                    int num_best_mappings, num_second_best_mappings;
                    F1R2_best_mappings.clear();
                    F2R1_best_mappings.clear();
                    if (split_alignment_) {
                      F1F2_best_mappings.clear();
                      R1R2_best_mappings.clear();
                    }
                    std::vector<std::vector<MappingRecord>>
                        &mappings_on_diff_ref_seqs =
                            mappings_on_diff_ref_seqs_for_diff_threads
                                [omp_get_thread_num()];
                    if (!split_alignment_) {
                      // GenerateBestMappingsForPairedEndRead assumes the
                      // mappings are sorted by coordinate for non split
                      // alignments In split alignment, we don't want to sort
                      // and this keeps mapping and split_sites vectors
                      // consistent.
                      std::sort(positive_mappings1.begin(),
                                positive_mappings1.end(),
                                [](const std::pair<int, uint64_t> &a,
                                   const std::pair<int, uint64_t> &b) {
                                  return a.second < b.second;
                                });
                      std::sort(positive_mappings2.begin(),
                                positive_mappings2.end(),
                                [](const std::pair<int, uint64_t> &a,
                                   const std::pair<int, uint64_t> &b) {
                                  return a.second < b.second;
                                });
                      std::sort(negative_mappings1.begin(),
                                negative_mappings1.end(),
                                [](const std::pair<int, uint64_t> &a,
                                   const std::pair<int, uint64_t> &b) {
                                  return a.second < b.second;
                                });
                      std::sort(negative_mappings2.begin(),
                                negative_mappings2.end(),
                                [](const std::pair<int, uint64_t> &a,
                                   const std::pair<int, uint64_t> &b) {
                                  return a.second < b.second;
                                });
                    }

                    int force_mapq = -1;
                    if (supplementCandidateResult != 0) {
                      force_mapq = 0;
                    }
                    GenerateBestMappingsForPairedEndRead(
                        pair_index, positive_candidates1.size(),
                        negative_candidates1.size(), repetitive_seed_length1,
                        min_num_errors1, num_best_mappings1,
                        second_min_num_errors1, num_second_best_mappings1,
                        read_batch1, positive_mappings1, positive_split_sites1,
                        negative_mappings1, negative_split_sites1,
                        positive_candidates2.size(),
                        negative_candidates2.size(), repetitive_seed_length2,
                        min_num_errors2, num_best_mappings2,
                        second_min_num_errors2, num_second_best_mappings2,
                        read_batch2, reference, barcode_batch,
                        positive_mappings2, positive_split_sites2,
                        negative_mappings2, negative_split_sites2,
                        &best_mapping_indices, &generator, &F1R2_best_mappings,
                        &F2R1_best_mappings, &F1F2_best_mappings,
                        &R1R2_best_mappings, &min_sum_errors,
                        &num_best_mappings, &second_min_sum_errors,
                        &num_second_best_mappings, force_mapq,
                        &mappings_on_diff_ref_seqs);
                    if (num_best_mappings == 1) {
                      ++thread_num_uniquely_mapped_reads;
                      ++thread_num_uniquely_mapped_reads;
                    }
                    thread_num_mappings +=
                        std::min(num_best_mappings, max_num_best_mappings_);
                    thread_num_mappings +=
                        std::min(num_best_mappings, max_num_best_mappings_);
                    if (num_best_mappings > 0) {
                      ++thread_num_mapped_reads;
                      ++thread_num_mapped_reads;
                    }
                  }
                }  // verify candidate
              }
            }
          }  // for pair_index
          // if (num_reads_ / 2 > initial_num_sample_barcodes_) {
          //  if (!is_bulk_data_) {
          //    if (!barcode_whitelist_file_path_.empty()) {
          //      UpdateBarcodeAbundance(num_loaded_pairs, barcode_batch);
          //    }
          //  }
          //}
          // Update cache
#pragma omp taskwait
          for (uint32_t pair_index = 0; pair_index < num_loaded_pairs;
               ++pair_index) {
            if (num_reads_ > 5000000 &&
                pair_index >= num_loaded_pairs / num_threads_) {
              break;
            }
            if (mm_history1[pair_index].timestamp != num_reads_) continue;

            mm_to_candidates_cache.Update(
                mm_history1[pair_index].minimizers,
                mm_history1[pair_index].positive_candidates,
                mm_history1[pair_index].negative_candidates,
                mm_history1[pair_index].repetitive_seed_length);
            mm_to_candidates_cache.Update(
                mm_history2[pair_index].minimizers,
                mm_history2[pair_index].positive_candidates,
                mm_history2[pair_index].negative_candidates,
                mm_history2[pair_index].repetitive_seed_length);

            if (mm_history1[pair_index].positive_candidates.size() > 50) {
              std::vector<Candidate>().swap(
                  mm_history1[pair_index].positive_candidates);
            }
            if (mm_history1[pair_index].negative_candidates.size() > 50) {
              std::vector<Candidate>().swap(
                  mm_history1[pair_index].negative_candidates);
            }
            if (mm_history2[pair_index].positive_candidates.size() > 50) {
              std::vector<Candidate>().swap(
                  mm_history2[pair_index].positive_candidates);
            }
            if (mm_history2[pair_index].negative_candidates.size() > 50) {
              std::vector<Candidate>().swap(
                  mm_history2[pair_index].negative_candidates);
            }
          }
          std::cerr << "Mapped " << num_loaded_pairs << " read pairs in "
                    << Chromap<>::GetRealTime() - real_batch_start_time
                    << "s.\n";
          real_batch_start_time = Chromap<>::GetRealTime();
          // Swap to next batch
          num_loaded_pairs = num_loaded_pairs_for_loading;
          read_batch1_for_loading.SwapSequenceBatch(read_batch1);
          read_batch2_for_loading.SwapSequenceBatch(read_batch2);
          barcode_batch_for_loading.SwapSequenceBatch(barcode_batch);
          mappings_on_diff_ref_seqs_for_diff_threads.swap(
              mappings_on_diff_ref_seqs_for_diff_threads_for_saving);
#pragma omp task
          {
            // Handle output
            num_mappings_in_mem += MoveMappingsInBuffersToMappingContainer(
                num_reference_sequences,
                &mappings_on_diff_ref_seqs_for_diff_threads_for_saving);
            if (low_memory_mode_ &&
                num_mappings_in_mem > max_num_mappings_in_mem) {
              TempMappingFileHandle<MappingRecord> temp_mapping_file_handle;
              temp_mapping_file_handle.file_path =
                  mapping_output_file_path_ + ".temp" +
                  std::to_string(temp_mapping_file_handles_.size());
              temp_mapping_file_handles_.emplace_back(temp_mapping_file_handle);
              SortOutputMappings(num_reference_sequences,
                                 &mappings_on_diff_ref_seqs_);
              output_tools_.OutputTempMapping(
                  temp_mapping_file_handle.file_path, num_reference_sequences,
                  mappings_on_diff_ref_seqs_);
              num_mappings_in_mem = 0;
              for (uint32_t i = 0; i < num_reference_sequences; ++i) {
                mappings_on_diff_ref_seqs_[i].clear();
              }
            }
          }
        }
      }  // end of openmp single
      num_barcode_in_whitelist_ += thread_num_barcode_in_whitelist;
      num_corrected_barcode_ += thread_num_corrected_barcode;
      num_candidates_ += thread_num_candidates;
      num_mappings_ += thread_num_mappings;
      num_mapped_reads_ += thread_num_mapped_reads;
      num_uniquely_mapped_reads_ += thread_num_uniquely_mapped_reads;
    }  // end of openmp parallel region
    read_batch1_for_loading.FinalizeLoading();
    read_batch2_for_loading.FinalizeLoading();
    if (!is_bulk_data_) {
      barcode_batch_for_loading.FinalizeLoading();
    }
  }
  std::cerr << "Mapped all reads in "
            << Chromap<>::GetRealTime() - real_start_mapping_time << "s.\n";
  delete[] mm_history1;
  delete[] mm_history2;
  OutputMappingStatistics();
  if (!is_bulk_data_) {
    OutputBarcodeStatistics();
  }
  index.Destroy();
  if (low_memory_mode_) {
    PostProcessingInLowMemory(num_mappings_in_mem, num_reference_sequences,
                              reference);
  } else {
    // OutputMappingStatistics(num_reference_sequences,
    // mappings_on_diff_ref_seqs_, mappings_on_diff_ref_seqs_);
    if (Tn5_shift_) {
      ApplyTn5ShiftOnPairedEndMapping(num_reference_sequences,
                                      &mappings_on_diff_ref_seqs_);
    }
    if (remove_pcr_duplicates_) {
      RemovePCRDuplicate(num_reference_sequences);
      std::cerr << "After removing PCR duplications, ";
      OutputMappingStatistics(num_reference_sequences,
                              deduped_mappings_on_diff_ref_seqs_,
                              deduped_mappings_on_diff_ref_seqs_);
    } else {
      SortOutputMappings(num_reference_sequences, &mappings_on_diff_ref_seqs_);
    }
    if (allocate_multi_mappings_) {
      AllocateMultiMappings(num_reference_sequences);
      std::cerr << "After allocating multi-mappings, ";
      OutputMappingStatistics(num_reference_sequences,
                              allocated_mappings_on_diff_ref_seqs_,
                              allocated_mappings_on_diff_ref_seqs_);
      SortOutputMappings(num_reference_sequences,
                         &allocated_mappings_on_diff_ref_seqs_);
      OutputMappings(num_reference_sequences, reference,
                     allocated_mappings_on_diff_ref_seqs_);
    } else {
      std::vector<std::vector<MappingRecord>> &mappings =
          remove_pcr_duplicates_ ? deduped_mappings_on_diff_ref_seqs_
                                 : mappings_on_diff_ref_seqs_;
      OutputMappings(num_reference_sequences, reference, mappings);
    }
    if (!is_bulk_data_ && !matrix_output_prefix_.empty()) {
      output_tools_.InitializeMatrixOutput(matrix_output_prefix_);
      OutputFeatureMatrix(num_reference_sequences, reference);
      output_tools_.FinalizeMatrixOutput();
    }
  }
  output_tools_.FinalizeMappingOutput();
  reference.FinalizeLoading();
  std::cerr << "Total time: " << Chromap<>::GetRealTime() - real_start_time
            << "s.\n";
}

template <typename MappingRecord>
void Chromap<MappingRecord>::OutputMappingsInVector(
    uint8_t mapq_threshold, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const std::vector<std::vector<MappingRecord>> &mappings) {
  uint64_t num_mappings_passing_filters = 0;
  for (uint32_t ri = 0; ri < num_reference_sequences; ++ri) {
    for (auto it = mappings[ri].begin(); it != mappings[ri].end(); ++it) {
      uint8_t mapq = (it->mapq);
      // uint8_t is_unique = (it->is_unique);
      if (mapq >= mapq_threshold) {
        // if (allocate_multi_mappings_ || (only_output_unique_mappings_ &&
        // is_unique == 1)) {
        output_tools_.AppendMapping(ri, reference, *it);
        ++num_mappings_passing_filters;
        //}
      }
    }
  }
  std::cerr << "Number of output mappings (passed filters): "
            << num_mappings_passing_filters << "\n";
}

template <typename MappingRecord>
void Chromap<MappingRecord>::OutputMappings(
    uint32_t num_reference_sequences, const SequenceBatch &reference,
    const std::vector<std::vector<MappingRecord>> &mappings) {
  // if (only_output_unique_mappings_ && mapq_threshold_ < 4)
  //  mapq_threshold_ = 4;
  OutputMappingsInVector(mapq_threshold_, num_reference_sequences, reference,
                         mappings);
}

template <typename MappingRecord>
void Chromap<MappingRecord>::ReduceCandidatesForPairedEndReadOnOneDirection(
    const std::vector<Candidate> &candidates1,
    const std::vector<Candidate> &candidates2,
    std::vector<Candidate> *filtered_candidates1,
    std::vector<Candidate> *filtered_candidates2) {
  uint32_t i1 = 0;
  uint32_t i2 = 0;
  uint32_t mapping_positions_distance = max_insert_size_;
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
        filtered_candidates2->emplace_back(candidates2[i2]);
        ++num_unpaired_candidate2;
      }
      ++i2;
    } else if (candidates2[i2].position >
               candidates1[i1].position + mapping_positions_distance) {
      if (num_unpaired_candidate1 < num_unpaired_candidate_threshold &&
          (candidates1[i1].position >> 32) ==
              (candidates2[i2].position >> 32) &&
          candidates1[i1].count >= max_candidate_count1) {
        filtered_candidates1->emplace_back(candidates1[i1]);
        ++num_unpaired_candidate1;
      }
      ++i1;
    } else {
      // ok, find a pair, we store current ni2 somewhere and keep looking until
      // we go out of the range, then we go back and then move to next pi1 and
      // keep doing the similar thing.
      filtered_candidates1->emplace_back(candidates1[i1]);
      if (candidates1[i1].count > max_candidate_count1) {
        max_candidate_count1 = candidates1[i1].count;
      }
      uint32_t current_i2 = i2;
      while (current_i2 < candidates2.size() &&
             candidates2[current_i2].position <=
                 candidates1[i1].position + mapping_positions_distance) {
        if (current_i2 >= previous_end_i2) {
          filtered_candidates2->emplace_back(candidates2[current_i2]);
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

template <typename MappingRecord>
void Chromap<MappingRecord>::ReduceCandidatesForPairedEndRead(
    const std::vector<Candidate> &positive_candidates1,
    const std::vector<Candidate> &negative_candidates1,
    const std::vector<Candidate> &positive_candidates2,
    const std::vector<Candidate> &negative_candidates2,
    std::vector<Candidate> *filtered_positive_candidates1,
    std::vector<Candidate> *filtered_negative_candidates1,
    std::vector<Candidate> *filtered_positive_candidates2,
    std::vector<Candidate> *filtered_negative_candidates2) {
  ReduceCandidatesForPairedEndReadOnOneDirection(
      positive_candidates1, negative_candidates2, filtered_positive_candidates1,
      filtered_negative_candidates2);
  ReduceCandidatesForPairedEndReadOnOneDirection(
      negative_candidates1, positive_candidates2, filtered_negative_candidates1,
      filtered_positive_candidates2);
}

template <typename MappingRecord>
void Chromap<MappingRecord>::
    RecalibrateBestMappingsForPairedEndReadOnOneDirection(
        Direction first_read_direction, uint32_t pair_index, int min_sum_errors,
        int second_min_sum_errors, int min_num_errors1, int num_best_mappings1,
        int second_min_num_errors1, int num_second_best_mappings1,
        const SequenceBatch &read_batch1,
        const std::vector<std::pair<int, uint64_t>> &mappings1,
        int min_num_errors2, int num_best_mappings2, int second_min_num_errors2,
        int num_second_best_mappings2, const SequenceBatch &read_batch2,
        const SequenceBatch &reference,
        const std::vector<std::pair<int, uint64_t>> &mappings2,
        const std::vector<std::pair<uint32_t, uint32_t>> &edit_best_mappings,
        std::vector<std::pair<uint32_t, uint32_t>> *best_mappings,
        int *best_alignment_score, int *num_best_mappings,
        int *second_best_alignment_score, int *num_second_best_mappings) {
  int8_t mat[25];
  int i, j, k;
  for (i = k = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j)
      mat[k++] = i == j ? match_score_ : -mismatch_penalty_;
    mat[k++] = 0;  // ambiguous base
  }
  for (j = 0; j < 5; ++j) mat[k++] = 0;
  best_mappings->reserve(*num_best_mappings);
  const char *read1 = read_batch1.GetSequenceAt(pair_index);
  const char *read2 = read_batch2.GetSequenceAt(pair_index);
  uint32_t read1_length = read_batch1.GetSequenceLengthAt(pair_index);
  uint32_t read2_length = read_batch2.GetSequenceLengthAt(pair_index);
  const std::string &negative_read1 =
      read_batch1.GetNegativeSequenceAt(pair_index);
  const std::string &negative_read2 =
      read_batch2.GetNegativeSequenceAt(pair_index);
  // uint32_t read_id = read_batch1.GetSequenceIdAt(pair_index);
  for (uint32_t mi = 0; mi < edit_best_mappings.size(); ++mi) {
    uint32_t i1 = edit_best_mappings[mi].first;
    uint32_t i2 = edit_best_mappings[mi].second;
    int current_sum_errors = mappings1[i1].first + mappings2[i2].first;
    if (current_sum_errors == min_sum_errors) {
      uint32_t rid1 = mappings1[i1].second >> 32;
      uint32_t position1 = mappings1[i1].second;
      uint32_t verification_window_start_position1 =
          position1 + 1 > (read1_length + error_threshold_)
              ? position1 + 1 - read1_length - error_threshold_
              : 0;
      if (position1 >= reference.GetSequenceLengthAt(rid1)) {
        verification_window_start_position1 =
            reference.GetSequenceLengthAt(rid1) - error_threshold_ -
            read1_length;
      }
      // int mapping_start_position1;
      uint32_t rid2 = mappings2[i2].second >> 32;
      uint32_t position2 = mappings2[i2].second;
      uint32_t verification_window_start_position2 =
          position2 + 1 > (read2_length + error_threshold_)
              ? position2 + 1 - read2_length - error_threshold_
              : 0;
      if (position2 >= reference.GetSequenceLengthAt(rid2)) {
        verification_window_start_position2 =
            reference.GetSequenceLengthAt(rid2) - error_threshold_ -
            read2_length;
      }
      int current_alignment_score1, current_alignment_score2,
          current_alignment_score;
      if (first_read_direction == kPositive) {
        current_alignment_score1 = ksw_semi_global2(
            read1_length + 2 * error_threshold_,
            reference.GetSequenceAt(rid1) + verification_window_start_position1,
            read1_length, read1, 5, mat, gap_open_penalties_[0],
            gap_extension_penalties_[0], gap_open_penalties_[1],
            gap_extension_penalties_[1], error_threshold_ * 2 + 1, NULL, NULL);
        current_alignment_score2 = ksw_semi_global2(
            read2_length + 2 * error_threshold_,
            reference.GetSequenceAt(rid2) + verification_window_start_position2,
            read2_length, negative_read2.data(), 5, mat, gap_open_penalties_[0],
            gap_extension_penalties_[0], gap_open_penalties_[1],
            gap_extension_penalties_[1], error_threshold_ * 2 + 1, NULL, NULL);
      } else {
        current_alignment_score1 = ksw_semi_global2(
            read1_length + 2 * error_threshold_,
            reference.GetSequenceAt(rid1) + verification_window_start_position1,
            read1_length, negative_read1.data(), 5, mat, gap_open_penalties_[0],
            gap_extension_penalties_[0], gap_open_penalties_[1],
            gap_extension_penalties_[1], error_threshold_ * 2 + 1, NULL, NULL);
        current_alignment_score2 = ksw_semi_global2(
            read1_length + 2 * error_threshold_,
            reference.GetSequenceAt(rid2) + verification_window_start_position2,
            read2_length, read2, 5, mat, gap_open_penalties_[0],
            gap_extension_penalties_[0], gap_open_penalties_[1],
            gap_extension_penalties_[1], error_threshold_ * 2 + 1, NULL, NULL);
      }
      current_alignment_score =
          current_alignment_score1 + current_alignment_score2;
      // std::cerr << current_alignment_score1 << " " <<
      // current_alignment_score2 << " " << current_alignment_score <<"\n";
      if (current_alignment_score > *best_alignment_score) {
        *second_best_alignment_score = *best_alignment_score;
        *num_second_best_mappings = *num_best_mappings;
        *best_alignment_score = current_alignment_score;
        *num_best_mappings = 1;
        best_mappings->clear();
        best_mappings->emplace_back(i1, i2);
      } else if (current_alignment_score == *best_alignment_score) {
        (*num_best_mappings)++;
        best_mappings->emplace_back(i1, i2);
      } else if (current_alignment_score == *second_best_alignment_score) {
        (*num_second_best_mappings)++;
      }
    } else if (current_sum_errors == second_min_sum_errors) {
    }
  }
}

template <typename MappingRecord>
void Chromap<MappingRecord>::GenerateBestMappingsForPairedEndReadOnOneDirection(
    Direction first_read_direction, uint32_t pair_index, int num_candidates1,
    int min_num_errors1, int num_best_mappings1, int second_min_num_errors1,
    int num_second_best_mappings1, const SequenceBatch &read_batch1,
    const std::vector<std::pair<int, uint64_t>> &mappings1, int num_candidates2,
    int min_num_errors2, int num_best_mappings2, int second_min_num_errors2,
    int num_second_best_mappings2, const SequenceBatch &read_batch2,
    const SequenceBatch &reference,
    const std::vector<std::pair<int, uint64_t>> &mappings2,
    std::vector<std::pair<uint32_t, uint32_t>> *best_mappings,
    int *min_sum_errors, int *num_best_mappings, int *second_min_sum_errors,
    int *num_second_best_mappings) {
  uint32_t i1 = 0;
  uint32_t i2 = 0;
  uint32_t min_overlap_length = min_read_length_;
  uint32_t read1_length = read_batch1.GetSequenceLengthAt(pair_index);
  uint32_t read2_length = read_batch2.GetSequenceLengthAt(pair_index);
#ifdef LI_DEBUG
  for (int i = 0; i < mappings1.size(); ++i)
    printf("mappings1 %d %d:%d\n", i, (int)(mappings1[i].second >> 32),
           (int)mappings1[i].second);
  for (int i = 0; i < mappings1.size(); ++i)
    printf("mappings2 %d %d:%d\n", i, (int)(mappings2[i].second >> 32),
           (int)mappings2[i].second);
#endif

  if (split_alignment_) {
    if (mappings1.size() == 0 || mappings2.size() == 0) {
      return;
    }
    // For split alignment, selecting the pairs whose both single-end are the
    // best.
    for (i1 = 0; i1 < mappings1.size(); ++i1) {
      if (mappings1[i1].first != min_num_errors1) {
        continue;
      }
      for (i2 = 0; i2 < mappings2.size(); ++i2) {
        if (mappings2[i2].first != min_num_errors2) {
          continue;
        }
        best_mappings->emplace_back(i1, i2);
        *min_sum_errors = min_num_errors1 + min_num_errors2;
        //*second_min_sum_errors = min_num_errors1 + min_num_errors2 + 1;
        (*num_best_mappings)++;
      }
    }

    return;
  }

  while (i1 < mappings1.size() && i2 < mappings2.size()) {
    if ((first_read_direction == kNegative &&
         mappings1[i1].second >
             mappings2[i2].second + max_insert_size_ - read1_length) ||
        (first_read_direction == kPositive &&
         mappings1[i1].second >
             mappings2[i2].second + read2_length - min_overlap_length)) {
      ++i2;
    } else if ((first_read_direction == kPositive &&
                mappings2[i2].second >
                    mappings1[i1].second + max_insert_size_ - read2_length) ||
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
                   mappings1[i1].second + max_insert_size_ - read2_length) ||
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
        if (current_sum_errors < *min_sum_errors) {
          *second_min_sum_errors = *min_sum_errors;
          *num_second_best_mappings = *num_best_mappings;
          *min_sum_errors = current_sum_errors;
          *num_best_mappings = 1;
          best_mappings->clear();
          best_mappings->emplace_back(i1, current_i2);
        } else if (current_sum_errors == *min_sum_errors) {
          (*num_best_mappings)++;
          best_mappings->emplace_back(i1, current_i2);
        } else if (current_sum_errors == *second_min_sum_errors) {
          (*num_second_best_mappings)++;
        } else if (current_sum_errors < *second_min_sum_errors) {
          *second_min_sum_errors = current_sum_errors;
          *num_second_best_mappings = 1;
        }
        ++current_i2;
      }
      ++i1;
    }
  }
}

template <typename MappingRecord>
void Chromap<MappingRecord>::EmplaceBackMappingRecord(
    uint32_t read_id, uint64_t barcode, uint32_t fragment_start_position,
    uint16_t fragment_length, uint8_t mapq, uint8_t direction,
    uint8_t is_unique, uint8_t num_dups, uint16_t positive_alignment_length,
    uint16_t negative_alignment_length,
    std::vector<MappingRecord> *mappings_on_diff_ref_seqs) {}

template <>
void Chromap<PairedEndMappingWithoutBarcode>::EmplaceBackMappingRecord(
    uint32_t read_id, uint64_t barcode, uint32_t fragment_start_position,
    uint16_t fragment_length, uint8_t mapq, uint8_t direction,
    uint8_t is_unique, uint8_t num_dups, uint16_t positive_alignment_length,
    uint16_t negative_alignment_length,
    std::vector<PairedEndMappingWithoutBarcode> *mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs->emplace_back(PairedEndMappingWithoutBarcode{
      read_id, fragment_start_position, fragment_length, mapq, direction,
      is_unique, num_dups, positive_alignment_length,
      negative_alignment_length});
}

template <>
void Chromap<PairedEndMappingWithBarcode>::EmplaceBackMappingRecord(
    uint32_t read_id, uint64_t barcode, uint32_t fragment_start_position,
    uint16_t fragment_length, uint8_t mapq, uint8_t direction,
    uint8_t is_unique, uint8_t num_dups, uint16_t positive_alignment_length,
    uint16_t negative_alignment_length,
    std::vector<PairedEndMappingWithBarcode> *mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs->emplace_back(PairedEndMappingWithBarcode{
      read_id, barcode, fragment_start_position, fragment_length, mapq,
      direction, is_unique, num_dups, positive_alignment_length,
      negative_alignment_length});
}

template <typename MappingRecord>
void Chromap<MappingRecord>::EmplaceBackMappingRecord(
    uint32_t read_id, const char *read1_name, const char *read2_name,
    uint16_t read1_length, uint16_t read2_length, uint64_t barcode,
    uint32_t fragment_start_position, uint16_t fragment_length, uint8_t mapq1,
    uint8_t mapq2, uint8_t direction, uint8_t is_unique, uint8_t num_dups,
    uint16_t positive_alignment_length, uint16_t negative_alignment_length,
    std::vector<MappingRecord> *mappings_on_diff_ref_seqs) {}

template <>
void Chromap<PairedPAFMapping>::EmplaceBackMappingRecord(
    uint32_t read_id, const char *read1_name, const char *read2_name,
    uint16_t read1_length, uint16_t read2_length, uint64_t barcode,
    uint32_t fragment_start_position, uint16_t fragment_length, uint8_t mapq1,
    uint8_t mapq2, uint8_t direction, uint8_t is_unique, uint8_t num_dups,
    uint16_t positive_alignment_length, uint16_t negative_alignment_length,
    std::vector<PairedPAFMapping> *mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs->emplace_back(PairedPAFMapping{
      read_id, std::string(read1_name), std::string(read2_name), read1_length,
      read2_length, fragment_start_position, fragment_length,
      positive_alignment_length, negative_alignment_length,
      mapq1 < mapq2 ? mapq1 : mapq2, mapq1, mapq2, direction, is_unique,
      num_dups});
}

template <typename MappingRecord>
void Chromap<MappingRecord>::EmplaceBackMappingRecord(
    uint32_t read_id, const char *read_name, uint64_t cell_barcode, int rid1,
    int rid2, uint32_t pos1, uint32_t pos2, int direction1, int direction2,
    uint8_t mapq, uint8_t is_unique, uint8_t num_dups,
    std::vector<MappingRecord> *mappings_on_diff_ref_seqs) {}

template <>
void Chromap<PairsMapping>::EmplaceBackMappingRecord(
    uint32_t read_id, const char *read_name, uint64_t cell_barcode, int rid1,
    int rid2, uint32_t pos1, uint32_t pos2, int direction1, int direction2,
    uint8_t mapq, uint8_t is_unique, uint8_t num_dups,
    std::vector<PairsMapping> *mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs->emplace_back(PairsMapping{
      read_id, std::string(read_name), cell_barcode, rid1, rid2, pos1, pos2,
      direction1, direction2, mapq, is_unique, num_dups});
}

template <typename MappingRecord>
void Chromap<MappingRecord>::ProcessBestMappingsForPairedEndReadOnOneDirection(
    Direction first_read_direction, Direction second_read_direction,
    uint32_t pair_index, uint8_t mapq, int num_candidates1,
    uint32_t repetitive_seed_length1, int min_num_errors1,
    int num_best_mappings1, int second_min_num_errors1,
    int num_second_best_mappings1, const SequenceBatch &read_batch1,
    const std::vector<std::pair<int, uint64_t>> &mappings1,
    const std::vector<int> &split_sites1, int num_candidates2,
    uint32_t repetitive_seed_length2, int min_num_errors2,
    int num_best_mappings2, int second_min_num_errors2,
    int num_second_best_mappings2, const SequenceBatch &read_batch2,
    const SequenceBatch &reference, const SequenceBatch &barcode_batch,
    const std::vector<int> &best_mapping_indices,
    const std::vector<std::pair<int, uint64_t>> &mappings2,
    const std::vector<int> &split_sites2,
    const std::vector<std::pair<uint32_t, uint32_t>> &best_mappings,
    int min_sum_errors, int num_best_mappings, int second_min_sum_errors,
    int num_second_best_mappings, int *best_mapping_index,
    int *num_best_mappings_reported, int force_mapq,
    std::vector<std::vector<MappingRecord>> *mappings_on_diff_ref_seqs) {
  const char *read1 = read_batch1.GetSequenceAt(pair_index);
  const char *read2 = read_batch2.GetSequenceAt(pair_index);
  uint32_t read1_length = read_batch1.GetSequenceLengthAt(pair_index);
  uint32_t read2_length = read_batch2.GetSequenceLengthAt(pair_index);
  const char *read1_name = read_batch1.GetSequenceNameAt(pair_index);
  const char *read2_name = read_batch2.GetSequenceNameAt(pair_index);
  const std::string &negative_read1 =
      read_batch1.GetNegativeSequenceAt(pair_index);
  const std::string &negative_read2 =
      read_batch2.GetNegativeSequenceAt(pair_index);
  uint32_t read_id = read_batch1.GetSequenceIdAt(pair_index);
  uint8_t is_unique = (num_best_mappings == 1 || num_best_mappings1 == 1 ||
                       num_best_mappings2 == 1)
                          ? 1
                          : 0;
  uint64_t barcode_key = 0;
  if (!is_bulk_data_) {
    barcode_key = barcode_batch.GenerateSeedFromSequenceAt(
        pair_index, 0, barcode_batch.GetSequenceLengthAt(pair_index));
  }

  // uint8_t is_unique = num_best_mappings == 1 ? 1 : 0;
  for (uint32_t mi = 0; mi < best_mappings.size(); ++mi) {
    uint32_t i1 = best_mappings[mi].first;
    uint32_t i2 = best_mappings[mi].second;
    int current_sum_errors = mappings1[i1].first + mappings2[i2].first;
    if (current_sum_errors == min_sum_errors) {
      if (*best_mapping_index ==
          best_mapping_indices[*num_best_mappings_reported]) {
        uint32_t rid1 = mappings1[i1].second >> 32;
        uint32_t rid2 = mappings2[i2].second >> 32;
        uint32_t ref_start_position1, ref_end_position1;
        uint32_t ref_start_position2, ref_end_position2;
        int split_site1 = 0;
        int split_site2 = 0;
        const char *effect_read1 = read1;
        const char *effect_read2 = read2;
        uint32_t *cigar1, *cigar2;
        int n_cigar1 = 0, n_cigar2 = 0;
        int NM1 = 0, NM2 = 0;
        std::string MD_tag1 = "", MD_tag2 = "";
        uint8_t mapq1 = 0;
        uint8_t mapq2 = 0;
        if (first_read_direction == kNegative) {
          effect_read1 = negative_read1.data();
        }
        if (second_read_direction == kNegative) {
          effect_read2 = negative_read2.data();
        }
        if (split_alignment_) {
          split_site1 = split_sites1[i1];
          split_site2 = split_sites2[i2];
        }
        GetRefStartEndPositionForReadFromMapping(
            first_read_direction, mappings1[i1], effect_read1, read1_length,
            split_site1, reference, &ref_start_position1, &ref_end_position1,
            &n_cigar1, &cigar1, &NM1, MD_tag1);
        GetRefStartEndPositionForReadFromMapping(
            second_read_direction, mappings2[i2], effect_read2, read2_length,
            split_site2, reference, &ref_start_position2, &ref_end_position2,
            &n_cigar2, &cigar2, &NM2, MD_tag2);
        /*if (!strcmp(read_batch1.GetSequenceNameAt(pair_index),
        "SIM3C:17:3C:1:1:1:100013")) { printf("%d\n", best_mappings.size());
          printf("%d. %d %d\n", min_sum_errors, i1, i2);
          printf("%d %d %d. %d %d\n", first_read_direction, mappings1[i1].first,
        (int)mappings1[i1].second, (split_site1 >> 16)&0xff, split_site1 &
        0xffff); printf("%d %d %d. %d %d\n", second_read_direction,
        mappings1[i2].first, (int)mappings2[i2].second, (split_site2 >>
        16)&0xff, split_site2 & 0xffff); printf("%d %d\n", ref_start_position1,
        ref_end_position1); printf("%d %d\n", ref_start_position2,
        ref_end_position2);
        }*/
        mapq = GetMAPQForPairedEndRead(
            num_candidates1, num_candidates2, repetitive_seed_length1,
            repetitive_seed_length2,
            ref_end_position1 - ref_start_position1 + 1,
            ref_end_position2 - ref_start_position2 + 1, min_sum_errors,
            num_best_mappings, second_min_sum_errors, num_second_best_mappings,
            mappings1[i1].first, mappings2[i2].first, min_num_errors1,
            min_num_errors2, num_best_mappings1, num_best_mappings2,
            second_min_num_errors1, second_min_num_errors2,
            num_second_best_mappings1, num_second_best_mappings2, read1_length,
            read2_length, force_mapq, mapq1, mapq2);
        uint8_t direction = 1;
        if (first_read_direction == kNegative) {
          direction = 0;
        }
        // mapq |= (uint8_t)1;
        if (mapping_output_format_ == MAPPINGFORMAT_SAM) {
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
          if (*num_best_mappings_reported >= 1) {
            flag1 |= BAM_FSECONDARY;
            flag2 |= BAM_FSECONDARY;
          }
          // printf("%d %d\n", ref_start_position1, ref_end_position1);
          // printf("%d %d\n", ref_start_position2, ref_end_position2);
          EmplaceBackMappingRecord(
              read_id, read1_name, barcode_key, 1, ref_start_position1, rid1,
              flag1, first_read_direction == kPositive ? 1 : 0, is_unique, mapq,
              NM1, n_cigar1, cigar1, MD_tag1, read1,
              read_batch1.GetSequenceQualAt(pair_index),
              &((*mappings_on_diff_ref_seqs)[rid1]));
          EmplaceBackMappingRecord(
              read_id, read2_name, barcode_key, 1, ref_start_position2, rid2,
              flag2, second_read_direction == kPositive ? 1 : 0, is_unique,
              mapq, NM2, n_cigar2, cigar2, MD_tag2, read2,
              read_batch2.GetSequenceQualAt(pair_index),
              &((*mappings_on_diff_ref_seqs)[rid2]));
        } else if (mapping_output_format_ == MAPPINGFORMAT_PAIRS) {
          int position1 = ref_start_position1;
          int position2 = ref_start_position2;

          uint8_t direction2 = 1;
          if (second_read_direction == kNegative) {
            direction2 = 0;
            position2 = ref_end_position2;
          }
          if (direction == 0) {
            position1 = ref_end_position1;
          }
          int rid1_rank = pairs_custom_rid_rank_[rid1];
          int rid2_rank = pairs_custom_rid_rank_[rid2];
          if (rid1_rank < rid2_rank ||
              (rid1 == rid2 && position1 < position2)) {
            EmplaceBackMappingRecord(read_id, read1_name, barcode_key, rid1,
                                     rid2, position1, position2, direction,
                                     direction2, mapq, is_unique, 1,
                                     &((*mappings_on_diff_ref_seqs)[rid1]));
          } else {
            EmplaceBackMappingRecord(read_id, read1_name, barcode_key, rid2,
                                     rid1, position2, position1, direction2,
                                     direction, mapq, is_unique, 1,
                                     &((*mappings_on_diff_ref_seqs)[rid2]));
          }
        } else if (mapping_output_format_ == MAPPINGFORMAT_PAF) {
          uint32_t fragment_start_position = ref_start_position1;
          uint16_t fragment_length =
              ref_end_position2 - ref_start_position1 + 1;
          ;
          uint16_t positive_alignment_length =
              ref_end_position1 - ref_start_position1 + 1;
          uint16_t negative_alignment_length =
              ref_end_position2 - ref_start_position2 + 1;
          if (direction == 0) {
            fragment_start_position = ref_start_position2;
            fragment_length = ref_end_position1 - ref_start_position2 + 1;
            positive_alignment_length =
                ref_end_position2 - ref_start_position2 + 1;
            negative_alignment_length =
                ref_end_position1 - ref_start_position1 + 1;
          }
          EmplaceBackMappingRecord(
              read_id, read1_name, read2_name,
              (uint16_t)read_batch1.GetSequenceLengthAt(pair_index),
              (uint16_t)read_batch2.GetSequenceLengthAt(pair_index),
              barcode_key, fragment_start_position, fragment_length, mapq1,
              mapq2, direction, is_unique, 1, positive_alignment_length,
              negative_alignment_length, &((*mappings_on_diff_ref_seqs)[rid1]));
        } else {
          uint32_t fragment_start_position = ref_start_position1;
          uint16_t fragment_length =
              ref_end_position2 - ref_start_position1 + 1;
          ;
          uint16_t positive_alignment_length =
              ref_end_position1 - ref_start_position1 + 1;
          uint16_t negative_alignment_length =
              ref_end_position2 - ref_start_position2 + 1;
          if (direction == 0) {
            fragment_start_position = ref_start_position2;
            fragment_length = ref_end_position1 - ref_start_position2 + 1;
            positive_alignment_length =
                ref_end_position2 - ref_start_position2 + 1;
            negative_alignment_length =
                ref_end_position1 - ref_start_position1 + 1;
          }
          EmplaceBackMappingRecord(
              read_id, barcode_key, fragment_start_position, fragment_length,
              mapq, direction, is_unique, 1, positive_alignment_length,
              negative_alignment_length, &((*mappings_on_diff_ref_seqs)[rid1]));
        }
        (*num_best_mappings_reported)++;
        if (*num_best_mappings_reported ==
            std::min(max_num_best_mappings_, num_best_mappings)) {
          break;
        }
      }
      (*best_mapping_index)++;
    }
  }
}

template <typename MappingRecord>
void Chromap<MappingRecord>::GenerateBestMappingsForPairedEndRead(
    uint32_t pair_index, int num_positive_candidates1,
    int num_negative_candidates1, uint32_t repetitive_seed_length1,
    int min_num_errors1, int num_best_mappings1, int second_min_num_errors1,
    int num_second_best_mappings1, const SequenceBatch &read_batch1,
    const std::vector<std::pair<int, uint64_t>> &positive_mappings1,
    const std::vector<int> &positive_split_sites1,
    const std::vector<std::pair<int, uint64_t>> &negative_mappings1,
    const std::vector<int> &negative_split_sites1, int num_positive_candidates2,
    int num_negative_candidates2, uint32_t repetitive_seed_length2,
    int min_num_errors2, int num_best_mappings2, int second_min_num_errors2,
    int num_second_best_mappings2, const SequenceBatch &read_batch2,
    const SequenceBatch &reference, const SequenceBatch &barcode_batch,
    const std::vector<std::pair<int, uint64_t>> &positive_mappings2,
    const std::vector<int> &positive_split_sites2,
    const std::vector<std::pair<int, uint64_t>> &negative_mappings2,
    const std::vector<int> &negative_split_sites2,
    std::vector<int> *best_mapping_indices, std::mt19937 *generator,
    std::vector<std::pair<uint32_t, uint32_t>> *F1R2_best_mappings,
    std::vector<std::pair<uint32_t, uint32_t>> *F2R1_best_mappings,
    std::vector<std::pair<uint32_t, uint32_t>> *F1F2_best_mappings,
    std::vector<std::pair<uint32_t, uint32_t>> *R1R2_best_mappings,
    int *min_sum_errors, int *num_best_mappings, int *second_min_sum_errors,
    int *num_second_best_mappings, int force_mapq,
    std::vector<std::vector<MappingRecord>> *mappings_on_diff_ref_seqs) {
  *min_sum_errors = 2 * error_threshold_ + 1;
  *num_best_mappings = 0;
  *second_min_sum_errors = *min_sum_errors;
  *num_second_best_mappings = 0;
  GenerateBestMappingsForPairedEndReadOnOneDirection(
      kPositive, pair_index, num_positive_candidates1, min_num_errors1,
      num_best_mappings1, second_min_num_errors1, num_second_best_mappings1,
      read_batch1, positive_mappings1, num_negative_candidates2,
      min_num_errors2, num_best_mappings2, second_min_num_errors2,
      num_second_best_mappings2, read_batch2, reference, negative_mappings2,
      F1R2_best_mappings, min_sum_errors, num_best_mappings,
      second_min_sum_errors, num_second_best_mappings);
  GenerateBestMappingsForPairedEndReadOnOneDirection(
      kNegative, pair_index, num_negative_candidates1, min_num_errors1,
      num_best_mappings1, second_min_num_errors1, num_second_best_mappings1,
      read_batch1, negative_mappings1, num_positive_candidates2,
      min_num_errors2, num_best_mappings2, second_min_num_errors2,
      num_second_best_mappings2, read_batch2, reference, positive_mappings2,
      F2R1_best_mappings, min_sum_errors, num_best_mappings,
      second_min_sum_errors, num_second_best_mappings);
  if (split_alignment_) {
    GenerateBestMappingsForPairedEndReadOnOneDirection(
        kPositive, pair_index, num_positive_candidates1, min_num_errors1,
        num_best_mappings1, second_min_num_errors1, num_second_best_mappings1,
        read_batch1, positive_mappings1, num_positive_candidates2,
        min_num_errors2, num_best_mappings2, second_min_num_errors2,
        num_second_best_mappings2, read_batch2, reference, positive_mappings2,
        F1F2_best_mappings, min_sum_errors, num_best_mappings,
        second_min_sum_errors, num_second_best_mappings);
    GenerateBestMappingsForPairedEndReadOnOneDirection(
        kNegative, pair_index, num_negative_candidates1, min_num_errors1,
        num_best_mappings1, second_min_num_errors1, num_second_best_mappings1,
        read_batch1, negative_mappings1, num_negative_candidates2,
        min_num_errors2, num_best_mappings2, second_min_num_errors2,
        num_second_best_mappings2, read_batch2, reference, negative_mappings2,
        R1R2_best_mappings, min_sum_errors, num_best_mappings,
        second_min_sum_errors, num_second_best_mappings);
  }
  // int best_alignment_score, second_best_alignment_score;
  // best_alignment_score = -1;
  //*num_best_mappings = 0;
  // second_best_alignment_score = best_alignment_score;
  //*num_second_best_mappings = 0;
  // std::vector<std::pair<uint32_t, uint32_t> > edit_F1R2_best_mappings,
  // edit_F2R1_best_mappings; edit_F1R2_best_mappings.swap(*F1R2_best_mappings);
  // edit_F2R1_best_mappings.swap(*F2R1_best_mappings);
  // RecalibrateBestMappingsForPairedEndReadOnOneDirection(kPositive,
  // pair_index, *min_sum_errors, *second_min_sum_errors, min_num_errors1,
  // num_best_mappings1, second_min_num_errors1, num_second_best_mappings1,
  // read_batch1, positive_mappings1, min_num_errors2, num_best_mappings2,
  // second_min_num_errors2, num_second_best_mappings2, read_batch2, reference,
  // negative_mappings2, edit_F1R2_best_mappings, F1R2_best_mappings,
  // &best_alignment_score, num_best_mappings, &second_best_alignment_score,
  // num_second_best_mappings);
  // RecalibrateBestMappingsForPairedEndReadOnOneDirection(kNegative,
  // pair_index, *min_sum_errors, *second_min_sum_errors, min_num_errors1,
  // num_best_mappings1, second_min_num_errors1, num_second_best_mappings1,
  // read_batch1, negative_mappings1, min_num_errors2, num_best_mappings2,
  // second_min_num_errors2, num_second_best_mappings2, read_batch2, reference,
  // positive_mappings2, edit_F2R1_best_mappings, F2R1_best_mappings,
  // &best_alignment_score, num_best_mappings, &second_best_alignment_score,
  // num_second_best_mappings); uint8_t mapq = GetMAPQ(*num_best_mappings,
  // *num_second_best_mappings);
  uint8_t mapq = 0;
  if (*num_best_mappings <= drop_repetitive_reads_) {
    // we will use reservoir sampling
    // std::vector<int> best_mapping_indices(max_num_best_mappings_);
    std::iota(best_mapping_indices->begin(), best_mapping_indices->end(), 0);
    if (*num_best_mappings > max_num_best_mappings_) {
      // std::mt19937 tmp_generator(11);
      for (int i = max_num_best_mappings_; i < *num_best_mappings; ++i) {
        std::uniform_int_distribution<int> distribution(
            0, i);  // important: inclusive range
        int j = distribution(*generator);
        // int j = distribution(tmp_generator);
        if (j < max_num_best_mappings_) {
          (*best_mapping_indices)[j] = i;
        }
      }
      std::sort(best_mapping_indices->begin(), best_mapping_indices->end());
    }
    int best_mapping_index = 0;
    int num_best_mappings_reported = 0;
    ProcessBestMappingsForPairedEndReadOnOneDirection(
        kPositive, kNegative, pair_index, mapq, num_positive_candidates1,
        repetitive_seed_length1, min_num_errors1, num_best_mappings1,
        second_min_num_errors1, num_second_best_mappings1, read_batch1,
        positive_mappings1, positive_split_sites1, num_negative_candidates2,
        repetitive_seed_length2, min_num_errors2, num_best_mappings2,
        second_min_num_errors2, num_second_best_mappings2, read_batch2,
        reference, barcode_batch, *best_mapping_indices, negative_mappings2,
        negative_split_sites2, *F1R2_best_mappings, *min_sum_errors,
        *num_best_mappings, *second_min_sum_errors, *num_second_best_mappings,
        &best_mapping_index, &num_best_mappings_reported, force_mapq,
        mappings_on_diff_ref_seqs);
    if (num_best_mappings_reported !=
        std::min(max_num_best_mappings_, *num_best_mappings)) {
      ProcessBestMappingsForPairedEndReadOnOneDirection(
          kNegative, kPositive, pair_index, mapq, num_negative_candidates1,
          repetitive_seed_length1, min_num_errors1, num_best_mappings1,
          second_min_num_errors1, num_second_best_mappings1, read_batch1,
          negative_mappings1, negative_split_sites1, num_positive_candidates2,
          repetitive_seed_length2, min_num_errors2, num_best_mappings2,
          second_min_num_errors2, num_second_best_mappings2, read_batch2,
          reference, barcode_batch, *best_mapping_indices, positive_mappings2,
          positive_split_sites2, *F2R1_best_mappings, *min_sum_errors,
          *num_best_mappings, *second_min_sum_errors, *num_second_best_mappings,
          &best_mapping_index, &num_best_mappings_reported, force_mapq,
          mappings_on_diff_ref_seqs);
    }

    if (split_alignment_ &&
        num_best_mappings_reported !=
            std::min(max_num_best_mappings_, *num_best_mappings)) {
      ProcessBestMappingsForPairedEndReadOnOneDirection(
          kPositive, kPositive, pair_index, mapq, num_positive_candidates1,
          repetitive_seed_length1, min_num_errors1, num_best_mappings1,
          second_min_num_errors1, num_second_best_mappings1, read_batch1,
          positive_mappings1, positive_split_sites1, num_positive_candidates2,
          repetitive_seed_length2, min_num_errors2, num_best_mappings2,
          second_min_num_errors2, num_second_best_mappings2, read_batch2,
          reference, barcode_batch, *best_mapping_indices, positive_mappings2,
          positive_split_sites2, *F1F2_best_mappings, *min_sum_errors,
          *num_best_mappings, *second_min_sum_errors, *num_second_best_mappings,
          &best_mapping_index, &num_best_mappings_reported, force_mapq,
          mappings_on_diff_ref_seqs);
    }
    if (split_alignment_ &&
        num_best_mappings_reported !=
            std::min(max_num_best_mappings_, *num_best_mappings)) {
      ProcessBestMappingsForPairedEndReadOnOneDirection(
          kNegative, kNegative, pair_index, mapq, num_negative_candidates1,
          repetitive_seed_length1, min_num_errors1, num_best_mappings1,
          second_min_num_errors1, num_second_best_mappings1, read_batch1,
          negative_mappings1, negative_split_sites1, num_positive_candidates2,
          repetitive_seed_length2, min_num_errors2, num_best_mappings2,
          second_min_num_errors2, num_second_best_mappings2, read_batch2,
          reference, barcode_batch, *best_mapping_indices, negative_mappings2,
          negative_split_sites2, *R1R2_best_mappings, *min_sum_errors,
          *num_best_mappings, *second_min_sum_errors, *num_second_best_mappings,
          &best_mapping_index, &num_best_mappings_reported, force_mapq,
          mappings_on_diff_ref_seqs);
    }
  }
}

template <typename MappingRecord>
void Chromap<MappingRecord>::ApplyTn5ShiftOnPairedEndMapping(
    uint32_t num_reference_sequences,
    std::vector<std::vector<MappingRecord>> *mappings) {
  uint64_t num_shifted_mappings = 0;
  for (auto &mappings_on_one_ref_seq : *mappings) {
    for (auto &mapping : mappings_on_one_ref_seq) {
      // mapping.fragment_start_position += 4;
      // mapping.positive_alignment_length -= 4;
      // mapping.fragment_length -= 9;
      // mapping.negative_alignment_length -= 5;
      mapping.Tn5Shift();
      ++num_shifted_mappings;
    }
  }
  std::cerr << "# shifted mappings: " << num_shifted_mappings << ".\n";
}

template <typename MappingRecord>
void Chromap<MappingRecord>::MapSingleEndReads() {
  double real_start_time = Chromap<>::GetRealTime();
  SequenceBatch reference;
  reference.InitializeLoading(reference_file_path_);
  uint32_t num_reference_sequences = reference.LoadAllSequences();
  if (custom_rid_order_path_.length() > 0) {
    GetRidRank(custom_rid_order_path_, reference, custom_rid_rank_);
    reference.ReorderSequences(custom_rid_rank_);
  }
  Index index(min_num_seeds_required_for_mapping_, max_seed_frequencies_,
              index_file_path_);
  index.Load();
  kmer_size_ = index.GetKmerSize();
  window_size_ = index.GetWindowSize();
  // index.Statistics(num_sequences, reference);
  SequenceBatch read_batch(read_batch_size_);
  SequenceBatch read_batch_for_loading(read_batch_size_);
  SequenceBatch barcode_batch(read_batch_size_);
  SequenceBatch barcode_batch_for_loading(read_batch_size_);
  read_batch.SetSeqEffectiveRange(read1_format_[0], read1_format_[1]);
  read_batch_for_loading.SetSeqEffectiveRange(read1_format_[0],
                                              read1_format_[1]);
  barcode_batch.SetSeqEffectiveRange(barcode_format_[0], barcode_format_[1]);
  barcode_batch_for_loading.SetSeqEffectiveRange(barcode_format_[0],
                                                 barcode_format_[1]);
  mappings_on_diff_ref_seqs_.reserve(num_reference_sequences);
  deduped_mappings_on_diff_ref_seqs_.reserve(num_reference_sequences);
  for (uint32_t i = 0; i < num_reference_sequences; ++i) {
    mappings_on_diff_ref_seqs_.emplace_back(std::vector<MappingRecord>());
    deduped_mappings_on_diff_ref_seqs_.emplace_back(
        std::vector<MappingRecord>());
  }
  // Preprocess barcodes for single cell data
  if (!is_bulk_data_) {
    if (!barcode_whitelist_file_path_.empty()) {
      LoadBarcodeWhitelist();
      ComputeBarcodeAbundance(initial_num_sample_barcodes_);
    }
  }

  output_tools_.InitializeMappingOutput(
      barcode_length_, mapping_output_file_path_, mapping_output_format_);
  output_tools_.OutputHeader(num_reference_sequences, reference);

  mm_cache mm_to_candidates_cache(2000003);
  mm_to_candidates_cache.SetKmerLength(kmer_size_);
  struct _mm_history *mm_history = new struct _mm_history[read_batch_size_];
  static uint64_t thread_num_candidates = 0;
  static uint64_t thread_num_mappings = 0;
  static uint64_t thread_num_mapped_reads = 0;
  static uint64_t thread_num_uniquely_mapped_reads = 0;
#pragma omp threadprivate(thread_num_candidates, thread_num_mappings, \
                          thread_num_mapped_reads,                    \
                          thread_num_uniquely_mapped_reads)
  double real_start_mapping_time = Chromap<>::GetRealTime();
  for (size_t read_file_index = 0; read_file_index < read_file1_paths_.size();
       ++read_file_index) {
    read_batch_for_loading.InitializeLoading(
        read_file1_paths_[read_file_index]);
    if (!is_bulk_data_) {
      barcode_batch_for_loading.InitializeLoading(
          barcode_file_paths_[read_file_index]);
    }
    uint32_t num_loaded_reads_for_loading = 0;
    uint32_t num_loaded_reads = LoadSingleEndReadsWithBarcodes(
        &read_batch_for_loading, &barcode_batch_for_loading);
    read_batch_for_loading.SwapSequenceBatch(read_batch);
    barcode_batch_for_loading.SwapSequenceBatch(barcode_batch);
    std::vector<std::vector<std::vector<MappingRecord>>>
        mappings_on_diff_ref_seqs_for_diff_threads;
    std::vector<std::vector<std::vector<MappingRecord>>>
        mappings_on_diff_ref_seqs_for_diff_threads_for_saving;
    mappings_on_diff_ref_seqs_for_diff_threads.reserve(num_threads_);
    mappings_on_diff_ref_seqs_for_diff_threads_for_saving.reserve(num_threads_);
    for (int ti = 0; ti < num_threads_; ++ti) {
      mappings_on_diff_ref_seqs_for_diff_threads.emplace_back(
          std::vector<std::vector<MappingRecord>>(num_reference_sequences));
      mappings_on_diff_ref_seqs_for_diff_threads_for_saving.emplace_back(
          std::vector<std::vector<MappingRecord>>(num_reference_sequences));
      for (uint32_t i = 0; i < num_reference_sequences; ++i) {
        mappings_on_diff_ref_seqs_for_diff_threads[ti][i].reserve(
            (num_loaded_reads +
             num_loaded_reads / 1000 * max_num_best_mappings_) /
            num_threads_ / num_reference_sequences);
        mappings_on_diff_ref_seqs_for_diff_threads_for_saving[ti][i].reserve(
            (num_loaded_reads +
             num_loaded_reads / 1000 * max_num_best_mappings_) /
            num_threads_ / num_reference_sequences);
      }
    }
#pragma omp parallel default(none) shared(window_size_, custom_rid_order_path_, error_threshold_, max_num_best_mappings_, num_threads_, num_reads_, mm_history, reference, index, read_batch, barcode_batch, read_batch_for_loading, barcode_batch_for_loading, std::cerr, num_loaded_reads_for_loading, num_loaded_reads, num_reference_sequences, mappings_on_diff_ref_seqs_for_diff_threads, mappings_on_diff_ref_seqs_for_diff_threads_for_saving, mm_to_candidates_cache) num_threads(num_threads_) reduction(+:num_candidates_, num_mappings_, num_mapped_reads_, num_uniquely_mapped_reads_)
    {
      thread_num_candidates = 0;
      thread_num_mappings = 0;
      thread_num_mapped_reads = 0;
      thread_num_uniquely_mapped_reads = 0;
      std::vector<std::pair<uint64_t, uint64_t>> minimizers;
      std::vector<uint64_t> positive_hits;
      std::vector<uint64_t> negative_hits;
      positive_hits.reserve(max_seed_frequencies_[0]);
      negative_hits.reserve(max_seed_frequencies_[0]);
      std::vector<Candidate> positive_candidates;
      std::vector<Candidate> negative_candidates;
      positive_candidates.reserve(max_seed_frequencies_[0]);
      negative_candidates.reserve(max_seed_frequencies_[0]);
      std::vector<std::pair<int, uint64_t>> positive_mappings;
      std::vector<std::pair<int, uint64_t>> negative_mappings;
      positive_mappings.reserve(max_seed_frequencies_[0]);
      negative_mappings.reserve(max_seed_frequencies_[0]);
      std::vector<int> positive_split_sites;
      std::vector<int> negative_split_sites;
      positive_split_sites.reserve(max_seed_frequencies_[0]);
      negative_split_sites.reserve(max_seed_frequencies_[0]);
#pragma omp single
      {
        while (num_loaded_reads > 0) {
          double real_batch_start_time = Chromap<>::GetRealTime();
          num_reads_ += num_loaded_reads;
#pragma omp task
          {
            num_loaded_reads_for_loading = LoadSingleEndReadsWithBarcodes(
                &read_batch_for_loading, &barcode_batch_for_loading);
          }  // end of openmp loading task
             // int grain_size = 10000;
//#pragma omp taskloop grainsize(grain_size) //num_tasks(num_threads_* 50)
#pragma omp taskloop num_tasks(num_threads_ *num_threads_)
          for (uint32_t read_index = 0; read_index < num_loaded_reads;
               ++read_index) {
            read_batch.PrepareNegativeSequenceAt(read_index);
            minimizers.clear();
            minimizers.reserve(read_batch.GetSequenceLengthAt(read_index) /
                               window_size_ * 2);
            index.GenerateMinimizerSketch(read_batch, read_index, &minimizers);
            if (minimizers.size() > 0) {
              if (custom_rid_order_path_.length() > 0) {
                RerankCandidatesRid(positive_candidates);
                RerankCandidatesRid(negative_candidates);
              }
              positive_hits.clear();
              negative_hits.clear();
              positive_candidates.clear();
              negative_candidates.clear();
              uint32_t repetitive_seed_length = 0;
              if (mm_to_candidates_cache.Query(
                      minimizers, positive_candidates, negative_candidates,
                      repetitive_seed_length,
                      read_batch.GetSequenceLengthAt(read_index)) == -1) {
                index.GenerateCandidates(
                    error_threshold_, minimizers, &repetitive_seed_length,
                    &positive_hits, &negative_hits, &positive_candidates,
                    &negative_candidates);
              }
              if (read_index < num_loaded_reads &&
                  (read_index < num_loaded_reads / num_threads_ ||
                   num_reads_ <= 2500000)) {
                mm_history[read_index].timestamp = num_reads_;
                mm_history[read_index].minimizers = minimizers;
                mm_history[read_index].positive_candidates =
                    positive_candidates;
                mm_history[read_index].negative_candidates =
                    negative_candidates;
              }
              uint32_t current_num_candidates =
                  positive_candidates.size() + negative_candidates.size();
              if (current_num_candidates > 0) {
                thread_num_candidates += current_num_candidates;
                positive_mappings.clear();
                negative_mappings.clear();
                positive_split_sites.clear();
                negative_split_sites.clear();
                int min_num_errors, second_min_num_errors;
                int num_best_mappings, num_second_best_mappings;
                VerifyCandidates(read_batch, read_index, reference, minimizers,
                                 positive_candidates, negative_candidates,
                                 &positive_mappings, &positive_split_sites,
                                 &negative_mappings, &negative_split_sites,
                                 &min_num_errors, &num_best_mappings,
                                 &second_min_num_errors,
                                 &num_second_best_mappings);
                uint32_t current_num_mappings =
                    positive_mappings.size() + negative_mappings.size();
                if (current_num_mappings > 0) {
                  std::vector<std::vector<MappingRecord>>
                      &mappings_on_diff_ref_seqs =
                          mappings_on_diff_ref_seqs_for_diff_threads
                              [omp_get_thread_num()];
                  GenerateBestMappingsForSingleEndRead(
                      positive_candidates.size(), negative_candidates.size(),
                      repetitive_seed_length, min_num_errors, num_best_mappings,
                      second_min_num_errors, num_second_best_mappings,
                      read_batch, read_index, reference, barcode_batch,
                      positive_mappings, positive_split_sites,
                      negative_mappings, negative_split_sites,
                      &mappings_on_diff_ref_seqs);
                  thread_num_mappings +=
                      std::min(num_best_mappings, max_num_best_mappings_);
                  ++thread_num_mapped_reads;
                  if (num_best_mappings == 1) {
                    ++thread_num_uniquely_mapped_reads;
                  }
                }
              }
            }
          }
#pragma omp taskwait
          for (uint32_t read_index = 0; read_index < num_loaded_reads;
               ++read_index) {
            if (num_reads_ > 2500000 &&
                read_index >= num_loaded_reads / num_threads_) {
              break;
            }
            if (mm_history[read_index].timestamp != num_reads_) continue;
            mm_to_candidates_cache.Update(
                mm_history[read_index].minimizers,
                mm_history[read_index].positive_candidates,
                mm_history[read_index].negative_candidates,
                mm_history[read_index].repetitive_seed_length);
            if (mm_history[read_index].positive_candidates.size() <
                mm_history[read_index].positive_candidates.capacity() / 2) {
              std::vector<Candidate>().swap(
                  mm_history[read_index].positive_candidates);
            }
            if (mm_history[read_index].negative_candidates.size() <
                mm_history[read_index].negative_candidates.capacity() / 2) {
              std::vector<Candidate>().swap(
                  mm_history[read_index].negative_candidates);
            }
          }
          // std::cerr<<"cache memusage: " <<
          // mm_to_candidates_cache.GetMemoryBytes() <<"\n" ;
          num_loaded_reads = num_loaded_reads_for_loading;
          read_batch_for_loading.SwapSequenceBatch(read_batch);
          barcode_batch_for_loading.SwapSequenceBatch(barcode_batch);
          mappings_on_diff_ref_seqs_for_diff_threads.swap(
              mappings_on_diff_ref_seqs_for_diff_threads_for_saving);
#pragma omp task
          {
            MoveMappingsInBuffersToMappingContainer(
                num_reference_sequences,
                &mappings_on_diff_ref_seqs_for_diff_threads_for_saving);
          }
          std::cerr << "Mapped in "
                    << Chromap<>::GetRealTime() - real_batch_start_time
                    << "s.\n";
        }
      }  // end of openmp single
      {
        num_candidates_ += thread_num_candidates;
        num_mappings_ += thread_num_mappings;
        num_mapped_reads_ += thread_num_mapped_reads;
        num_uniquely_mapped_reads_ += thread_num_uniquely_mapped_reads;
      }  // end of updating shared mapping stats
    }    // end of openmp parallel region
    read_batch_for_loading.FinalizeLoading();
    if (!is_bulk_data_) {
      barcode_batch_for_loading.FinalizeLoading();
    }
  }
  delete[] mm_history;
  OutputMappingStatistics();
  std::cerr << "Mapped all reads in "
            << Chromap<>::GetRealTime() - real_start_mapping_time << "s.\n";
  OutputMappingStatistics(num_reference_sequences, mappings_on_diff_ref_seqs_,
                          mappings_on_diff_ref_seqs_);
  if (Tn5_shift_) {
    ApplyTn5ShiftOnSingleEndMapping(num_reference_sequences,
                                    &mappings_on_diff_ref_seqs_);
  }
  if (remove_pcr_duplicates_) {
    RemovePCRDuplicate(num_reference_sequences);
    std::cerr << "After removing PCR duplications, ";
    OutputMappingStatistics(num_reference_sequences,
                            deduped_mappings_on_diff_ref_seqs_,
                            deduped_mappings_on_diff_ref_seqs_);
  } else {
    SortOutputMappings(num_reference_sequences, &mappings_on_diff_ref_seqs_);
  }
  if (allocate_multi_mappings_) {
    AllocateMultiMappings(num_reference_sequences);
    std::cerr << "After allocating multi-mappings, ";
    OutputMappingStatistics(num_reference_sequences,
                            allocated_mappings_on_diff_ref_seqs_,
                            allocated_mappings_on_diff_ref_seqs_);
    SortOutputMappings(num_reference_sequences,
                       &allocated_mappings_on_diff_ref_seqs_);
    OutputMappings(num_reference_sequences, reference,
                   allocated_mappings_on_diff_ref_seqs_);
  } else {
    std::vector<std::vector<MappingRecord>> &mappings =
        remove_pcr_duplicates_ ? deduped_mappings_on_diff_ref_seqs_
                               : mappings_on_diff_ref_seqs_;
    OutputMappings(num_reference_sequences, reference, mappings);
  }
  output_tools_.FinalizeMappingOutput();
  reference.FinalizeLoading();
  std::cerr << "Total time: " << Chromap<>::GetRealTime() - real_start_time
            << "s.\n";
}

template <typename MappingRecord>
void Chromap<MappingRecord>::EmplaceBackMappingRecord(
    uint32_t read_id, uint32_t barcode, uint32_t fragment_start_position,
    uint16_t fragment_length, uint8_t mapq, uint8_t direction,
    uint8_t is_unique, uint8_t num_dups,
    std::vector<MappingRecord> *mappings_on_diff_ref_seqs) {}

template <>
void Chromap<MappingWithoutBarcode>::EmplaceBackMappingRecord(
    uint32_t read_id, uint32_t barcode, uint32_t fragment_start_position,
    uint16_t fragment_length, uint8_t mapq, uint8_t direction,
    uint8_t is_unique, uint8_t num_dups,
    std::vector<MappingWithoutBarcode> *mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs->emplace_back(
      MappingWithoutBarcode{read_id, fragment_start_position, fragment_length,
                            mapq, direction, is_unique, num_dups});
}

template <>
void Chromap<MappingWithBarcode>::EmplaceBackMappingRecord(
    uint32_t read_id, uint32_t barcode, uint32_t fragment_start_position,
    uint16_t fragment_length, uint8_t mapq, uint8_t direction,
    uint8_t is_unique, uint8_t num_dups,
    std::vector<MappingWithBarcode> *mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs->emplace_back(MappingWithBarcode{
      read_id, barcode, fragment_start_position, fragment_length, mapq,
      direction, is_unique, num_dups});
}

template <typename MappingRecord>
void Chromap<MappingRecord>::EmplaceBackMappingRecord(
    uint32_t read_id, const char *read_name, uint16_t read_length,
    uint32_t barcode, uint32_t fragment_start_position,
    uint16_t fragment_length, uint8_t mapq, uint8_t direction,
    uint8_t is_unique, uint8_t num_dups,
    std::vector<MappingRecord> *mappings_on_diff_ref_seqs) {}

template <>
void Chromap<PAFMapping>::EmplaceBackMappingRecord(
    uint32_t read_id, const char *read_name, uint16_t read_length,
    uint32_t barcode, uint32_t fragment_start_position,
    uint16_t fragment_length, uint8_t mapq, uint8_t direction,
    uint8_t is_unique, uint8_t num_dups,
    std::vector<PAFMapping> *mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs->emplace_back(PAFMapping{
      read_id, std::string(read_name), read_length, fragment_start_position,
      fragment_length, mapq, direction, is_unique, num_dups});
}

template <typename MappingRecord>
void Chromap<MappingRecord>::EmplaceBackMappingRecord(
    uint32_t read_id, const char *read_name, uint32_t cell_barcode,
    uint8_t num_dups, int64_t position, int rid, int flag, uint8_t direction,
    uint8_t is_unique, uint8_t mapq, uint32_t NM, int n_cigar, uint32_t *cigar,
    std::string &MD_tag, const char *read, const char *read_qual,
    std::vector<MappingRecord> *mappings_on_diff_ref_seqs) {}

template <>
void Chromap<SAMMapping>::EmplaceBackMappingRecord(
    uint32_t read_id, const char *read_name, uint32_t cell_barcode,
    uint8_t num_dups, int64_t position, int rid, int flag, uint8_t direction,
    uint8_t is_unique, uint8_t mapq, uint32_t NM, int n_cigar, uint32_t *cigar,
    std::string &MD_tag, const char *read, const char *read_qual,
    std::vector<SAMMapping> *mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs->emplace_back(SAMMapping{
      read_id, std::string(read_name), cell_barcode, num_dups, position, rid,
      flag, direction, 0, is_unique, mapq, NM, n_cigar, cigar, MD_tag,
      std::string(read), std::string(read_qual)});
}

template <typename MappingRecord>
void Chromap<MappingRecord>::GenerateMDTag(const char *pattern,
                                           const char *text,
                                           int mapping_start_position,
                                           int n_cigar, const uint32_t *cigar,
                                           int &NM, std::string &MD_tag) {
  int num_matches = 0;
  const char *read = text;
  const char *reference = pattern + mapping_start_position;
  int read_position = 0;
  int reference_position = 0;
  for (int ci = 0; ci < n_cigar; ++ci) {
    uint32_t current_cigar_uint = cigar[ci];
    uint8_t cigar_operation = bam_cigar_op(current_cigar_uint);
    int num_cigar_operations = bam_cigar_oplen(current_cigar_uint);
    if (cigar_operation == BAM_CMATCH) {
      for (int opi = 0; opi < num_cigar_operations; ++opi) {
        if (reference[reference_position] == read[read_position] ||
            reference[reference_position] - 'a' + 'A' == read[read_position]) {
          // a match
          ++num_matches;
        } else {
          // a mismatch
          ++NM;
          if (num_matches != 0) {
            MD_tag.append(std::to_string(num_matches));
            num_matches = 0;
          }
          MD_tag.push_back(reference[reference_position]);
        }
        ++reference_position;
        ++read_position;
      }
    } else if (cigar_operation == BAM_CINS) {
      NM += num_cigar_operations;
      read_position += num_cigar_operations;
    } else if (cigar_operation == BAM_CDEL) {
      NM += num_cigar_operations;
      if (num_matches != 0) {
        MD_tag.append(std::to_string(num_matches));
        num_matches = 0;
      }
      MD_tag.push_back('^');
      for (int opi = 0; opi < num_cigar_operations; ++opi) {
        MD_tag.push_back(reference[reference_position]);
        ++reference_position;
      }
    } else {
      std::cerr << "Unexpected cigar op: " << (int)cigar_operation << "\n";
    }
  }
  if (num_matches != 0) {
    MD_tag.append(std::to_string(num_matches));
  }
}

// return: newly adjust reference start/end position (kPositive for start,
// kNegative for end)
template <typename MappingRecord>
int Chromap<MappingRecord>::AdjustGapBeginning(
    Direction mapping_direction, const char *ref, const char *read,
    int *gap_beginning, int read_end, int ref_start_position,
    int ref_end_position, int *n_cigar, uint32_t **cigar) {
  int i, j;
  if (mapping_direction == kPositive) {
    if (*gap_beginning <= 0) {
      return ref_start_position;
    }
    // printf("%d\n", *gap_beginning);
    for (i = *gap_beginning - 1, j = ref_start_position - 1; i >= 0 && j >= 0;
         --i, --j) {
      // printf("%c %c\n", read[i], ref[j]);
      if (read[i] != ref[j] && read[i] != ref[j] - 'a' + 'A') {
        break;
      }
    }
    *gap_beginning = i + 1;
    // TODO: add soft clip in cigar
    if (n_cigar && *n_cigar > 0) {
      if (((*cigar)[0] & 0xf) == BAM_CMATCH) {
        (*cigar)[0] += (ref_start_position - 1 - j) << 4;
      }
    }
    return j + 1;
  } else {
    if (*gap_beginning <= 0) {
      return ref_end_position;
    }
    // printf("%d\n", *gap_beginning);
    /*char *tmp = new char[255] ;
    strncpy(tmp, ref + ref_start_position, ref_end_position - ref_start_position
    + 1 + 10) ; printf("%s %d. %d %d\n", tmp, strlen(tmp), ref_end_position -
    ref_start_position + 1 + 10, strlen(ref)) ; delete[] tmp;*/
    for (i = read_end + 1, j = ref_end_position + 1; read[i] && ref[j];
         ++i, ++j) {
      // printf("%c %c %c %c %c %c\n", read[i], ref[j - 1], ref[j], ref[j + 1],
      // ref[j + 2], ref[j + 3]);
      if (read[i] != ref[j] && read[i] != ref[j] - 'a' + 'A') {
        break;
      }
    }
    *gap_beginning = *gap_beginning + i - (read_end + 1);
    if (n_cigar && *n_cigar > 0) {
      if (((*cigar)[*n_cigar - 1] & 0xf) == BAM_CMATCH) {
        (*cigar)[*n_cigar - 1] += (j - (ref_end_position + 1)) << 4;
      }
    }

    return j - 1;
  }
}

// The returned coordinate is left closed and right closed, and is without chrom
// id.
template <typename MappingRecord>
void Chromap<MappingRecord>::GetRefStartEndPositionForReadFromMapping(
    Direction mapping_direction, const std::pair<int, uint64_t> &mapping,
    const char *read, int read_length, int in_split_site,
    const SequenceBatch &reference, uint32_t *ref_start_position,
    uint32_t *ref_end_position, int *n_cigar, uint32_t **cigar, int *NM,
    std::string &MD_tag) {
  int8_t mat[25];
  // if (output_mapping_in_SAM_) {
  int i, j, k;
  for (i = k = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j)
      mat[k++] = i == j ? match_score_ : -mismatch_penalty_;
    mat[k++] = 0;  // ambiguous base
  }
  for (j = 0; j < 5; ++j) mat[k++] = 0;
  //}

  uint32_t rid = mapping.second >> 32;
  uint32_t position = mapping.second;
  int full_read_length = read_length;
  int min_num_errors = mapping.first;
  int split_site = mapping_direction == kPositive ? 0 : read_length;
  int gap_beginning = 0;
  int actual_num_errors = 0;
  if (split_alignment_) {
    split_site = in_split_site & 0xffff;
    gap_beginning =
        (in_split_site >> 16) & 0xff;  // beginning means the 5' end of the read
    actual_num_errors =
        (in_split_site >> 24) &
        0xff;  // in split alignment, -num_errors is the number of matches.
    read_length = split_site - gap_beginning;
  }
  uint32_t verification_window_start_position =
      position + 1 > (uint32_t)(read_length + error_threshold_)
          ? position + 1 - read_length - error_threshold_
          : 0;
  // printf("ne4: %d %d. %d %d %d %d.\n", position,
  // verification_window_start_position, read_length, split_site, gap_beginning,
  // actual_num_errors);
  if (position + error_threshold_ >= reference.GetSequenceLengthAt(rid)) {
    verification_window_start_position =
        reference.GetSequenceLengthAt(rid) - error_threshold_ - read_length;
  }
  if (verification_window_start_position < 0) {
    verification_window_start_position = 0;
  }
  if (split_alignment_) {
    if (split_site < full_read_length &&
        mapping_output_format_ == MAPPINGFORMAT_SAM &&
        split_site > 3 * error_threshold_) {
      split_site -= 3 * error_threshold_;
    }
    read_length = split_site - gap_beginning;
  }
  int mapping_start_position;
  if (mapping_direction == kPositive) {
    if (mapping_output_format_ == MAPPINGFORMAT_SAM) {
      *n_cigar = 0;
      int mapping_end_position;
      ksw_semi_global3(
          read_length + 2 * error_threshold_,
          reference.GetSequenceAt(rid) + verification_window_start_position,
          read_length, read + gap_beginning, 5, mat, gap_open_penalties_[0],
          gap_extension_penalties_[0], gap_open_penalties_[1],
          gap_extension_penalties_[1], error_threshold_ * 2 + 1, n_cigar, cigar,
          &mapping_start_position, &mapping_end_position);
      // std::cerr << verification_window_start_position << " " << read_length
      // << " " << split_site << " " << mapping_start_position << " " <<
      // mapping_end_position << "\n";
      if (gap_beginning > 0) {
        int new_ref_start_position = AdjustGapBeginning(
            mapping_direction, reference.GetSequenceAt(rid), read,
            &gap_beginning, read_length - 1,
            verification_window_start_position + mapping_start_position,
            verification_window_start_position + mapping_end_position - 1,
            n_cigar, cigar);
        mapping_start_position =
            new_ref_start_position - verification_window_start_position;
      }
      GenerateMDTag(reference.GetSequenceAt(rid), read + gap_beginning,
                    verification_window_start_position + mapping_start_position,
                    *n_cigar, *cigar, *NM, MD_tag);

      *ref_start_position =
          verification_window_start_position + mapping_start_position;
      *ref_end_position =
          verification_window_start_position + mapping_end_position - 1;
    } else {
      if (!split_alignment_) {
        BandedTraceback(
            min_num_errors,
            reference.GetSequenceAt(rid) + verification_window_start_position,
            read, read_length, &mapping_start_position);
      } else {
        BandedTraceback(
            actual_num_errors,
            reference.GetSequenceAt(rid) + verification_window_start_position,
            read + gap_beginning, read_length, &mapping_start_position);
      }
      if (gap_beginning > 0) {
        int new_ref_start_position = AdjustGapBeginning(
            mapping_direction, reference.GetSequenceAt(rid), read,
            &gap_beginning, read_length - 1,
            verification_window_start_position + mapping_start_position,
            position, n_cigar, cigar);
        mapping_start_position =
            new_ref_start_position - verification_window_start_position;
      }

      *ref_start_position =
          verification_window_start_position + mapping_start_position;
      *ref_end_position = position;
    }
  } else {  // reverse strand
    int read_start_site = full_read_length - split_site;
    if (mapping_output_format_ == MAPPINGFORMAT_SAM) {
      *n_cigar = 0;
      int mapping_end_position;

      //  reversed read looks like:
      //
      //      veri_start_pos       position
      //  ref   --|-------------------|------------------->
      //  read     <-|--read_length---|--gap_beginning--
      //          split_site
      //

      // Is the map start/end position left close right open as in bed format?
      /*if (read_start_site + read_length > full_read_length || read_start_site
      < 0
        || verification_window_start_position + read_start_site + read_length +
      2 * error_threshold_ > reference.GetSequenceLengthAt(rid)) {
        printf("ERROR! %d %d %d %d\n", full_read_length, split_site,
      gap_beginning, read_length);
      }*/
      ksw_semi_global3(read_length + 2 * error_threshold_,
                       reference.GetSequenceAt(rid) +
                           verification_window_start_position + read_start_site,
                       read_length, read + read_start_site, 5, mat,
                       gap_open_penalties_[0], gap_extension_penalties_[0],
                       gap_open_penalties_[1], gap_extension_penalties_[1],
                       error_threshold_ * 2 + 1, n_cigar, cigar,
                       &mapping_start_position, &mapping_end_position);
      if (gap_beginning > 0) {
        // printf("before adjust %d %d %d. %d %d\n", split_site,
        // read_start_site, read_length, verification_window_start_position +
        // mapping_start_position, verification_window_start_position +
        // mapping_end_position);
        int new_ref_end_position = AdjustGapBeginning(
            mapping_direction, reference.GetSequenceAt(rid),
            read + read_start_site, &gap_beginning, read_length - 1,
            verification_window_start_position + mapping_start_position,
            verification_window_start_position + mapping_end_position - 1,
            n_cigar, cigar);
        // The returned position is right-closed, so need to plus one to match
        // bed convention
        mapping_end_position = new_ref_end_position + 1 -
                               verification_window_start_position -
                               read_start_site;
        read_length = split_site - gap_beginning;
      }
      GenerateMDTag(reference.GetSequenceAt(rid), read + read_start_site,
                    verification_window_start_position + read_start_site +
                        mapping_start_position,
                    *n_cigar, *cigar, *NM, MD_tag);
      *ref_start_position = verification_window_start_position +
                            read_start_site + mapping_start_position;
      *ref_end_position = verification_window_start_position + read_start_site +
                          mapping_end_position - 1;
    } else {
      // int n_cigar = 0;
      // uint32_t *cigar;
      int mapping_end_position =
          position - verification_window_start_position + 1;
      mapping_start_position = error_threshold_;
      // ksw_semi_global3(read_length + 2 * error_threshold_,
      // reference.GetSequenceAt(rid) + verification_window_start_position,
      // read_length, negative_read.data() + split_sites[mi], 5, mat,
      // gap_open_penalties_[0], gap_extension_penalties_[0],
      // gap_open_penalties_[1], gap_extension_penalties_[1], error_threshold_ *
      // 2 + 1, &n_cigar, &cigar, &mapping_start_position,
      // &mapping_end_position); mapq =
      // GetMAPQForSingleEndRead(error_threshold_, 0, 0, mapping_end_position -
      // mapping_start_position + 1, min_num_errors, num_best_mappings,
      // second_min_num_errors, num_second_best_mappings); uint32_t
      // fragment_start_position = verification_window_start_position +
      // mapping_start_position; uint16_t fragment_length = mapping_end_position
      // - mapping_start_position + 1;
      if (!split_alignment_) {
        BandedTraceback(
            min_num_errors,
            reference.GetSequenceAt(rid) + verification_window_start_position,
            read + read_start_site, read_length, &mapping_start_position);
      } else {
        // BandedTracebackToEnd(actual_num_errors, reference.GetSequenceAt(rid)
        // + verification_window_start_position, read + read_start_site,
        // read_length, &mapping_end_position);
        BandedAlignPatternToText(
            reference.GetSequenceAt(rid) + verification_window_start_position,
            read + read_start_site, read_length, &mapping_end_position);
        mapping_end_position +=
            1;  // seems banded align's mapping end position is included?
      }
      // int mapping_end_position = (int)position;
      if (gap_beginning > 0) {
        // printf("before adjust %d %d %d. %d %d\n", split_site,
        // read_start_site, read_length, verification_window_start_position +
        // mapping_start_position, verification_window_start_position +
        // mapping_end_position);
        int new_ref_end_position = AdjustGapBeginning(
            mapping_direction, reference.GetSequenceAt(rid),
            read + read_start_site, &gap_beginning, read_length - 1,
            verification_window_start_position + mapping_start_position,
            verification_window_start_position + mapping_end_position, n_cigar,
            cigar);
        // The returned position is right-closed, so need to plus one to match
        // bed convention
        mapping_end_position =
            new_ref_end_position - verification_window_start_position + 1;
        read_length = split_site - gap_beginning;
        // position = new_ref_end_position;
      }
      *ref_start_position =
          verification_window_start_position + mapping_start_position;
      //*ref_end_position = position;
      *ref_end_position =
          verification_window_start_position + mapping_end_position - 1;
    }
  }
}

template <typename MappingRecord>
void Chromap<MappingRecord>::ProcessBestMappingsForSingleEndRead(
    Direction mapping_direction, uint8_t mapq, int num_candidates,
    uint32_t repetitive_seed_length, int min_num_errors, int num_best_mappings,
    int second_min_num_errors, int num_second_best_mappings,
    const SequenceBatch &read_batch, uint32_t read_index,
    const SequenceBatch &reference, const SequenceBatch &barcode_batch,
    const std::vector<int> &best_mapping_indices,
    const std::vector<std::pair<int, uint64_t>> &mappings,
    const std::vector<int> &split_sites, int *best_mapping_index,
    int *num_best_mappings_reported,
    std::vector<std::vector<MappingRecord>> *mappings_on_diff_ref_seqs) {
  const char *read = read_batch.GetSequenceAt(read_index);
  uint32_t read_id = read_batch.GetSequenceIdAt(read_index);
  const char *read_name = read_batch.GetSequenceNameAt(read_index);
  uint32_t read_length = read_batch.GetSequenceLengthAt(read_index);
  const std::string &negative_read =
      read_batch.GetNegativeSequenceAt(read_index);
  uint8_t is_unique = num_best_mappings == 1 ? 1 : 0;
  uint64_t barcode_key = 0;
  if (!is_bulk_data_) {
    barcode_key = barcode_batch.GenerateSeedFromSequenceAt(
        read_index, 0, barcode_batch.GetSequenceLengthAt(read_index));
  }
  for (uint32_t mi = 0; mi < mappings.size(); ++mi) {
    if (mappings[mi].first == min_num_errors) {
      if (*best_mapping_index ==
          best_mapping_indices[*num_best_mappings_reported]) {
        read_length = read_batch.GetSequenceLengthAt(read_index);
        uint32_t ref_start_position;
        uint32_t ref_end_position;
        uint8_t direction = 1;

        uint32_t *cigar;
        int n_cigar = 0;
        int NM = 0;
        std::string MD_tag = "";
        const char *effect_read = read;
        if (mapping_direction == kNegative) {
          direction = 0;
          effect_read = negative_read.data();
        }
        uint32_t rid = mappings[mi].second >> 32;

        int split_site = 0;
        if (split_alignment_) {
          split_site = split_sites[mi];
        }
        // printf("%d %d\n", split_site, read_length);
        GetRefStartEndPositionForReadFromMapping(
            mapping_direction, mappings[mi], effect_read, read_length,
            split_site, reference, &ref_start_position, &ref_end_position,
            &n_cigar, &cigar, &NM, MD_tag);
        mapq = GetMAPQForSingleEndRead(
            error_threshold_, num_candidates, repetitive_seed_length,
            ref_end_position - ref_start_position + 1, min_num_errors,
            num_best_mappings, second_min_num_errors, num_second_best_mappings,
            error_threshold_, read_length);

        if (mapping_output_format_ == MAPPINGFORMAT_SAM) {
          uint16_t flag = mapping_direction == kPositive ? 0 : BAM_FREVERSE;
          if (*num_best_mappings_reported >= 1) {
            flag |= BAM_FSECONDARY;
          }
          EmplaceBackMappingRecord(read_id, read_name, barcode_key, 1,
                                   ref_start_position, rid, flag, 0, is_unique,
                                   mapq, NM, n_cigar, cigar, MD_tag, read,
                                   read_batch.GetSequenceQualAt(read_index),
                                   &((*mappings_on_diff_ref_seqs)[rid]));
        } else if (mapping_output_format_ == MAPPINGFORMAT_PAF) {
          EmplaceBackMappingRecord(
              read_id, read_name, read_length, barcode_key, ref_start_position,
              ref_end_position - ref_start_position + 1, mapq, direction,
              is_unique, 1, &((*mappings_on_diff_ref_seqs)[rid]));
        } else {
          EmplaceBackMappingRecord(read_id, barcode_key, ref_start_position,
                                   ref_end_position - ref_start_position + 1,
                                   mapq, direction, is_unique, 1,
                                   &((*mappings_on_diff_ref_seqs)[rid]));
        }
        (*num_best_mappings_reported)++;
        if (*num_best_mappings_reported ==
            std::min(max_num_best_mappings_, num_best_mappings)) {
          break;
        }
      }
      (*best_mapping_index)++;
    }
  }
}

template <typename MappingRecord>
void Chromap<MappingRecord>::GenerateBestMappingsForSingleEndRead(
    int num_positive_candidates, int num_negative_candidates,
    uint32_t repetitive_seed_length, int min_num_errors, int num_best_mappings,
    int second_min_num_errors, int num_second_best_mappings,
    const SequenceBatch &read_batch, uint32_t read_index,
    const SequenceBatch &reference, const SequenceBatch &barcode_batch,
    const std::vector<std::pair<int, uint64_t>> &positive_mappings,
    const std::vector<int> &positive_split_sites,
    const std::vector<std::pair<int, uint64_t>> &negative_mappings,
    const std::vector<int> &negative_split_sites,
    std::vector<std::vector<MappingRecord>> *mappings_on_diff_ref_seqs) {
  // uint8_t mapq = GetMAPQ(num_best_mappings, num_second_best_mappings);
  uint8_t mapq = 0;
  // we will use reservoir sampling
  std::vector<int> best_mapping_indices(max_num_best_mappings_);
  std::iota(best_mapping_indices.begin(), best_mapping_indices.end(), 0);
  if (num_best_mappings > max_num_best_mappings_) {
    std::mt19937 generator(11);
    for (int i = max_num_best_mappings_; i < num_best_mappings; ++i) {
      std::uniform_int_distribution<int> distribution(
          0, i);  // important: inclusive range
      int j = distribution(generator);
      if (j < max_num_best_mappings_) {
        best_mapping_indices[j] = i;
      }
    }
    std::sort(best_mapping_indices.begin(), best_mapping_indices.end());
  }
  int best_mapping_index = 0;
  int num_best_mappings_reported = 0;
  ProcessBestMappingsForSingleEndRead(
      kPositive, mapq, num_positive_candidates, repetitive_seed_length,
      min_num_errors, num_best_mappings, second_min_num_errors,
      num_second_best_mappings, read_batch, read_index, reference,
      barcode_batch, best_mapping_indices, positive_mappings,
      positive_split_sites, &best_mapping_index, &num_best_mappings_reported,
      mappings_on_diff_ref_seqs);
  if (num_best_mappings_reported !=
      std::min(num_best_mappings, max_num_best_mappings_)) {
    ProcessBestMappingsForSingleEndRead(
        kNegative, num_negative_candidates, repetitive_seed_length, mapq,
        min_num_errors, num_best_mappings, second_min_num_errors,
        num_second_best_mappings, read_batch, read_index, reference,
        barcode_batch, best_mapping_indices, negative_mappings,
        negative_split_sites, &best_mapping_index, &num_best_mappings_reported,
        mappings_on_diff_ref_seqs);
  }
}

template <typename MappingRecord>
void Chromap<MappingRecord>::ApplyTn5ShiftOnSingleEndMapping(
    uint32_t num_reference_sequences,
    std::vector<std::vector<MappingRecord>> *mappings) {
  uint64_t num_shifted_mappings = 0;
  for (auto &mappings_on_one_ref_seq : *mappings) {
    for (auto &mapping : mappings_on_one_ref_seq) {
      mapping.Tn5Shift();
      // uint8_t strand = mapping.direction & 1;
      // if (strand == 1) {
      //  mapping.fragment_start_position += 4;
      //  mapping.fragment_length -= 4;
      //} else {
      //  mapping.fragment_length -= 5;
      //}
      ++num_shifted_mappings;
    }
  }
  std::cerr << "# shifted mappings: " << num_shifted_mappings << ".\n";
}

template <typename MappingRecord>
uint32_t Chromap<MappingRecord>::LoadSingleEndReadsWithBarcodes(
    SequenceBatch *read_batch, SequenceBatch *barcode_batch) {
  double real_start_time = Chromap<>::GetRealTime();
  uint32_t num_loaded_reads = 0;
  while (num_loaded_reads < read_batch_size_) {
    bool no_more_read = read_batch->LoadOneSequenceAndSaveAt(num_loaded_reads);
    bool no_more_barcode = no_more_read;
    if (!is_bulk_data_) {
      no_more_barcode =
          barcode_batch->LoadOneSequenceAndSaveAt(num_loaded_reads);
    }
    if ((!no_more_read) && (!no_more_barcode)) {
      if (read_batch->GetSequenceLengthAt(num_loaded_reads) <
          (uint32_t)min_read_length_) {
        continue;  // reads are too short, just drop.
      }
      // if (PairedEndReadWithBarcodeIsDuplicate(num_loaded_pairs,
      // (*barcode_batch), (*read_batch1), (*read_batch2))) {
      //  num_duplicated_reads_ += 2;
      //  continue;
      //}
    } else if (no_more_read && no_more_barcode) {
      break;
    } else {
      Chromap<>::ExitWithMessage("Numbers of reads and barcodes don't match!");
    }
    ++num_loaded_reads;
  }
  if (num_loaded_reads > 0) {
    std::cerr << "Loaded " << num_loaded_reads << " reads in "
              << Chromap<>::GetRealTime() - real_start_time << "s.\n";
  } else {
    std::cerr << "No more reads.\n";
  }
  return num_loaded_reads;
}

template <typename MappingRecord>
void Chromap<MappingRecord>::ConstructIndex() {
  // TODO(Haowen): Need a faster algorithm
  // Load all sequences in the reference into one batch
  SequenceBatch reference;
  reference.InitializeLoading(reference_file_path_);
  uint32_t num_sequences = reference.LoadAllSequences();
  Index index(kmer_size_, window_size_, num_threads_, index_file_path_);
  index.Construct(num_sequences, reference);
  index.Statistics(num_sequences, reference);
  index.Save();
  reference.FinalizeLoading();
}

template <typename MappingRecord>
uint32_t Chromap<MappingRecord>::MoveMappingsInBuffersToMappingContainer(
    uint32_t num_reference_sequences,
    std::vector<std::vector<std::vector<MappingRecord>>>
        *mappings_on_diff_ref_seqs_for_diff_threads_for_saving) {
  // double real_start_time = Chromap<>::GetRealTime();
  uint32_t num_moved_mappings = 0;
  for (int ti = 0; ti < num_threads_; ++ti) {
    for (uint32_t i = 0; i < num_reference_sequences; ++i) {
      num_moved_mappings +=
          (*mappings_on_diff_ref_seqs_for_diff_threads_for_saving)[ti][i]
              .size();
      mappings_on_diff_ref_seqs_[i].insert(
          mappings_on_diff_ref_seqs_[i].end(),
          std::make_move_iterator(
              (*mappings_on_diff_ref_seqs_for_diff_threads_for_saving)[ti][i]
                  .begin()),
          std::make_move_iterator(
              (*mappings_on_diff_ref_seqs_for_diff_threads_for_saving)[ti][i]
                  .end()));
      (*mappings_on_diff_ref_seqs_for_diff_threads_for_saving)[ti][i].clear();
    }
  }
  // std::cerr << "Moved mappings in " << Chromap<>::GetRealTime() -
  // real_start_time << "s.\n";
  return num_moved_mappings;
}

template <typename MappingRecord>
void Chromap<MappingRecord>::SortOutputMappings(
    uint32_t num_reference_sequences,
    std::vector<std::vector<MappingRecord>> *mappings) {
  // double real_dedupe_start_time = Chromap<>::GetRealTime();
  uint32_t num_mappings = 0;
  for (uint32_t ri = 0; ri < num_reference_sequences; ++ri) {
    std::sort((*mappings)[ri].begin(), (*mappings)[ri].end());
    num_mappings += (*mappings)[ri].size();
  }
  // std::cerr << "Sorted " << num_mappings << " elements in " <<
  // Chromap<>::GetRealTime() - real_dedupe_start_time << "s.\n";
}

template <typename MappingRecord>
void Chromap<MappingRecord>::RemovePCRDuplicate(
    uint32_t num_reference_sequences) {
  uint32_t num_mappings = 0;
  double real_dedupe_start_time = Chromap<>::GetRealTime();
  for (uint32_t ri = 0; ri < num_reference_sequences; ++ri) {
    // double real_start_time = Chromap<>::GetRealTime();
    // radix_sort_with_barcode(mappings_on_diff_ref_seqs_[ri].data(),
    // mappings_on_diff_ref_seqs_[ri].data() +
    // mappings_on_diff_ref_seqs_[ri].size());
    std::sort(mappings_on_diff_ref_seqs_[ri].begin(),
              mappings_on_diff_ref_seqs_[ri].end());
    num_mappings += mappings_on_diff_ref_seqs_[ri].size();
    // std::cerr << "Sorted " << mappings_on_diff_ref_seqs_[ri].size() << "
    // elements on " << ri << " in " << Chromap<>::GetRealTime() -
    // real_start_time << "s.\n";
  }
  std::cerr << "Sorted " << num_mappings << " elements in "
            << Chromap<>::GetRealTime() - real_dedupe_start_time << "s.\n";
  num_mappings = 0;
  for (uint32_t ri = 0; ri < num_reference_sequences; ++ri) {
    if (mappings_on_diff_ref_seqs_[ri].size() != 0) {
      deduped_mappings_on_diff_ref_seqs_[ri].emplace_back(
          mappings_on_diff_ref_seqs_[ri]
              .front());  // ideally I should output the last of the dups of
                          // first mappings.
      // std::vector<MappingRecord>::iterator last_it =
      // mappings_on_diff_ref_seqs_[ri].begin();
      auto last_it = mappings_on_diff_ref_seqs_[ri].begin();
      uint32_t last_dup_count = 1;
      // for (std::vector<MappingRecord>::iterator it =
      // ++(mappings_on_diff_ref_seqs_[ri].begin()); it !=
      // mappings_on_diff_ref_seqs_[ri].end(); ++it) {
      for (auto it = ++(mappings_on_diff_ref_seqs_[ri].begin());
           it != mappings_on_diff_ref_seqs_[ri].end(); ++it) {
        if (!((*it) == (*last_it))) {
          // last_it->num_dups = last_dup_count;
          deduped_mappings_on_diff_ref_seqs_[ri].back().num_dups = std::min(
              (uint32_t)std::numeric_limits<uint8_t>::max(), last_dup_count);
          last_dup_count = 1;
          deduped_mappings_on_diff_ref_seqs_[ri].emplace_back((*it));
          last_it = it;
        } else {
          ++last_dup_count;
        }
      }
      deduped_mappings_on_diff_ref_seqs_[ri].back().num_dups = std::min(
          (uint32_t)std::numeric_limits<uint8_t>::max(), last_dup_count);
      std::vector<MappingRecord>().swap(mappings_on_diff_ref_seqs_[ri]);
      num_mappings += deduped_mappings_on_diff_ref_seqs_[ri].size();
    }
  }
  std::cerr << num_mappings << " mappings left after dedupe in "
            << Chromap<>::GetRealTime() - real_dedupe_start_time << "s.\n";
}

template <typename MappingRecord>
void Chromap<MappingRecord>::BuildAugmentedTree(uint32_t ref_id) {
  // std::sort(mappings.begin(), mappings.end(), IntervalLess());
  int max_level = 0;
  size_t i, last_i = 0;  // last_i points to the rightmost node in the tree
  uint32_t last = 0;     // last is the max value at node last_i
  int k;
  std::vector<MappingRecord> &mappings =
      allocated_mappings_on_diff_ref_seqs_[ref_id];
  std::vector<uint32_t> &extras = tree_extras_on_diff_ref_seqs_[ref_id];
  if (mappings.size() == 0) {
    max_level = -1;
  }
  for (i = 0; i < mappings.size(); i += 2) {
    last_i = i;
    // last = mappings[i].max = mappings[i].en; // leaves (i.e. at level 0)
    last = extras[i] =
        mappings[i].GetEndPosition();  // leaves (i.e. at level 0)
  }
  for (k = 1; 1LL << k <= (int64_t)mappings.size();
       ++k) {  // process internal nodes in the bottom-up order
    size_t x = 1LL << (k - 1);
    size_t i0 = (x << 1) - 1;
    size_t step = x << 2;  // i0 is the first node
    for (i = i0; i < mappings.size();
         i += step) {               // traverse all nodes at level k
      uint32_t el = extras[i - x];  // max value of the left child
      uint32_t er =
          i + x < mappings.size() ? extras[i + x] : last;  // of the right child
      uint32_t e = mappings[i].GetEndPosition();
      e = e > el ? e : el;
      e = e > er ? e : er;
      extras[i] = e;  // set the max value for node i
    }
    last_i =
        last_i >> k & 1
            ? last_i - x
            : last_i +
                  x;  // last_i now points to the parent of the original last_i
    if (last_i < mappings.size() &&
        extras[last_i] > last)  // update last accordingly
      last = extras[last_i];
  }
  max_level = k - 1;
  tree_info_on_diff_ref_seqs_.emplace_back(max_level, mappings.size());
}

template <typename MappingRecord>
uint32_t Chromap<MappingRecord>::GetNumOverlappedMappings(
    uint32_t ref_id, const MappingRecord &mapping) {
  int t = 0;
  StackCell stack[64];
  // out.clear();
  int num_overlapped_mappings = 0;
  int max_level = tree_info_on_diff_ref_seqs_[ref_id].first;
  uint32_t num_tree_nodes = tree_info_on_diff_ref_seqs_[ref_id].second;
  std::vector<MappingRecord> &mappings =
      allocated_mappings_on_diff_ref_seqs_[ref_id];
  std::vector<uint32_t> &extras = tree_extras_on_diff_ref_seqs_[ref_id];
  // uint32_t interval_start = mapping.fragment_start_position;
  uint32_t interval_start =
      mapping.GetStartPosition() > (uint32_t)multi_mapping_allocation_distance_
          ? mapping.GetStartPosition() - multi_mapping_allocation_distance_
          : 0;
  uint32_t interval_end =
      mapping.GetEndPosition() + (uint32_t)multi_mapping_allocation_distance_;
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
      for (i = i0; i < i1 && mappings[i].GetStartPosition() < interval_end;
           ++i) {
        if (interval_start <
            mappings[i].GetEndPosition()) {  // if overlap, append to out[]
          // out.push_back(i);
          ++num_overlapped_mappings;
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
               mappings[z.x].GetStartPosition() <
                   interval_end) {  // need to push the right child
      if (interval_start < mappings[z.x].GetEndPosition()) {
        // out.push_back(z.x); // test if z.x overlaps the query; if yes, append
        // to out[]
        ++num_overlapped_mappings;
      }
      stack[t++] = StackCell(z.k - 1, z.x + (1LL << (z.k - 1)),
                             0);  // push the right child
    }
  }
  return num_overlapped_mappings;
}

template <typename MappingRecord>
void Chromap<MappingRecord>::AllocateMultiMappings(
    uint32_t num_reference_sequences) {
  double real_start_time = Chromap<>::GetRealTime();
  std::vector<std::vector<MappingRecord>> &mappings =
      remove_pcr_duplicates_ ? deduped_mappings_on_diff_ref_seqs_
                             : mappings_on_diff_ref_seqs_;
  multi_mappings_.reserve((num_mapped_reads_ - num_uniquely_mapped_reads_) / 2);
  allocated_mappings_on_diff_ref_seqs_.reserve(num_reference_sequences);
  tree_extras_on_diff_ref_seqs_.reserve(num_reference_sequences);
  tree_info_on_diff_ref_seqs_.reserve(num_reference_sequences);
  // two passes, one for memory pre-allocation, another to move the mappings.
  for (uint32_t ri = 0; ri < num_reference_sequences; ++ri) {
    allocated_mappings_on_diff_ref_seqs_.emplace_back(
        std::vector<MappingRecord>());
    tree_extras_on_diff_ref_seqs_.emplace_back(std::vector<uint32_t>());
    uint32_t num_uni_mappings = 0;
    uint32_t num_multi_mappings = 0;
    for (uint32_t mi = 0; mi < mappings[ri].size(); ++mi) {
      MappingRecord &mapping = mappings[ri][mi];
      if ((mapping.mapq) <
          min_unique_mapping_mapq_) {  // we have to ensure that the mapq is
                                       // lower than this if and only if it is a
                                       // multi-read.
        ++num_multi_mappings;
      } else {
        ++num_uni_mappings;
      }
    }
    allocated_mappings_on_diff_ref_seqs_[ri].reserve(num_uni_mappings);
    tree_extras_on_diff_ref_seqs_[ri].reserve(num_uni_mappings);
    for (uint32_t mi = 0; mi < mappings[ri].size(); ++mi) {
      MappingRecord &mapping = mappings[ri][mi];
      if ((mapping.mapq) < min_unique_mapping_mapq_) {
        multi_mappings_.emplace_back(ri, mapping);
      } else {
        allocated_mappings_on_diff_ref_seqs_[ri].emplace_back(mapping);
        tree_extras_on_diff_ref_seqs_[ri].emplace_back(0);
      }
    }
    std::vector<MappingRecord>().swap(mappings[ri]);
    BuildAugmentedTree(ri);
  }
  std::cerr << "Got all " << multi_mappings_.size() << " multi-mappings!\n";
  std::stable_sort(multi_mappings_.begin(), multi_mappings_.end(),
                   ReadIdLess<MappingRecord>);
  std::vector<uint32_t> weights;
  weights.reserve(max_num_best_mappings_);
  uint32_t sum_weight = 0;
  assert(multi_mappings_.size() > 0);
  uint32_t previous_read_id = multi_mappings_[0].second.read_id;
  uint32_t start_mapping_index = 0;
  // add a fake mapping at the end and make sure its id is different from the
  // last one
  assert(multi_mappings_.size() != UINT32_MAX);
  std::pair<uint32_t, MappingRecord> foo_mapping = multi_mappings_.back();
  foo_mapping.second.read_id = UINT32_MAX;
  multi_mappings_.emplace_back(foo_mapping);
  std::mt19937 generator(multi_mapping_allocation_seed_);
  uint32_t current_read_id;  //, reference_id, mapping_index;
  // uint32_t allocated_read_id, allocated_reference_id,
  // allocated_mapping_index;
  uint32_t num_allocated_multi_mappings = 0;
  uint32_t num_multi_mappings_without_overlapping_unique_mappings = 0;
  for (uint32_t mi = 0; mi < multi_mappings_.size(); ++mi) {
    std::pair<uint32_t, MappingRecord> &current_multi_mapping =
        multi_mappings_[mi];  // mappings[reference_id][mapping_index];
    current_read_id = current_multi_mapping.second.read_id;
    uint32_t num_overlaps = GetNumOverlappedMappings(
        current_multi_mapping.first, current_multi_mapping.second);
    // std::cerr << mi << " " << current_read_id << " " << previous_read_id << "
    // " << reference_id << " " << mapping_index << " " << interval_start << " "
    // << num_overlaps << " " << sum_weight << "\n";
    if (current_read_id == previous_read_id) {
      weights.emplace_back(num_overlaps);
      sum_weight += num_overlaps;
    } else {
      // deal with the previous one.
      if (sum_weight == 0) {
        ++num_multi_mappings_without_overlapping_unique_mappings;
        // assert(weights.size() > 1); // After PCR dedupe, some multi-reads may
        // become uni-reads. For now, we just assign it to that unique mapping
        // positions. std::fill(weights.begin(), weights.end(), 1); // We drop
        // the multi-mappings that have no overlap with uni-mappings.
      } else {
        std::discrete_distribution<uint32_t> distribution(weights.begin(),
                                                          weights.end());
        uint32_t randomly_assigned_mapping_index = distribution(generator);
        allocated_mappings_on_diff_ref_seqs_
            [multi_mappings_[start_mapping_index +
                             randomly_assigned_mapping_index]
                 .first]
                .emplace_back(multi_mappings_[start_mapping_index +
                                              randomly_assigned_mapping_index]
                                  .second);
        ++num_allocated_multi_mappings;
      }
      // update current
      weights.clear();
      weights.emplace_back(num_overlaps);
      sum_weight = num_overlaps;
      start_mapping_index = mi;
      previous_read_id = current_read_id;
    }
  }
  std::cerr << "Allocated " << num_allocated_multi_mappings
            << " multi-mappings in "
            << Chromap<>::GetRealTime() - real_start_time << "s.\n";
  std::cerr << "# multi-mappings that have no uni-mapping overlaps: "
            << num_multi_mappings_without_overlapping_unique_mappings << ".\n";
}

template <typename MappingRecord>
void Chromap<MappingRecord>::MergeCandidates(std::vector<Candidate> &c1,
                                             std::vector<Candidate> &c2,
                                             std::vector<Candidate> &buffer) {
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
          c1[i].position > buffer.back().position + error_threshold_) {
        if (c1[i].count > c2[j].count) {
          buffer.push_back(c1[i]);
        } else {
          buffer.push_back(c2[j]);
        }
      }
      ++i, ++j;
    } else if (c1[i].position < c2[j].position) {
      if (buffer.empty() ||
          c1[i].position > buffer.back().position + error_threshold_) {
        buffer.push_back(c1[i]);
      }
      ++i;
    } else {
      if (buffer.empty() ||
          c2[j].position > buffer.back().position + error_threshold_) {
        buffer.push_back(c2[j]);
      }
      ++j;
    }
  }
  while (i < size1) {
    if (buffer.empty() ||
        c1[i].position > buffer.back().position + error_threshold_) {
      buffer.push_back(c1[i]);
    }
    ++i;
  }
  while (j < size2) {
    if (buffer.empty() ||
        c2[j].position > buffer.back().position + error_threshold_) {
      buffer.push_back(c2[j]);
    }
    ++j;
  }
  // c1 = buffer;
  c1.swap(buffer);
}

template <typename MappingRecord>
void Chromap<MappingRecord>::VerifyCandidatesOnOneDirectionUsingSIMD(
    Direction candidate_direction, const SequenceBatch &read_batch,
    uint32_t read_index, const SequenceBatch &reference,
    const std::vector<Candidate> &candidates,
    std::vector<std::pair<int, uint64_t>> *mappings, int *min_num_errors,
    int *num_best_mappings, int *second_min_num_errors,
    int *num_second_best_mappings) {
  const char *read = read_batch.GetSequenceAt(read_index);
  uint32_t read_length = read_batch.GetSequenceLengthAt(read_index);
  const std::string &negative_read =
      read_batch.GetNegativeSequenceAt(read_index);

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
    if (position < (uint32_t)error_threshold_ ||
        position >= reference.GetSequenceLengthAt(rid) ||
        position + read_length + error_threshold_ >=
            reference.GetSequenceLengthAt(rid)) {
      // not a valid candidate
      ++candidate_index;
      continue;
    } else {
      valid_candidates[valid_candidate_index] =
          candidates[candidate_index];  // reference.GetSequenceAt(rid) +
                                        // position - error_threshold_;
      valid_candidate_starts[valid_candidate_index] =
          reference.GetSequenceAt(rid) + position - error_threshold_;
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
          BandedAlign8PatternsToText(valid_candidate_starts, read, read_length,
                                     mapping_edit_distances,
                                     mapping_end_positions);
        } else {
          BandedAlign8PatternsToText(
              valid_candidate_starts, negative_read.data(), read_length,
              mapping_edit_distances, mapping_end_positions);
        }
        for (int mi = 0; mi < NUM_VPU_LANES_; ++mi) {
          if (mapping_edit_distances[mi] <= error_threshold_) {
            if (mapping_edit_distances[mi] < *min_num_errors) {
              *second_min_num_errors = *min_num_errors;
              *num_second_best_mappings = *num_best_mappings;
              *min_num_errors = mapping_edit_distances[mi];
              *num_best_mappings = 1;
            } else if (mapping_edit_distances[mi] == *min_num_errors) {
              (*num_best_mappings)++;
            } else if (mapping_edit_distances[mi] == *second_min_num_errors) {
              (*num_second_best_mappings)++;
            } else if (mapping_edit_distances[mi] < *second_min_num_errors) {
              *num_second_best_mappings = 1;
              *second_min_num_errors = mapping_edit_distances[mi];
            }
            if (candidate_direction == kPositive) {
              mappings->emplace_back((uint8_t)mapping_edit_distances[mi],
                                     valid_candidates[mi].position -
                                         error_threshold_ +
                                         mapping_end_positions[mi]);
            } else {
              mappings->emplace_back((uint8_t)mapping_edit_distances[mi],
                                     valid_candidates[mi].position -
                                         read_length + 1 - error_threshold_ +
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
          BandedAlign4PatternsToText(valid_candidate_starts, read, read_length,
                                     mapping_edit_distances,
                                     mapping_end_positions);
        } else {
          BandedAlign4PatternsToText(
              valid_candidate_starts, negative_read.data(), read_length,
              mapping_edit_distances, mapping_end_positions);
        }
        for (int mi = 0; mi < NUM_VPU_LANES_; ++mi) {
          if (mapping_edit_distances[mi] <= error_threshold_) {
            if (mapping_edit_distances[mi] < *min_num_errors) {
              *second_min_num_errors = *min_num_errors;
              *num_second_best_mappings = *num_best_mappings;
              *min_num_errors = mapping_edit_distances[mi];
              *num_best_mappings = 1;
            } else if (mapping_edit_distances[mi] == *min_num_errors) {
              (*num_best_mappings)++;
            } else if (mapping_edit_distances[mi] == *second_min_num_errors) {
              (*num_second_best_mappings)++;
            } else if (mapping_edit_distances[mi] < *second_min_num_errors) {
              *num_second_best_mappings = 1;
              *second_min_num_errors = mapping_edit_distances[mi];
            }
            if (candidate_direction == kPositive) {
              mappings->emplace_back((uint8_t)mapping_edit_distances[mi],
                                     valid_candidates[mi].position -
                                         error_threshold_ +
                                         mapping_end_positions[mi]);
            } else {
              mappings->emplace_back((uint8_t)mapping_edit_distances[mi],
                                     valid_candidates[mi].position -
                                         read_length + 1 - error_threshold_ +
                                         mapping_end_positions[mi]);
            }
          } else {
            candidate_count_threshold = valid_candidates[mi].count;
          }
        }
      }
      valid_candidate_index = 0;
      // Check whether we should stop early. Assuming the candidates are sorted
      // if (GetMAPQForSingleEndRead(error_threshold_, num_candidates, 0,
      // read_length + error_threshold_, *min_num_errors, *num_best_mappings,
      // *second_min_num_errors, *num_second_best_mappings) == 0 &&
      // candidate_count_threshold + 1 < candidates[candidate_index].count)
      //  candidate_count_threshold = candidates[candidate_index].count - 1 ;
    }
    ++candidate_index;
  }

  for (uint32_t ci = 0; ci < valid_candidate_index; ++ci) {
    uint32_t rid = valid_candidates[ci].position >> 32;
    uint32_t position = valid_candidates[ci].position;
    if (candidate_direction == kNegative) {
      position = position - read_length + 1;
    }
    if (position < (uint32_t)error_threshold_ ||
        position >= reference.GetSequenceLengthAt(rid) ||
        position + read_length + error_threshold_ >=
            reference.GetSequenceLengthAt(rid)) {
      continue;
    }
    int mapping_end_position;
    int num_errors;
    if (candidate_direction == kPositive) {
      num_errors = BandedAlignPatternToText(
          reference.GetSequenceAt(rid) + position - error_threshold_, read,
          read_length, &mapping_end_position);
    } else {
      num_errors = BandedAlignPatternToText(
          reference.GetSequenceAt(rid) + position - error_threshold_,
          negative_read.data(), read_length, &mapping_end_position);
    }
    if (num_errors <= error_threshold_) {
      if (num_errors < *min_num_errors) {
        *second_min_num_errors = *min_num_errors;
        *num_second_best_mappings = *num_best_mappings;
        *min_num_errors = num_errors;
        *num_best_mappings = 1;
      } else if (num_errors == *min_num_errors) {
        (*num_best_mappings)++;
      } else if (num_errors == *second_min_num_errors) {
        (*num_second_best_mappings)++;
      } else if (num_errors < *second_min_num_errors) {
        *num_second_best_mappings = 1;
        *second_min_num_errors = num_errors;
      }
      if (candidate_direction == kPositive) {
        mappings->emplace_back(num_errors, valid_candidates[ci].position -
                                               error_threshold_ +
                                               mapping_end_position);
      } else {
        mappings->emplace_back(num_errors,
                               valid_candidates[ci].position - read_length + 1 -
                                   error_threshold_ + mapping_end_position);
      }
    }
  }
}

template <typename MappingRecord>
void Chromap<MappingRecord>::VerifyCandidatesOnOneDirection(
    Direction candidate_direction, const SequenceBatch &read_batch,
    uint32_t read_index, const SequenceBatch &reference,
    const std::vector<Candidate> &candidates,
    std::vector<std::pair<int, uint64_t>> *mappings,
    std::vector<int> *split_sites, int *min_num_errors, int *num_best_mappings,
    int *second_min_num_errors, int *num_second_best_mappings) {
  const char *read = read_batch.GetSequenceAt(read_index);
  uint32_t read_length = read_batch.GetSequenceLengthAt(read_index);
  const std::string &negative_read =
      read_batch.GetNegativeSequenceAt(read_index);
  uint32_t candidate_count_threshold = 0;

  for (uint32_t ci = 0; ci < candidates.size(); ++ci) {
    if (candidates[ci].count < candidate_count_threshold) break;
    uint32_t rid = candidates[ci].position >> 32;
    uint32_t position = candidates[ci].position;
    if (candidate_direction == kNegative) {
      position = position - read_length + 1;
    }
    if (position < (uint32_t)error_threshold_ ||
        position >= reference.GetSequenceLengthAt(rid) ||
        position + read_length + error_threshold_ >=
            reference.GetSequenceLengthAt(rid)) {
      continue;
    }
    int mapping_end_position = read_length;
    int gap_beginning = 0;
    int num_errors;
    int allow_gap_beginning_ = 20;
    int mapping_length_threshold = 30;
    int allow_gap_beginning = allow_gap_beginning_ - error_threshold_;
    int actual_num_errors = 0;
    int read_mapping_length = 0;
    int best_mapping_longest_match = 0;
    int longest_match = 0;

    if (split_alignment_) {
      if (candidate_direction == kPositive) {
        num_errors = BandedAlignPatternToTextWithDropOff(
            reference.GetSequenceAt(rid) + position - error_threshold_, read,
            read_length, &mapping_end_position, &read_mapping_length);
        if (mapping_end_position < 0 && allow_gap_beginning > 0) {
          int backup_num_errors = num_errors;
          int backup_mapping_end_position = -mapping_end_position;
          int backup_read_mapping_length = read_mapping_length;
          num_errors = BandedAlignPatternToTextWithDropOff(
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
            mapping_end_position +=
                gap_beginning;  // realign the mapping end position as it is the
                                // alignment from the whole read
            // I use this adjustment since "position" is based on the whole
            // read, and it will be more consistent with no gap beginning case
            read_mapping_length += gap_beginning;
          }
        }
      } else {
        num_errors = BandedAlignPatternToTextWithDropOffFrom3End(
            reference.GetSequenceAt(rid) + position - error_threshold_,
            negative_read.data(), read_length, &mapping_end_position,
            &read_mapping_length);
        if (mapping_end_position < 0 && allow_gap_beginning > 0) {
          int backup_num_errors = num_errors;
          int backup_mapping_end_position = -mapping_end_position;
          int backup_read_mapping_length = read_mapping_length;
          num_errors = BandedAlignPatternToTextWithDropOffFrom3End(
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
      // std::cerr << "ne1: " << num_errors << " " << mapping_end_position <<
      // "\n"; if (num_errors > 2 * error_threshold_) {
      //  if (mapping_end_position - error_threshold_ - num_errors >=
      //  mapping_length_threshold) {
      //    mapping_end_position -= num_errors;
      //    num_errors = -(mapping_end_position - error_threshold_);
      //  }
      //} else {
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
      //}
      // std::cerr << "ne2: " << num_errors << " " << mapping_end_position <<
      // "\n";
    } else {
      if (candidate_direction == kPositive) {
        num_errors = BandedAlignPatternToText(
            reference.GetSequenceAt(rid) + position - error_threshold_, read,
            read_length, &mapping_end_position);
      } else {
        num_errors = BandedAlignPatternToText(
            reference.GetSequenceAt(rid) + position - error_threshold_,
            negative_read.data(), read_length, &mapping_end_position);
      }
    }

    // std::cerr << "ne3: " << num_errors << " " << mapping_end_position << " "
    // << actual_num_errors << " "<< reference.GetSequenceNameAt(rid) <<" " <<
    // (int)candidates[ci].position << " " << position<<"\n";
    if (num_errors <= error_threshold_) {
      if (num_errors < *min_num_errors) {
        *second_min_num_errors = *min_num_errors;
        *num_second_best_mappings = *num_best_mappings;
        *min_num_errors = num_errors;
        *num_best_mappings = 1;
        if (split_alignment_) {
          if (candidates.size() > 50) {
            candidate_count_threshold = candidates[ci].count;
          } else {
            candidate_count_threshold = candidates[ci].count / 2;
          }
          if (*second_min_num_errors < *min_num_errors + error_threshold_ / 2 &&
              best_mapping_longest_match > longest_match &&
              candidates.size() > 200) {
            *second_min_num_errors = *min_num_errors;
          }
        }
        best_mapping_longest_match = longest_match;
      } else if (num_errors == *min_num_errors) {
        (*num_best_mappings)++;
        /*if (split_alignment_ && candidates.size() > 50) {
                candidate_count_threshold = candidates[ci].count + 1;
        }*/
      } else if (num_errors == *second_min_num_errors) {
        (*num_second_best_mappings)++;
      } else if (num_errors < *second_min_num_errors) {
        *num_second_best_mappings = 1;
        *second_min_num_errors = num_errors;
      }
      if (candidate_direction == kPositive) {
        mappings->emplace_back(
            num_errors,
            candidates[ci].position - error_threshold_ + mapping_end_position);
      } else {
        if (split_alignment_ && mapping_output_format_ != MAPPINGFORMAT_SAM) {
          // mappings->emplace_back(num_errors, candidates[ci].position +
          // error_threshold_ - 1 - mapping_end_position
          //					+ read_mapping_length - 1 -
          // gap_beginning);
          mappings->emplace_back(num_errors,
                                 candidates[ci].position - gap_beginning);
        } else {
          // Need to minus gap_beginning because mapping_end_position is
          // adjusted by it, but read_length is not.
          // printf("%d %d %d\n", candidates[ci].position, mapping_end_position,
          // gap_beginning);
          mappings->emplace_back(num_errors,
                                 candidates[ci].position - read_length + 1 -
                                     error_threshold_ + mapping_end_position);
        }
      }
      if (split_alignment_) {
        /*if (mapping_end_position - error_threshold_ < 0 ||
           mapping_end_position - error_threshold_ > 200 || mapping_end_position
           - error_threshold_ < 20) { printf("ERROR! %d %d %d %d %d\n",
           mapping_end_position, error_threshold_,
           read_length,(int)candidates[ci].position, gap_beginning);
                }*/
        /*if (num_errors < *min_num_errors + error_threshold_ / 2 && num_errors
        > *min_num_errors
                        && longest_match > best_mapping_longest_match &&
        candidates.size() > 200) {
                (*num_second_best_mappings)++;
                *second_min_num_errors = *min_num_errors;
        }*/
        split_sites->emplace_back(((actual_num_errors & 0xff) << 24) |
                                  ((gap_beginning & 0xff) << 16) |
                                  (read_mapping_length & 0xffff));
      }
    }
  }
}

template <typename MappingRecord>
void Chromap<MappingRecord>::VerifyCandidates(
    const SequenceBatch &read_batch, uint32_t read_index,
    const SequenceBatch &reference,
    const std::vector<std::pair<uint64_t, uint64_t>> &minimizers,
    const std::vector<Candidate> &positive_candidates,
    const std::vector<Candidate> &negative_candidates,
    std::vector<std::pair<int, uint64_t>> *positive_mappings,
    std::vector<int> *positive_split_sites,
    std::vector<std::pair<int, uint64_t>> *negative_mappings,
    std::vector<int> *negative_split_sites, int *min_num_errors,
    int *num_best_mappings, int *second_min_num_errors,
    int *num_second_best_mappings) {
  *min_num_errors = error_threshold_ + 1;
  *num_best_mappings = 0;
  *second_min_num_errors = error_threshold_ + 1;
  *num_second_best_mappings = 0;

  if (!split_alignment_) {
    // Directly obtain the mapping in ideal case.
    uint32_t i;
    int maxCnt = 0;
    int maxTag = 0;
    // printf("LI_DEBUG: %u %u\n", positive_candidates.size() +
    // negative_candidates.size(), minimizers.size()) ;
    for (i = 0; i < positive_candidates.size(); ++i) {
#ifdef LI_DEBUG
      printf("%s + %u %u %d:%d\n", __func__, i, positive_candidates[i].count,
             (int)(positive_candidates[i].position >> 32),
             (int)positive_candidates[i].position);
#endif
      if (positive_candidates[i].count == minimizers.size()) {
        maxTag = i << 1;
        ++maxCnt;
      }
    }
    for (i = 0; i < negative_candidates.size(); ++i) {
#ifdef LI_DEBUG
      printf("%s - %u %u %d:%d\n", __func__, i, negative_candidates[i].count,
             (int)(negative_candidates[i].position >> 32),
             (int)negative_candidates[i].position);
#endif
      if (negative_candidates[i].count == minimizers.size()) {
        maxTag = (i << 1) | 1;
        ++maxCnt;
      }
    }
    if (maxCnt == 1 &&
        positive_candidates.size() + negative_candidates.size() == 1) {
      Direction candidate_direction = (maxTag & 1) ? kNegative : kPositive;
      uint32_t ci = maxTag >> 1;
      *num_best_mappings = 1;
      *num_second_best_mappings = 0;
      *min_num_errors = 0;

      uint32_t rid = 0;
      uint32_t position = 0;
      uint32_t read_length = read_batch.GetSequenceLengthAt(read_index);
      if (candidate_direction == kPositive) {
        rid = positive_candidates[ci].position >> 32;
        position = positive_candidates[ci].position;
      } else {
        rid = negative_candidates[ci].position >> 32;
        position = (uint32_t)negative_candidates[ci].position - read_length + 1;
      }
      bool flag = true;
      if (position < (uint32_t)error_threshold_ ||
          position >= reference.GetSequenceLengthAt(rid) ||
          position + read_length + error_threshold_ >=
              reference.GetSequenceLengthAt(rid)) {
        flag = false;
      }
      if (flag) {
        if (candidate_direction == kPositive) {
          positive_mappings->emplace_back(
              0, positive_candidates[ci].position + read_length - 1);
        } else {
          negative_mappings->emplace_back(0, negative_candidates[ci].position);
        }
        // fprintf(stderr, "Saved %d\n", positive_candidates.size() +
        // negative_candidates.size() ) ;
        return;
      }
    }
  }
  // printf("Notsaved %d\n", positive_candidates.size() +
  // negative_candidates.size()) ;

  // Use more sophicated approach to obtain the mapping
  if (split_alignment_) {
    std::vector<Candidate> sorted_candidates(positive_candidates);
    std::sort(sorted_candidates.begin(), sorted_candidates.end());
    VerifyCandidatesOnOneDirection(
        kPositive, read_batch, read_index, reference, sorted_candidates,
        positive_mappings, positive_split_sites, min_num_errors,
        num_best_mappings, second_min_num_errors, num_second_best_mappings);

    sorted_candidates = negative_candidates;
    std::sort(sorted_candidates.begin(), sorted_candidates.end());
    VerifyCandidatesOnOneDirection(
        kNegative, read_batch, read_index, reference, sorted_candidates,
        negative_mappings, negative_split_sites, min_num_errors,
        num_best_mappings, second_min_num_errors, num_second_best_mappings);
  } else {
    if (positive_candidates.size() < (size_t)NUM_VPU_LANES_) {
      VerifyCandidatesOnOneDirection(
          kPositive, read_batch, read_index, reference, positive_candidates,
          positive_mappings, positive_split_sites, min_num_errors,
          num_best_mappings, second_min_num_errors, num_second_best_mappings);
    } else {
      std::vector<Candidate> sorted_candidates(positive_candidates);
      std::sort(sorted_candidates.begin(), sorted_candidates.end());
      // std::cerr << "best: " << sorted_candidates[0].count << " " << "second
      // best: " << sorted_candidates[1].count << "\n";
      VerifyCandidatesOnOneDirectionUsingSIMD(
          kPositive, read_batch, read_index, reference, sorted_candidates,
          positive_mappings, min_num_errors, num_best_mappings,
          second_min_num_errors, num_second_best_mappings);
      // VerifyCandidatesOnOneDirectionUsingSIMD(kPositive, read_batch,
      // read_index, reference, positive_candidates, positive_mappings,
      // min_num_errors, num_best_mappings, second_min_num_errors,
      // num_second_best_mappings);
    }
    if (negative_candidates.size() < (size_t)NUM_VPU_LANES_) {
      VerifyCandidatesOnOneDirection(
          kNegative, read_batch, read_index, reference, negative_candidates,
          negative_mappings, negative_split_sites, min_num_errors,
          num_best_mappings, second_min_num_errors, num_second_best_mappings);
    } else {
      std::vector<Candidate> sorted_candidates(negative_candidates);
      std::sort(sorted_candidates.begin(), sorted_candidates.end());
      // std::cerr << "best: " << sorted_candidates[0].count << " " << "second
      // best: " << sorted_candidates[1].count << "\n";
      VerifyCandidatesOnOneDirectionUsingSIMD(
          kNegative, read_batch, read_index, reference, sorted_candidates,
          negative_mappings, min_num_errors, num_best_mappings,
          second_min_num_errors, num_second_best_mappings);
      // VerifyCandidatesOnOneDirectionUsingSIMD(kNegative, read_batch,
      // read_index, reference, negative_candidates, negative_mappings,
      // min_num_errors, num_best_mappings, second_min_num_errors,
      // num_second_best_mappings);
    }
  }
}

template <typename MappingRecord>
int Chromap<MappingRecord>::GetLongestMatchLength(const char *pattern,
                                                  const char *text,
                                                  const int read_length) {
  int max_match = 0;
  int tmp = 0;
  for (int i = 0; i < read_length; ++i) {
    if (SequenceBatch::CharToUint8(pattern[i]) ==
        SequenceBatch::CharToUint8(text[i])) {
      ++tmp;
    } else if (tmp > max_match) {
      max_match = tmp;
    }
  }
  if (tmp > max_match) {
    max_match = tmp;
  }
  return max_match;
}

template <typename MappingRecord>
int Chromap<MappingRecord>::BandedAlignPatternToText(
    const char *pattern, const char *text, const int read_length,
    int *mapping_end_position) {
  // int error_count = 0;
  // for (int i = 0; i < read_length; ++i) {
  //  if (pattern[i + error_threshold_] != text[i]) {
  //    ++error_count;
  //    if (error_count > 1) {
  //      break;
  //    }
  //  }
  //}
  // if (error_count <= 1) {
  //  *mapping_end_position = read_length - 1 + error_threshold_;
  //  return error_count;
  //}
  uint32_t Peq[5] = {0, 0, 0, 0, 0};
  for (int i = 0; i < 2 * error_threshold_; i++) {
    uint8_t base = SequenceBatch::CharToUint8(pattern[i]);
    Peq[base] = Peq[base] | (1 << i);
  }
  uint32_t highest_bit_in_band_mask = 1 << (2 * error_threshold_);
  uint32_t lowest_bit_in_band_mask = 1;
  uint32_t VP = 0;
  uint32_t VN = 0;
  uint32_t X = 0;
  uint32_t D0 = 0;
  uint32_t HN = 0;
  uint32_t HP = 0;
  int num_errors_at_band_start_position = 0;
  for (int i = 0; i < read_length; i++) {
    uint8_t pattern_base =
        SequenceBatch::CharToUint8(pattern[i + 2 * error_threshold_]);
    Peq[pattern_base] = Peq[pattern_base] | highest_bit_in_band_mask;
    X = Peq[SequenceBatch::CharToUint8(text[i])] | VN;
    D0 = ((VP + (X & VP)) ^ VP) | X;
    HN = VP & D0;
    HP = VN | ~(VP | D0);
    X = D0 >> 1;
    VN = X & HP;
    VP = HN | ~(X | HP);
    num_errors_at_band_start_position += 1 - (D0 & lowest_bit_in_band_mask);
    if (num_errors_at_band_start_position > 3 * error_threshold_) {
      return error_threshold_ + 1;
    }
    for (int ai = 0; ai < 5; ai++) {
      Peq[ai] >>= 1;
    }
  }
  int band_start_position = read_length - 1;
  int min_num_errors = num_errors_at_band_start_position;
  *mapping_end_position = band_start_position;
  for (int i = 0; i < 2 * error_threshold_; i++) {
    num_errors_at_band_start_position =
        num_errors_at_band_start_position + ((VP >> i) & (uint32_t)1);
    num_errors_at_band_start_position =
        num_errors_at_band_start_position - ((VN >> i) & (uint32_t)1);
    if (num_errors_at_band_start_position < min_num_errors ||
        (num_errors_at_band_start_position == min_num_errors &&
         i + 1 == error_threshold_)) {
      min_num_errors = num_errors_at_band_start_position;
      *mapping_end_position = band_start_position + 1 + i;
    }
  }
  return min_num_errors;
}

// Return negative number if the termination are deemed at the beginning of the
// read mappping_end_position is relative to pattern (reference)
// read_mapping_length is for text (read)
template <typename MappingRecord>
int Chromap<MappingRecord>::BandedAlignPatternToTextWithDropOff(
    const char *pattern, const char *text, const int read_length,
    int *mapping_end_position, int *read_mapping_length) {
  uint32_t Peq[5] = {0, 0, 0, 0, 0};
  for (int i = 0; i < 2 * error_threshold_; i++) {
    uint8_t base = SequenceBatch::CharToUint8(pattern[i]);
    Peq[base] = Peq[base] | (1 << i);
  }
  uint32_t highest_bit_in_band_mask = 1 << (2 * error_threshold_);
  uint32_t lowest_bit_in_band_mask = 1;
  uint32_t VP = 0;
  uint32_t VN = 0;
  uint32_t X = 0;
  uint32_t D0 = 0;
  uint32_t HN = 0;
  uint32_t HP = 0;
  uint32_t prev_VP = 0;
  uint32_t prev_VN = 0;
  int num_errors_at_band_start_position = 0;
  int i = 0;
  int fail_beginning = 0;  // the alignment failed at the beginning part
  int prev_num_errors_at_band_start_position = 0;
  for (; i < read_length; i++) {
    uint8_t pattern_base =
        SequenceBatch::CharToUint8(pattern[i + 2 * error_threshold_]);
    Peq[pattern_base] = Peq[pattern_base] | highest_bit_in_band_mask;
    X = Peq[SequenceBatch::CharToUint8(text[i])] | VN;
    D0 = ((VP + (X & VP)) ^ VP) | X;
    HN = VP & D0;
    HP = VN | ~(VP | D0);
    X = D0 >> 1;
    prev_VN = VN;
    prev_VP = VP;
    VN = X & HP;
    VP = HN | ~(X | HP);
    prev_num_errors_at_band_start_position = num_errors_at_band_start_position;
    num_errors_at_band_start_position += 1 - (D0 & lowest_bit_in_band_mask);
    if (num_errors_at_band_start_position > 2 * error_threshold_) {
      // return error_threshold_ + 1;
      // the min error in this band could be still less than the
      // error_threshold, and could but this should be fine since it does not
      // affect the 5' end of the read.
      if (i < 4 * error_threshold_ && i < read_length / 2) {
        fail_beginning = 1;
      }
      break;
    }
    for (int ai = 0; ai < 5; ai++) {
      Peq[ai] >>= 1;
    }
  }

  /*char tmp[255] ;
  strncpy(tmp, pattern, read_length + 2 * error_threshold_);
  printf("%s\n%s\n", tmp, text);
  printf("%d\n", i) ;
  fflush(stdout);*/
  if (i < read_length) {
    num_errors_at_band_start_position = prev_num_errors_at_band_start_position;
    VN = prev_VN;
    VP = prev_VP;
  }
  int band_start_position = i - 1;
  int min_num_errors = num_errors_at_band_start_position;
  *read_mapping_length = i;
  *mapping_end_position = band_start_position;

  for (i = 0; i < 2 * error_threshold_; i++) {
    num_errors_at_band_start_position =
        num_errors_at_band_start_position + ((VP >> i) & (uint32_t)1);
    num_errors_at_band_start_position =
        num_errors_at_band_start_position - ((VN >> i) & (uint32_t)1);
    if (num_errors_at_band_start_position < min_num_errors ||
        (num_errors_at_band_start_position == min_num_errors &&
         i + 1 == error_threshold_)) {
      min_num_errors = num_errors_at_band_start_position;
      *mapping_end_position = band_start_position + 1 + i;
    }
  }
  if (fail_beginning ||
      (read_length > 60 &&
       *mapping_end_position + 1 - error_threshold_ - min_num_errors < 30)) {
    *mapping_end_position = -*mapping_end_position;
  }
  return min_num_errors;
}

template <typename MappingRecord>
int Chromap<MappingRecord>::BandedAlignPatternToTextWithDropOffFrom3End(
    const char *pattern, const char *text, const int read_length,
    int *mapping_end_position, int *read_mapping_length) {
  uint32_t Peq[5] = {0, 0, 0, 0, 0};
  for (int i = 0; i < 2 * error_threshold_; i++) {
    uint8_t base = SequenceBatch::CharToUint8(
        pattern[read_length + 2 * error_threshold_ - 1 - i]);
    Peq[base] = Peq[base] | (1 << i);
  }
  uint32_t highest_bit_in_band_mask = 1 << (2 * error_threshold_);
  uint32_t lowest_bit_in_band_mask = 1;
  uint32_t VP = 0;
  uint32_t VN = 0;
  uint32_t X = 0;
  uint32_t D0 = 0;
  uint32_t HN = 0;
  uint32_t HP = 0;
  uint32_t prev_VP = 0;
  uint32_t prev_VN = 0;
  int num_errors_at_band_start_position = 0;
  int i = 0;
  int fail_beginning = 0;  // the alignment failed at the beginning part
  int prev_num_errors_at_band_start_position = 0;
  for (; i < read_length; i++) {
    // printf("%c %c %d\n", pattern[read_length - 1 - i], pattern[read_length -
    // 1 - i + error_threshold_], text[read_length - 1 - i]);
    uint8_t pattern_base =
        SequenceBatch::CharToUint8(pattern[read_length - 1 - i]);
    Peq[pattern_base] = Peq[pattern_base] | highest_bit_in_band_mask;
    X = Peq[SequenceBatch::CharToUint8(text[read_length - 1 - i])] | VN;
    D0 = ((VP + (X & VP)) ^ VP) | X;
    HN = VP & D0;
    HP = VN | ~(VP | D0);
    X = D0 >> 1;
    prev_VN = VN;
    prev_VP = VP;
    VN = X & HP;
    VP = HN | ~(X | HP);
    prev_num_errors_at_band_start_position = num_errors_at_band_start_position;
    num_errors_at_band_start_position += 1 - (D0 & lowest_bit_in_band_mask);
    /*printf("->%d %d %c %c", i, num_errors_at_band_start_position,
    pattern[read_length - 1 - i], text[read_length - 1 - i]) ; int tmp =
    num_errors_at_band_start_position; for (int j = 0; j < 2 * error_threshold_;
    j++) { tmp = tmp + ((VP >> j) & (uint32_t) 1); tmp = tmp - ((VN >> j) &
    (uint32_t) 1); printf(" %d", tmp);
    }
    printf("\n");*/
    if (num_errors_at_band_start_position > 2 * error_threshold_) {
      // return error_threshold_ + 1;
      if (i < 4 * error_threshold_ && i < read_length / 2) {
        fail_beginning = 1;
      }
      break;
    }
    for (int ai = 0; ai < 5; ai++) {
      Peq[ai] >>= 1;
    }
  }
  // printf("li %d: %d %d %d\n", fail_beginning, i, error_threshold_,
  // read_length);
  if (i < read_length) {
    num_errors_at_band_start_position = prev_num_errors_at_band_start_position;
    VN = prev_VN;
    VP = prev_VP;
  }
  int band_start_position = i - 1;
  int min_num_errors = num_errors_at_band_start_position;
  *read_mapping_length = i;
  *mapping_end_position = band_start_position;
  // printf("-1: %d\n", num_errors_at_band_start_position);
  for (i = 0; i < 2 * error_threshold_; i++) {
    num_errors_at_band_start_position =
        num_errors_at_band_start_position + ((VP >> i) & (uint32_t)1);
    num_errors_at_band_start_position =
        num_errors_at_band_start_position - ((VN >> i) & (uint32_t)1);
    // printf("%d: %d\n", i, num_errors_at_band_start_position);
    if (num_errors_at_band_start_position < min_num_errors ||
        (num_errors_at_band_start_position == min_num_errors &&
         i + 1 == error_threshold_)) {
      min_num_errors = num_errors_at_band_start_position;
      *mapping_end_position = band_start_position + (1 + i);
    }
  }
  if (fail_beginning ||
      (read_length > 60 &&
       *mapping_end_position + 1 - error_threshold_ - min_num_errors < 30)) {
    *mapping_end_position = -*mapping_end_position;
  }
  return min_num_errors;
}

template <typename MappingRecord>
void Chromap<MappingRecord>::BandedAlign4PatternsToText(
    const char **patterns, const char *text, int read_length,
    int32_t *mapping_edit_distances, int32_t *mapping_end_positions) {
  int ALPHABET_SIZE = 5;
  const char *reference_sequence0 = patterns[0];
  const char *reference_sequence1 = patterns[1];
  const char *reference_sequence2 = patterns[2];
  const char *reference_sequence3 = patterns[3];
  uint32_t highest_bit_in_band_mask = 1 << (2 * error_threshold_);
  __m128i highest_bit_in_band_mask_vpu0 =
      _mm_set_epi32(0, 0, 0, highest_bit_in_band_mask);
  __m128i highest_bit_in_band_mask_vpu1 =
      _mm_set_epi32(0, 0, highest_bit_in_band_mask, 0);
  __m128i highest_bit_in_band_mask_vpu2 =
      _mm_set_epi32(0, highest_bit_in_band_mask, 0, 0);
  __m128i highest_bit_in_band_mask_vpu3 =
      _mm_set_epi32(highest_bit_in_band_mask, 0, 0, 0);
  // Init Peq
  __m128i Peq[ALPHABET_SIZE];
  for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
    Peq[ai] = _mm_setzero_si128();
  }
  for (int i = 0; i < 2 * error_threshold_; i++) {
    uint8_t base0 = SequenceBatch::CharToUint8(reference_sequence0[i]);
    uint8_t base1 = SequenceBatch::CharToUint8(reference_sequence1[i]);
    uint8_t base2 = SequenceBatch::CharToUint8(reference_sequence2[i]);
    uint8_t base3 = SequenceBatch::CharToUint8(reference_sequence3[i]);
    Peq[base0] = _mm_or_si128(highest_bit_in_band_mask_vpu0, Peq[base0]);
    Peq[base1] = _mm_or_si128(highest_bit_in_band_mask_vpu1, Peq[base1]);
    Peq[base2] = _mm_or_si128(highest_bit_in_band_mask_vpu2, Peq[base2]);
    Peq[base3] = _mm_or_si128(highest_bit_in_band_mask_vpu3, Peq[base3]);
    for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
      Peq[ai] = _mm_srli_epi32(Peq[ai], 1);
    }
  }

  uint32_t lowest_bit_in_band_mask = 1;
  __m128i lowest_bit_in_band_mask_vpu = _mm_set1_epi32(lowest_bit_in_band_mask);
  __m128i VP = _mm_setzero_si128();
  __m128i VN = _mm_setzero_si128();
  __m128i X = _mm_setzero_si128();
  __m128i D0 = _mm_setzero_si128();
  __m128i HN = _mm_setzero_si128();
  __m128i HP = _mm_setzero_si128();
  __m128i max_mask_vpu = _mm_set1_epi32(0xffffffff);
  __m128i num_errors_at_band_start_position_vpu = _mm_setzero_si128();
  __m128i early_stop_threshold_vpu = _mm_set1_epi32(error_threshold_ * 3);
  for (int i = 0; i < read_length; i++) {
    uint8_t base0 = SequenceBatch::CharToUint8(
        reference_sequence0[i + 2 * error_threshold_]);
    uint8_t base1 = SequenceBatch::CharToUint8(
        reference_sequence1[i + 2 * error_threshold_]);
    uint8_t base2 = SequenceBatch::CharToUint8(
        reference_sequence2[i + 2 * error_threshold_]);
    uint8_t base3 = SequenceBatch::CharToUint8(
        reference_sequence3[i + 2 * error_threshold_]);
    Peq[base0] = _mm_or_si128(highest_bit_in_band_mask_vpu0, Peq[base0]);
    Peq[base1] = _mm_or_si128(highest_bit_in_band_mask_vpu1, Peq[base1]);
    Peq[base2] = _mm_or_si128(highest_bit_in_band_mask_vpu2, Peq[base2]);
    Peq[base3] = _mm_or_si128(highest_bit_in_band_mask_vpu3, Peq[base3]);
    X = _mm_or_si128(Peq[SequenceBatch::CharToUint8(text[i])], VN);
    D0 = _mm_and_si128(X, VP);
    D0 = _mm_add_epi32(D0, VP);
    D0 = _mm_xor_si128(D0, VP);
    D0 = _mm_or_si128(D0, X);
    HN = _mm_and_si128(VP, D0);
    HP = _mm_or_si128(VP, D0);
    HP = _mm_xor_si128(HP, max_mask_vpu);
    HP = _mm_or_si128(HP, VN);
    X = _mm_srli_epi32(D0, 1);
    VN = _mm_and_si128(X, HP);
    VP = _mm_or_si128(X, HP);
    VP = _mm_xor_si128(VP, max_mask_vpu);
    VP = _mm_or_si128(VP, HN);
    __m128i E = _mm_and_si128(D0, lowest_bit_in_band_mask_vpu);
    E = _mm_xor_si128(E, lowest_bit_in_band_mask_vpu);
    num_errors_at_band_start_position_vpu =
        _mm_add_epi32(num_errors_at_band_start_position_vpu, E);
    __m128i early_stop = _mm_cmpgt_epi32(num_errors_at_band_start_position_vpu,
                                         early_stop_threshold_vpu);
    int tmp = _mm_movemask_epi8(early_stop);
    if (tmp == 0xffff) {
      _mm_store_si128((__m128i *)mapping_edit_distances,
                      num_errors_at_band_start_position_vpu);
      return;
    }
    for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
      Peq[ai] = _mm_srli_epi32(Peq[ai], 1);
    }
  }
  int band_start_position = read_length - 1;
  __m128i min_num_errors_vpu = num_errors_at_band_start_position_vpu;
  for (int i = 0; i < 2 * error_threshold_; i++) {
    __m128i lowest_bit_in_VP_vpu =
        _mm_and_si128(VP, lowest_bit_in_band_mask_vpu);
    __m128i lowest_bit_in_VN_vpu =
        _mm_and_si128(VN, lowest_bit_in_band_mask_vpu);
    num_errors_at_band_start_position_vpu = _mm_add_epi32(
        num_errors_at_band_start_position_vpu, lowest_bit_in_VP_vpu);
    num_errors_at_band_start_position_vpu = _mm_sub_epi32(
        num_errors_at_band_start_position_vpu, lowest_bit_in_VN_vpu);
    __m128i mapping_end_positions_update_mask_vpu = _mm_cmplt_epi32(
        num_errors_at_band_start_position_vpu, min_num_errors_vpu);
    __m128i mapping_end_positions_update_mask_vpu1 = _mm_cmpeq_epi32(
        num_errors_at_band_start_position_vpu, min_num_errors_vpu);
    int mapping_end_positions_update_mask =
        _mm_movemask_epi8(mapping_end_positions_update_mask_vpu);
    int mapping_end_positions_update_mask1 =
        _mm_movemask_epi8(mapping_end_positions_update_mask_vpu1);
    for (int li = 0; li < 4; ++li) {
      if ((mapping_end_positions_update_mask & 1) == 1 ||
          ((mapping_end_positions_update_mask1 & 1) == 1 &&
           i + 1 == error_threshold_)) {
        mapping_end_positions[li] = band_start_position + 1 + i;
      }
      mapping_end_positions_update_mask =
          mapping_end_positions_update_mask >> 4;
      mapping_end_positions_update_mask1 =
          mapping_end_positions_update_mask1 >> 4;
    }
    min_num_errors_vpu = _mm_min_epi32(min_num_errors_vpu,
                                       num_errors_at_band_start_position_vpu);
    VP = _mm_srli_epi32(VP, 1);
    VN = _mm_srli_epi32(VN, 1);
  }
  _mm_store_si128((__m128i *)mapping_edit_distances, min_num_errors_vpu);
}

template <typename MappingRecord>
void Chromap<MappingRecord>::BandedAlign8PatternsToText(
    const char **patterns, const char *text, int read_length,
    int16_t *mapping_edit_distances, int16_t *mapping_end_positions) {
  int ALPHABET_SIZE = 5;
  const char *reference_sequence0 = patterns[0];
  const char *reference_sequence1 = patterns[1];
  const char *reference_sequence2 = patterns[2];
  const char *reference_sequence3 = patterns[3];
  const char *reference_sequence4 = patterns[4];
  const char *reference_sequence5 = patterns[5];
  const char *reference_sequence6 = patterns[6];
  const char *reference_sequence7 = patterns[7];
  uint16_t highest_bit_in_band_mask = 1 << (2 * error_threshold_);
  __m128i highest_bit_in_band_mask_vpu0 =
      _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, highest_bit_in_band_mask);
  __m128i highest_bit_in_band_mask_vpu1 =
      _mm_set_epi16(0, 0, 0, 0, 0, 0, highest_bit_in_band_mask, 0);
  __m128i highest_bit_in_band_mask_vpu2 =
      _mm_set_epi16(0, 0, 0, 0, 0, highest_bit_in_band_mask, 0, 0);
  __m128i highest_bit_in_band_mask_vpu3 =
      _mm_set_epi16(0, 0, 0, 0, highest_bit_in_band_mask, 0, 0, 0);
  __m128i highest_bit_in_band_mask_vpu4 =
      _mm_set_epi16(0, 0, 0, highest_bit_in_band_mask, 0, 0, 0, 0);
  __m128i highest_bit_in_band_mask_vpu5 =
      _mm_set_epi16(0, 0, highest_bit_in_band_mask, 0, 0, 0, 0, 0);
  __m128i highest_bit_in_band_mask_vpu6 =
      _mm_set_epi16(0, highest_bit_in_band_mask, 0, 0, 0, 0, 0, 0);
  __m128i highest_bit_in_band_mask_vpu7 =
      _mm_set_epi16(highest_bit_in_band_mask, 0, 0, 0, 0, 0, 0, 0);
  // Init Peq
  __m128i Peq[ALPHABET_SIZE];
  for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
    Peq[ai] = _mm_setzero_si128();
  }
  for (int i = 0; i < 2 * error_threshold_; i++) {
    uint8_t base0 = SequenceBatch::CharToUint8(reference_sequence0[i]);
    uint8_t base1 = SequenceBatch::CharToUint8(reference_sequence1[i]);
    uint8_t base2 = SequenceBatch::CharToUint8(reference_sequence2[i]);
    uint8_t base3 = SequenceBatch::CharToUint8(reference_sequence3[i]);
    uint8_t base4 = SequenceBatch::CharToUint8(reference_sequence4[i]);
    uint8_t base5 = SequenceBatch::CharToUint8(reference_sequence5[i]);
    uint8_t base6 = SequenceBatch::CharToUint8(reference_sequence6[i]);
    uint8_t base7 = SequenceBatch::CharToUint8(reference_sequence7[i]);
    Peq[base0] = _mm_or_si128(highest_bit_in_band_mask_vpu0, Peq[base0]);
    Peq[base1] = _mm_or_si128(highest_bit_in_band_mask_vpu1, Peq[base1]);
    Peq[base2] = _mm_or_si128(highest_bit_in_band_mask_vpu2, Peq[base2]);
    Peq[base3] = _mm_or_si128(highest_bit_in_band_mask_vpu3, Peq[base3]);
    Peq[base4] = _mm_or_si128(highest_bit_in_band_mask_vpu4, Peq[base4]);
    Peq[base5] = _mm_or_si128(highest_bit_in_band_mask_vpu5, Peq[base5]);
    Peq[base6] = _mm_or_si128(highest_bit_in_band_mask_vpu6, Peq[base6]);
    Peq[base7] = _mm_or_si128(highest_bit_in_band_mask_vpu7, Peq[base7]);
    for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
      Peq[ai] = _mm_srli_epi16(Peq[ai], 1);
    }
  }

  uint16_t lowest_bit_in_band_mask = 1;
  __m128i lowest_bit_in_band_mask_vpu = _mm_set1_epi16(lowest_bit_in_band_mask);
  __m128i VP = _mm_setzero_si128();
  __m128i VN = _mm_setzero_si128();
  __m128i X = _mm_setzero_si128();
  __m128i D0 = _mm_setzero_si128();
  __m128i HN = _mm_setzero_si128();
  __m128i HP = _mm_setzero_si128();
  __m128i max_mask_vpu = _mm_set1_epi16(0xffff);
  __m128i num_errors_at_band_start_position_vpu = _mm_setzero_si128();
  __m128i early_stop_threshold_vpu = _mm_set1_epi16(error_threshold_ * 3);
  for (int i = 0; i < read_length; i++) {
    uint8_t base0 = SequenceBatch::CharToUint8(
        reference_sequence0[i + 2 * error_threshold_]);
    uint8_t base1 = SequenceBatch::CharToUint8(
        reference_sequence1[i + 2 * error_threshold_]);
    uint8_t base2 = SequenceBatch::CharToUint8(
        reference_sequence2[i + 2 * error_threshold_]);
    uint8_t base3 = SequenceBatch::CharToUint8(
        reference_sequence3[i + 2 * error_threshold_]);
    uint8_t base4 = SequenceBatch::CharToUint8(
        reference_sequence4[i + 2 * error_threshold_]);
    uint8_t base5 = SequenceBatch::CharToUint8(
        reference_sequence5[i + 2 * error_threshold_]);
    uint8_t base6 = SequenceBatch::CharToUint8(
        reference_sequence6[i + 2 * error_threshold_]);
    uint8_t base7 = SequenceBatch::CharToUint8(
        reference_sequence7[i + 2 * error_threshold_]);
    Peq[base0] = _mm_or_si128(highest_bit_in_band_mask_vpu0, Peq[base0]);
    Peq[base1] = _mm_or_si128(highest_bit_in_band_mask_vpu1, Peq[base1]);
    Peq[base2] = _mm_or_si128(highest_bit_in_band_mask_vpu2, Peq[base2]);
    Peq[base3] = _mm_or_si128(highest_bit_in_band_mask_vpu3, Peq[base3]);
    Peq[base4] = _mm_or_si128(highest_bit_in_band_mask_vpu4, Peq[base4]);
    Peq[base5] = _mm_or_si128(highest_bit_in_band_mask_vpu5, Peq[base5]);
    Peq[base6] = _mm_or_si128(highest_bit_in_band_mask_vpu6, Peq[base6]);
    Peq[base7] = _mm_or_si128(highest_bit_in_band_mask_vpu7, Peq[base7]);
    X = _mm_or_si128(Peq[SequenceBatch::CharToUint8(text[i])], VN);
    D0 = _mm_and_si128(X, VP);
    D0 = _mm_add_epi16(D0, VP);
    D0 = _mm_xor_si128(D0, VP);
    D0 = _mm_or_si128(D0, X);
    HN = _mm_and_si128(VP, D0);
    HP = _mm_or_si128(VP, D0);
    HP = _mm_xor_si128(HP, max_mask_vpu);
    HP = _mm_or_si128(HP, VN);
    X = _mm_srli_epi16(D0, 1);
    VN = _mm_and_si128(X, HP);
    VP = _mm_or_si128(X, HP);
    VP = _mm_xor_si128(VP, max_mask_vpu);
    VP = _mm_or_si128(VP, HN);
    __m128i E = _mm_and_si128(D0, lowest_bit_in_band_mask_vpu);
    E = _mm_xor_si128(E, lowest_bit_in_band_mask_vpu);
    num_errors_at_band_start_position_vpu =
        _mm_add_epi16(num_errors_at_band_start_position_vpu, E);
    __m128i early_stop = _mm_cmpgt_epi16(num_errors_at_band_start_position_vpu,
                                         early_stop_threshold_vpu);
    int tmp = _mm_movemask_epi8(early_stop);
    if (tmp == 0xffff) {
      _mm_store_si128((__m128i *)mapping_edit_distances,
                      num_errors_at_band_start_position_vpu);
      return;
    }
    for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
      Peq[ai] = _mm_srli_epi16(Peq[ai], 1);
    }
  }
  int band_start_position = read_length - 1;
  __m128i min_num_errors_vpu = num_errors_at_band_start_position_vpu;
  for (int i = 0; i < 2 * error_threshold_; i++) {
    __m128i lowest_bit_in_VP_vpu =
        _mm_and_si128(VP, lowest_bit_in_band_mask_vpu);
    __m128i lowest_bit_in_VN_vpu =
        _mm_and_si128(VN, lowest_bit_in_band_mask_vpu);
    num_errors_at_band_start_position_vpu = _mm_add_epi16(
        num_errors_at_band_start_position_vpu, lowest_bit_in_VP_vpu);
    num_errors_at_band_start_position_vpu = _mm_sub_epi16(
        num_errors_at_band_start_position_vpu, lowest_bit_in_VN_vpu);
    __m128i mapping_end_positions_update_mask_vpu = _mm_cmplt_epi16(
        num_errors_at_band_start_position_vpu, min_num_errors_vpu);
    __m128i mapping_end_positions_update_mask_vpu1 = _mm_cmpeq_epi16(
        num_errors_at_band_start_position_vpu, min_num_errors_vpu);
    int mapping_end_positions_update_mask =
        _mm_movemask_epi8(mapping_end_positions_update_mask_vpu);
    int mapping_end_positions_update_mask1 =
        _mm_movemask_epi8(mapping_end_positions_update_mask_vpu1);
    for (int li = 0; li < 8; ++li) {
      if ((mapping_end_positions_update_mask & 1) == 1 ||
          ((mapping_end_positions_update_mask1 & 1) == 1 &&
           i + 1 == error_threshold_)) {
        mapping_end_positions[li] = band_start_position + 1 + i;
      }
      mapping_end_positions_update_mask =
          mapping_end_positions_update_mask >> 2;
      mapping_end_positions_update_mask1 =
          mapping_end_positions_update_mask1 >> 2;
    }
    min_num_errors_vpu = _mm_min_epi16(min_num_errors_vpu,
                                       num_errors_at_band_start_position_vpu);
    VP = _mm_srli_epi16(VP, 1);
    VN = _mm_srli_epi16(VN, 1);
  }
  _mm_store_si128((__m128i *)mapping_edit_distances, min_num_errors_vpu);
}

template <typename MappingRecord>
void Chromap<MappingRecord>::BandedTraceback(int min_num_errors,
                                             const char *pattern,
                                             const char *text,
                                             const int read_length,
                                             int *mapping_start_position) {
  // fisrt calculate the hamming distance and see whether it's equal to # errors
  if (min_num_errors == 0) {
    *mapping_start_position = error_threshold_;
    return;
  }
  int error_count = 0;
  for (int i = 0; i < read_length; ++i) {
    if (pattern[i + error_threshold_] != text[i]) {
      ++error_count;
    }
  }
  if (error_count == min_num_errors) {
    *mapping_start_position = error_threshold_;
    return;
  }
  // if not then there are gaps so that we have to traceback with edit distance.
  uint32_t Peq[5] = {0, 0, 0, 0, 0};
  for (int i = 0; i < 2 * error_threshold_; i++) {
    uint8_t base = SequenceBatch::CharToUint8(
        pattern[read_length - 1 + 2 * error_threshold_ - i]);
    Peq[base] = Peq[base] | (1 << i);
  }
  uint32_t highest_bit_in_band_mask = 1 << (2 * error_threshold_);
  uint32_t lowest_bit_in_band_mask = 1;
  uint32_t VP = 0;
  uint32_t VN = 0;
  uint32_t X = 0;
  uint32_t D0 = 0;
  uint32_t HN = 0;
  uint32_t HP = 0;
  int num_errors_at_band_start_position = 0;
  for (int i = 0; i < read_length; i++) {
    uint8_t pattern_base =
        SequenceBatch::CharToUint8(pattern[read_length - 1 - i]);
    Peq[pattern_base] = Peq[pattern_base] | highest_bit_in_band_mask;
    X = Peq[SequenceBatch::CharToUint8(text[read_length - 1 - i])] | VN;
    D0 = ((VP + (X & VP)) ^ VP) | X;
    HN = VP & D0;
    HP = VN | ~(VP | D0);
    X = D0 >> 1;
    VN = X & HP;
    VP = HN | ~(X | HP);
    num_errors_at_band_start_position += 1 - (D0 & lowest_bit_in_band_mask);
    for (int ai = 0; ai < 5; ai++) {
      Peq[ai] >>= 1;
    }
  }
  *mapping_start_position = 2 * error_threshold_;
  for (int i = 0; i < 2 * error_threshold_; i++) {
    num_errors_at_band_start_position =
        num_errors_at_band_start_position + ((VP >> i) & (uint32_t)1);
    num_errors_at_band_start_position =
        num_errors_at_band_start_position - ((VN >> i) & (uint32_t)1);
    if (num_errors_at_band_start_position == min_num_errors) {
      *mapping_start_position = 2 * error_threshold_ - (1 + i);
      if (i + 1 == error_threshold_) {
        return;
      }
    }
  }
}

template <typename MappingRecord>
void Chromap<MappingRecord>::BandedTracebackToEnd(int min_num_errors,
                                                  const char *pattern,
                                                  const char *text,
                                                  const int read_length,
                                                  int *mapping_end_position) {
  // fisrt calculate the hamming distance and see whether it's equal to # errors
  if (min_num_errors == 0) {
    *mapping_end_position = read_length + error_threshold_;
    return;
  }
  int error_count = 0;
  for (int i = 0; i < read_length; ++i) {
    if (pattern[i + error_threshold_] != text[i]) {
      ++error_count;
    }
  }
  if (error_count == min_num_errors) {
    *mapping_end_position = read_length + error_threshold_;
    return;
  }
  // if not then there are gaps so that we have to traceback with edit distance.
  uint32_t Peq[5] = {0, 0, 0, 0, 0};
  for (int i = 0; i < 2 * error_threshold_; i++) {
    uint8_t base = SequenceBatch::CharToUint8(pattern[i]);
    Peq[base] = Peq[base] | (1 << i);
  }
  uint32_t highest_bit_in_band_mask = 1 << (2 * error_threshold_);
  uint32_t lowest_bit_in_band_mask = 1;
  uint32_t VP = 0;
  uint32_t VN = 0;
  uint32_t X = 0;
  uint32_t D0 = 0;
  uint32_t HN = 0;
  uint32_t HP = 0;
  int num_errors_at_band_start_position = 0;
  for (int i = 0; i < read_length; i++) {
    // printf("=>%d %d %c %c\n", i, num_errors_at_band_start_position, pattern[i
    // + 2 * error_threshold_], text[i]) ;
    uint8_t pattern_base =
        SequenceBatch::CharToUint8(pattern[i + 2 * error_threshold_]);
    Peq[pattern_base] = Peq[pattern_base] | highest_bit_in_band_mask;
    X = Peq[SequenceBatch::CharToUint8(text[i])] | VN;
    D0 = ((VP + (X & VP)) ^ VP) | X;
    HN = VP & D0;
    HP = VN | ~(VP | D0);
    X = D0 >> 1;
    VN = X & HP;
    VP = HN | ~(X | HP);
    num_errors_at_band_start_position += 1 - (D0 & lowest_bit_in_band_mask);
    for (int ai = 0; ai < 5; ai++) {
      Peq[ai] >>= 1;
    }
  }
  int band_start_position = read_length;
  *mapping_end_position = band_start_position + 1;
  for (int i = 0; i < 2 * error_threshold_; i++) {
    num_errors_at_band_start_position =
        num_errors_at_band_start_position + ((VP >> i) & (uint32_t)1);
    num_errors_at_band_start_position =
        num_errors_at_band_start_position - ((VN >> i) & (uint32_t)1);
    if (num_errors_at_band_start_position == min_num_errors) {
      *mapping_end_position = band_start_position + (i + 1);
      if (i + 1 == error_threshold_) {
        return;
      }
    }
  }
}

template <typename MappingRecord>
void Chromap<MappingRecord>::OutputBarcodeStatistics() {
  std::cerr << "Number of barcodes in whitelist: " << num_barcode_in_whitelist_
            << ".\n";
  std::cerr << "Number of corrected barcodes: " << num_corrected_barcode_
            << ".\n";
}

template <typename MappingRecord>
void Chromap<MappingRecord>::OutputMappingStatistics() {
  std::cerr << "Number of reads: " << num_reads_ << ".\n";
  // std::cerr << "Number of duplicated reads: " << num_duplicated_reads_ <<
  // ".\n";
  std::cerr << "Number of mapped reads: " << num_mapped_reads_ << ".\n";
  std::cerr << "Number of uniquely mapped reads: " << num_uniquely_mapped_reads_
            << ".\n";
  std::cerr << "Number of reads have multi-mappings: "
            << num_mapped_reads_ - num_uniquely_mapped_reads_ << ".\n";
  std::cerr << "Number of candidates: " << num_candidates_ << ".\n";
  std::cerr << "Number of mappings: " << num_mappings_ << ".\n";
  std::cerr << "Number of uni-mappings: " << num_uniquely_mapped_reads_
            << ".\n";
  std::cerr << "Number of multi-mappings: "
            << num_mappings_ - num_uniquely_mapped_reads_ << ".\n";
}

template <typename MappingRecord>
void Chromap<MappingRecord>::OutputMappingStatistics(
    uint32_t num_reference_sequences,
    const std::vector<std::vector<MappingRecord>> &uni_mappings,
    const std::vector<std::vector<MappingRecord>> &multi_mappings) {
  uint64_t num_uni_mappings = 0;
  uint64_t num_multi_mappings = 0;
  for (auto &uni_mappings_on_one_ref_seq : uni_mappings) {
    for (auto &uni_mapping : uni_mappings_on_one_ref_seq) {
      if ((uni_mapping.is_unique) == 1) {
        ++num_uni_mappings;
      }
    }
  }
  for (auto &multi_mappings_on_one_ref_seq : multi_mappings) {
    for (auto &multi_mapping : multi_mappings_on_one_ref_seq) {
      if ((multi_mapping.is_unique) != 1) {
        ++num_multi_mappings;
      }
    }
  }
  std::cerr << "# uni-mappings: " << num_uni_mappings
            << ", # multi-mappings: " << num_multi_mappings
            << ", total: " << num_uni_mappings + num_multi_mappings << ".\n";
}

template <typename MappingRecord>
void Chromap<MappingRecord>::LoadBarcodeWhitelist() {
  double real_start_time = GetRealTime();
  int num_barcodes = 0;
  std::ifstream barcode_whitelist_file_stream(barcode_whitelist_file_path_);
  std::string barcode_whitelist_file_line;
  // bool first_line = true;
  while (getline(barcode_whitelist_file_stream, barcode_whitelist_file_line)) {
    std::stringstream barcode_whitelist_file_line_string_stream(
        barcode_whitelist_file_line);
    //// skip the header
    // if (barcode_whitelist_file_line[0] == '#' ||
    // barcode_whitelist_file_line.find("kmer") == 0) {
    //  continue;
    //}
    std::string barcode;
    barcode_whitelist_file_line_string_stream >> barcode;
    size_t barcode_length = barcode.length();
    if (num_barcodes == 0) {
      barcode_length_ = barcode_length;
    } else {
      assert(barcode_length == barcode_length_);
    }
    // if (first_line) {
    //  //size_t barcode_length = kmer.length();
    //  // Allocate memory to save pore model parameters
    //  //size_t num_pore_models = 1 << (kmer_size_ * 2);
    //  //pore_models_.assign(num_pore_models, PoreModelParameters());
    //  //first_line = false;
    //}
    // assert(kmer.length() == (size_t)kmer_size_);
    uint64_t barcode_key = SequenceBatch::GenerateSeedFromSequence(
        barcode.data(), barcode_length, 0, barcode_length);
    // PoreModelParameters &pore_model_parameters =
    // pore_models_[kmer_hash_value]; barcode_whitelist_file_line_string_stream
    // >> pore_model_parameters.level_mean >> pore_model_parameters.level_stdv
    // >> pore_model_parameters.sd_mean >> pore_model_parameters.sd_stdv;
    int khash_return_code;
    khiter_t barcode_whitelist_lookup_table_iterator =
        kh_put(k64_seq, barcode_whitelist_lookup_table_, barcode_key,
               &khash_return_code);
    kh_value(barcode_whitelist_lookup_table_,
             barcode_whitelist_lookup_table_iterator) = 0;
    assert(khash_return_code != -1 && khash_return_code != 0);
    ++num_barcodes;
  }
  barcode_whitelist_file_stream.close();
  std::cerr << "Loaded " << num_barcodes << " barcodes in "
            << GetRealTime() - real_start_time << "s.\n";
}

template <typename MappingRecord>
uint8_t Chromap<MappingRecord>::GetMAPQForSingleEndRead(
    int error_threshold, int num_candidates, uint32_t repetitive_seed_length,
    uint16_t alignment_length, int min_num_errors, int num_best_mappings,
    int second_min_num_errors, int num_second_best_mappings,
    int max_num_error_difference, int read_length) {
  // int mapq_coef = 60;
  // int mapq_coef_length = 45;
  int mapq_coef_length = 50;
  // double mapq_coef_fraction = log(mapq_coef_length);
  int mapq_coef_fraction = log(mapq_coef_length);
  if (!split_alignment_) {
    alignment_length =
        alignment_length > read_length ? alignment_length : read_length;
  }
  double alignment_identity = 1 - (double)min_num_errors / alignment_length;
  if (split_alignment_) {
    alignment_identity = (double)(-min_num_errors) / alignment_length;
    if (alignment_identity > 1) alignment_identity = 1;
  }
  int mapq = 0;
  if (num_best_mappings > 1) {
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
    if (second_min_num_errors > min_num_errors + max_num_error_difference) {
      second_min_num_errors = min_num_errors + max_num_error_difference;
    }
    // mapq = (int)(mapq_coef * ((double)(second_min_num_errors -
    // min_num_errors) / second_min_num_errors) + .499); double tmp =
    // alignment_identity * alignment_identity; tmp = tmp * tmp; tmp = tmp *
    // tmp; mapq = alignment_identity < 0.98 ? (int)(mapq * tmp + .499) : mapq;
    double tmp = alignment_length < mapq_coef_length
                     ? 1.0
                     : mapq_coef_fraction / log(alignment_length);
    tmp *= alignment_identity * alignment_identity;
    if (!split_alignment_) {
      // mapq = 6 * 6.02 * (second_min_num_errors - min_num_errors) * tmp * tmp
      // + 0.499 + 10;
      mapq = 5 * 6.02 * (second_min_num_errors - min_num_errors) * tmp * tmp +
             0.499;
      // std::cerr << "sne: " << second_min_num_errors << " min_e: " <<
      // min_num_errors << " aln_len: " << alignment_length << "id: " <<
      // alignment_identity << " tmp: " << tmp << " 1: mapq:" << (int)mapq <<
      // "\n";
    } else {
      mapq = 5 * 6.02 * (second_min_num_errors - min_num_errors) * tmp * tmp +
             0.499;
      // if (second_min_num_errors - min_num_errors < error_threshold_ + 1) {
      //  mapq = 6 * 6.02 * (second_min_num_errors - min_num_errors) * tmp * tmp
      //  + 0.499 ;
      //} else {
      //  mapq = 6 * 6.02 * (error_threshold_ + 1) * tmp * tmp + 0.499 ;
      //}
    }
    // mapq = 30 - 34.0 / error_threshold + 34.0 / error_threshold *
    // (second_min_num_errors - min_num_errors) * tmp * tmp + 0.499;
  }
  // printf("%d %d %d %d. %d\n", alignment_length, min_num_errors,
  // second_min_num_errors, mapq, read_length);
  if (num_second_best_mappings > 0) {
    mapq -= (int)(4.343 * log(num_second_best_mappings + 1) + 0.499);
    // std::cerr << " 2: mapq:" << (int)mapq << "\n";
  }
  // if (split_alignment_ && num_candidates > 1) {
  //  mapq -= (int)(4.343 * log(num_candidates) + 0.499);
  //}
  if (mapq > 60) {
    mapq = 60;
  }
  if (mapq < 0) {
    mapq = 0;
  }
  // printf("%d\n", repetitive_seed_length);
  if (repetitive_seed_length > 0) {
    // double frac_rep = (repetitive_seed_length) / (double)alignment_length;
    double frac_rep = (repetitive_seed_length) / (double)read_length;
    if (repetitive_seed_length >= (uint32_t)read_length) {
      frac_rep = 0.999;
    }
    // mapq = mapq * (1 - frac_rep / 2) + 0.499;
    // if (alignment_identity <= 0.95 && second_min_num_errors >
    // error_threshold_) {
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
  if (split_alignment_ && alignment_length < read_length - error_threshold_ &&
      second_min_num_errors != min_num_errors) {
    if (repetitive_seed_length >= alignment_length &&
        repetitive_seed_length < (uint32_t)read_length &&
        alignment_length < read_length / 3) {
      mapq = 0;
    }
    int diff = second_min_num_errors - min_num_errors;
    if (second_min_num_errors - min_num_errors <= error_threshold_ * 3 / 4 &&
        num_candidates >= 5) {
      mapq -= (num_candidates / 5 / diff);
    }
    if (mapq < 0) {
      mapq = 0;
    }
    if (num_second_best_mappings > 0 &&
        second_min_num_errors - min_num_errors <= error_threshold_ * 3 / 4) {
      mapq /= (num_second_best_mappings / diff + 1);
    }
  }
  return (uint8_t)mapq;
}

#define raw_mapq(diff, a) ((int)(5 * 6.02 * (diff) / (a) + .499))

template <typename MappingRecord>
uint8_t Chromap<MappingRecord>::GetMAPQForPairedEndRead(
    int num_positive_candidates, int num_negative_candidates,
    uint32_t repetitive_seed_length1, uint32_t repetitive_seed_length2,
    uint16_t positive_alignment_length, uint16_t negative_alignment_length,
    int min_sum_errors, int num_best_mappings, int second_min_sum_errors,
    int num_second_best_mappings, int num_errors1, int num_errors2,
    int min_num_errors1, int min_num_errors2, int num_best_mappings1,
    int num_best_mappings2, int second_min_num_errors1,
    int second_min_num_errors2, int num_second_best_mappings1,
    int num_second_best_mappings2, int read1_length, int read2_length,
    int force_mapq, uint8_t &mapq1, uint8_t &mapq2) {
  // std::cerr << " rl1:" << (int)repetitive_seed_length1 << " rl2:" <<
  // (int)repetitive_seed_length2 << " pal:" << (int)positive_alignment_length
  // << " nal:" << (int)negative_alignment_length << " me:" << min_sum_errors <<
  // " #bm:" << num_best_mappings << " sme:" << second_min_sum_errors << "
  // #sbm:"
  // << num_second_best_mappings << " ne1:" << num_errors1 << " ne2:" <<
  // num_errors2 << " me1:" << min_num_errors1 << " me2:" << min_num_errors2 <<
  // " #bm1:" << num_best_mappings1 << " #bm2:" << num_best_mappings2 << "
  // sme1:"
  // << second_min_num_errors1 << " sme2:" << second_min_num_errors2 << "
  // #sbm1:"
  // << num_second_best_mappings1 << " #sbm2:" << num_second_best_mappings2 <<
  // "\n";
  uint8_t mapq_pe = 0;
  int min_num_unpaired_sum_errors = min_num_errors1 + min_num_errors2 + 3;
  // bool is_paired = (min_num_errors1 == num_errors1 && min_num_errors2 ==
  // num_errors2);
  if (num_best_mappings <= 1) {
    int adjusted_second_min_sum_errors =
        second_min_sum_errors < min_num_unpaired_sum_errors
            ? second_min_sum_errors
            : min_num_unpaired_sum_errors;
    // mapq_pe = GetMAPQForSingleEndRead(2 * error_threshold_,
    // num_positive_candidates + num_negative_candidates,
    // repetitive_seed_length1
    // + repetitive_seed_length2, positive_alignment_length +
    // negative_alignment_length, min_sum_errors, num_best_mappings,
    // second_min_sum_errors, num_second_best_mappings, read1_length +
    // read2_length);
    mapq_pe = raw_mapq(adjusted_second_min_sum_errors - min_sum_errors, 1);
    // std::cerr << "mapqpe: " << (int)mapq_pe << "\n";
    if (num_second_best_mappings > 0) {
      mapq_pe -= (int)(4.343 * log(num_second_best_mappings + 1) + 0.499);
    }
    // if (num_positive_candidates > 10 && num_negative_candidates > 10 &&
    // second_min_sum_errors > 2 * error_threshold_) {
    //  mapq_pe -= (int)(4.343 * log(num_positive_candidates +
    //  num_negative_candidates) + 0.499);
    //  //mapq_pe *= 0.8;
    //}
    // if (num_positive_candidates > 10 && num_errors1 + 2 >=
    // second_min_num_errors1 && second_min_num_errors1 > error_threshold_) {
    //  mapq_pe -= (int)(4.343 * log(num_positive_candidates) + 0.499);
    //}
    // if (num_negative_candidates > 10 && num_errors2 + 2 >=
    // second_min_num_errors2 && second_min_num_errors2 > error_threshold_) {
    //  mapq_pe -= (int)(4.343 * log(num_negative_candidates) + 0.499);
    //}
    if (mapq_pe > 60) {
      mapq_pe = 60;
    }
    if (mapq_pe < 0) {
      mapq_pe = 0;
    }
    // std::cerr << "mapqpe: " << (int)mapq_pe << "\n";
    int repetitive_seed_length =
        repetitive_seed_length1 + repetitive_seed_length2;
    // int total_alignment_length = positive_alignment_length +
    // negative_alignment_length;
    if (repetitive_seed_length > 0) {
      double total_read_length = read1_length + read2_length;
      double frac_rep = (double)repetitive_seed_length / total_read_length;
      // double frac_rep1 = (double)repetitive_seed_length1 / read1_length;
      // double frac_rep2 = (double)repetitive_seed_length2 / read2_length;
      // frac_rep = frac_rep1 > frac_rep2 ? frac_rep1 : frac_rep2;
      if (repetitive_seed_length >= total_read_length) {
        frac_rep = 0.999;
      }
      // mapq_pe = mapq_pe * (1 - frac_rep / 2) + 0.499;
      // double alignment_identity = 1 - (double)min_sum_errors /
      // (total_read_length > total_alignment_length ? total_read_length :
      // total_alignment_length);
      double alignment_identity1 =
          1 - (double)num_errors1 / (read1_length > positive_alignment_length
                                         ? read1_length
                                         : positive_alignment_length);
      double alignment_identity2 =
          1 - (double)num_errors2 / (read2_length > negative_alignment_length
                                         ? read2_length
                                         : negative_alignment_length);
      double alignment_identity = alignment_identity1 < alignment_identity2
                                      ? alignment_identity1
                                      : alignment_identity2;
      // if (alignment_identity <= 0.95 && second_min_sum_errors > 2 *
      // error_threshold_) {
      if (alignment_identity <= 0.95) {
        mapq_pe = mapq_pe * (1 - sqrt(frac_rep)) + 0.499;
      } else if (alignment_identity <= 0.97) {
        mapq_pe = mapq_pe * (1 - frac_rep) + 0.499;
      } else if (alignment_identity >= 0.999) {
        mapq_pe =
            mapq_pe * (1 - frac_rep * frac_rep * frac_rep * frac_rep) + 0.499;
        // std::cerr << "mapqpe: " << (int)mapq_pe << "\n";
      } else {
        mapq_pe = mapq_pe * (1 - frac_rep * frac_rep) + 0.499;
      }
      // if (frac_rep > 0.8) {
      //} else {
      // mapq_pe = mapq_pe * (1 - frac_rep * frac_rep) + 0.499;
      //}
      // std::cerr <<"id1,2: " << alignment_identity1 << ", " <<
      // alignment_identity2 << " frep: " << frac_rep << " frep1: " << frac_rep1
      // << " frep2: " << frac_rep2 << " id: " << (1 - (double)min_sum_errors /
      // (total_read_length > total_alignment_length ? total_read_length :
      // total_alignment_length)) << " " << alignment_identity << " mapqpe: " <<
      // (int)mapq_pe << "\n";
    }
  }
  // mapq1 = GetMAPQForSingleEndRead(error_threshold_, num_positive_candidates,
  // repetitive_seed_length1, positive_alignment_length, min_num_errors1,
  // num_best_mappings1, second_min_num_errors1, num_second_best_mappings1,
  // read1_length);
  mapq1 = GetMAPQForSingleEndRead(
      error_threshold_, num_positive_candidates, repetitive_seed_length1,
      positive_alignment_length, num_errors1, num_best_mappings1,
      second_min_num_errors1, num_second_best_mappings1, 2, read1_length);
  // mapq2 = GetMAPQForSingleEndRead(error_threshold_, num_negative_candidates,
  // repetitive_seed_length2, negative_alignment_length, min_num_errors2,
  // num_best_mappings2, second_min_num_errors2, num_second_best_mappings2,
  // read2_length);
  mapq2 = GetMAPQForSingleEndRead(
      error_threshold_, num_negative_candidates, repetitive_seed_length2,
      negative_alignment_length, num_errors2, num_best_mappings2,
      second_min_num_errors2, num_second_best_mappings2, 2, read2_length);
  // std::cerr << " 1:" << (int)mapq1 << " 2:" << (int)mapq2 << " mapq_pe:" <<
  // (int)mapq_pe << "\n";
  if (!split_alignment_) {
    // mapq1 = mapq1 > mapq_pe ? mapq1 : mapq_pe < mapq1 + 40? mapq_pe : mapq1 +
    // 40; mapq2 = mapq2 > mapq_pe ? mapq2 : mapq_pe < mapq2 + 40? mapq_pe :
    // mapq2 + 40;
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
  // std::cerr << " 1:" << (int)mapq1 << " 2:" << (int)mapq2 << "\n\n";
  // if (second_min_num_errors1 > error_threshold_) {
  //  second_min_num_errors1 = 2 * error_threshold_ + 1;
  //}
  // if (second_min_num_errors2 > error_threshold_) {
  //  second_min_num_errors2 = 2 * error_threshold_ + 1;
  //}
  // int mapq_coef = 60;
  // uint8_t max_mapq1 = (int)(mapq_coef * ((double)(second_min_num_errors1 -
  // min_num_errors1) / second_min_num_errors1) + .499); uint8_t max_mapq2 =
  // (int)(mapq_coef * ((double)(second_min_num_errors2 - min_num_errors2) /
  // second_min_num_errors2) + .499); uint8_t max_mapq1 = 30 - 34.0 /
  // error_threshold_ + 34.0 / error_threshold_ * (second_min_num_errors1 -
  // min_num_errors1); uint8_t max_mapq2 = 30 - 34.0 / error_threshold_ + 34.0 /
  // error_threshold_ * (second_min_num_errors2 - min_num_errors2); mapq1 =
  // mapq1 < max_mapq1 ? mapq1 : max_mapq1; mapq2 = mapq2 < max_mapq2 ? mapq2 :
  // max_mapq2; mapq1 = mapq1 < raw_mapq(second_min_num_errors1 -
  // min_num_errors1, 1) ? mapq1 : raw_mapq(second_min_num_errors1 -
  // min_num_errors1, 1); if (second_min_num_errors1 < num_errors1) {
  //  second_min_num_errors1 = num_errors1;
  //}
  // mapq1 = mapq1 < raw_mapq(second_min_num_errors1 - num_errors1, 1) ? mapq1 :
  // raw_mapq(second_min_num_errors1 - num_errors1, 1);
  ////mapq2 = mapq2 < raw_mapq(second_min_num_errors2 - min_num_errors2, 1) ?
  /// mapq2 : raw_mapq(second_min_num_errors2 - min_num_errors2, 1);
  // if (second_min_num_errors2 < num_errors2) {
  //  second_min_num_errors2 = num_errors2;
  //}
  // mapq2 = mapq2 < raw_mapq(second_min_num_errors2 - num_errors2, 1) ? mapq2 :
  // raw_mapq(second_min_num_errors2 - num_errors2, 1);
  uint8_t mapq = mapq1 < mapq2 ? mapq1 : mapq2;
  // uint8_t mapq = (mapq1 + mapq2) / 2;
  // uint8_t mapq = (mapq1 + mapq2) / 2 + 0.499;
  // if ((mapq1 < mapq2 && mapq1 + 30 >= mapq2) || (mapq2 < mapq1 && mapq2 + 30
  // >= mapq1)) {
  //  uint8_t min_mapq = mapq1 < mapq2 ? mapq1 : mapq2;
  //  mapq = min_mapq;
  //  //mapq = (mapq + 1 * min_mapq) / 2 + 0.499;
  //} else {
  //  uint8_t max_mapq = mapq1 > mapq2 ? mapq1 : mapq2;
  //  mapq = (mapq + 1 * max_mapq) / 2 + 0.499;
  //}
  // std::cerr << " 1:" << (int)mapq1 << " 2:" << (int)mapq2 << " mapq_pe:" <<
  // (int)mapq_pe << " mapq:" << (int)mapq << "\n";

  if (mapq < 60 && force_mapq >= 0 && force_mapq < mapq) {
    mapq = force_mapq;
  }

  return (uint8_t)mapq;
}

template <typename MappingRecord>
void Chromap<MappingRecord>::ParseReadFormat(const std::string &read_format) {
  uint32_t i;
  int j = 0;
  int k = 0;  // for read1, read2, or barcode
  read1_format_[0] = 0;
  read1_format_[1] = -1;
  read2_format_[0] = 0;
  read2_format_[1] = -1;
  barcode_format_[0] = 0;
  barcode_format_[1] = -1;
  int fields[2] = {0, -1};
  char buffer[20];
  int blen = 0;
  for (i = 0; i < read_format.size(); ++i) {
    if (read_format[i] == ',' || i == 0) {
      if (i > 0) {
        buffer[blen] = '\0';
        fields[j] = atoi(buffer);

        if (k == 0)
          memcpy(read1_format_, fields, sizeof(fields));
        else if (k == 1)
          memcpy(read2_format_, fields, sizeof(fields));
        else
          memcpy(barcode_format_, fields, sizeof(fields));

        ++i;
      }
      if (read_format[i] == 'r' && read_format[i + 1] == '1')
        k = 0;
      else if (read_format[i] == 'r' && read_format[i + 1] == '2')
        k = 1;
      else if (read_format[i] == 'b' && read_format[i + 1] == 'c')
        k = 2;
      else
        ExitWithMessage("Unknown read format: " + read_format + "\n");
      j = 0;
      fields[0] = fields[1] = 0;
      i += 2;
      blen = 0;
    } else {
      if (read_format[i] != ':') {
        if (j < 2) {
          buffer[blen] = read_format[i];
          ++blen;
        } else {
          ExitWithMessage("Unknown read format: " + read_format + "\n");
        }
      } else {
        buffer[blen] = '\0';
        fields[j] = atoi(buffer);
        ++j;
        blen = 0;
      }
    }
  }
  buffer[blen] = '\0';
  fields[j] = atoi(buffer);
  // By initialization, it is fine even if there is no read_format specified
  if (k == 0)
    memcpy(read1_format_, fields, sizeof(fields));
  else if (k == 1)
    memcpy(read2_format_, fields, sizeof(fields));
  else
    memcpy(barcode_format_, fields, sizeof(fields));
}

void ChromapDriver::ParseArgsAndRun(int argc, char *argv[]) {
  cxxopts::Options options(
      "chromap", "Fast alignment and preprocessing of chromatin profiles");
  options.add_options("Indexing")("i,build-index", "Build index")(
      "min-frag-length",
      "Min fragment length for choosing k and w automatically [30]",
      cxxopts::value<int>(),
      "INT")("k,kmer", "Kmer length [17]", cxxopts::value<int>(), "INT")(
      "w,window", "Window size [7]", cxxopts::value<int>(), "INT");
  options.set_width(120).add_options("Mapping")
      //("m,map", "Map reads")
      ("preset",
       "Preset parameters for mapping reads (always applied before other "
       "options) []\natac: mapping ATAC-seq/scATAC-seq reads\nchip: mapping "
       "ChIP-seq reads\nhic: mapping Hi-C reads",
       cxxopts::value<std::string>(),
       "STR")("split-alignment", "Allow split alignments")(
          "e,error-threshold", "Max # errors allowed to map a read [8]",
          cxxopts::value<int>(), "INT")
      //("A,match-score", "Match score [1]", cxxopts::value<int>(), "INT")
      //("B,mismatch-penalty", "Mismatch penalty [4]", cxxopts::value<int>(),
      //"INT")
      //("O,gap-open-penalties", "Gap open penalty [6,6]",
      // cxxopts::value<std::vector<int>>(), "INT[,INT]")
      //("E,gap-extension-penalties", "Gap extension penalty [1,1]",
      // cxxopts::value<std::vector<int>>(), "INT[,INT]")
      ("s,min-num-seeds", "Min # seeds to try to map a read [2]",
       cxxopts::value<int>(),
       "INT")("f,max-seed-frequencies",
              "Max seed frequencies for a seed to be selected [500,1000]",
              cxxopts::value<std::vector<int>>(), "INT[,INT]")
      //("n,max-num-best-mappings", "Only report n best mappings [1]",
      // cxxopts::value<int>(), "INT")
      ("l,max-insert-size",
       "Max insert size, only for paired-end read mapping [1000]",
       cxxopts::value<int>(),
       "INT")("q,MAPQ-threshold",
              "Min MAPQ in range [0, 60] for mappings to be output [30]",
              cxxopts::value<uint8_t>(),
              "INT")("min-read-length", "Min read length [30]",
                     cxxopts::value<int>(), "INT")
      //("multi-mapping-allocation-distance", "Uni-mappings within this distance
      // from any end of multi-mappings are used for allocation [0]",
      // cxxopts::value<int>(), "INT")
      //("multi-mapping-allocation-seed", "Seed for random number generator in
      // multi-mapping allocation [11]", cxxopts::value<int>(), "INT")
      //("drop-repetitive-reads", "Drop reads with too many best mappings
      //[500000]", cxxopts::value<int>(), "INT")
      ("trim-adapters", "Try to trim adapters on 3'")("remove-pcr-duplicates",
                                                      "Remove PCR duplicates")(
          "remove-pcr-duplicates-at-bulk-level",
          "Remove PCR duplicates at bulk level for single cell data")(
          "remove-pcr-duplicates-at-cell-level",
          "Remove PCR duplicates at cell level for single cell data")
      //("allocate-multi-mappings", "Allocate multi-mappings")
      ("Tn5-shift", "Perform Tn5 shift")("low-mem", "Use low memory mode")(
          "bc-error-threshold",
          "Max Hamming distance allowed to correct a barcode [1]",
          cxxopts::value<int>(),
          "INT")("bc-probability-threshold",
                 "Min probability to correct a barcode [0.9]",
                 cxxopts::value<double>(),
                 "FLT")("t,num-threads", "# threads for mapping [1]",
                        cxxopts::value<int>(), "INT");
  // options.add_options("Peak")
  //  ("cell-by-bin", "Generate cell-by-bin matrix")
  //  ("bin-size", "Bin size to generate cell-by-bin matrix [5000]",
  //  cxxopts::value<int>(), "INT")
  //  ("depth-cutoff", "Depth cutoff for peak calling [3]",
  //  cxxopts::value<int>(), "INT")
  //  ("peak-min-length", "Min length of peaks to report [30]",
  //  cxxopts::value<int>(), "INT")
  //  ("peak-merge-max-length", "Peaks within this length will be merged [30]",
  //  cxxopts::value<int>(), "INT");
  options.add_options("Input")("r,ref", "Reference file",
                               cxxopts::value<std::string>(), "FILE")(
      "x,index", "Index file", cxxopts::value<std::string>(), "FILE")(
      "1,read1", "Single-end read files or paired-end read files 1",
      cxxopts::value<std::vector<std::string>>(),
      "FILE")("2,read2", "Paired-end read files 2",
              cxxopts::value<std::vector<std::string>>(),
              "FILE")("b,barcode", "Cell barcode files",
                      cxxopts::value<std::vector<std::string>>(), "FILE")(
      "barcode-whitelist", "Cell barcode whitelist file",
      cxxopts::value<std::string>(),
      "FILE")("read-format",
              "Format for read files and barcode files  [\"r1:0:-1,bc:0:-1\" "
              "as 10x Genomics single-end format]",
              cxxopts::value<std::string>(), "STR");
  options.add_options("Output")("o,output", "Output file",
                                cxxopts::value<std::string>(), "FILE")
      //("p,matrix-output-prefix", "Prefix of matrix output files",
      // cxxopts::value<std::string>(), "FILE")
      ("output-mappings-not-in-whitelist",
       "Output mappings with barcode not in the whitelist")
      ("chr-order", "custom chromsome order", cxxopts::value<std::string>(),
          "FILE")("BED", "Output mappings in BED/BEDPE format")(
          "TagAlign", "Output mappings in TagAlign/PairedTagAlign format")(
          "SAM", "Output mappings in SAM format")(
          "pairs",
          "Output mappings in pairs format (defined by 4DN for HiC data)")(
          "pairs-natural-chr-order",
          "natural chromosome order for pairs flipping",
          cxxopts::value<std::string>(), "FILE");
  //("PAF", "Output mappings in PAF format (only for test)");
  options.add_options()("v,version", "Print version")("h,help", "Print help");
  options.add_options("Development options")("A,match-score", "Match score [1]",
                                             cxxopts::value<int>(), "INT")(
      "B,mismatch-penalty", "Mismatch penalty [4]", cxxopts::value<int>(),
      "INT")("O,gap-open-penalties", "Gap open penalty [6,6]",
             cxxopts::value<std::vector<int>>(), "INT[,INT]")(
      "E,gap-extension-penalties", "Gap extension penalty [1,1]",
      cxxopts::value<std::vector<int>>(),
      "INT[,INT]")("n,max-num-best-mappings", "Only report n best mappings [1]",
                   cxxopts::value<int>(),
                   "INT")("multi-mapping-allocation-distance",
                          "Uni-mappings within this distance from any end of "
                          "multi-mappings are used for allocation [0]",
                          cxxopts::value<int>(), "INT")(
      "multi-mapping-allocation-seed",
      "Seed for random number generator in multi-mapping allocation [11]",
      cxxopts::value<int>(), "INT")(
      "drop-repetitive-reads",
      "Drop reads with too many best mappings [500000]", cxxopts::value<int>(),
      "INT")("allocate-multi-mappings", "Allocate multi-mappings")(
      "PAF", "Output mappings in PAF format (only for test)")
      ("skip-barcode-check", "Do not check whether too few barcodes are in the whitelist")
        ;

  auto result = options.parse(argc, argv);
  if (result.count("h")) {
    std::cerr << options.help(
        {"", "Indexing", "Mapping", "Peak", "Input", "Output"});
    return;
  }
  if (result.count("v")) {
    std::cerr << CHROMAP_VERSION << "\n";
    return;
  }
  // Parameters and their default
  IndexParameters index_parameters;
  MappingParameters mapping_parameters;

  if (result.count("preset")) {
    std::string read_type = result["preset"].as<std::string>();
    if (read_type == "atac") {
      std::cerr << "Preset parameters for ATAC-seq/scATAC-seq are used.\n";
      mapping_parameters.max_insert_size = 2000;
      mapping_parameters.trim_adapters = true;
      mapping_parameters.remove_pcr_duplicates = true;
      mapping_parameters.remove_pcr_duplicates_at_bulk_level = false;
      mapping_parameters.Tn5_shift = true;
      mapping_parameters.mapping_output_format = MAPPINGFORMAT_BED;
      mapping_parameters.low_memory_mode = true;
    } else if (read_type == "chip") {
      std::cerr << "Preset parameters for ChIP-seq are used.\n";
      mapping_parameters.max_insert_size = 2000;
      mapping_parameters.remove_pcr_duplicates = true;
      mapping_parameters.low_memory_mode = true;
      mapping_parameters.mapping_output_format = MAPPINGFORMAT_BED;
    } else if (read_type == "hic") {
      std::cerr << "Preset parameters for Hi-C are used.\n";
      mapping_parameters.error_threshold = 4;
      mapping_parameters.mapq_threshold = 1;
      mapping_parameters.split_alignment = true;
      mapping_parameters.low_memory_mode = true;
      mapping_parameters.mapping_output_format = MAPPINGFORMAT_PAIRS;
    } else {
      chromap::Chromap<>::ExitWithMessage("Unrecognized preset parameters " +
                                          read_type + "\n");
    }
  }
  // Optional parameters
  if (result.count("min-frag-length")) {
    int min_fragment_length = result["min-frag-length"].as<int>();
    if (min_fragment_length <= 60) {
      index_parameters.kmer_size = 17;
      index_parameters.window_size = 7;
    } else if (min_fragment_length <= 80) {
      index_parameters.kmer_size = 19;
      index_parameters.window_size = 10;
    } else {
      index_parameters.kmer_size = 23;
      index_parameters.window_size = 11;
    }
  }
  if (result.count("k")) {
    index_parameters.kmer_size = result["kmer"].as<int>();
  }
  if (result.count("w")) {
    index_parameters.window_size = result["window"].as<int>();
  }
  if (result.count("e")) {
    mapping_parameters.error_threshold = result["error-threshold"].as<int>();
  }
  if (result.count("A")) {
    mapping_parameters.match_score = result["match-score"].as<int>();
  }
  if (result.count("B")) {
    mapping_parameters.mismatch_penalty = result["mismatch-penalty"].as<int>();
  }
  if (result.count("O")) {
    mapping_parameters.gap_open_penalties =
        result["gap-open-penalties"].as<std::vector<int>>();
  }
  if (result.count("E")) {
    mapping_parameters.gap_extension_penalties =
        result["gap-extension-penalties"].as<std::vector<int>>();
  }
  if (result.count("s")) {
    mapping_parameters.min_num_seeds_required_for_mapping =
        result["min-num-seeds"].as<int>();
  }
  if (result.count("f")) {
    mapping_parameters.max_seed_frequencies =
        result["max-seed-frequencies"].as<std::vector<int>>();
  }
  if (result.count("n")) {
    mapping_parameters.max_num_best_mappings =
        result["max-num-best-mappings"].as<int>();
  }
  if (result.count("l")) {
    mapping_parameters.max_insert_size = result["max-insert-size"].as<int>();
  }
  if (result.count("q")) {
    mapping_parameters.mapq_threshold = result["MAPQ-threshold"].as<uint8_t>();
  }
  if (result.count("t")) {
    mapping_parameters.num_threads = result["num-threads"].as<int>();
  }
  if (result.count("min-read-length")) {
    mapping_parameters.min_read_length = result["min-read-length"].as<int>();
  }
  if (result.count("bc-error-threshold")) {
    mapping_parameters.barcode_correction_error_threshold =
        result["bc-error-threshold"].as<int>();
  }
  if (result.count("bc-probability-threshold")) {
    mapping_parameters.barcode_correction_probability_threshold =
        result["bc-probability-threshold"].as<double>();
  }
  if (result.count("multi-mapping-allocation-distance")) {
    mapping_parameters.multi_mapping_allocation_distance =
        result["multi-mapping-allocation-distance"].as<int>();
  }
  if (result.count("multi-mapping-allocation-seed")) {
    mapping_parameters.multi_mapping_allocation_seed =
        result["multi-mapping-allocation-seed"].as<int>();
  }
  if (result.count("drop-repetitive-reads")) {
    mapping_parameters.drop_repetitive_reads =
        result["drop-repetitive-reads"].as<int>();
  }
  if (result.count("trim-adapters")) {
    mapping_parameters.trim_adapters = true;
  }
  if (result.count("remove-pcr-duplicates")) {
    mapping_parameters.remove_pcr_duplicates = true;
  }
  if (result.count("remove-pcr-duplicates-at-bulk-level")) {
    mapping_parameters.remove_pcr_duplicates_at_bulk_level = true;
  }
  if (result.count("remove-pcr-duplicates-at-cell-level")) {
    mapping_parameters.remove_pcr_duplicates_at_bulk_level = false;
  }
  if (result.count("allocate-multi-mappings")) {
    mapping_parameters.allocate_multi_mappings = true;
    mapping_parameters.only_output_unique_mappings = false;
  }
  if (result.count("Tn5-shift")) {
    mapping_parameters.Tn5_shift = true;
  }
  if (result.count("split-alignment")) {
    mapping_parameters.split_alignment = true;
  }
  if (result.count("output-mappings-not-in-whitelist")) {
    mapping_parameters.output_mappings_not_in_whitelist = true;
  }
  if (result.count("BED")) {
    mapping_parameters.mapping_output_format = MAPPINGFORMAT_BED;
  }
  if (result.count("TagAlign")) {
    mapping_parameters.mapping_output_format = MAPPINGFORMAT_TAGALIGN;
  }
  if (result.count("PAF")) {
    mapping_parameters.mapping_output_format = MAPPINGFORMAT_PAF;
  }
  if (result.count("pairs")) {
    mapping_parameters.mapping_output_format = MAPPINGFORMAT_PAIRS;
  }
  if (result.count("SAM")) {
    mapping_parameters.mapping_output_format = MAPPINGFORMAT_SAM;
  }
  if (result.count("low-mem")) {
    mapping_parameters.low_memory_mode = true;
  }
  if (result.count("cell-by-bin")) {
    mapping_parameters.cell_by_bin = true;
  }
  if (result.count("bin-size")) {
    mapping_parameters.bin_size = result["bin-size"].as<int>();
  }
  if (result.count("depth-cutoff")) {
    mapping_parameters.depth_cutoff_to_call_peak =
        result["depth-cutoff"].as<uint16_t>();
  }
  if (result.count("peak-min-length")) {
    mapping_parameters.peak_min_length = result["peak-min-length"].as<int>();
  }
  if (result.count("peak-merge-max-length")) {
    mapping_parameters.peak_merge_max_length =
        result["peak-merge-max-length"].as<int>();
  }

  std::cerr << std::setprecision(2) << std::fixed;
  if (result.count("i")) {
    if (result.count("r")) {
      index_parameters.reference_file_path = result["ref"].as<std::string>();
    } else {
      chromap::Chromap<>::ExitWithMessage("No reference specified!");
    }
    if (result.count("o")) {
      index_parameters.index_output_file_path =
          result["output"].as<std::string>();
    } else {
      chromap::Chromap<>::ExitWithMessage("No output file specified!");
    }
    std::cerr << "Build index for the reference.\n";
    std::cerr << "Kmer length: " << index_parameters.kmer_size
              << ", window size: " << index_parameters.window_size << "\n";
    std::cerr << "Reference file: " << index_parameters.reference_file_path
              << "\n";
    std::cerr << "Output file: " << index_parameters.index_output_file_path
              << "\n";
    chromap::Chromap<> chromap_for_indexing(index_parameters);
    chromap_for_indexing.ConstructIndex();
  } else if (result.count("1")) {
    std::cerr << "Start to map reads.\n";
    if (result.count("r")) {
      mapping_parameters.reference_file_path = result["ref"].as<std::string>();
    } else {
      chromap::Chromap<>::ExitWithMessage("No reference specified!");
    }
    if (result.count("o")) {
      mapping_parameters.mapping_output_file_path =
          result["output"].as<std::string>();
    } else {
      chromap::Chromap<>::ExitWithMessage("No output file specified!");
    }
    if (result.count("x")) {
      mapping_parameters.index_file_path = result["index"].as<std::string>();
    } else {
      chromap::Chromap<>::ExitWithMessage("No index file specified!");
    }
    if (result.count("1")) {
      mapping_parameters.read_file1_paths =
          result["read1"].as<std::vector<std::string>>();
    } else {
      chromap::Chromap<>::ExitWithMessage("No read file specified!");
    }
    if (result.count("2")) {
      mapping_parameters.read_file2_paths =
          result["read2"].as<std::vector<std::string>>();
    }
    if (result.count("b")) {
      mapping_parameters.is_bulk_data = false;
      mapping_parameters.barcode_file_paths =
          result["barcode"].as<std::vector<std::string>>();
      if (result.count("barcode-whitelist") == 0) {
        chromap::Chromap<>::ExitWithMessage(
            "There are input barcode files but a barcode whitelist file is "
            "missing!");
      }
    }
    if (result.count("barcode-whitelist")) {
      if (mapping_parameters.is_bulk_data) {
        chromap::Chromap<>::ExitWithMessage(
            "No barcode file specified but the barcode whitelist file is "
            "given!");
      }
      mapping_parameters.barcode_whitelist_file_path =
          result["barcode-whitelist"].as<std::string>();
    }
    if (result.count("p")) {
      mapping_parameters.matrix_output_prefix =
          result["matrix-output-prefix"].as<std::string>();
      if (mapping_parameters.is_bulk_data) {
        chromap::Chromap<>::ExitWithMessage(
            "No barcode file specified but asked to output matrix files!");
      }
    }
    if (result.count("read-format")) {
      mapping_parameters.read_format = result["read-format"].as<std::string>();
    }

    if (result.count("chr-order")) {
      mapping_parameters.custom_rid_order_path =
          result["chr-order"].as<std::string>();
    }

    if (result.count("pairs-natural-chr-order")) {
      mapping_parameters.pairs_custom_rid_order_path =
          result["pairs-natural-chr-order"].as<std::string>();
    }

    if (result.count("skip-barcode-check")) {
      mapping_parameters.skip_barcode_check = true;
    }

    // std::cerr << "Parameters: error threshold: " << error_threshold << ",
    // match score: " << match_score << ", mismatch_penalty: " <<
    // mismatch_penalty << ", gap open penalties for deletions and insertions: "
    // << gap_open_penalties[0] << "," << gap_open_penalties[1] << ", gap
    // extension penalties for deletions and insertions: " <<
    // gap_extension_penalties[0] << "," << gap_extension_penalties[1] << ",
    // min-num-seeds: " << min_num_seeds_required_for_mapping << ",
    // max-seed-frequency: " << max_seed_frequencies[0] << "," <<
    // max_seed_frequencies[1] << ", max-num-best-mappings: " <<
    // max_num_best_mappings << ", max-insert-size: " << max_insert_size << ",
    // MAPQ-threshold: " << (int)mapq_threshold << ", min-read-length: " <<
    // min_read_length << ", multi-mapping-allocation-distance: " <<
    // multi_mapping_allocation_distance << ", multi-mapping-allocation-seed: "
    // << multi_mapping_allocation_seed << ", drop-repetitive-reads: " <<
    // drop_repetitive_reads << "\n";
    std::cerr << "Parameters: error threshold: "
              << mapping_parameters.error_threshold << ", min-num-seeds: "
              << mapping_parameters.min_num_seeds_required_for_mapping
              << ", max-seed-frequency: "
              << mapping_parameters.max_seed_frequencies[0] << ","
              << mapping_parameters.max_seed_frequencies[1]
              << ", max-num-best-mappings: "
              << mapping_parameters.max_num_best_mappings
              << ", max-insert-size: " << mapping_parameters.max_insert_size
              << ", MAPQ-threshold: " << (int)mapping_parameters.mapq_threshold
              << ", min-read-length: " << mapping_parameters.min_read_length
              << ", bc-error-threshold: "
              << mapping_parameters.barcode_correction_error_threshold
              << ", bc-probability-threshold: "
              << mapping_parameters.barcode_correction_probability_threshold
              << "\n";
    std::cerr << "Number of threads: " << mapping_parameters.num_threads
              << "\n";
    if (mapping_parameters.is_bulk_data) {
      std::cerr << "Analyze bulk data.\n";
    } else {
      std::cerr << "Analyze single-cell data.\n";
    }
    if (mapping_parameters.trim_adapters) {
      std::cerr << "Will try to remove adapters on 3'.\n";
    } else {
      std::cerr << "Won't try to remove adapters on 3'.\n";
    }
    if (mapping_parameters.remove_pcr_duplicates) {
      std::cerr << "Will remove PCR duplicates after mapping.\n";
    } else {
      std::cerr << "Won't remove PCR duplicates after mapping.\n";
    }
    if (mapping_parameters.remove_pcr_duplicates_at_bulk_level) {
      std::cerr << "Will remove PCR duplicates at bulk level.\n";
    } else {
      std::cerr << "Will remove PCR duplicates at cell level.\n";
    }
    if (mapping_parameters.allocate_multi_mappings) {
      std::cerr << "Will allocate multi-mappings after mapping.\n";
    } else {
      std::cerr << "Won't allocate multi-mappings after mapping.\n";
    }
    if (mapping_parameters.only_output_unique_mappings) {
      std::cerr << "Only output unique mappings after mapping.\n";
    }
    if (!mapping_parameters.output_mappings_not_in_whitelist) {
      std::cerr << "Only output mappings of which barcodes are in whitelist.\n";
    } else {
      std::cerr << "No filtering of mappings based on whether their barcodes "
                   "are in whitelist.\n";
    }
    // if (allocate_multi_mappings && only_output_unique_mappings) {
    //  std::cerr << "WARNING: you want to output unique mappings only but you
    //  ask to allocate multi-mappings! In this case, it won't allocate
    //  multi-mappings and will only output unique mappings.\n";
    //  allocate_multi_mappings = false;
    //}
    if (mapping_parameters.max_num_best_mappings >
        mapping_parameters.drop_repetitive_reads) {
      std::cerr << "WARNING: you want to drop mapped reads with more than "
                << mapping_parameters.drop_repetitive_reads
                << " mappings. But you want to output top "
                << mapping_parameters.max_num_best_mappings
                << " best mappings. In this case, only reads with <="
                << mapping_parameters.drop_repetitive_reads
                << " best mappings will be output.\n";
      mapping_parameters.max_num_best_mappings =
          mapping_parameters.drop_repetitive_reads;
    }
    if (mapping_parameters.Tn5_shift) {
      std::cerr << "Perform Tn5 shift.\n";
    }
    if (mapping_parameters.split_alignment) {
      std::cerr << "Allow split alignment.\n";
    }

    switch (mapping_parameters.mapping_output_format) {
      case MAPPINGFORMAT_BED:
        std::cerr << "Output mappings in BED/BEDPE format.\n";
        break;
      case MAPPINGFORMAT_TAGALIGN:
        std::cerr << "Output mappings in TagAlign/PairedTagAlign format.\n";
        break;
      case MAPPINGFORMAT_PAF:
        std::cerr << "Output mappings in PAF format.\n";
        break;
      case MAPPINGFORMAT_SAM:
        std::cerr << "Output mappings in SAM format.\n";
        break;
      case MAPPINGFORMAT_PAIRS:
        std::cerr << "Output mappings in pairs format.\n";
        break;
      default:
        chromap::Chromap<>::ExitWithMessage("Unknown mapping output format!");
        break;
    }

    std::cerr << "Reference file: " << mapping_parameters.reference_file_path
              << "\n";
    std::cerr << "Index file: " << mapping_parameters.index_file_path << "\n";
    for (size_t i = 0; i < mapping_parameters.read_file1_paths.size(); ++i) {
      std::cerr << i + 1
                << "th read 1 file: " << mapping_parameters.read_file1_paths[i]
                << "\n";
    }
    if (result.count("2") != 0) {
      for (size_t i = 0; i < mapping_parameters.read_file2_paths.size(); ++i) {
        std::cerr << i + 1 << "th read 2 file: "
                  << mapping_parameters.read_file2_paths[i] << "\n";
      }
    }
    if (result.count("b") != 0) {
      for (size_t i = 0; i < mapping_parameters.barcode_file_paths.size();
           ++i) {
        std::cerr << i + 1 << "th cell barcode file: "
                  << mapping_parameters.barcode_file_paths[i] << "\n";
      }
    }
    if (result.count("barcode-whitelist") != 0) {
      std::cerr << "Cell barcode whitelist file: "
                << mapping_parameters.barcode_whitelist_file_path << "\n";
    }
    std::cerr << "Output file: " << mapping_parameters.mapping_output_file_path
              << "\n";
    if (result.count("matrix-output-prefix") != 0) {
      std::cerr << "Matrix output prefix: "
                << mapping_parameters.matrix_output_prefix << "\n";
    }

    if (result.count("2") == 0) {
      // Single-end reads.
      switch (mapping_parameters.mapping_output_format) {
        case MAPPINGFORMAT_PAF: {
          chromap::Chromap<chromap::PAFMapping> chromap_for_mapping(
              mapping_parameters);
          chromap_for_mapping.MapSingleEndReads();
          break;
        }
        case MAPPINGFORMAT_SAM: {
          chromap::Chromap<chromap::SAMMapping> chromap_for_mapping(
              mapping_parameters);
          chromap_for_mapping.MapSingleEndReads();
          break;
        }
        case MAPPINGFORMAT_PAIRS:
          chromap::Chromap<>::ExitWithMessage(
              "No support for single-end HiC yet!");
          break;
        case MAPPINGFORMAT_BED:
        case MAPPINGFORMAT_TAGALIGN:
          if (result.count("b") != 0) {
            chromap::Chromap<chromap::MappingWithBarcode> chromap_for_mapping(
                mapping_parameters);
            chromap_for_mapping.MapSingleEndReads();
          } else {
            chromap::Chromap<chromap::MappingWithoutBarcode>
                chromap_for_mapping(mapping_parameters);
            chromap_for_mapping.MapSingleEndReads();
          }
          break;
        default:
          chromap::Chromap<>::ExitWithMessage("Unknown mapping output format!");
          break;
      }
    } else {
      // Paired-end reads.
      switch (mapping_parameters.mapping_output_format) {
        case MAPPINGFORMAT_PAF: {
          chromap::Chromap<chromap::PairedPAFMapping> chromap_for_mapping(
              mapping_parameters);
          chromap_for_mapping.MapPairedEndReads();
          break;
        }
        case MAPPINGFORMAT_SAM: {
          chromap::Chromap<chromap::SAMMapping> chromap_for_mapping(
              mapping_parameters);
          chromap_for_mapping.MapPairedEndReads();
          break;
        }
        case MAPPINGFORMAT_PAIRS: {
          chromap::Chromap<chromap::PairsMapping> chromap_for_mapping(
              mapping_parameters);
          chromap_for_mapping.MapPairedEndReads();
          break;
        }
        case MAPPINGFORMAT_BED:
        case MAPPINGFORMAT_TAGALIGN:
          if (result.count("b") != 0) {
            chromap::Chromap<chromap::PairedEndMappingWithBarcode>
                chromap_for_mapping(mapping_parameters);
            chromap_for_mapping.MapPairedEndReads();
          } else {
            chromap::Chromap<chromap::PairedEndMappingWithoutBarcode>
                chromap_for_mapping(mapping_parameters);
            chromap_for_mapping.MapPairedEndReads();
          }
          break;
        default:
          chromap::Chromap<>::ExitWithMessage("Unknown mapping output format!");
          break;
      }
    }
  } else {
    std::cerr << options.help(
        {"", "Indexing", "Mapping", "Peak", "Input", "Output"});
  }
}
}  // namespace chromap

int main(int argc, char *argv[]) {
  chromap::ChromapDriver chromap_driver;
  chromap_driver.ParseArgsAndRun(argc, argv);
  return 0;
}
