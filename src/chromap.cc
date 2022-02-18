#include "chromap.h"

#include <assert.h>
#include <math.h>
#include <omp.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <type_traits>

#include "alignment.h"
#include "candidate_processor.h"
#include "cxxopts.hpp"
#include "mapping_generator.h"
//#include "ksw.h"
#include "mmcache.hpp"
#include "utils.h"

namespace chromap {
template <typename MappingRecord>
void Chromap<MappingRecord>::TrimAdapterForPairedEndRead(
    uint32_t pair_index, SequenceBatch &read_batch1,
    SequenceBatch &read_batch2) {
  const char *read1 = read_batch1.GetSequenceAt(pair_index);
  uint32_t read2_length = read_batch2.GetSequenceLengthAt(pair_index);
  const std::string &negative_read2 =
      read_batch2.GetNegativeSequenceAt(pair_index);
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
        read_batch1.TrimSequenceAt(pair_index, overlap_length);
        read_batch2.TrimSequenceAt(pair_index, overlap_length);
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
    uint32_t barcode_index, SequenceBatch &barcode_batch,
    uint64_t &num_barcode_in_whitelist, uint64_t &num_corrected_barcode) {
  const uint32_t barcode_length =
      barcode_batch.GetSequenceLengthAt(barcode_index);
  const uint64_t barcode_key = barcode_batch.GenerateSeedFromSequenceAt(
      barcode_index, 0, barcode_length);
  khiter_t barcode_whitelist_lookup_table_iterator =
      kh_get(k64_seq, barcode_whitelist_lookup_table_, barcode_key);
  if (barcode_whitelist_lookup_table_iterator !=
      kh_end(barcode_whitelist_lookup_table_)) {
    // Correct barcode
    ++num_barcode_in_whitelist;
    return true;
  } else if (barcode_correction_error_threshold_ > 0) {
    // Need to correct this barcode
    // const char *barcode = barcode_batch->GetSequenceAt(barcode_index);
    // std::cerr << barcode_index << " barcode " << barcode << " needs
    // correction\n";
    const char *barcode_qual = barcode_batch.GetSequenceQualAt(barcode_index);
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
      barcode_batch.CorrectBaseAt(
          barcode_index, corrected_barcodes_with_quals[0].corrected_base_index1,
          corrected_barcodes_with_quals[0].correct_base1);
      if (corrected_barcodes_with_quals[0].correct_base2 != 0) {
        barcode_batch.CorrectBaseAt(
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
      ++num_corrected_barcode;
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
        barcode_batch.CorrectBaseAt(
            barcode_index,
            corrected_barcodes_with_quals[best_corrected_barcode_index]
                .corrected_base_index1,
            corrected_barcodes_with_quals[best_corrected_barcode_index]
                .correct_base1);
        if (corrected_barcodes_with_quals[best_corrected_barcode_index]
                .correct_base2 != 0) {
          barcode_batch.CorrectBaseAt(
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
        ++num_corrected_barcode;
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
uint32_t Chromap<MappingRecord>::LoadPairedEndReadsWithBarcodes(
    SequenceBatch &read_batch1, SequenceBatch &read_batch2,
    SequenceBatch &barcode_batch) {
  // double real_start_time = Chromap<>::GetRealTime();
  uint32_t num_loaded_pairs = 0;
  while (num_loaded_pairs < read_batch_size_) {
    bool no_more_read1 = read_batch1.LoadOneSequenceAndSaveAt(num_loaded_pairs);
    bool no_more_read2 = read_batch2.LoadOneSequenceAndSaveAt(num_loaded_pairs);
    bool no_more_barcode = no_more_read2;
    if (!is_bulk_data_) {
      no_more_barcode =
          barcode_batch.LoadOneSequenceAndSaveAt(num_loaded_pairs);
    }
    if ((!no_more_read1) && (!no_more_read2) && (!no_more_barcode)) {
      if (read_batch1.GetSequenceLengthAt(num_loaded_pairs) <
              (uint32_t)min_read_length_ ||
          read_batch2.GetSequenceLengthAt(num_loaded_pairs) <
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
      ExitWithMessage("Numbers of reads and barcodes don't match!");
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
  double real_start_time = GetRealTime();
  SequenceBatch barcode_batch(read_batch_size_);
  barcode_batch.SetSeqEffectiveRange(barcode_format_[0], barcode_format_[1],
                                     barcode_format_[2]);
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
      if (!skip_barcode_check_ &&
          num_sample_barcodes_ * 20 < num_loaded_barcodes) {
        // Since num_loaded_pairs is a constant, this if is actuaclly only
        // effective in the first iteration
        ExitWithMessage(
            "Less than 5\% barcodes can be found or corrected based on the "
            "barcode whitelist.\nPlease check whether the barcode whitelist "
            "matches the data, e.g. length, reverse-complement. If this is a "
            "false warning, please run Chromap with the option "
            "--skip-barcode-check.");
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
            << " in " << GetRealTime() - real_start_time << "s.\n";
}

template <typename MappingRecord>
void Chromap<MappingRecord>::UpdateBarcodeAbundance(
    uint32_t num_loaded_barcodes, const SequenceBatch &barcode_batch) {
  double real_start_time = GetRealTime();
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
            << " in " << GetRealTime() - real_start_time << "s.\n";
}

template <typename MappingRecord>
size_t Chromap<MappingRecord>::FindBestMappingIndexFromDuplicates(
    const std::vector<MappingRecord> &duplicates) {
  // Find the best barcode, break ties first by the number of the
  // barcodes in the dups, then by the barcode abundance.
  size_t best_mapping_index = 0;

  khiter_t barcode_whitelist_lookup_table_iterator =
      kh_get(k64_seq, barcode_whitelist_lookup_table_,
             duplicates[best_mapping_index].GetBarcode());

  double best_mapping_barcode_abundance =
      kh_value(barcode_whitelist_lookup_table_,
               barcode_whitelist_lookup_table_iterator) /
      (double)num_sample_barcodes_;

  for (size_t bulk_dup_i = 1; bulk_dup_i < duplicates.size(); ++bulk_dup_i) {
    barcode_whitelist_lookup_table_iterator =
        kh_get(k64_seq, barcode_whitelist_lookup_table_,
               duplicates[bulk_dup_i].GetBarcode());

    const double current_mapping_barcode_abundance =
        kh_value(barcode_whitelist_lookup_table_,
                 barcode_whitelist_lookup_table_iterator) /
        (double)num_sample_barcodes_;

    const bool same_num_dups_with_higer_barcode_abundance =
        duplicates[bulk_dup_i].num_dups_ ==
            duplicates[best_mapping_index].num_dups_ &&
        current_mapping_barcode_abundance > best_mapping_barcode_abundance;

    if (duplicates[bulk_dup_i].num_dups_ >
            duplicates[best_mapping_index].num_dups_ ||
        same_num_dups_with_higer_barcode_abundance) {
      best_mapping_index = bulk_dup_i;
      best_mapping_barcode_abundance = current_mapping_barcode_abundance;
    }
  }
  return best_mapping_index;
}

template <typename MappingRecord>
void Chromap<MappingRecord>::ProcessAndOutputMappingsInLowMemory(
    uint32_t num_mappings_in_mem, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const MappingProcessor<MappingRecord> &mapping_processor) {
  // First, process the mappings in the memory and save them on disk.
  if (num_mappings_in_mem > 0) {
    OutputTempMappings(num_reference_sequences, mapping_processor);
    num_mappings_in_mem = 0;
  }

  if (temp_mapping_file_handles_.size() == 0) {
    return;
  }

  double sort_and_dedupe_start_time = GetRealTime();

  // Calculate block size and initialize
  uint64_t max_mem_size = 10 * ((uint64_t)1 << 30);
  if (mapping_output_format_ == MAPPINGFORMAT_SAM ||
      mapping_output_format_ == MAPPINGFORMAT_PAIRS ||
      mapping_output_format_ == MAPPINGFORMAT_PAF) {
    max_mem_size = (uint64_t)1 << 30;
  }
  for (size_t hi = 0; hi < temp_mapping_file_handles_.size(); ++hi) {
    const uint32_t temp_mapping_block_size = max_mem_size /
                                             temp_mapping_file_handles_.size() /
                                             sizeof(MappingRecord);

    temp_mapping_file_handles_[hi].InitializeTempMappingLoading(
        temp_mapping_block_size);
    temp_mapping_file_handles_[hi].LoadTempMappingBlock(
        num_reference_sequences);
  }

  // Merge and dedupe.
  uint32_t last_rid = std::numeric_limits<uint32_t>::max();
  MappingRecord last_mapping = MappingRecord();
  uint32_t num_last_mapping_dups = 0;
  uint64_t num_uni_mappings = 0;
  uint64_t num_multi_mappings = 0;
  uint64_t num_mappings_passing_filters = 0;
  uint64_t num_total_mappings = 0;
  std::vector<MappingRecord> temp_dups_for_bulk_level_dedup;
  temp_dups_for_bulk_level_dedup.reserve(255);

  const bool deduplicate_at_bulk_level_for_single_cell_data =
      remove_pcr_duplicates_ && !is_bulk_data_ &&
      remove_pcr_duplicates_at_bulk_level_;

  while (true) {
    // Merge, dedupe and output.
    // Find min first (sorted by rid and then barcode and then positions).
    size_t min_handle_index = temp_mapping_file_handles_.size();
    uint32_t min_rid = std::numeric_limits<uint32_t>::max();

    for (size_t hi = 0; hi < temp_mapping_file_handles_.size(); ++hi) {
      const TempMappingFileHandle<MappingRecord> &current_handle =
          temp_mapping_file_handles_[hi];
      if (current_handle.HasMappings()) {
        const bool rid_is_smaller = current_handle.current_rid < min_rid;
        const bool same_rid_smaller_mapping =
            current_handle.current_rid == min_rid &&
            current_handle.GetCurrentMapping() <
                temp_mapping_file_handles_[min_handle_index]
                    .GetCurrentMapping();

        if (rid_is_smaller || same_rid_smaller_mapping) {
          min_handle_index = hi;
          min_rid = current_handle.current_rid;
        }
      }
    }

    // All mappings are merged. We only have to handle the case when the last
    // mapping is a duplicate.
    if (min_handle_index == temp_mapping_file_handles_.size()) {
      break;
    }

    ++num_total_mappings;

    // Output the current min mapping if it is not a duplicate.
    const MappingRecord &current_min_mapping =
        temp_mapping_file_handles_[min_handle_index].GetCurrentMapping();

    const bool is_first_iteration = num_total_mappings == 1;
    const bool current_mapping_is_duplicated_at_cell_level =
        !is_first_iteration && current_min_mapping == last_mapping;
    const bool current_mapping_is_duplicated_at_bulk_level =
        !is_first_iteration && deduplicate_at_bulk_level_for_single_cell_data &&
        current_min_mapping.IsSamePosition(last_mapping);
    const bool current_mapping_is_duplicated =
        last_rid == min_rid && (current_mapping_is_duplicated_at_cell_level ||
                                current_mapping_is_duplicated_at_bulk_level);

    if (remove_pcr_duplicates_ && current_mapping_is_duplicated) {
      ++num_last_mapping_dups;
      if (deduplicate_at_bulk_level_for_single_cell_data) {
        if (!temp_dups_for_bulk_level_dedup.empty() &&
            current_min_mapping == temp_dups_for_bulk_level_dedup.back()) {
          // Merge if their barcodes are the same. Be careful of "==" here!
          temp_dups_for_bulk_level_dedup.back() = current_min_mapping;
          temp_dups_for_bulk_level_dedup.back().num_dups_ += 1;
        } else {
          temp_dups_for_bulk_level_dedup.push_back(current_min_mapping);
          temp_dups_for_bulk_level_dedup.back().num_dups_ = 1;
        }
      }
    } else {
      if (!is_first_iteration) {
        if (deduplicate_at_bulk_level_for_single_cell_data) {
          size_t best_mapping_index = FindBestMappingIndexFromDuplicates(
              temp_dups_for_bulk_level_dedup);
          last_mapping = temp_dups_for_bulk_level_dedup[best_mapping_index];

          temp_dups_for_bulk_level_dedup.clear();
        }

        if (last_mapping.mapq_ >= mapq_threshold_) {
          last_mapping.num_dups_ =
              std::min((uint32_t)std::numeric_limits<uint8_t>::max(),
                       num_last_mapping_dups);
          if (Tn5_shift_) {
            last_mapping.Tn5Shift();
          }
          mapping_writer_.AppendMapping(last_rid, reference, last_mapping);
          ++num_mappings_passing_filters;
        }

        if (last_mapping.is_unique_ == 1) {
          ++num_uni_mappings;
        } else {
          ++num_multi_mappings;
        }
      }

      last_mapping = current_min_mapping;
      last_rid = min_rid;
      num_last_mapping_dups = 1;

      if (deduplicate_at_bulk_level_for_single_cell_data) {
        temp_dups_for_bulk_level_dedup.push_back(current_min_mapping);
        temp_dups_for_bulk_level_dedup.back().num_dups_ = 1;
      }
    }

    temp_mapping_file_handles_[min_handle_index].Next(num_reference_sequences);
  }

  if (last_mapping.mapq_ >= mapq_threshold_) {
    if (deduplicate_at_bulk_level_for_single_cell_data) {
      size_t best_mapping_index =
          FindBestMappingIndexFromDuplicates(temp_dups_for_bulk_level_dedup);
      last_mapping = temp_dups_for_bulk_level_dedup[best_mapping_index];

      temp_dups_for_bulk_level_dedup.clear();
    }

    last_mapping.num_dups_ = std::min(
        (uint32_t)std::numeric_limits<uint8_t>::max(), num_last_mapping_dups);
    if (Tn5_shift_) {
      last_mapping.Tn5Shift();
    }
    mapping_writer_.AppendMapping(last_rid, reference, last_mapping);
    ++num_mappings_passing_filters;
  }

  if (last_mapping.is_unique_ == 1) {
    ++num_uni_mappings;
  } else {
    ++num_multi_mappings;
  }

  // Delete temp files.
  for (size_t hi = 0; hi < temp_mapping_file_handles_.size(); ++hi) {
    temp_mapping_file_handles_[hi].FinalizeTempMappingLoading();
    remove(temp_mapping_file_handles_[hi].file_path.c_str());
  }

  if (remove_pcr_duplicates_) {
    std::cerr << "Sorted, deduped and outputed mappings in "
              << GetRealTime() - sort_and_dedupe_start_time << "s.\n";
  } else {
    std::cerr << "Sorted and outputed mappings in "
              << GetRealTime() - sort_and_dedupe_start_time << "s.\n";
  }
  std::cerr << "# uni-mappings: " << num_uni_mappings
            << ", # multi-mappings: " << num_multi_mappings
            << ", total: " << num_uni_mappings + num_multi_mappings << ".\n";
  std::cerr << "Number of output mappings (passed filters): "
            << num_mappings_passing_filters << "\n";
}

template <typename MappingRecord>
void Chromap<MappingRecord>::GenerateCustomizedRidRank(
    const std::string rid_order_path, const SequenceBatch &reference,
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
    chromap::ExitWithMessage(
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
  double real_start_time = GetRealTime();

  // Load reference
  SequenceBatch reference;
  reference.InitializeLoading(reference_file_path_);
  uint32_t num_reference_sequences = reference.LoadAllSequences();
  if (custom_rid_order_path_.length() > 0) {
    GenerateCustomizedRidRank(custom_rid_order_path_, reference,
                              custom_rid_rank_);
    reference.ReorderSequences(custom_rid_rank_);
  }
  if (mapping_output_format_ == MAPPINGFORMAT_PAIRS) {
    GenerateCustomizedRidRank(pairs_custom_rid_order_path_, reference,
                              pairs_custom_rid_rank_);
  }

  // Load index
  Index index(index_file_path_);
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

  read_batch1.SetSeqEffectiveRange(read1_format_[0], read1_format_[1], 1);
  read_batch2.SetSeqEffectiveRange(read2_format_[0], read2_format_[1], 1);
  barcode_batch.SetSeqEffectiveRange(barcode_format_[0], barcode_format_[1],
                                     barcode_format_[2]);
  read_batch1_for_loading.SetSeqEffectiveRange(read1_format_[0],
                                               read1_format_[1], 1);
  read_batch2_for_loading.SetSeqEffectiveRange(read2_format_[0],
                                               read2_format_[1], 1);
  barcode_batch_for_loading.SetSeqEffectiveRange(
      barcode_format_[0], barcode_format_[1], barcode_format_[2]);

  // Initialize cache
  mm_cache mm_to_candidates_cache(2000003);
  mm_to_candidates_cache.SetKmerLength(kmer_size_);
  struct _mm_history *mm_history1 = new struct _mm_history[read_batch_size_];
  struct _mm_history *mm_history2 = new struct _mm_history[read_batch_size_];

  // Initialize mapping container
  mappings_on_diff_ref_seqs_.reserve(num_reference_sequences);
  for (uint32_t i = 0; i < num_reference_sequences; ++i) {
    mappings_on_diff_ref_seqs_.emplace_back(std::vector<MappingRecord>());
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
    mapping_writer_.SetPairsCustomRidRank(pairs_custom_rid_rank_);
  }

  mapping_writer_.InitializeMappingOutput(
      barcode_length_, mapping_output_file_path_, mapping_output_format_);
  mapping_writer_.OutputHeader(num_reference_sequences, reference);

  CandidateProcessor candidate_processor(min_num_seeds_required_for_mapping_,
                                         max_seed_frequencies_);

  MappingProcessor<MappingRecord> mapping_processor(
      min_unique_mapping_mapq_, multi_mapping_allocation_seed_,
      multi_mapping_allocation_distance_, max_num_best_mappings_);

  MappingGenerator<MappingRecord> mapping_generator(
      error_threshold_, NUM_VPU_LANES_, match_score_, mismatch_penalty_,
      gap_open_penalties_, gap_extension_penalties_, max_num_best_mappings_,
      max_insert_size_, min_read_length_, drop_repetitive_reads_, is_bulk_data_,
      split_alignment_, mapping_output_format_, pairs_custom_rid_rank_);

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
  double real_start_mapping_time = GetRealTime();
  for (size_t read_file_index = 0; read_file_index < read_file1_paths_.size();
       ++read_file_index) {
    // Set read batches to the current read files.
    read_batch1_for_loading.InitializeLoading(
        read_file1_paths_[read_file_index]);
    read_batch2_for_loading.InitializeLoading(
        read_file2_paths_[read_file_index]);
    if (!is_bulk_data_) {
      barcode_batch_for_loading.InitializeLoading(
          barcode_file_paths_[read_file_index]);
    }

    // Load the first batches.
    uint32_t num_loaded_pairs_for_loading = 0;
    uint32_t num_loaded_pairs = LoadPairedEndReadsWithBarcodes(
        read_batch1_for_loading, read_batch2_for_loading,
        barcode_batch_for_loading);
    read_batch1_for_loading.SwapSequenceBatch(read_batch1);
    read_batch2_for_loading.SwapSequenceBatch(read_batch2);
    if (!is_bulk_data_) {
      barcode_batch_for_loading.SwapSequenceBatch(barcode_batch);
      // TODO(Haowen): simplify this condition check.
      if (num_loaded_pairs > 0) {
        mapping_writer_.SetBarcodeLength(barcode_batch.GetSequenceLengthAt(0));
      }
    }

    // Setup thread private vectors to save mapping results.
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

#pragma omp parallel default(none) shared(output_mappings_not_in_whitelist_, barcode_whitelist_file_path_, trim_adapters_, window_size_, custom_rid_order_path_, split_alignment_, error_threshold_, max_seed_frequencies_, max_num_best_mappings_, max_insert_size_, num_threads_, num_reads_, low_memory_mode_, num_reference_sequences, reference, index, read_batch1, read_batch2, barcode_batch, read_batch1_for_loading, read_batch2_for_loading, barcode_batch_for_loading, candidate_processor, mapping_processor, mapping_generator, mapping_writer_, std::cerr, num_loaded_pairs_for_loading, num_loaded_pairs, mappings_on_diff_ref_seqs_for_diff_threads, mappings_on_diff_ref_seqs_for_diff_threads_for_saving, mappings_on_diff_ref_seqs_, mapping_output_file_path_, num_mappings_in_mem, max_num_mappings_in_mem, temp_mapping_file_handles_, mm_to_candidates_cache, mm_history1, mm_history2) num_threads(num_threads_) reduction(+:num_candidates_, num_mappings_, num_mapped_reads_, num_uniquely_mapped_reads_, num_barcode_in_whitelist_, num_corrected_barcode_)
    {
      thread_num_candidates = 0;
      thread_num_mappings = 0;
      thread_num_mapped_reads = 0;
      thread_num_uniquely_mapped_reads = 0;
      thread_num_barcode_in_whitelist = 0;
      thread_num_corrected_barcode = 0;
      PairedEndMappingMetadata paired_end_mapping_metadata;

      std::vector<int> best_mapping_indices(max_num_best_mappings_);
      std::mt19937 generator(11);
#pragma omp single
      {
        double real_batch_start_time = GetRealTime();
        while (num_loaded_pairs > 0) {
          num_reads_ += num_loaded_pairs;
          num_reads_ += num_loaded_pairs;

#pragma omp task
          {
            num_loaded_pairs_for_loading = LoadPairedEndReadsWithBarcodes(
                read_batch1_for_loading, read_batch2_for_loading,
                barcode_batch_for_loading);
          }  // end of openmp loading task

          int grain_size = 5000;
#pragma omp taskloop grainsize(grain_size)
          for (uint32_t pair_index = 0; pair_index < num_loaded_pairs;
               ++pair_index) {
            bool current_barcode_is_whitelisted = true;
            if (!barcode_whitelist_file_path_.empty()) {
              current_barcode_is_whitelisted = CorrectBarcodeAt(
                  pair_index, barcode_batch, thread_num_barcode_in_whitelist,
                  thread_num_corrected_barcode);
            }

            if (current_barcode_is_whitelisted ||
                output_mappings_not_in_whitelist_) {
              read_batch1.PrepareNegativeSequenceAt(pair_index);
              read_batch2.PrepareNegativeSequenceAt(pair_index);

              if (trim_adapters_) {
                TrimAdapterForPairedEndRead(pair_index, read_batch1,
                                            read_batch2);
              }

              paired_end_mapping_metadata.PreparedForMappingNextReadPair(
                  max_seed_frequencies_[0]);

              index.GenerateMinimizerSketch(
                  read_batch1, pair_index,
                  paired_end_mapping_metadata.mapping_metadata1_.minimizers_);
              index.GenerateMinimizerSketch(
                  read_batch2, pair_index,
                  paired_end_mapping_metadata.mapping_metadata2_.minimizers_);

              if (paired_end_mapping_metadata.BothEndsHaveMinimizers()) {
                // Generate candidates
                if (mm_to_candidates_cache.Query(
                        paired_end_mapping_metadata.mapping_metadata1_,
                        read_batch1.GetSequenceLengthAt(pair_index)) == -1) {
                  candidate_processor.GenerateCandidates(
                      error_threshold_, index,
                      paired_end_mapping_metadata.mapping_metadata1_);
                }

                size_t current_num_candidates1 =
                    paired_end_mapping_metadata.mapping_metadata1_
                        .GetNumCandidates();

                if (mm_to_candidates_cache.Query(
                        paired_end_mapping_metadata.mapping_metadata2_,
                        read_batch2.GetSequenceLengthAt(pair_index)) == -1) {
                  candidate_processor.GenerateCandidates(
                      error_threshold_, index,
                      paired_end_mapping_metadata.mapping_metadata2_);
                }

                size_t current_num_candidates2 =
                    paired_end_mapping_metadata.mapping_metadata2_
                        .GetNumCandidates();

                if (pair_index < num_loaded_pairs &&
                    (pair_index < num_loaded_pairs / num_threads_ ||
                     num_reads_ <= 5000000)) {
                  mm_history1[pair_index].timestamp =
                      mm_history2[pair_index].timestamp = num_reads_;
                  mm_history1[pair_index].minimizers =
                      paired_end_mapping_metadata.mapping_metadata1_
                          .minimizers_;
                  mm_history1[pair_index].positive_candidates =
                      paired_end_mapping_metadata.mapping_metadata1_
                          .positive_candidates_;
                  mm_history1[pair_index].negative_candidates =
                      paired_end_mapping_metadata.mapping_metadata1_
                          .negative_candidates_;
                  mm_history1[pair_index].repetitive_seed_length =
                      paired_end_mapping_metadata.mapping_metadata1_
                          .repetitive_seed_length_;
                  mm_history2[pair_index].minimizers =
                      paired_end_mapping_metadata.mapping_metadata2_
                          .minimizers_;
                  mm_history2[pair_index].positive_candidates =
                      paired_end_mapping_metadata.mapping_metadata2_
                          .positive_candidates_;
                  mm_history2[pair_index].negative_candidates =
                      paired_end_mapping_metadata.mapping_metadata2_
                          .negative_candidates_;
                  mm_history2[pair_index].repetitive_seed_length =
                      paired_end_mapping_metadata.mapping_metadata2_
                          .repetitive_seed_length_;
                }

                // Test whether we need to augment the candidate list with mate
                // information.
                int supplementCandidateResult = 0;
                if (!split_alignment_) {
                  supplementCandidateResult =
                      candidate_processor.SupplementCandidates(
                          error_threshold_,
                          /*search_range=*/2 * max_insert_size_, index,
                          paired_end_mapping_metadata);
                  current_num_candidates1 =
                      paired_end_mapping_metadata.mapping_metadata1_
                          .GetNumCandidates();
                  current_num_candidates2 =
                      paired_end_mapping_metadata.mapping_metadata2_
                          .GetNumCandidates();
                }

                if (current_num_candidates1 > 0 &&
                    current_num_candidates2 > 0 && !split_alignment_) {
                  paired_end_mapping_metadata.MoveCandidiatesToBuffer();

                  // Paired-end filter
                  candidate_processor.ReduceCandidatesForPairedEndRead(
                      max_insert_size_, paired_end_mapping_metadata);

                  current_num_candidates1 =
                      paired_end_mapping_metadata.mapping_metadata1_
                          .GetNumCandidates();
                  current_num_candidates2 =
                      paired_end_mapping_metadata.mapping_metadata2_
                          .GetNumCandidates();
                }

                // Verify candidates
                if (current_num_candidates1 > 0 &&
                    current_num_candidates2 > 0) {
                  thread_num_candidates +=
                      current_num_candidates1 + current_num_candidates2;

                  if (custom_rid_order_path_.length() > 0) {
                    RerankCandidatesRid(
                        paired_end_mapping_metadata.mapping_metadata1_
                            .positive_candidates_);
                    RerankCandidatesRid(
                        paired_end_mapping_metadata.mapping_metadata1_
                            .negative_candidates_);
                    RerankCandidatesRid(
                        paired_end_mapping_metadata.mapping_metadata2_
                            .positive_candidates_);
                    RerankCandidatesRid(
                        paired_end_mapping_metadata.mapping_metadata2_
                            .negative_candidates_);
                  }

                  VerifyCandidates(
                      read_batch1, pair_index, reference,
                      paired_end_mapping_metadata.mapping_metadata1_);

                  size_t current_num_mappings1 =
                      paired_end_mapping_metadata.mapping_metadata1_
                          .GetNumMappings();

                  VerifyCandidates(
                      read_batch2, pair_index, reference,
                      paired_end_mapping_metadata.mapping_metadata2_);

                  size_t current_num_mappings2 =
                      paired_end_mapping_metadata.mapping_metadata2_
                          .GetNumMappings();

                  if (current_num_mappings1 > 0 && current_num_mappings2 > 0) {
                    std::vector<std::vector<MappingRecord>>
                        &mappings_on_diff_ref_seqs =
                            mappings_on_diff_ref_seqs_for_diff_threads
                                [omp_get_thread_num()];

                    if (!split_alignment_) {
                      // GenerateBestMappingsForPairedEndRead assumes the
                      // mappings are sorted by coordinate for non split
                      // alignments. In split alignment, we don't want to sort
                      // and this keeps mapping and split_sites vectors
                      // consistent.
                      paired_end_mapping_metadata.SortMappingsByPositions();
                    }

                    int force_mapq = -1;
                    if (supplementCandidateResult != 0) {
                      force_mapq = 0;
                    }

                    mapping_generator.GenerateBestMappingsForPairedEndRead(
                        pair_index, read_batch1, read_batch2, barcode_batch,
                        reference, best_mapping_indices, generator, force_mapq,
                        paired_end_mapping_metadata, mappings_on_diff_ref_seqs);

                    if (paired_end_mapping_metadata.GetNumBestMappings() == 1) {
                      ++thread_num_uniquely_mapped_reads;
                      ++thread_num_uniquely_mapped_reads;
                    }

                    thread_num_mappings += std::min(
                        paired_end_mapping_metadata.GetNumBestMappings(),
                        max_num_best_mappings_);
                    thread_num_mappings += std::min(
                        paired_end_mapping_metadata.GetNumBestMappings(),
                        max_num_best_mappings_);
                    if (paired_end_mapping_metadata.GetNumBestMappings() > 0) {
                      ++thread_num_mapped_reads;
                      ++thread_num_mapped_reads;
                    }
                  }
                }  // verify candidate
              }
            }
          }  // end of for pair_index

          // if (num_reads_ / 2 > initial_num_sample_barcodes_) {
          //  if (!is_bulk_data_) {
          //    if (!barcode_whitelist_file_path_.empty()) {
          //      UpdateBarcodeAbundance(num_loaded_pairs, barcode_batch);
          //    }
          //  }
          //}
#pragma omp taskwait
          // Update cache
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
                    << GetRealTime() - real_batch_start_time << "s.\n";
          real_batch_start_time = GetRealTime();

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
                mappings_on_diff_ref_seqs_for_diff_threads_for_saving);
            if (low_memory_mode_ &&
                num_mappings_in_mem > max_num_mappings_in_mem) {
              OutputTempMappings(num_reference_sequences, mapping_processor);
              num_mappings_in_mem = 0;
            }
          }  // end of omp task to handle output
        }    // end of while num_loaded_pairs
      }      // end of openmp single

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
  }  // end of for read_file_index

  std::cerr << "Mapped all reads in " << GetRealTime() - real_start_mapping_time
            << "s.\n";

  delete[] mm_history1;
  delete[] mm_history2;

  OutputMappingStatistics();
  if (!is_bulk_data_) {
    OutputBarcodeStatistics();
  }

  index.Destroy();

  if (low_memory_mode_) {
    ProcessAndOutputMappingsInLowMemory(num_mappings_in_mem,
                                        num_reference_sequences, reference,
                                        mapping_processor);
  } else {
    // OutputMappingStatistics(num_reference_sequences,
    // mappings_on_diff_ref_seqs_, mappings_on_diff_ref_seqs_);
    if (Tn5_shift_) {
      ApplyTn5ShiftOnPairedEndMapping(num_reference_sequences,
                                      mappings_on_diff_ref_seqs_);
    }

    if (remove_pcr_duplicates_) {
      mapping_processor.RemovePCRDuplicate(num_reference_sequences,
                                           mappings_on_diff_ref_seqs_);
      std::cerr << "After removing PCR duplications, ";
      OutputMappingStatistics(num_reference_sequences,
                              mappings_on_diff_ref_seqs_,
                              mappings_on_diff_ref_seqs_);
    } else {
      mapping_processor.SortOutputMappings(num_reference_sequences,
                                           mappings_on_diff_ref_seqs_);
    }

    if (allocate_multi_mappings_) {
      const uint64_t num_multi_mappings =
          num_mapped_reads_ - num_uniquely_mapped_reads_;
      mapping_processor.AllocateMultiMappings(
          num_reference_sequences, num_multi_mappings,
          multi_mapping_allocation_distance_, mappings_on_diff_ref_seqs_);
      std::cerr << "After allocating multi-mappings, ";
      OutputMappingStatistics(num_reference_sequences,
                              mappings_on_diff_ref_seqs_,
                              mappings_on_diff_ref_seqs_);
      mapping_processor.SortOutputMappings(num_reference_sequences,
                                           mappings_on_diff_ref_seqs_);
    }
    OutputMappings(num_reference_sequences, reference,
                   mappings_on_diff_ref_seqs_);

    // Temporarily disable feature matrix output. Do not delete the following
    // commented code.
    // if (!is_bulk_data_ && !matrix_output_prefix_.empty()) {
    //   if constexpr (std::is_same<MappingRecord,
    //                             PairedEndMappingWithBarcode>::value) {
    //    FeatureBarcodeMatrix feature_barcode_matrix(
    //        cell_by_bin_, bin_size_, multi_mapping_allocation_distance_,
    //        depth_cutoff_to_call_peak_);
    //    std::vector<std::vector<PairedEndMappingWithBarcode>> &mappings =
    //        allocate_multi_mappings_
    //            ? allocated_mappings_on_diff_ref_seqs_
    //            : (remove_pcr_duplicates_ ? deduped_mappings_on_diff_ref_seqs_
    //                                      : mappings_on_diff_ref_seqs_);

    //    feature_barcode_matrix.OutputFeatureMatrix(num_reference_sequences,
    //                                               reference, mappings,
    //                                               matrix_output_prefix_);
    //  }
    //}
  }

  mapping_writer_.FinalizeMappingOutput();

  reference.FinalizeLoading();

  std::cerr << "Total time: " << GetRealTime() - real_start_time << "s.\n";
}

template <typename MappingRecord>
void Chromap<MappingRecord>::OutputTempMappings(
    uint32_t num_reference_sequences,
    const MappingProcessor<MappingRecord> &mapping_processor) {
  TempMappingFileHandle<MappingRecord> temp_mapping_file_handle;
  temp_mapping_file_handle.file_path =
      mapping_output_file_path_ + ".temp" +
      std::to_string(temp_mapping_file_handles_.size());
  temp_mapping_file_handles_.emplace_back(temp_mapping_file_handle);

  mapping_processor.SortOutputMappings(num_reference_sequences,
                                       mappings_on_diff_ref_seqs_);

  mapping_writer_.OutputTempMapping(temp_mapping_file_handle.file_path,
                                    num_reference_sequences,
                                    mappings_on_diff_ref_seqs_);

  for (uint32_t i = 0; i < num_reference_sequences; ++i) {
    mappings_on_diff_ref_seqs_[i].clear();
  }
}

template <typename MappingRecord>
void Chromap<MappingRecord>::OutputMappingsInVector(
    uint8_t mapq_threshold, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const std::vector<std::vector<MappingRecord>> &mappings) {
  uint64_t num_mappings_passing_filters = 0;
  for (uint32_t ri = 0; ri < num_reference_sequences; ++ri) {
    for (auto it = mappings[ri].begin(); it != mappings[ri].end(); ++it) {
      uint8_t mapq = (it->mapq_);
      // uint8_t is_unique = (it->is_unique);
      if (mapq >= mapq_threshold) {
        // if (allocate_multi_mappings_ || (only_output_unique_mappings_ &&
        // is_unique == 1)) {
        mapping_writer_.AppendMapping(ri, reference, *it);
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
void Chromap<MappingRecord>::ApplyTn5ShiftOnPairedEndMapping(
    uint32_t num_reference_sequences,
    std::vector<std::vector<MappingRecord>> &mappings) {
  uint64_t num_shifted_mappings = 0;
  for (auto &mappings_on_one_ref_seq : mappings) {
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
  double real_start_time = GetRealTime();

  SequenceBatch reference;
  reference.InitializeLoading(reference_file_path_);
  uint32_t num_reference_sequences = reference.LoadAllSequences();
  if (custom_rid_order_path_.length() > 0) {
    GenerateCustomizedRidRank(custom_rid_order_path_, reference,
                              custom_rid_rank_);
    reference.ReorderSequences(custom_rid_rank_);
  }

  Index index(index_file_path_);
  index.Load();
  kmer_size_ = index.GetKmerSize();
  window_size_ = index.GetWindowSize();
  // index.Statistics(num_sequences, reference);

  SequenceBatch read_batch(read_batch_size_);
  SequenceBatch read_batch_for_loading(read_batch_size_);
  SequenceBatch barcode_batch(read_batch_size_);
  SequenceBatch barcode_batch_for_loading(read_batch_size_);
  read_batch.SetSeqEffectiveRange(read1_format_[0], read1_format_[1], 1);
  read_batch_for_loading.SetSeqEffectiveRange(read1_format_[0],
                                              read1_format_[1], 1);
  barcode_batch.SetSeqEffectiveRange(barcode_format_[0], barcode_format_[1],
                                     barcode_format_[2]);
  barcode_batch_for_loading.SetSeqEffectiveRange(
      barcode_format_[0], barcode_format_[1], barcode_format_[2]);

  mappings_on_diff_ref_seqs_.reserve(num_reference_sequences);
  for (uint32_t i = 0; i < num_reference_sequences; ++i) {
    mappings_on_diff_ref_seqs_.emplace_back(std::vector<MappingRecord>());
  }

  // Preprocess barcodes for single cell data
  if (!is_bulk_data_) {
    if (!barcode_whitelist_file_path_.empty()) {
      LoadBarcodeWhitelist();
      ComputeBarcodeAbundance(initial_num_sample_barcodes_);
    }
  }

  CandidateProcessor candidate_processor(min_num_seeds_required_for_mapping_,
                                         max_seed_frequencies_);

  MappingProcessor<MappingRecord> mapping_processor(
      min_unique_mapping_mapq_, multi_mapping_allocation_seed_,
      multi_mapping_allocation_distance_, max_num_best_mappings_);

  MappingGenerator<MappingRecord> mapping_generator(
      error_threshold_, NUM_VPU_LANES_, match_score_, mismatch_penalty_,
      gap_open_penalties_, gap_extension_penalties_, max_num_best_mappings_,
      max_insert_size_, min_read_length_, drop_repetitive_reads_, is_bulk_data_,
      split_alignment_, mapping_output_format_, pairs_custom_rid_rank_);

  mapping_writer_.InitializeMappingOutput(
      barcode_length_, mapping_output_file_path_, mapping_output_format_);
  mapping_writer_.OutputHeader(num_reference_sequences, reference);

  uint32_t num_mappings_in_mem = 0;
  uint64_t max_num_mappings_in_mem =
      1 * ((uint64_t)1 << 30) / sizeof(MappingRecord);
  if (mapping_output_format_ == MAPPINGFORMAT_SAM ||
      mapping_output_format_ == MAPPINGFORMAT_PAF ||
      mapping_output_format_ == MAPPINGFORMAT_PAIRS) {
    max_num_mappings_in_mem = 1 * ((uint64_t)1 << 29) / sizeof(MappingRecord);
  }

  mm_cache mm_to_candidates_cache(2000003);
  mm_to_candidates_cache.SetKmerLength(kmer_size_);
  struct _mm_history *mm_history = new struct _mm_history[read_batch_size_];
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
  double real_start_mapping_time = GetRealTime();
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

    if (!is_bulk_data_) {
      barcode_batch_for_loading.SwapSequenceBatch(barcode_batch);
      // TODO(Haowen): simplify this condition check.
      if (num_loaded_reads > 0) {
        mapping_writer_.SetBarcodeLength(barcode_batch.GetSequenceLengthAt(0));
      }
    }

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
#pragma omp parallel default(none) shared(barcode_whitelist_file_path_, output_mappings_not_in_whitelist_, window_size_, custom_rid_order_path_, error_threshold_, max_num_best_mappings_, max_seed_frequencies_, num_threads_, num_reads_, mm_history, reference, index, read_batch, barcode_batch, read_batch_for_loading, barcode_batch_for_loading, std::cerr, num_loaded_reads_for_loading, num_loaded_reads, num_reference_sequences, mappings_on_diff_ref_seqs_for_diff_threads, mappings_on_diff_ref_seqs_for_diff_threads_for_saving, mm_to_candidates_cache, mapping_writer_, candidate_processor, mapping_processor, mapping_generator, num_mappings_in_mem, max_num_mappings_in_mem) num_threads(num_threads_) reduction(+:num_candidates_, num_mappings_, num_mapped_reads_, num_uniquely_mapped_reads_, num_barcode_in_whitelist_, num_corrected_barcode_)
    {
      thread_num_candidates = 0;
      thread_num_mappings = 0;
      thread_num_mapped_reads = 0;
      thread_num_uniquely_mapped_reads = 0;
      thread_num_barcode_in_whitelist = 0;
      thread_num_corrected_barcode = 0;
      MappingMetadata mapping_metadata;
#pragma omp single
      {
        while (num_loaded_reads > 0) {
          double real_batch_start_time = GetRealTime();
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
            bool current_barcode_is_whitelisted = true;
            if (!barcode_whitelist_file_path_.empty()) {
              current_barcode_is_whitelisted = CorrectBarcodeAt(
                  read_index, barcode_batch, thread_num_barcode_in_whitelist,
                  thread_num_corrected_barcode);
            }

            if (!(current_barcode_is_whitelisted ||
                  output_mappings_not_in_whitelist_))
              continue;

            read_batch.PrepareNegativeSequenceAt(read_index);

            mapping_metadata.PrepareForMappingNextRead(
                max_seed_frequencies_[0]);

            index.GenerateMinimizerSketch(read_batch, read_index,
                                          mapping_metadata.minimizers_);

            if (mapping_metadata.minimizers_.size() > 0) {
              if (custom_rid_order_path_.length() > 0) {
                RerankCandidatesRid(mapping_metadata.positive_candidates_);
                RerankCandidatesRid(mapping_metadata.negative_candidates_);
              }

              if (mm_to_candidates_cache.Query(
                      mapping_metadata,
                      read_batch.GetSequenceLengthAt(read_index)) == -1) {
                candidate_processor.GenerateCandidates(error_threshold_, index,
                                                       mapping_metadata);
              }

              if (read_index < num_loaded_reads &&
                  (read_index < num_loaded_reads / num_threads_ ||
                   num_reads_ <= 2500000)) {
                mm_history[read_index].timestamp = num_reads_;
                mm_history[read_index].minimizers =
                    mapping_metadata.minimizers_;
                mm_history[read_index].positive_candidates =
                    mapping_metadata.positive_candidates_;
                mm_history[read_index].negative_candidates =
                    mapping_metadata.negative_candidates_;
                mm_history[read_index].repetitive_seed_length =
                    mapping_metadata.repetitive_seed_length_;
              }

              size_t current_num_candidates =
                  mapping_metadata.GetNumCandidates();
              if (current_num_candidates > 0) {
                thread_num_candidates += current_num_candidates;
                VerifyCandidates(read_batch, read_index, reference,
                                 mapping_metadata);

                size_t current_num_mappings = mapping_metadata.GetNumMappings();
                if (current_num_mappings > 0) {
                  std::vector<std::vector<MappingRecord>>
                      &mappings_on_diff_ref_seqs =
                          mappings_on_diff_ref_seqs_for_diff_threads
                              [omp_get_thread_num()];

                  mapping_generator.GenerateBestMappingsForSingleEndRead(
                      read_batch, read_index, reference, barcode_batch,
                      mapping_metadata, mappings_on_diff_ref_seqs);

                  thread_num_mappings +=
                      std::min(mapping_metadata.GetNumBestMappings(),
                               max_num_best_mappings_);
                  ++thread_num_mapped_reads;

                  if (mapping_metadata.GetNumBestMappings() == 1) {
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
            num_mappings_in_mem += MoveMappingsInBuffersToMappingContainer(
                num_reference_sequences,
                mappings_on_diff_ref_seqs_for_diff_threads_for_saving);
            if (low_memory_mode_ &&
                num_mappings_in_mem > max_num_mappings_in_mem) {
              OutputTempMappings(num_reference_sequences, mapping_processor);
              num_mappings_in_mem = 0;
            }
          }
          std::cerr << "Mapped in " << GetRealTime() - real_batch_start_time
                    << "s.\n";
        }
      }  // end of openmp single
      {
        num_barcode_in_whitelist_ += thread_num_barcode_in_whitelist;
        num_corrected_barcode_ += thread_num_corrected_barcode;
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

  std::cerr << "Mapped all reads in " << GetRealTime() - real_start_mapping_time
            << "s.\n";

  delete[] mm_history;

  OutputMappingStatistics();
  if (!is_bulk_data_) {
    OutputBarcodeStatistics();
  }

  index.Destroy();

  if (low_memory_mode_) {
    ProcessAndOutputMappingsInLowMemory(num_mappings_in_mem,
                                        num_reference_sequences, reference,
                                        mapping_processor);
  } else {
    if (Tn5_shift_) {
      ApplyTn5ShiftOnSingleEndMapping(num_reference_sequences,
                                      &mappings_on_diff_ref_seqs_);
    }

    if (remove_pcr_duplicates_) {
      mapping_processor.RemovePCRDuplicate(num_reference_sequences,
                                           mappings_on_diff_ref_seqs_);
      std::cerr << "After removing PCR duplications, ";
      OutputMappingStatistics(num_reference_sequences,
                              mappings_on_diff_ref_seqs_,
                              mappings_on_diff_ref_seqs_);
    } else {
      mapping_processor.SortOutputMappings(num_reference_sequences,
                                           mappings_on_diff_ref_seqs_);
    }

    if (allocate_multi_mappings_) {
      const uint64_t num_multi_mappings =
          num_mapped_reads_ - num_uniquely_mapped_reads_;
      mapping_processor.AllocateMultiMappings(
          num_reference_sequences, num_multi_mappings,
          multi_mapping_allocation_distance_, mappings_on_diff_ref_seqs_);
      std::cerr << "After allocating multi-mappings, ";
      OutputMappingStatistics(num_reference_sequences,
                              mappings_on_diff_ref_seqs_,
                              mappings_on_diff_ref_seqs_);
      mapping_processor.SortOutputMappings(num_reference_sequences,
                                           mappings_on_diff_ref_seqs_);
    }
    OutputMappings(num_reference_sequences, reference,
                   mappings_on_diff_ref_seqs_);
  }

  mapping_writer_.FinalizeMappingOutput();
  reference.FinalizeLoading();
  std::cerr << "Total time: " << GetRealTime() - real_start_time << "s.\n";
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
  double real_start_time = GetRealTime();
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
      ExitWithMessage("Numbers of reads and barcodes don't match!");
    }
    ++num_loaded_reads;
  }
  if (num_loaded_reads > 0) {
    std::cerr << "Loaded " << num_loaded_reads << " reads in "
              << GetRealTime() - real_start_time << "s.\n";
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
        &mappings_on_diff_ref_seqs_for_diff_threads_for_saving) {
  // double real_start_time = Chromap<>::GetRealTime();
  uint32_t num_moved_mappings = 0;
  for (int ti = 0; ti < num_threads_; ++ti) {
    for (uint32_t i = 0; i < num_reference_sequences; ++i) {
      num_moved_mappings +=
          mappings_on_diff_ref_seqs_for_diff_threads_for_saving[ti][i].size();
      mappings_on_diff_ref_seqs_[i].insert(
          mappings_on_diff_ref_seqs_[i].end(),
          std::make_move_iterator(
              mappings_on_diff_ref_seqs_for_diff_threads_for_saving[ti][i]
                  .begin()),
          std::make_move_iterator(
              mappings_on_diff_ref_seqs_for_diff_threads_for_saving[ti][i]
                  .end()));
      mappings_on_diff_ref_seqs_for_diff_threads_for_saving[ti][i].clear();
    }
  }
  // std::cerr << "Moved mappings in " << Chromap<>::GetRealTime() -
  // real_start_time << "s.\n";
  return num_moved_mappings;
}

template <typename MappingRecord>
void Chromap<MappingRecord>::VerifyCandidatesOnOneDirectionUsingSIMD(
    Direction candidate_direction, const SequenceBatch &read_batch,
    uint32_t read_index, const SequenceBatch &reference,
    const std::vector<Candidate> &candidates,
    std::vector<std::pair<int, uint64_t>> &mappings, int &min_num_errors,
    int &num_best_mappings, int &second_min_num_errors,
    int &num_second_best_mappings) {
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
          BandedAlign8PatternsToText(error_threshold_, valid_candidate_starts,
                                     read, read_length, mapping_edit_distances,
                                     mapping_end_positions);
        } else {
          BandedAlign8PatternsToText(
              error_threshold_, valid_candidate_starts, negative_read.data(),
              read_length, mapping_edit_distances, mapping_end_positions);
        }
        for (int mi = 0; mi < NUM_VPU_LANES_; ++mi) {
          if (mapping_edit_distances[mi] <= error_threshold_) {
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
                                        error_threshold_ +
                                        mapping_end_positions[mi]);
            } else {
              mappings.emplace_back((uint8_t)mapping_edit_distances[mi],
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
          BandedAlign4PatternsToText(error_threshold_, valid_candidate_starts,
                                     read, read_length, mapping_edit_distances,
                                     mapping_end_positions);
        } else {
          BandedAlign4PatternsToText(
              error_threshold_, valid_candidate_starts, negative_read.data(),
              read_length, mapping_edit_distances, mapping_end_positions);
        }
        for (int mi = 0; mi < NUM_VPU_LANES_; ++mi) {
          if (mapping_edit_distances[mi] <= error_threshold_) {
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
                                        error_threshold_ +
                                        mapping_end_positions[mi]);
            } else {
              mappings.emplace_back((uint8_t)mapping_edit_distances[mi],
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
          error_threshold_,
          reference.GetSequenceAt(rid) + position - error_threshold_, read,
          read_length, &mapping_end_position);
    } else {
      num_errors = BandedAlignPatternToText(
          error_threshold_,
          reference.GetSequenceAt(rid) + position - error_threshold_,
          negative_read.data(), read_length, &mapping_end_position);
    }
    if (num_errors <= error_threshold_) {
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
        mappings.emplace_back(num_errors, valid_candidates[ci].position -
                                              error_threshold_ +
                                              mapping_end_position);
      } else {
        mappings.emplace_back(num_errors,
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
    std::vector<std::pair<int, uint64_t>> &mappings,
    std::vector<int> &split_sites, int &min_num_errors, int &num_best_mappings,
    int &second_min_num_errors, int &num_second_best_mappings) {
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
            error_threshold_,
            reference.GetSequenceAt(rid) + position - error_threshold_, read,
            read_length, &mapping_end_position, &read_mapping_length);
        if (mapping_end_position < 0 && allow_gap_beginning > 0) {
          int backup_num_errors = num_errors;
          int backup_mapping_end_position = -mapping_end_position;
          int backup_read_mapping_length = read_mapping_length;
          num_errors = BandedAlignPatternToTextWithDropOff(
              error_threshold_,
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
            error_threshold_,
            reference.GetSequenceAt(rid) + position - error_threshold_,
            negative_read.data(), read_length, &mapping_end_position,
            &read_mapping_length);
        if (mapping_end_position < 0 && allow_gap_beginning > 0) {
          int backup_num_errors = num_errors;
          int backup_mapping_end_position = -mapping_end_position;
          int backup_read_mapping_length = read_mapping_length;
          num_errors = BandedAlignPatternToTextWithDropOffFrom3End(
              error_threshold_,
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
            error_threshold_,
            reference.GetSequenceAt(rid) + position - error_threshold_, read,
            read_length, &mapping_end_position);
      } else {
        num_errors = BandedAlignPatternToText(
            error_threshold_,
            reference.GetSequenceAt(rid) + position - error_threshold_,
            negative_read.data(), read_length, &mapping_end_position);
      }
    }

    // std::cerr << "ne3: " << num_errors << " " << mapping_end_position << " "
    // << actual_num_errors << " "<< reference.GetSequenceNameAt(rid) <<" " <<
    // (int)candidates[ci].position << " " << position<<"\n";
    if (num_errors <= error_threshold_) {
      if (num_errors < min_num_errors) {
        second_min_num_errors = min_num_errors;
        num_second_best_mappings = num_best_mappings;
        min_num_errors = num_errors;
        num_best_mappings = 1;
        if (split_alignment_) {
          if (candidates.size() > 50) {
            candidate_count_threshold = candidates[ci].count;
          } else {
            candidate_count_threshold = candidates[ci].count / 2;
          }
          if (second_min_num_errors < min_num_errors + error_threshold_ / 2 &&
              best_mapping_longest_match > longest_match &&
              candidates.size() > 200) {
            second_min_num_errors = min_num_errors;
          }
        }
        best_mapping_longest_match = longest_match;
      } else if (num_errors == min_num_errors) {
        num_best_mappings++;
        /*if (split_alignment_ && candidates.size() > 50) {
                candidate_count_threshold = candidates[ci].count + 1;
        }*/
      } else if (num_errors == second_min_num_errors) {
        num_second_best_mappings++;
      } else if (num_errors < second_min_num_errors) {
        num_second_best_mappings = 1;
        second_min_num_errors = num_errors;
      }
      if (candidate_direction == kPositive) {
        mappings.emplace_back(
            num_errors,
            candidates[ci].position - error_threshold_ + mapping_end_position);
      } else {
        if (split_alignment_ && mapping_output_format_ != MAPPINGFORMAT_SAM) {
          // mappings->emplace_back(num_errors, candidates[ci].position +
          // error_threshold_ - 1 - mapping_end_position
          //					+ read_mapping_length - 1 -
          // gap_beginning);
          mappings.emplace_back(num_errors,
                                candidates[ci].position - gap_beginning);
        } else {
          // Need to minus gap_beginning because mapping_end_position is
          // adjusted by it, but read_length is not.
          // printf("%d %d %d\n", candidates[ci].position, mapping_end_position,
          // gap_beginning);
          mappings.emplace_back(num_errors,
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
        split_sites.emplace_back(((actual_num_errors & 0xff) << 24) |
                                 ((gap_beginning & 0xff) << 16) |
                                 (read_mapping_length & 0xffff));
      }
    }
  }
}

template <typename MappingRecord>
void Chromap<MappingRecord>::VerifyCandidates(
    const SequenceBatch &read_batch, uint32_t read_index,
    const SequenceBatch &reference, MappingMetadata &mapping_metadata) {
  const std::vector<std::pair<uint64_t, uint64_t>> &minimizers =
      mapping_metadata.minimizers_;
  const std::vector<Candidate> &positive_candidates =
      mapping_metadata.positive_candidates_;
  const std::vector<Candidate> &negative_candidates =
      mapping_metadata.negative_candidates_;
  std::vector<std::pair<int, uint64_t>> &positive_mappings =
      mapping_metadata.positive_mappings_;
  std::vector<std::pair<int, uint64_t>> &negative_mappings =
      mapping_metadata.negative_mappings_;
  std::vector<int> &positive_split_sites =
      mapping_metadata.positive_split_sites_;
  std::vector<int> &negative_split_sites =
      mapping_metadata.negative_split_sites_;
  int &min_num_errors = mapping_metadata.min_num_errors_;
  int &num_best_mappings = mapping_metadata.num_best_mappings_;
  int &second_min_num_errors = mapping_metadata.second_min_num_errors_;
  int &num_second_best_mappings = mapping_metadata.num_second_best_mappings_;

  min_num_errors = error_threshold_ + 1;
  num_best_mappings = 0;
  second_min_num_errors = error_threshold_ + 1;
  num_second_best_mappings = 0;

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
      num_best_mappings = 1;
      num_second_best_mappings = 0;
      min_num_errors = 0;

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
          positive_mappings.emplace_back(
              0, positive_candidates[ci].position + read_length - 1);
        } else {
          negative_mappings.emplace_back(0, negative_candidates[ci].position);
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
      if ((uni_mapping.is_unique_) == 1) {
        ++num_uni_mappings;
      }
    }
  }
  for (auto &multi_mappings_on_one_ref_seq : multi_mappings) {
    for (auto &multi_mapping : multi_mappings_on_one_ref_seq) {
      if ((multi_mapping.is_unique_) != 1) {
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
void Chromap<MappingRecord>::ParseReadFormat(const std::string &read_format) {
  uint32_t i;
  int j = 0;
  int k = 0;  // for read1, read2, or barcode
  read1_format_[0] = 0;
  read1_format_[1] = -1;
  read1_format_[2] = 1;
  read2_format_[0] = 0;
  read2_format_[1] = -1;
  read2_format_[2] = 1;
  barcode_format_[0] = 0;
  barcode_format_[1] = -1;
  barcode_format_[2] = 1;
  int fields[3] = {0, -1, 1};
  char buffer[20];
  int blen = 0;
  for (i = 0; i < read_format.size(); ++i) {
    if (read_format[i] == ',' || i == 0) {
      if (i > 0) {
        buffer[blen] = '\0';
        if (j <= 1) {
          fields[j] = atoi(buffer);
        } else {
          fields[j] = buffer[0] == '+' ? 1 : -1;
        }
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
        if (j < 3) {
          buffer[blen] = read_format[i];
          ++blen;
        } else {
          ExitWithMessage("Unknown read format: " + read_format + "\n");
        }
      } else {
        buffer[blen] = '\0';
        if (j <= 1) {
          fields[j] = atoi(buffer);
        } else {
          fields[j] = buffer[0] == '+' ? 1 : -1;
        }
        ++j;
        blen = 0;
      }
    }
  }
  buffer[blen] = '\0';
  if (j <= 1) {
    fields[j] = atoi(buffer);
  } else {
    fields[j] = buffer[0] == '+' ? 1 : -1;
  }
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
       "Output mappings with barcode not in the whitelist")(
          "chr-order", "custom chromsome order", cxxopts::value<std::string>(),
          "FILE")("BED", "Output mappings in BED/BEDPE format")(
          "TagAlign", "Output mappings in TagAlign/PairedTagAlign format")(
          "SAM", "Output mappings in SAM format")(
          "pairs",
          "Output mappings in pairs format (defined by 4DN for HiC data)")(
          "pairs-natural-chr-order",
          "natural chromosome order for pairs flipping",
          cxxopts::value<std::string>(),
          "FILE")("barcode-translate",
                  "Convert barcode to the specified sequences during output",
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
      "PAF", "Output mappings in PAF format (only for test)")(
      "skip-barcode-check",
      "Do not check whether too few barcodes are in the whitelist");

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
      chromap::ExitWithMessage("Unrecognized preset parameters " + read_type +
                               "\n");
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
      chromap::ExitWithMessage("No reference specified!");
    }
    if (result.count("o")) {
      index_parameters.index_output_file_path =
          result["output"].as<std::string>();
    } else {
      chromap::ExitWithMessage("No output file specified!");
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
      chromap::ExitWithMessage("No reference specified!");
    }
    if (result.count("o")) {
      mapping_parameters.mapping_output_file_path =
          result["output"].as<std::string>();
    } else {
      chromap::ExitWithMessage("No output file specified!");
    }
    if (result.count("x")) {
      mapping_parameters.index_file_path = result["index"].as<std::string>();
    } else {
      chromap::ExitWithMessage("No index file specified!");
    }
    if (result.count("1")) {
      mapping_parameters.read_file1_paths =
          result["read1"].as<std::vector<std::string>>();
    } else {
      chromap::ExitWithMessage("No read file specified!");
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
        std::cerr << "WARNING: there are input barcode files but a barcode "
                     "whitelist file is missing!\n";
      }
    }

    if (result.count("barcode-whitelist")) {
      if (mapping_parameters.is_bulk_data) {
        chromap::ExitWithMessage(
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
        chromap::ExitWithMessage(
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

    if (result.count("barcode-translate")) {
      mapping_parameters.barcode_translate_table_path =
          result["barcode-translate"].as<std::string>();
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
        chromap::ExitWithMessage("Unknown mapping output format!");
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
          chromap::ExitWithMessage("No support for single-end HiC yet!");
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
          chromap::ExitWithMessage("Unknown mapping output format!");
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
          chromap::ExitWithMessage("Unknown mapping output format!");
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
