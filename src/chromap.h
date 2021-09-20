#ifndef CHROMAP_H_
#define CHROMAP_H_

#include <sys/resource.h>
#include <sys/time.h>

#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "index.h"
#include "khash.h"
#include "ksort.h"
#include "output_tools.h"
#include "sequence_batch.h"

#define CHROMAP_VERSION "0.1.2-r244"

namespace chromap {
struct uint128_t {
  uint64_t first;
  uint64_t second;
};

struct StackCell {
  size_t x;  // node
  int k, w;  // k: level; w: 0 if left child hasn't been processed
  StackCell(){};
  StackCell(int k_, size_t x_, int w_) : x(x_), k(k_), w(w_){};
};

struct _mm_history {
  unsigned int timestamp;
  std::vector<std::pair<uint64_t, uint64_t> > minimizers;
  std::vector<Candidate> positive_candidates;
  std::vector<Candidate> negative_candidates;
  uint32_t repetitive_seed_length;
};

struct Peak {
  uint32_t start_position;
  uint16_t length;
  uint32_t index;
};

KHASH_MAP_INIT_INT64(k128, uint128_t);
KHASH_MAP_INIT_INT64(k64_seq, uint64_t);
KHASH_SET_INIT_INT(k32_set);
KHASH_MAP_INIT_INT64(kmatrix, uint32_t);

struct BarcodeWithQual {
  uint32_t corrected_base_index1;
  char correct_base1;
  uint32_t corrected_base_index2;
  char correct_base2;
  double score;
  bool operator>(const BarcodeWithQual &b) const {
    return std::tie(score, corrected_base_index1, correct_base1,
                    corrected_base_index2, correct_base2) >
           std::tie(b.score, b.corrected_base_index1, b.correct_base1,
                    b.corrected_base_index2, b.correct_base2);
  }
};

struct IndexParameters {
  int kmer_size = 17;
  int window_size = 7;
  int num_threads = 1;
  std::string reference_file_path;
  std::string index_output_file_path;
};

struct MappingParameters {
  int error_threshold = 8;
  int match_score = 1;
  int mismatch_penalty = 4;
  std::vector<int> gap_open_penalties = {6, 6};
  std::vector<int> gap_extension_penalties = {1, 1};
  int min_num_seeds_required_for_mapping = 2;
  std::vector<int> max_seed_frequencies = {500, 1000};
  int max_num_best_mappings = 1;
  int max_insert_size = 1000;
  uint8_t mapq_threshold = 30;
  int num_threads = 1;
  int min_read_length = 30;
  int barcode_correction_error_threshold = 1;
  double barcode_correction_probability_threshold = 0.9;
  int multi_mapping_allocation_distance = 0;
  int multi_mapping_allocation_seed = 11;
  int drop_repetitive_reads = 500000;
  bool trim_adapters = false;
  bool remove_pcr_duplicates = false;
  bool remove_pcr_duplicates_at_bulk_level = true;
  bool is_bulk_data = true;
  bool allocate_multi_mappings = false;
  bool only_output_unique_mappings = true;
  bool output_mappings_not_in_whitelist = false;
  bool Tn5_shift = false;
  bool split_alignment = false;
  MappingOutputFormat mapping_output_format = MAPPINGFORMAT_BED;
  bool low_memory_mode = false;
  bool cell_by_bin = false;
  int bin_size = 5000;
  uint16_t depth_cutoff_to_call_peak = 3;
  int peak_min_length = 30;
  int peak_merge_max_length = 30;
  std::string reference_file_path;
  std::string index_file_path;
  std::vector<std::string> read_file1_paths;
  std::vector<std::string> read_file2_paths;
  std::vector<std::string> barcode_file_paths;
  std::string barcode_whitelist_file_path;
  std::string read_format;
  std::string mapping_output_file_path;
  std::string matrix_output_prefix;
  std::string custom_rid_order_path;
  std::string pairs_custom_rid_order_path;
};

#define SortMappingWithoutBarcode(m)                                    \
  (((((m).fragment_start_position << 16) | (m).fragment_length) << 8) | \
   (m).mapq)
//#define SortMappingWithoutBarcode(m) (m)

class ChromapDriver {
 public:
  ChromapDriver() {}
  ~ChromapDriver() {}
  void ParseArgsAndRun(int argc, char *argv[]);
};

template <typename MappingRecord = MappingWithoutBarcode>
class Chromap {
 public:
  // For index construction
  Chromap(const IndexParameters &index_parameters)
      : kmer_size_(index_parameters.kmer_size),
        window_size_(index_parameters.window_size),
        num_threads_(index_parameters.num_threads),
        reference_file_path_(index_parameters.reference_file_path),
        index_file_path_(index_parameters.index_output_file_path) {
    barcode_lookup_table_ = NULL;
    barcode_whitelist_lookup_table_ = NULL;
    barcode_histogram_ = NULL;
    barcode_index_table_ = NULL;
  }

  // For mapping
  Chromap(const MappingParameters &mapping_parameters)
      : error_threshold_(mapping_parameters.error_threshold),
        match_score_(mapping_parameters.match_score),
        mismatch_penalty_(mapping_parameters.mismatch_penalty),
        gap_open_penalties_(mapping_parameters.gap_open_penalties),
        gap_extension_penalties_(mapping_parameters.gap_extension_penalties),
        min_num_seeds_required_for_mapping_(
            mapping_parameters.min_num_seeds_required_for_mapping),
        max_seed_frequencies_(mapping_parameters.max_seed_frequencies),
        max_num_best_mappings_(mapping_parameters.max_num_best_mappings),
        max_insert_size_(mapping_parameters.max_insert_size),
        mapq_threshold_(mapping_parameters.mapq_threshold),
        num_threads_(mapping_parameters.num_threads),
        min_read_length_(mapping_parameters.min_read_length),
        barcode_correction_error_threshold_(
            mapping_parameters.barcode_correction_error_threshold),
        barcode_correction_probability_threshold_(
            mapping_parameters.barcode_correction_probability_threshold),
        multi_mapping_allocation_distance_(
            mapping_parameters.multi_mapping_allocation_distance),
        multi_mapping_allocation_seed_(
            mapping_parameters.multi_mapping_allocation_seed),
        drop_repetitive_reads_(mapping_parameters.drop_repetitive_reads),
        trim_adapters_(mapping_parameters.trim_adapters),
        remove_pcr_duplicates_(mapping_parameters.remove_pcr_duplicates),
        remove_pcr_duplicates_at_bulk_level_(
            mapping_parameters.remove_pcr_duplicates_at_bulk_level),
        is_bulk_data_(mapping_parameters.is_bulk_data),
        allocate_multi_mappings_(mapping_parameters.allocate_multi_mappings),
        only_output_unique_mappings_(
            mapping_parameters.only_output_unique_mappings),
        output_mappings_not_in_whitelist_(
            mapping_parameters.output_mappings_not_in_whitelist),
        Tn5_shift_(mapping_parameters.Tn5_shift),
        split_alignment_(mapping_parameters.split_alignment),
        mapping_output_format_(mapping_parameters.mapping_output_format),
        low_memory_mode_(mapping_parameters.low_memory_mode),
        cell_by_bin_(mapping_parameters.cell_by_bin),
        bin_size_(mapping_parameters.bin_size),
        depth_cutoff_to_call_peak_(
            mapping_parameters.depth_cutoff_to_call_peak),
        peak_min_length_(mapping_parameters.peak_min_length),
        peak_merge_max_length_(mapping_parameters.peak_merge_max_length),
        reference_file_path_(mapping_parameters.reference_file_path),
        index_file_path_(mapping_parameters.index_file_path),
        read_file1_paths_(mapping_parameters.read_file1_paths),
        read_file2_paths_(mapping_parameters.read_file2_paths),
        barcode_file_paths_(mapping_parameters.barcode_file_paths),
        barcode_whitelist_file_path_(
            mapping_parameters.barcode_whitelist_file_path),
        mapping_output_file_path_(mapping_parameters.mapping_output_file_path),
        matrix_output_prefix_(mapping_parameters.matrix_output_prefix),
        custom_rid_order_path_(mapping_parameters.custom_rid_order_path),
        pairs_custom_rid_order_path_(
            mapping_parameters.pairs_custom_rid_order_path) {
    barcode_lookup_table_ = kh_init(k64_seq);
    barcode_whitelist_lookup_table_ = kh_init(k64_seq);
    barcode_histogram_ = kh_init(k64_seq);
    barcode_index_table_ = kh_init(k64_seq);

    NUM_VPU_LANES_ = 0;
    if (error_threshold_ < 8) {
      NUM_VPU_LANES_ = 8;
    } else if (error_threshold_ < 16) {
      NUM_VPU_LANES_ = 4;
    }

    ParseReadFormat(mapping_parameters.read_format);
  }

  ~Chromap() {
    if (barcode_whitelist_lookup_table_ != NULL) {
      kh_destroy(k64_seq, barcode_whitelist_lookup_table_);
    }
    if (barcode_histogram_ != NULL) {
      kh_destroy(k64_seq, barcode_histogram_);
    }
    if (barcode_index_table_ != NULL) {
      kh_destroy(k64_seq, barcode_index_table_);
    }
    if (barcode_lookup_table_ != NULL) {
      kh_destroy(k64_seq, barcode_lookup_table_);
    }
    if (read_lookup_tables_.size() > 0) {
      for (uint32_t i = 0; i < read_lookup_tables_.size(); ++i) {
        kh_destroy(k128, read_lookup_tables_[i]);
      }
    }
  }

  KRADIX_SORT_INIT(without_barcode, MappingWithoutBarcode,
                   SortMappingWithoutBarcode, 11);
  KRADIX_SORT_INIT(with_barcode, MappingWithBarcode, SortMappingWithoutBarcode,
                   15);

  // For paired-end read mapping
  void MapPairedEndReads();
  uint32_t LoadPairedEndReadsWithBarcodes(SequenceBatch *read_batch1,
                                          SequenceBatch *read_batch2,
                                          SequenceBatch *barcode_batch);
  void TrimAdapterForPairedEndRead(uint32_t pair_index,
                                   SequenceBatch *read_batch1,
                                   SequenceBatch *read_batch2);
  bool PairedEndReadWithBarcodeIsDuplicate(uint32_t pair_index,
                                           const SequenceBatch &barcode_batch,
                                           const SequenceBatch &read_batch1,
                                           const SequenceBatch &read_batch2);
  void ReduceCandidatesForPairedEndReadOnOneDirection(
      const std::vector<Candidate> &candidates1,
      const std::vector<Candidate> &candidates2,
      std::vector<Candidate> *filtered_candidates1,
      std::vector<Candidate> *filtered_candidates2);
  void ReduceCandidatesForPairedEndRead(
      const std::vector<Candidate> &positive_candidates1,
      const std::vector<Candidate> &negative_candidates1,
      const std::vector<Candidate> &positive_candidates2,
      const std::vector<Candidate> &negative_candidates2,
      std::vector<Candidate> *filtered_positive_candidates1,
      std::vector<Candidate> *filtered_negative_candidates1,
      std::vector<Candidate> *filtered_positive_candidates2,
      std::vector<Candidate> *filtered_negative_candidates2);
  void GenerateBestMappingsForPairedEndReadOnOneDirection(
      Direction first_read_direction, uint32_t pair_index, int num_candidates1,
      int min_num_errors1, int num_best_mappings1, int second_min_num_errors1,
      int num_second_best_mappings1, const SequenceBatch &read_batch1,
      const std::vector<std::pair<int, uint64_t> > &mappings1,
      int num_candidates2, int min_num_errors2, int num_best_mappings2,
      int second_min_num_errors2, int num_second_best_mappings2,
      const SequenceBatch &read_batch2, const SequenceBatch &reference,
      const std::vector<std::pair<int, uint64_t> > &mappings2,
      std::vector<std::pair<uint32_t, uint32_t> > *best_mappings,
      int *min_sum_errors, int *num_best_mappings, int *second_min_sum_errors,
      int *num_second_best_mappings);
  void RecalibrateBestMappingsForPairedEndReadOnOneDirection(
      Direction first_read_direction, uint32_t pair_index, int min_sum_errors,
      int second_min_sum_errors, int min_num_errors1, int num_best_mappings1,
      int second_min_num_errors1, int num_second_best_mappings1,
      const SequenceBatch &read_batch1,
      const std::vector<std::pair<int, uint64_t> > &mappings1,
      int min_num_errors2, int num_best_mappings2, int second_min_num_errors2,
      int num_second_best_mappings2, const SequenceBatch &read_batch2,
      const SequenceBatch &reference,
      const std::vector<std::pair<int, uint64_t> > &mappings2,
      const std::vector<std::pair<uint32_t, uint32_t> > &edit_best_mappings,
      std::vector<std::pair<uint32_t, uint32_t> > *best_mappings,
      int *best_alignment_score, int *num_best_mappings,
      int *second_best_alignment_score, int *num_second_best_mappings);
  void ProcessBestMappingsForPairedEndReadOnOneDirection(
      Direction first_read_direction, Direction second_read_direction,
      uint32_t pair_index, uint8_t mapq, int num_candidates1,
      uint32_t repetitive_seed_length1, int min_num_errors1,
      int num_best_mappings1, int second_min_num_errors1,
      int num_second_best_mappings1, const SequenceBatch &read_batch1,
      const std::vector<std::pair<int, uint64_t> > &mappings1,
      const std::vector<int> &split_sites1, int num_candidates2,
      uint32_t repetitive_seed_length2, int min_num_errors2,
      int num_best_mappings2, int second_min_num_errors2,
      int num_second_best_mappings2, const SequenceBatch &read_batch2,
      const SequenceBatch &reference, const SequenceBatch &barcode_batch,
      const std::vector<int> &best_mapping_indices,
      const std::vector<std::pair<int, uint64_t> > &mappings2,
      const std::vector<int> &split_sites2,
      const std::vector<std::pair<uint32_t, uint32_t> > &best_mappings,
      int min_sum_errors, int num_best_mappings, int second_min_sum_errors,
      int num_second_best_mappings, int *best_mapping_index,
      int *num_best_mappings_reported, int force_mapq,
      std::vector<std::vector<MappingRecord> > *mappings_on_diff_ref_seqs);
  void GenerateBestMappingsForPairedEndRead(
      uint32_t pair_index, int num_positive_candidates1,
      int num_negative_candidates1, uint32_t repetitive_seed_length1,
      int min_num_errors1, int num_best_mappings1, int second_min_num_errors1,
      int num_second_best_mappings1, const SequenceBatch &read_batch1,
      const std::vector<std::pair<int, uint64_t> > &positive_mappings1,
      const std::vector<int> &positive_split_sites1,
      const std::vector<std::pair<int, uint64_t> > &negative_mappings1,
      const std::vector<int> &negative_split_sites1,
      int num_positive_candidates2, int num_negative_candidates2,
      uint32_t repetitive_seed_length2, int min_num_errors2,
      int num_best_mappings2, int second_min_num_errors2,
      int num_second_best_mappings2, const SequenceBatch &read_batch2,
      const SequenceBatch &reference, const SequenceBatch &barcode_batch,
      const std::vector<std::pair<int, uint64_t> > &positive_mappings2,
      const std::vector<int> &positive_split_sites2,
      const std::vector<std::pair<int, uint64_t> > &negative_mappings2,
      const std::vector<int> &negative_split_sites2,
      std::vector<int> *best_mapping_indices, std::mt19937 *generator,
      std::vector<std::pair<uint32_t, uint32_t> > *F1R2_best_mappings,
      std::vector<std::pair<uint32_t, uint32_t> > *F2R1_best_mappings,
      std::vector<std::pair<uint32_t, uint32_t> > *F1F2_best_mappings,
      std::vector<std::pair<uint32_t, uint32_t> > *R1R2_best_mappings,
      int *min_sum_errors, int *num_best_mappings, int *second_min_sum_errors,
      int *num_second_best_mappings, int force_mapq,
      std::vector<std::vector<MappingRecord> > *mappings_on_diff_ref_seqs);
  void EmplaceBackMappingRecord(
      uint32_t read_id, uint64_t barcode, uint32_t fragment_start_position,
      uint16_t fragment_length, uint8_t mapq, uint8_t direction,
      uint8_t is_unique, uint8_t num_dups, uint16_t positive_alignment_length,
      uint16_t negative_alignment_length,
      std::vector<MappingRecord> *mappings_on_diff_ref_seqs);
  void EmplaceBackMappingRecord(
      uint32_t read_id, const char *read1_name, const char *read2_name,
      uint16_t read1_length, uint16_t read2_length, uint64_t barcode,
      uint32_t fragment_start_position, uint16_t fragment_length, uint8_t mapq1,
      uint8_t mapq2, uint8_t direction, uint8_t is_unique, uint8_t num_dups,
      uint16_t positive_alignment_length, uint16_t negative_alignment_length,
      std::vector<MappingRecord> *mappings_on_diff_ref_seqs);
  void EmplaceBackMappingRecord(
      uint32_t read_id, const char *read_name, uint64_t cell_barcode, int rid1,
      int rid2, uint32_t pos1, uint32_t pos2, int direction1, int direction2,
      uint8_t mapq, uint8_t is_unique, uint8_t num_dups,
      std::vector<MappingRecord> *mappings_on_diff_ref_seqs);
  void ApplyTn5ShiftOnPairedEndMapping(
      uint32_t num_reference_sequences,
      std::vector<std::vector<MappingRecord> > *mappings);

  // For single-end read mapping
  void MapSingleEndReads();
  void GenerateBestMappingsForSingleEndRead(
      int num_positive_candidates, int num_negative_candidates,
      uint32_t repetitive_seed_length, int min_num_errors,
      int num_best_mappings, int second_min_num_errors,
      int num_second_best_mappings, const SequenceBatch &read_batch,
      uint32_t read_index, const SequenceBatch &reference,
      const SequenceBatch &barcode_batch,
      const std::vector<std::pair<int, uint64_t> > &positive_mappings,
      const std::vector<int> &positive_split_sites,
      const std::vector<std::pair<int, uint64_t> > &negative_mappings,
      const std::vector<int> &negative_split_sites,
      std::vector<std::vector<MappingRecord> > *mappings_on_diff_ref_seqs);
  void ProcessBestMappingsForSingleEndRead(
      Direction mapping_direction, uint8_t mapq, int num_candidates,
      uint32_t repetitive_seed_length, int min_num_errors,
      int num_best_mappings, int second_min_num_errors,
      int num_second_best_mappings, const SequenceBatch &read_batch,
      uint32_t read_index, const SequenceBatch &reference,
      const SequenceBatch &barcode_batch,
      const std::vector<int> &best_mapping_indices,
      const std::vector<std::pair<int, uint64_t> > &mappings,
      const std::vector<int> &split_sites, int *best_mapping_index,
      int *num_best_mappings_reported,
      std::vector<std::vector<MappingRecord> > *mappings_on_diff_ref_seqs);
  void EmplaceBackMappingRecord(
      uint32_t read_id, uint32_t barcode, uint32_t fragment_start_position,
      uint16_t fragment_length, uint8_t mapq, uint8_t direction,
      uint8_t is_unique, uint8_t num_dups,
      std::vector<MappingRecord> *mappings_on_diff_ref_seqs);
  void EmplaceBackMappingRecord(
      uint32_t read_id, const char *read_name, uint16_t read_length,
      uint32_t barcode, uint32_t fragment_start_position,
      uint16_t fragment_length, uint8_t mapq, uint8_t direction,
      uint8_t is_unique, uint8_t num_dups,
      std::vector<MappingRecord> *mappings_on_diff_ref_seqs);
  void EmplaceBackMappingRecord(
      uint32_t read_id, const char *read_name, uint32_t cell_barcode,
      uint8_t num_dups, int64_t position, int rid, int flag, uint8_t direction,
      uint8_t is_unique, uint8_t mapq, uint32_t NM, int n_cigar,
      uint32_t *cigar, std::string &MD_tag, const char *read,
      const char *read_qual,
      std::vector<MappingRecord> *mappings_on_diff_ref_seqs);
  void ApplyTn5ShiftOnSingleEndMapping(
      uint32_t num_reference_sequences,
      std::vector<std::vector<MappingRecord> > *mappings);
  uint32_t LoadSingleEndReadsWithBarcodes(SequenceBatch *read_batch,
                                          SequenceBatch *barcode_batch);

  // Supportive functions
  void ConstructIndex();
  int BandedAlignPatternToText(const char *pattern, const char *text,
                               const int read_length,
                               int *mapping_end_location);
  int BandedAlignPatternToTextWithDropOff(const char *pattern, const char *text,
                                          const int read_length,
                                          int *mapping_end_location,
                                          int *read_mapping_length);
  int BandedAlignPatternToTextWithDropOffFrom3End(const char *pattern,
                                                  const char *text,
                                                  const int read_length,
                                                  int *mapping_end_position,
                                                  int *read_mapping_length);
  void BandedAlign4PatternsToText(const char **patterns, const char *text,
                                  int read_length,
                                  int32_t *mapping_edit_distances,
                                  int32_t *mapping_end_positions);
  void BandedAlign8PatternsToText(const char **patterns, const char *text,
                                  int read_length,
                                  int16_t *mapping_edit_distances,
                                  int16_t *mapping_end_positions);
  void BandedTraceback(int min_num_errors, const char *pattern,
                       const char *text, const int read_length,
                       int *mapping_start_position);
  void BandedTracebackToEnd(int min_num_errors, const char *pattern,
                            const char *text, const int read_length,
                            int *mapping_end_position);
  int GetLongestMatchLength(const char *pattern, const char *text,
                            const int read_length);
  void MergeCandidates(std::vector<Candidate> &c1, std::vector<Candidate> &c2,
                       std::vector<Candidate> &buffer);
  int SupplementCandidates(
      const Index &index, uint32_t repetitive_seed_length1,
      uint32_t repetitive_seed_length2,
      std::vector<std::pair<uint64_t, uint64_t> > &minimizers1,
      std::vector<std::pair<uint64_t, uint64_t> > &minimizers2,
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
      std::vector<Candidate> &negative_candidates2_buffer);
  void PostProcessingInLowMemory(uint32_t num_mappings_in_mem,
                                 uint32_t num_reference_sequences,
                                 const SequenceBatch &reference);
  void VerifyCandidatesOnOneDirectionUsingSIMD(
      Direction candidate_direction, const SequenceBatch &read_batch,
      uint32_t read_index, const SequenceBatch &reference,
      const std::vector<Candidate> &candidates,
      std::vector<std::pair<int, uint64_t> > *mappings, int *min_num_errors,
      int *num_best_mappings, int *second_min_num_errors,
      int *num_second_best_mappings);
  void VerifyCandidatesOnOneDirection(
      Direction candidate_direction, const SequenceBatch &read_batch,
      uint32_t read_index, const SequenceBatch &reference,
      const std::vector<Candidate> &candidates,
      std::vector<std::pair<int, uint64_t> > *mappings,
      std::vector<int> *split_sites, int *min_num_errors,
      int *num_best_mappings, int *second_min_num_errors,
      int *num_second_best_mappings);
  void VerifyCandidates(
      const SequenceBatch &read_batch, uint32_t read_index,
      const SequenceBatch &reference,
      const std::vector<std::pair<uint64_t, uint64_t> > &minimizers,
      const std::vector<Candidate> &positive_candidates,
      const std::vector<Candidate> &negative_candidates,
      std::vector<std::pair<int, uint64_t> > *positive_mappings,
      std::vector<int> *positive_split_sites,
      std::vector<std::pair<int, uint64_t> > *negative_mappings,
      std::vector<int> *negative_split_sites, int *min_num_errors,
      int *num_best_mappings, int *second_min_num_errors,
      int *num_second_best_mappings);
  void GenerateMDTag(const char *pattern, const char *text,
                     int mapping_start_position, int n_cigar,
                     const uint32_t *cigar, int &NM, std::string &MD_tag);
  void AllocateMultiMappings(uint32_t num_reference_sequences);
  void RemovePCRDuplicate(uint32_t num_reference_sequences);
  uint32_t MoveMappingsInBuffersToMappingContainer(
      uint32_t num_reference_sequences,
      std::vector<std::vector<std::vector<MappingRecord> > >
          *mappings_on_diff_ref_seqs_for_diff_threads_for_saving);
  void OutputBarcodeStatistics();
  void OutputMappingStatistics();
  void OutputMappingStatistics(
      uint32_t num_reference_sequences,
      const std::vector<std::vector<MappingRecord> > &uni_mappings,
      const std::vector<std::vector<MappingRecord> > &multi_mappings);
  uint8_t GetMAPQForSingleEndRead(
      int error_threshold, int num_candidates, uint32_t repetitive_seed_length,
      uint16_t alignment_length, int min_num_errors, int num_best_mappings,
      int second_min_num_errors, int num_second_best_mappings,
      int max_num_error_difference, int read_length);
  uint8_t GetMAPQForPairedEndRead(
      int num_positive_candidates, int num_negative_candidates,
      uint32_t repetitive_seed_length1, uint32_t repetitive_seed_length2,
      uint16_t positive_alignment_length, uint16_t negative_alignment_length,
      int min_sum_errors, int num_best_mappings, int second_min_sum_errors,
      int num_second_best_mappings, int num_errors1, int num_errors2,
      int min_num_errors1, int min_num_errors2, int num_best_mappings1,
      int num_best_mappings2, int second_min_num_errors1,
      int second_min_num_errors2, int num_second_best_mappings1,
      int num_second_best_mappings2, int read1_length, int read2_length,
      int force_mapq, uint8_t &mapq1, uint8_t &mapq2);
  void SortOutputMappings(uint32_t num_reference_sequences,
                          std::vector<std::vector<MappingRecord> > *mappings);
  void BuildAugmentedTree(uint32_t ref_id);
  uint32_t GetNumOverlappedMappings(uint32_t ref_id,
                                    const MappingRecord &mapping);
  void LoadBarcodeWhitelist();
  void ComputeBarcodeAbundance(uint64_t max_num_sample_barcodes);
  void UpdateBarcodeAbundance(uint32_t num_loaded_barcodes,
                              const SequenceBatch &barcode_batch);
  bool CorrectBarcodeAt(uint32_t barcode_index, SequenceBatch *barcode_batch,
                        uint64_t *num_barcode_in_whitelist,
                        uint64_t *num_corrected_barcode);
  uint32_t CallPeaks(uint16_t coverage_threshold,
                     uint32_t num_reference_sequences,
                     const SequenceBatch &reference);
  void OutputFeatureMatrix(uint32_t num_sequences,
                           const SequenceBatch &reference);
  void GetNumOverlappedBins(uint32_t rid, uint32_t start_position,
                            uint16_t mapping_length, uint32_t num_sequences,
                            const SequenceBatch &reference,
                            std::vector<uint32_t> &overlapped_peak_indices);
  uint32_t GetNumOverlappedPeaks(
      uint32_t ref_id, const MappingRecord &mapping,
      std::vector<uint32_t> &overlapped_peak_indices);
  void BuildAugmentedTreeForPeaks(uint32_t ref_id);
  void OutputMappingsInVector(
      uint8_t mapq_threshold, uint32_t num_reference_sequences,
      const SequenceBatch &reference,
      const std::vector<std::vector<MappingRecord> > &mappings);
  void OutputMappings(uint32_t num_reference_sequences,
                      const SequenceBatch &reference,
                      const std::vector<std::vector<MappingRecord> > &mappings);
  int AdjustGapBeginning(Direction mapping_direction, const char *ref,
                         const char *read, int *gap_beginning, int read_end,
                         int ref_start_position, int ref_end_position,
                         int *n_cigar, uint32_t **cigar);
  void GetRefStartEndPositionForReadFromMapping(
      Direction mapping_direction, const std::pair<int, uint64_t> &mapping,
      const char *read, int read_length, int in_split_site,
      const SequenceBatch &reference, uint32_t *ref_start_position,
      uint32_t *ref_end_position, int *n_cigar, uint32_t **cigar, int *NM,
      std::string &MD_TAG);
  void GetRidRank(const std::string rid_order_path,
                  const SequenceBatch &reference, std::vector<int> &rid_rank);
  void RerankCandidatesRid(std::vector<Candidate> &candidates);
  void ParseReadFormat(const std::string &read_format);

  inline static double GetRealTime() {
    struct timeval tp;
    struct timezone tzp;
    gettimeofday(&tp, &tzp);
    return tp.tv_sec + tp.tv_usec * 1e-6;
  }
  inline static double GetCPUTime() {
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    return r.ru_utime.tv_sec + r.ru_stime.tv_sec +
           1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
  }
  inline static void ExitWithMessage(const std::string &message) {
    std::cerr << message << std::endl;
    exit(-1);
  }

 private:
  // Parameters
  int kmer_size_;
  int window_size_;
  int error_threshold_;
  int NUM_VPU_LANES_;
  int match_score_;
  int mismatch_penalty_;
  std::vector<int> gap_open_penalties_;
  std::vector<int> gap_extension_penalties_;
  int min_num_seeds_required_for_mapping_;
  std::vector<int> max_seed_frequencies_;
  // Read with # best mappings greater than it will have this number of best
  // mappings reported.
  int max_num_best_mappings_;
  int max_insert_size_;
  uint8_t mapq_threshold_;
  int num_threads_;
  int min_read_length_;
  int barcode_correction_error_threshold_;
  double barcode_correction_probability_threshold_;
  int multi_mapping_allocation_distance_;
  int multi_mapping_allocation_seed_;
  // Read with more than this number of mappings will be dropped.
  int drop_repetitive_reads_;
  bool trim_adapters_;
  bool remove_pcr_duplicates_;
  bool remove_pcr_duplicates_at_bulk_level_;
  bool is_bulk_data_;
  bool allocate_multi_mappings_;
  bool only_output_unique_mappings_;
  bool output_mappings_not_in_whitelist_;
  bool Tn5_shift_;
  bool split_alignment_;
  MappingOutputFormat mapping_output_format_;
  // Default batch size, # reads for single-end reads, # read pairs for
  // paired-end reads.
  uint32_t read_batch_size_ = 500000;
  bool low_memory_mode_;
  bool cell_by_bin_;
  int bin_size_;
  uint16_t depth_cutoff_to_call_peak_;
  int peak_min_length_;
  int peak_merge_max_length_;
  std::string reference_file_path_;
  std::string index_file_path_;
  std::vector<std::string> read_file1_paths_;
  std::vector<std::string> read_file2_paths_;
  std::vector<std::string> barcode_file_paths_;
  std::string barcode_whitelist_file_path_;
  int barcode_format_[2];  // 0-start, 1-end (includsive)
  int read1_format_[2];
  int read2_format_[2];
  std::string mapping_output_file_path_;
  FILE *mapping_output_file_;
  std::string matrix_output_prefix_;
  // The order for general sorting.
  std::string custom_rid_order_path_;
  // The order for pairs format flipping.
  std::string pairs_custom_rid_order_path_;
  std::vector<int> custom_rid_rank_;
  std::vector<int> pairs_custom_rid_rank_;
  // khash_t(k32_set)* barcode_whitelist_lookup_table_;
  khash_t(k64_seq) * barcode_whitelist_lookup_table_;
  // For identical read dedupe
  int allocated_barcode_lookup_table_size_ = (1 << 10);
  khash_t(k64_seq) * barcode_lookup_table_;
  std::vector<khash_t(k128) *> read_lookup_tables_;
  // For mapping
  int min_unique_mapping_mapq_ = 4;
  std::vector<TempMappingFileHandle<MappingRecord> > temp_mapping_file_handles_;
  std::vector<std::vector<MappingRecord> > mappings_on_diff_ref_seqs_;
  std::vector<std::vector<MappingRecord> > deduped_mappings_on_diff_ref_seqs_;
  std::vector<std::pair<uint32_t, MappingRecord> > multi_mappings_;
  std::vector<std::vector<MappingRecord> > allocated_mappings_on_diff_ref_seqs_;
  std::vector<std::vector<uint32_t> > tree_extras_on_diff_ref_seqs_;  // max
  // (max_level, # nodes)
  std::vector<std::pair<int, uint32_t> > tree_info_on_diff_ref_seqs_;
  OutputTools<MappingRecord> output_tools_;
  // For mapping stats.
  uint64_t num_candidates_ = 0;
  uint64_t num_mappings_ = 0;
  uint64_t num_mapped_reads_ = 0;
  uint64_t num_uniquely_mapped_reads_ = 0;
  uint64_t num_reads_ = 0;
  uint64_t num_duplicated_reads_ = 0;  // # identical reads
  // For barcode stats
  uint64_t initial_num_sample_barcodes_ = 20000000;
  uint64_t num_sample_barcodes_ = 0;
  uint64_t num_barcode_in_whitelist_ = 0;
  uint64_t num_corrected_barcode_ = 0;
  uint32_t barcode_length_ = 0;
  khash_t(k64_seq) * barcode_histogram_;
  khash_t(k64_seq) * barcode_index_table_;
  // For peak calling
  std::vector<std::vector<uint16_t> > pileup_on_diff_ref_seqs_;
  std::vector<std::vector<Peak> > peaks_on_diff_ref_seqs_;
};
}  // namespace chromap

#endif  // CHROMAP_H_
