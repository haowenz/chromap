#ifndef MAPPING_GENERATOR_H_
#define MAPPING_GENERATOR_H_

#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>
#include <tuple>
#include <vector>

#include "alignment.h"
#include "bed_mapping.h"
#include "ksw.h"
#include "mapping.h"
#include "mapping_metadata.h"
#include "mapping_parameters.h"
#include "paf_mapping.h"
#include "paired_end_mapping_metadata.h"
#include "pairs_mapping.h"
#include "sam_mapping.h"
#include "sequence_batch.h"
#include "utils.h"

namespace chromap {

// Class to generate mappings from candidates. It supports multi-threadidng as
// only the parameters are owned by the class.
template <typename MappingRecord>
class MappingGenerator {
 public:
  MappingGenerator() = delete;
  MappingGenerator(int error_threshold, int NUM_VPU_LANES, int match_score,
                   int mismatch_penalty,
                   const std::vector<int> &gap_open_penalties,
                   const std::vector<int> &gap_extension_penalties,
                   int max_num_best_mappings, int max_insert_size,
                   int min_read_length, int drop_repetitive_reads,
                   bool is_bulk_data, bool split_alignment,
                   const MappingOutputFormat &mapping_output_format,
                   const std::vector<int> &pairs_custom_rid_rank)
      : error_threshold_(error_threshold),
        NUM_VPU_LANES_(NUM_VPU_LANES),
        match_score_(match_score),
        mismatch_penalty_(mismatch_penalty),
        gap_open_penalties_(gap_open_penalties),
        gap_extension_penalties_(gap_extension_penalties),
        max_num_best_mappings_(max_num_best_mappings),
        max_insert_size_(max_insert_size),
        min_read_length_(min_read_length),
        drop_repetitive_reads_(drop_repetitive_reads),
        is_bulk_data_(is_bulk_data),
        split_alignment_(split_alignment),
        mapping_output_format_(mapping_output_format),
        pairs_custom_rid_rank_(pairs_custom_rid_rank) {}

  ~MappingGenerator() = default;

  void GenerateBestMappingsForSingleEndRead(
      const SequenceBatch &read_batch, uint32_t read_index,
      const SequenceBatch &reference, const SequenceBatch &barcode_batch,
      MappingMetadata &mapping_metadata,
      std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs);

  void GenerateBestMappingsForPairedEndRead(
      uint32_t pair_index, const SequenceBatch &read_batch1,
      const SequenceBatch &read_batch2, const SequenceBatch &barcode_batch,
      const SequenceBatch &reference, std::vector<int> &best_mapping_indices,
      std::mt19937 &generator, int force_mapq,
      PairedEndMappingMetadata &paired_end_mapping_metadata,
      std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs);

 private:
  void ProcessBestMappingsForSingleEndRead(
      Direction mapping_direction, uint8_t mapq, int num_candidates,
      uint32_t repetitive_seed_length, int min_num_errors,
      int num_best_mappings, int second_min_num_errors,
      int num_second_best_mappings, const SequenceBatch &read_batch,
      uint32_t read_index, const SequenceBatch &reference,
      const SequenceBatch &barcode_batch,
      const std::vector<int> &best_mapping_indices,
      const std::vector<std::pair<int, uint64_t>> &mappings,
      const std::vector<int> &split_sites, int &best_mapping_index,
      int &num_best_mappings_reported,
      std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs);

  void GenerateBestMappingsForPairedEndReadOnOneDirection(
      Direction first_read_direction, uint32_t pair_index, int num_candidates1,
      int min_num_errors1, int num_best_mappings1, int second_min_num_errors1,
      int num_second_best_mappings1, const SequenceBatch &read_batch1,
      const std::vector<std::pair<int, uint64_t>> &mappings1,
      int num_candidates2, int min_num_errors2, int num_best_mappings2,
      int second_min_num_errors2, int num_second_best_mappings2,
      const SequenceBatch &read_batch2, const SequenceBatch &reference,
      const std::vector<std::pair<int, uint64_t>> &mappings2,
      std::vector<std::pair<uint32_t, uint32_t>> &best_mappings,
      int &min_sum_errors, int &num_best_mappings, int &second_min_sum_errors,
      int &num_second_best_mappings);

  void RecalibrateBestMappingsForPairedEndReadOnOneDirection(
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
      int *second_best_alignment_score, int *num_second_best_mappings);

  void ProcessBestMappingsForPairedEndReadOnOneDirection(
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
      int num_second_best_mappings, int &best_mapping_index,
      int &num_best_mappings_reported, int force_mapq,
      std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs);

  void GetRefStartEndPositionForReadFromMapping(
      Direction mapping_direction, const std::pair<int, uint64_t> &mapping,
      const char *read, int read_length, int in_split_site,
      const SequenceBatch &reference, uint32_t *ref_start_position,
      uint32_t *ref_end_position, int *n_cigar, uint32_t **cigar, int *NM,
      std::string &MD_TAG);

  void EmplaceBackMappingRecord(
      uint32_t read_id, uint64_t barcode, uint32_t fragment_start_position,
      uint16_t fragment_length, uint8_t mapq, uint8_t direction,
      uint8_t is_unique, uint8_t num_dups,
      std::vector<MappingRecord> *mappings_on_diff_ref_seqs);
  void EmplaceBackMappingRecord(
      uint32_t read_id, const char *read_name, uint16_t read_length,
      uint64_t barcode, uint32_t fragment_start_position,
      uint16_t fragment_length, uint8_t mapq, uint8_t direction,
      uint8_t is_unique, uint8_t num_dups,
      std::vector<MappingRecord> *mappings_on_diff_ref_seqs);
  void EmplaceBackMappingRecord(
      uint32_t read_id, const char *read_name, uint64_t cell_barcode,
      uint8_t num_dups, int64_t position, int rid, int flag, uint8_t direction,
      uint8_t is_unique, uint8_t mapq, uint32_t NM, int n_cigar,
      uint32_t *cigar, std::string &MD_tag, const char *read,
      const char *read_qual,
      std::vector<MappingRecord> *mappings_on_diff_ref_seqs);

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

  const int error_threshold_;
  const int NUM_VPU_LANES_;
  const int match_score_;
  const int mismatch_penalty_;
  const std::vector<int> gap_open_penalties_;
  const std::vector<int> gap_extension_penalties_;
  const int max_num_best_mappings_;
  const int max_insert_size_;
  const int min_read_length_;
  const int drop_repetitive_reads_;
  const bool is_bulk_data_;
  const bool split_alignment_;
  const MappingOutputFormat mapping_output_format_;
  const std::vector<int> pairs_custom_rid_rank_;
};

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::GenerateBestMappingsForSingleEndRead(
    const SequenceBatch &read_batch, uint32_t read_index,
    const SequenceBatch &reference, const SequenceBatch &barcode_batch,
    MappingMetadata &mapping_metadata,
    std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs) {
  int num_positive_candidates = mapping_metadata.positive_candidates_.size();
  int num_negative_candidates = mapping_metadata.negative_candidates_.size();

  uint32_t repetitive_seed_length = mapping_metadata.repetitive_seed_length_;

  int min_num_errors = mapping_metadata.min_num_errors_;
  int num_best_mappings = mapping_metadata.num_best_mappings_;
  int second_min_num_errors = mapping_metadata.second_min_num_errors_;
  int num_second_best_mappings = mapping_metadata.num_second_best_mappings_;
  const std::vector<std::pair<int, uint64_t>> &positive_mappings =
      mapping_metadata.positive_mappings_;
  const std::vector<std::pair<int, uint64_t>> &negative_mappings =
      mapping_metadata.negative_mappings_;
  const std::vector<int> &positive_split_sites =
      mapping_metadata.positive_split_sites_;
  const std::vector<int> &negative_split_sites =
      mapping_metadata.negative_split_sites_;

  // we will use reservoir sampling
  std::vector<int> best_mapping_indices(max_num_best_mappings_);
  std::iota(best_mapping_indices.begin(), best_mapping_indices.end(), 0);
  if (num_best_mappings > max_num_best_mappings_) {
    std::mt19937 generator(11);
    for (int i = max_num_best_mappings_; i < num_best_mappings; ++i) {
      // important: inclusive range
      std::uniform_int_distribution<int> distribution(0, i);
      int j = distribution(generator);
      if (j < max_num_best_mappings_) {
        best_mapping_indices[j] = i;
      }
    }
    std::sort(best_mapping_indices.begin(), best_mapping_indices.end());
  }

  uint8_t mapq = 0;
  int best_mapping_index = 0;
  int num_best_mappings_reported = 0;
  ProcessBestMappingsForSingleEndRead(
      kPositive, mapq, num_positive_candidates, repetitive_seed_length,
      min_num_errors, num_best_mappings, second_min_num_errors,
      num_second_best_mappings, read_batch, read_index, reference,
      barcode_batch, best_mapping_indices, positive_mappings,
      positive_split_sites, best_mapping_index, num_best_mappings_reported,
      mappings_on_diff_ref_seqs);

  if (num_best_mappings_reported !=
      std::min(num_best_mappings, max_num_best_mappings_)) {
    ProcessBestMappingsForSingleEndRead(
        kNegative, num_negative_candidates, repetitive_seed_length, mapq,
        min_num_errors, num_best_mappings, second_min_num_errors,
        num_second_best_mappings, read_batch, read_index, reference,
        barcode_batch, best_mapping_indices, negative_mappings,
        negative_split_sites, best_mapping_index, num_best_mappings_reported,
        mappings_on_diff_ref_seqs);
  }
}

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::GenerateBestMappingsForPairedEndRead(
    uint32_t pair_index, const SequenceBatch &read_batch1,
    const SequenceBatch &read_batch2, const SequenceBatch &barcode_batch,
    const SequenceBatch &reference, std::vector<int> &best_mapping_indices,
    std::mt19937 &generator, int force_mapq,
    PairedEndMappingMetadata &paired_end_mapping_metadata,
    std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs) {
  const int num_positive_candidates1 =
      paired_end_mapping_metadata.mapping_metadata1_.positive_candidates_
          .size();
  const int num_negative_candidates1 =
      paired_end_mapping_metadata.mapping_metadata1_.negative_candidates_
          .size();
  const uint32_t repetitive_seed_length1 =
      paired_end_mapping_metadata.mapping_metadata1_.repetitive_seed_length_;
  const int min_num_errors1 =
      paired_end_mapping_metadata.mapping_metadata1_.min_num_errors_;
  const int num_best_mappings1 =
      paired_end_mapping_metadata.mapping_metadata1_.num_best_mappings_;
  const int second_min_num_errors1 =
      paired_end_mapping_metadata.mapping_metadata1_.second_min_num_errors_;
  const int num_second_best_mappings1 =
      paired_end_mapping_metadata.mapping_metadata1_.num_second_best_mappings_;

  const int num_positive_candidates2 =
      paired_end_mapping_metadata.mapping_metadata2_.positive_candidates_
          .size();
  const int num_negative_candidates2 =
      paired_end_mapping_metadata.mapping_metadata2_.negative_candidates_
          .size();
  const uint32_t repetitive_seed_length2 =
      paired_end_mapping_metadata.mapping_metadata2_.repetitive_seed_length_;
  const int min_num_errors2 =
      paired_end_mapping_metadata.mapping_metadata2_.min_num_errors_;
  const int num_best_mappings2 =
      paired_end_mapping_metadata.mapping_metadata2_.num_best_mappings_;
  const int second_min_num_errors2 =
      paired_end_mapping_metadata.mapping_metadata2_.second_min_num_errors_;
  const int num_second_best_mappings2 =
      paired_end_mapping_metadata.mapping_metadata2_.num_second_best_mappings_;

  int &min_sum_errors = paired_end_mapping_metadata.min_sum_errors_;
  int &num_best_mappings = paired_end_mapping_metadata.num_best_mappings_;
  int &second_min_sum_errors =
      paired_end_mapping_metadata.second_min_sum_errors_;
  int &num_second_best_mappings =
      paired_end_mapping_metadata.num_second_best_mappings_;

  const std::vector<std::pair<int, uint64_t>> &positive_mappings1 =
      paired_end_mapping_metadata.mapping_metadata1_.positive_mappings_;
  const std::vector<std::pair<int, uint64_t>> &negative_mappings1 =
      paired_end_mapping_metadata.mapping_metadata1_.negative_mappings_;
  const std::vector<int> &positive_split_sites1 =
      paired_end_mapping_metadata.mapping_metadata1_.positive_split_sites_;
  const std::vector<int> &negative_split_sites1 =
      paired_end_mapping_metadata.mapping_metadata1_.negative_split_sites_;

  const std::vector<std::pair<int, uint64_t>> &positive_mappings2 =
      paired_end_mapping_metadata.mapping_metadata2_.positive_mappings_;
  const std::vector<std::pair<int, uint64_t>> &negative_mappings2 =
      paired_end_mapping_metadata.mapping_metadata2_.negative_mappings_;
  const std::vector<int> &positive_split_sites2 =
      paired_end_mapping_metadata.mapping_metadata2_.positive_split_sites_;
  const std::vector<int> &negative_split_sites2 =
      paired_end_mapping_metadata.mapping_metadata2_.negative_split_sites_;

  std::vector<std::pair<uint32_t, uint32_t>> &F1R2_best_mappings =
      paired_end_mapping_metadata.F1R2_best_mappings_;
  std::vector<std::pair<uint32_t, uint32_t>> &F2R1_best_mappings =
      paired_end_mapping_metadata.F2R1_best_mappings_;
  std::vector<std::pair<uint32_t, uint32_t>> &F1F2_best_mappings =
      paired_end_mapping_metadata.F1F2_best_mappings_;
  std::vector<std::pair<uint32_t, uint32_t>> &R1R2_best_mappings =
      paired_end_mapping_metadata.R1R2_best_mappings_;

  min_sum_errors = 2 * error_threshold_ + 1;
  num_best_mappings = 0;
  second_min_sum_errors = min_sum_errors;
  num_second_best_mappings = 0;

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
  if (num_best_mappings <= drop_repetitive_reads_) {
    // we will use reservoir sampling
    // std::vector<int> best_mapping_indices(max_num_best_mappings_);
    std::iota(best_mapping_indices.begin(), best_mapping_indices.end(), 0);
    if (num_best_mappings > max_num_best_mappings_) {
      // std::mt19937 generator(11);
      for (int i = max_num_best_mappings_; i < num_best_mappings; ++i) {
        std::uniform_int_distribution<int> distribution(
            0, i);  // important: inclusive range
        int j = distribution(generator);
        // int j = distribution(tmp_generator);
        if (j < max_num_best_mappings_) {
          best_mapping_indices[j] = i;
        }
      }
      std::sort(best_mapping_indices.begin(), best_mapping_indices.end());
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
        reference, barcode_batch, best_mapping_indices, negative_mappings2,
        negative_split_sites2, F1R2_best_mappings, min_sum_errors,
        num_best_mappings, second_min_sum_errors, num_second_best_mappings,
        best_mapping_index, num_best_mappings_reported, force_mapq,
        mappings_on_diff_ref_seqs);

    if (num_best_mappings_reported !=
        std::min(max_num_best_mappings_, num_best_mappings)) {
      ProcessBestMappingsForPairedEndReadOnOneDirection(
          kNegative, kPositive, pair_index, mapq, num_negative_candidates1,
          repetitive_seed_length1, min_num_errors1, num_best_mappings1,
          second_min_num_errors1, num_second_best_mappings1, read_batch1,
          negative_mappings1, negative_split_sites1, num_positive_candidates2,
          repetitive_seed_length2, min_num_errors2, num_best_mappings2,
          second_min_num_errors2, num_second_best_mappings2, read_batch2,
          reference, barcode_batch, best_mapping_indices, positive_mappings2,
          positive_split_sites2, F2R1_best_mappings, min_sum_errors,
          num_best_mappings, second_min_sum_errors, num_second_best_mappings,
          best_mapping_index, num_best_mappings_reported, force_mapq,
          mappings_on_diff_ref_seqs);
    }

    if (split_alignment_ &&
        num_best_mappings_reported !=
            std::min(max_num_best_mappings_, num_best_mappings)) {
      ProcessBestMappingsForPairedEndReadOnOneDirection(
          kPositive, kPositive, pair_index, mapq, num_positive_candidates1,
          repetitive_seed_length1, min_num_errors1, num_best_mappings1,
          second_min_num_errors1, num_second_best_mappings1, read_batch1,
          positive_mappings1, positive_split_sites1, num_positive_candidates2,
          repetitive_seed_length2, min_num_errors2, num_best_mappings2,
          second_min_num_errors2, num_second_best_mappings2, read_batch2,
          reference, barcode_batch, best_mapping_indices, positive_mappings2,
          positive_split_sites2, F1F2_best_mappings, min_sum_errors,
          num_best_mappings, second_min_sum_errors, num_second_best_mappings,
          best_mapping_index, num_best_mappings_reported, force_mapq,
          mappings_on_diff_ref_seqs);
    }

    if (split_alignment_ &&
        num_best_mappings_reported !=
            std::min(max_num_best_mappings_, num_best_mappings)) {
      ProcessBestMappingsForPairedEndReadOnOneDirection(
          kNegative, kNegative, pair_index, mapq, num_negative_candidates1,
          repetitive_seed_length1, min_num_errors1, num_best_mappings1,
          second_min_num_errors1, num_second_best_mappings1, read_batch1,
          negative_mappings1, negative_split_sites1, num_positive_candidates2,
          repetitive_seed_length2, min_num_errors2, num_best_mappings2,
          second_min_num_errors2, num_second_best_mappings2, read_batch2,
          reference, barcode_batch, best_mapping_indices, negative_mappings2,
          negative_split_sites2, R1R2_best_mappings, min_sum_errors,
          num_best_mappings, second_min_sum_errors, num_second_best_mappings,
          best_mapping_index, num_best_mappings_reported, force_mapq,
          mappings_on_diff_ref_seqs);
    }
  }
}

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::ProcessBestMappingsForSingleEndRead(
    Direction mapping_direction, uint8_t mapq, int num_candidates,
    uint32_t repetitive_seed_length, int min_num_errors, int num_best_mappings,
    int second_min_num_errors, int num_second_best_mappings,
    const SequenceBatch &read_batch, uint32_t read_index,
    const SequenceBatch &reference, const SequenceBatch &barcode_batch,
    const std::vector<int> &best_mapping_indices,
    const std::vector<std::pair<int, uint64_t>> &mappings,
    const std::vector<int> &split_sites, int &best_mapping_index,
    int &num_best_mappings_reported,
    std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs) {
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
      if (best_mapping_index ==
          best_mapping_indices[num_best_mappings_reported]) {
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
          if (num_best_mappings_reported >= 1) {
            flag |= BAM_FSECONDARY;
          }
          EmplaceBackMappingRecord(
              read_id, read_name, barcode_key, 1, ref_start_position, rid, flag,
              0, is_unique, mapq, NM, n_cigar, cigar, MD_tag, effect_read,
              read_batch.GetSequenceQualAt(read_index),
              &(mappings_on_diff_ref_seqs[rid]));
        } else if (mapping_output_format_ == MAPPINGFORMAT_PAF) {
          EmplaceBackMappingRecord(
              read_id, read_name, read_length, barcode_key, ref_start_position,
              ref_end_position - ref_start_position + 1, mapq, direction,
              is_unique, 1, &(mappings_on_diff_ref_seqs[rid]));
        } else {
          EmplaceBackMappingRecord(read_id, barcode_key, ref_start_position,
                                   ref_end_position - ref_start_position + 1,
                                   mapq, direction, is_unique, 1,
                                   &(mappings_on_diff_ref_seqs[rid]));
        }
        num_best_mappings_reported++;
        if (num_best_mappings_reported ==
            std::min(max_num_best_mappings_, num_best_mappings)) {
          break;
        }
      }
      best_mapping_index++;
    }
  }
}

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::
    GenerateBestMappingsForPairedEndReadOnOneDirection(
        Direction first_read_direction, uint32_t pair_index,
        int num_candidates1, int min_num_errors1, int num_best_mappings1,
        int second_min_num_errors1, int num_second_best_mappings1,
        const SequenceBatch &read_batch1,
        const std::vector<std::pair<int, uint64_t>> &mappings1,
        int num_candidates2, int min_num_errors2, int num_best_mappings2,
        int second_min_num_errors2, int num_second_best_mappings2,
        const SequenceBatch &read_batch2, const SequenceBatch &reference,
        const std::vector<std::pair<int, uint64_t>> &mappings2,
        std::vector<std::pair<uint32_t, uint32_t>> &best_mappings,
        int &min_sum_errors, int &num_best_mappings, int &second_min_sum_errors,
        int &num_second_best_mappings) {
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
        best_mappings.emplace_back(i1, i2);
        min_sum_errors = min_num_errors1 + min_num_errors2;
        //*second_min_sum_errors = min_num_errors1 + min_num_errors2 + 1;
        num_best_mappings++;
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
        if (current_sum_errors < min_sum_errors) {
          second_min_sum_errors = min_sum_errors;
          num_second_best_mappings = num_best_mappings;
          min_sum_errors = current_sum_errors;
          num_best_mappings = 1;
          best_mappings.clear();
          best_mappings.emplace_back(i1, current_i2);
        } else if (current_sum_errors == min_sum_errors) {
          num_best_mappings++;
          best_mappings.emplace_back(i1, current_i2);
        } else if (current_sum_errors == second_min_sum_errors) {
          num_second_best_mappings++;
        } else if (current_sum_errors < second_min_sum_errors) {
          second_min_sum_errors = current_sum_errors;
          num_second_best_mappings = 1;
        }
        ++current_i2;
      }
      ++i1;
    }
  }
}

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::
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
void MappingGenerator<MappingRecord>::
    ProcessBestMappingsForPairedEndReadOnOneDirection(
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
        int num_second_best_mappings, int &best_mapping_index,
        int &num_best_mappings_reported, int force_mapq,
        std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs) {
  const char *read1 = read_batch1.GetSequenceAt(pair_index);
  const char *read2 = read_batch2.GetSequenceAt(pair_index);
  const uint32_t read1_length = read_batch1.GetSequenceLengthAt(pair_index);
  const uint32_t read2_length = read_batch2.GetSequenceLengthAt(pair_index);
  const char *read1_name = read_batch1.GetSequenceNameAt(pair_index);
  const char *read2_name = read_batch2.GetSequenceNameAt(pair_index);
  const std::string &negative_read1 =
      read_batch1.GetNegativeSequenceAt(pair_index);
  const std::string &negative_read2 =
      read_batch2.GetNegativeSequenceAt(pair_index);
  const uint32_t read_id = read_batch1.GetSequenceIdAt(pair_index);

  const uint8_t is_unique = (num_best_mappings == 1 ||
                             num_best_mappings1 == 1 || num_best_mappings2 == 1)
                                ? 1
                                : 0;
  uint64_t barcode_key = 0;
  if (!is_bulk_data_) {
    barcode_key = barcode_batch.GenerateSeedFromSequenceAt(
        pair_index, 0, barcode_batch.GetSequenceLengthAt(pair_index));
  }

  for (uint32_t mi = 0; mi < best_mappings.size(); ++mi) {
    const uint32_t i1 = best_mappings[mi].first;
    const uint32_t i2 = best_mappings[mi].second;
    const int current_sum_errors = mappings1[i1].first + mappings2[i2].first;
    if (current_sum_errors == min_sum_errors) {
      if (best_mapping_index ==
          best_mapping_indices[num_best_mappings_reported]) {
        const uint32_t rid1 = mappings1[i1].second >> 32;
        const uint32_t rid2 = mappings2[i2].second >> 32;
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
          if (num_best_mappings_reported >= 1) {
            flag1 |= BAM_FSECONDARY;
            flag2 |= BAM_FSECONDARY;
          }
          // printf("%d %d\n", ref_start_position1, ref_end_position1);
          // printf("%d %d\n", ref_start_position2, ref_end_position2);
          EmplaceBackMappingRecord(
              read_id, read1_name, barcode_key, 1, ref_start_position1, rid1,
              flag1, first_read_direction == kPositive ? 1 : 0, is_unique, mapq,
              NM1, n_cigar1, cigar1, MD_tag1, effect_read1,
              read_batch1.GetSequenceQualAt(pair_index),
              &(mappings_on_diff_ref_seqs[rid1]));
          EmplaceBackMappingRecord(
              read_id, read2_name, barcode_key, 1, ref_start_position2, rid2,
              flag2, second_read_direction == kPositive ? 1 : 0, is_unique,
              mapq, NM2, n_cigar2, cigar2, MD_tag2, effect_read2,
              read_batch2.GetSequenceQualAt(pair_index),
              &(mappings_on_diff_ref_seqs[rid2]));
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
                                     &(mappings_on_diff_ref_seqs[rid1]));
          } else {
            EmplaceBackMappingRecord(read_id, read1_name, barcode_key, rid2,
                                     rid1, position2, position1, direction2,
                                     direction, mapq, is_unique, 1,
                                     &(mappings_on_diff_ref_seqs[rid2]));
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
              negative_alignment_length, &(mappings_on_diff_ref_seqs[rid1]));
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
              negative_alignment_length, &(mappings_on_diff_ref_seqs[rid1]));
        }
        num_best_mappings_reported++;
        if (num_best_mappings_reported ==
            std::min(max_num_best_mappings_, num_best_mappings)) {
          break;
        }
      }
      best_mapping_index++;
    }
  }
}

// The returned coordinate is left closed and right closed, and is without chrom
// id.
template <typename MappingRecord>
void MappingGenerator<MappingRecord>::GetRefStartEndPositionForReadFromMapping(
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
            error_threshold_, min_num_errors,
            reference.GetSequenceAt(rid) + verification_window_start_position,
            read, read_length, &mapping_start_position);
      } else {
        BandedTraceback(
            error_threshold_, actual_num_errors,
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
            error_threshold_, min_num_errors,
            reference.GetSequenceAt(rid) + verification_window_start_position,
            read + read_start_site, read_length, &mapping_start_position);
      } else {
        // BandedTracebackToEnd(error_threshold_,actual_num_errors,
        // reference.GetSequenceAt(rid)
        // + verification_window_start_position, read + read_start_site,
        // read_length, &mapping_end_position);
        BandedAlignPatternToText(
            error_threshold_,
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
void MappingGenerator<MappingRecord>::EmplaceBackMappingRecord(
    uint32_t read_id, uint64_t barcode, uint32_t fragment_start_position,
    uint16_t fragment_length, uint8_t mapq, uint8_t direction,
    uint8_t is_unique, uint8_t num_dups, uint16_t positive_alignment_length,
    uint16_t negative_alignment_length,
    std::vector<MappingRecord> *mappings_on_diff_ref_seqs) {}

template <>
void MappingGenerator<PairedEndMappingWithoutBarcode>::EmplaceBackMappingRecord(
    uint32_t read_id, uint64_t barcode, uint32_t fragment_start_position,
    uint16_t fragment_length, uint8_t mapq, uint8_t direction,
    uint8_t is_unique, uint8_t num_dups, uint16_t positive_alignment_length,
    uint16_t negative_alignment_length,
    std::vector<PairedEndMappingWithoutBarcode> *mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs->emplace_back(PairedEndMappingWithoutBarcode(
      read_id, fragment_start_position, fragment_length, mapq, direction,
      is_unique, num_dups, positive_alignment_length,
      negative_alignment_length));
}

template <>
void MappingGenerator<PairedEndMappingWithBarcode>::EmplaceBackMappingRecord(
    uint32_t read_id, uint64_t barcode, uint32_t fragment_start_position,
    uint16_t fragment_length, uint8_t mapq, uint8_t direction,
    uint8_t is_unique, uint8_t num_dups, uint16_t positive_alignment_length,
    uint16_t negative_alignment_length,
    std::vector<PairedEndMappingWithBarcode> *mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs->emplace_back(PairedEndMappingWithBarcode(
      read_id, barcode, fragment_start_position, fragment_length, mapq,
      direction, is_unique, num_dups, positive_alignment_length,
      negative_alignment_length));
}

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::EmplaceBackMappingRecord(
    uint32_t read_id, const char *read1_name, const char *read2_name,
    uint16_t read1_length, uint16_t read2_length, uint64_t barcode,
    uint32_t fragment_start_position, uint16_t fragment_length, uint8_t mapq1,
    uint8_t mapq2, uint8_t direction, uint8_t is_unique, uint8_t num_dups,
    uint16_t positive_alignment_length, uint16_t negative_alignment_length,
    std::vector<MappingRecord> *mappings_on_diff_ref_seqs) {}

template <>
void MappingGenerator<PairedPAFMapping>::EmplaceBackMappingRecord(
    uint32_t read_id, const char *read1_name, const char *read2_name,
    uint16_t read1_length, uint16_t read2_length, uint64_t barcode,
    uint32_t fragment_start_position, uint16_t fragment_length, uint8_t mapq1,
    uint8_t mapq2, uint8_t direction, uint8_t is_unique, uint8_t num_dups,
    uint16_t positive_alignment_length, uint16_t negative_alignment_length,
    std::vector<PairedPAFMapping> *mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs->emplace_back(PairedPAFMapping(
      read_id, std::string(read1_name), std::string(read2_name), read1_length,
      read2_length, fragment_start_position, fragment_length,
      positive_alignment_length, negative_alignment_length,
      mapq1 < mapq2 ? mapq1 : mapq2, mapq1, mapq2, direction, is_unique,
      num_dups));
}

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::EmplaceBackMappingRecord(
    uint32_t read_id, const char *read_name, uint64_t cell_barcode, int rid1,
    int rid2, uint32_t pos1, uint32_t pos2, int direction1, int direction2,
    uint8_t mapq, uint8_t is_unique, uint8_t num_dups,
    std::vector<MappingRecord> *mappings_on_diff_ref_seqs) {}

template <>
void MappingGenerator<PairsMapping>::EmplaceBackMappingRecord(
    uint32_t read_id, const char *read_name, uint64_t cell_barcode, int rid1,
    int rid2, uint32_t pos1, uint32_t pos2, int direction1, int direction2,
    uint8_t mapq, uint8_t is_unique, uint8_t num_dups,
    std::vector<PairsMapping> *mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs->emplace_back(PairsMapping{
      read_id, std::string(read_name), cell_barcode, rid1, rid2, pos1, pos2,
      direction1, direction2, mapq, is_unique, num_dups});
}

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::EmplaceBackMappingRecord(
    uint32_t read_id, uint64_t barcode, uint32_t fragment_start_position,
    uint16_t fragment_length, uint8_t mapq, uint8_t direction,
    uint8_t is_unique, uint8_t num_dups,
    std::vector<MappingRecord> *mappings_on_diff_ref_seqs) {}

template <>
void MappingGenerator<MappingWithoutBarcode>::EmplaceBackMappingRecord(
    uint32_t read_id, uint64_t barcode, uint32_t fragment_start_position,
    uint16_t fragment_length, uint8_t mapq, uint8_t direction,
    uint8_t is_unique, uint8_t num_dups,
    std::vector<MappingWithoutBarcode> *mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs->emplace_back(
      MappingWithoutBarcode{read_id, fragment_start_position, fragment_length,
                            mapq, direction, is_unique, num_dups});
}

template <>
void MappingGenerator<MappingWithBarcode>::EmplaceBackMappingRecord(
    uint32_t read_id, uint64_t barcode, uint32_t fragment_start_position,
    uint16_t fragment_length, uint8_t mapq, uint8_t direction,
    uint8_t is_unique, uint8_t num_dups,
    std::vector<MappingWithBarcode> *mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs->emplace_back(MappingWithBarcode{
      read_id, barcode, fragment_start_position, fragment_length, mapq,
      direction, is_unique, num_dups});
}

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::EmplaceBackMappingRecord(
    uint32_t read_id, const char *read_name, uint16_t read_length,
    uint64_t barcode, uint32_t fragment_start_position,
    uint16_t fragment_length, uint8_t mapq, uint8_t direction,
    uint8_t is_unique, uint8_t num_dups,
    std::vector<MappingRecord> *mappings_on_diff_ref_seqs) {}

template <>
void MappingGenerator<PAFMapping>::EmplaceBackMappingRecord(
    uint32_t read_id, const char *read_name, uint16_t read_length,
    uint64_t barcode, uint32_t fragment_start_position,
    uint16_t fragment_length, uint8_t mapq, uint8_t direction,
    uint8_t is_unique, uint8_t num_dups,
    std::vector<PAFMapping> *mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs->emplace_back(PAFMapping{
      read_id, std::string(read_name), read_length, fragment_start_position,
      fragment_length, mapq, direction, is_unique, num_dups});
}

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::EmplaceBackMappingRecord(
    uint32_t read_id, const char *read_name, uint64_t cell_barcode,
    uint8_t num_dups, int64_t position, int rid, int flag, uint8_t direction,
    uint8_t is_unique, uint8_t mapq, uint32_t NM, int n_cigar, uint32_t *cigar,
    std::string &MD_tag, const char *read, const char *read_qual,
    std::vector<MappingRecord> *mappings_on_diff_ref_seqs) {}

template <>
void MappingGenerator<SAMMapping>::EmplaceBackMappingRecord(
    uint32_t read_id, const char *read_name, uint64_t cell_barcode,
    uint8_t num_dups, int64_t position, int rid, int flag, uint8_t direction,
    uint8_t is_unique, uint8_t mapq, uint32_t NM, int n_cigar, uint32_t *cigar,
    std::string &MD_tag, const char *read, const char *read_qual,
    std::vector<SAMMapping> *mappings_on_diff_ref_seqs) {
  mappings_on_diff_ref_seqs->emplace_back(
      read_id, std::string(read_name), cell_barcode, num_dups, position, rid,
      flag, direction, 0, is_unique, mapq, NM, n_cigar, cigar, MD_tag,
      std::string(read), std::string(read_qual));
}

template <typename MappingRecord>
uint8_t MappingGenerator<MappingRecord>::GetMAPQForSingleEndRead(
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
uint8_t MappingGenerator<MappingRecord>::GetMAPQForPairedEndRead(
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

}  // namespace chromap

#endif  // MAPPING_GENERATOR_H_
