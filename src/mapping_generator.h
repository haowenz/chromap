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
#include "mapping_in_memory.h"
#include "mapping_metadata.h"
#include "mapping_parameters.h"
#include "paf_mapping.h"
#include "paired_end_mapping_metadata.h"
#include "pairs_mapping.h"
#include "sam_mapping.h"
#include "sequence_batch.h"
#include "utils.h"

namespace chromap {

// Class to process draft mappings and generate best mappings and alignments. It
// supports multi-threadidng as only the parameters are owned by the class.
template <typename MappingRecord>
class MappingGenerator {
 public:
  MappingGenerator() = delete;

  MappingGenerator(const MappingParameters &mapping_parameters,
                   const std::vector<int> &pairs_custom_rid_rank)
      : mapping_parameters_(mapping_parameters),
        pairs_custom_rid_rank_(pairs_custom_rid_rank) {}

  ~MappingGenerator() = default;

  void GenerateBestMappingsForSingleEndRead(
      const SequenceBatch &read_batch, uint32_t read_index,
      const SequenceBatch &reference, const SequenceBatch &barcode_batch,
      MappingMetadata &mapping_metadata,
      std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs);

  // When the number of supplemented candidates is greater than 0, the
  // force_mapq will be 0, and thereby setting the mapq to 0 (not mapq1 or
  // mapq2). Split alignment won't run candidate supplement.
  void GenerateBestMappingsForPairedEndRead(
      uint32_t pair_index, const SequenceBatch &read_batch1,
      const SequenceBatch &read_batch2, const SequenceBatch &barcode_batch,
      const SequenceBatch &reference, std::vector<int> &best_mapping_indices,
      std::mt19937 &generator, int force_mapq,
      PairedEndMappingMetadata &paired_end_mapping_metadata,
      std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs);

 private:
  void ProcessBestMappingsForSingleEndRead(
      const Strand mapping_strand, uint32_t read_index,
      const SequenceBatch &read_batch, const SequenceBatch &barcode_batch,
      const SequenceBatch &reference, const MappingMetadata &mapping_metadata,
      const std::vector<int> &best_mapping_indices, int &best_mapping_index,
      int &num_best_mappings_reported,
      std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs);

  void GenerateBestMappingsForPairedEndReadOnOneDirection(
      const Strand first_read_strand, const Strand second_read_strand,
      uint32_t pair_index, const SequenceBatch &read_batch1,
      const SequenceBatch &read_batch2, const SequenceBatch &reference,
      PairedEndMappingMetadata &paired_end_mapping_metadata);

  void ProcessBestMappingsForPairedEndReadOnOneDirection(
      const Strand first_read_strand, const Strand second_read_strand,
      uint32_t pair_index, const SequenceBatch &read_batch1,
      const SequenceBatch &read_batch2, const SequenceBatch &barcode_batch,
      const SequenceBatch &reference,
      const std::vector<int> &best_mapping_indices, int &best_mapping_index,
      int &num_best_mappings_reported, int force_mapq,
      const PairedEndMappingMetadata &paired_end_mapping_metadata,
      std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs);

  void GetRefStartEndPositionForReadFromMapping(
      const DraftMapping &mapping, const SequenceBatch &reference,
      MappingInMemory &mapping_in_memory);

  // For single-end. It should be fully specialized.
  void EmplaceBackSingleEndMappingRecord(
      MappingInMemory &mapping_in_memory,
      std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs);

  // For paired-end. It should be fully specialized.
  void EmplaceBackPairedEndMappingRecord(
      PairedEndMappingInMemory &paired_mapping_in_memory,
      std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs);

  uint8_t GetMAPQForSingleEndRead(const Strand strand, int num_errors,
                                  uint16_t alignment_length, int read_length,
                                  int max_num_error_difference,
                                  const MappingMetadata &mapping_metadata);

  uint8_t GetMAPQForPairedEndRead(
      const Strand first_read_strand, const Strand second_read_strand,
      int read1_num_errors, int read2_num_errors,
      uint16_t read1_alignment_length, uint16_t read2_alignment_length,
      int read1_length, int read2_length, int force_mapq,
      const PairedEndMappingMetadata &paired_end_mapping_metadata,
      uint8_t &mapq1, uint8_t &mapq2);

  const MappingParameters mapping_parameters_;
  const std::vector<int> pairs_custom_rid_rank_;
};

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::GenerateBestMappingsForSingleEndRead(
    const SequenceBatch &read_batch, uint32_t read_index,
    const SequenceBatch &reference, const SequenceBatch &barcode_batch,
    MappingMetadata &mapping_metadata,
    std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs) {
  const int num_best_mappings = mapping_metadata.num_best_mappings_;

  // We use reservoir sampling when the number of best mappings exceeds the
  // threshold.
  std::vector<int> best_mapping_indices(
      mapping_parameters_.max_num_best_mappings);
  std::iota(best_mapping_indices.begin(), best_mapping_indices.end(), 0);
  if (num_best_mappings > mapping_parameters_.max_num_best_mappings) {
    std::mt19937 generator(11);
    for (int i = mapping_parameters_.max_num_best_mappings;
         i < num_best_mappings; ++i) {
      // Important: inclusive range.
      std::uniform_int_distribution<int> distribution(0, i);
      int j = distribution(generator);
      if (j < mapping_parameters_.max_num_best_mappings) {
        best_mapping_indices[j] = i;
      }
    }
    std::sort(best_mapping_indices.begin(), best_mapping_indices.end());
  }

  int best_mapping_index = 0;
  int num_best_mappings_reported = 0;
  const int num_best_mappings_to_report =
      std::min(num_best_mappings, mapping_parameters_.max_num_best_mappings);

  ProcessBestMappingsForSingleEndRead(
      kPositive, read_index, read_batch, barcode_batch, reference,
      mapping_metadata, best_mapping_indices, best_mapping_index,
      num_best_mappings_reported, mappings_on_diff_ref_seqs);

  if (num_best_mappings_reported != num_best_mappings_to_report) {
    ProcessBestMappingsForSingleEndRead(
        kNegative, read_index, read_batch, barcode_batch, reference,
        mapping_metadata, best_mapping_indices, best_mapping_index,
        num_best_mappings_reported, mappings_on_diff_ref_seqs);
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
  paired_end_mapping_metadata.SetMinSumErrors(
      2 * mapping_parameters_.error_threshold + 1);
  paired_end_mapping_metadata.SetNumBestMappings(0);
  paired_end_mapping_metadata.SetSecondMinSumErrors(
      2 * mapping_parameters_.error_threshold + 1);
  paired_end_mapping_metadata.SetNumSecondBestMappings(0);

  GenerateBestMappingsForPairedEndReadOnOneDirection(
      kPositive, kNegative, pair_index, read_batch1, read_batch2, reference,
      paired_end_mapping_metadata);
  GenerateBestMappingsForPairedEndReadOnOneDirection(
      kNegative, kPositive, pair_index, read_batch1, read_batch2, reference,
      paired_end_mapping_metadata);

  if (mapping_parameters_.split_alignment) {
    GenerateBestMappingsForPairedEndReadOnOneDirection(
        kPositive, kPositive, pair_index, read_batch1, read_batch2, reference,
        paired_end_mapping_metadata);
    GenerateBestMappingsForPairedEndReadOnOneDirection(
        kNegative, kNegative, pair_index, read_batch1, read_batch2, reference,
        paired_end_mapping_metadata);
  }

  if (paired_end_mapping_metadata.GetNumBestMappings() >
      mapping_parameters_.drop_repetitive_reads) {
    return;
  }

  // We use reservoir sampling when the number of best mappings exceeds the
  // threshold.
  // std::vector<int>
  // best_mapping_indices(mapping_parameters_.max_num_best_mappings);
  std::iota(best_mapping_indices.begin(), best_mapping_indices.end(), 0);
  if (paired_end_mapping_metadata.GetNumBestMappings() >
      mapping_parameters_.max_num_best_mappings) {
    // std::mt19937 generator(11);
    for (int i = mapping_parameters_.max_num_best_mappings;
         i < paired_end_mapping_metadata.GetNumBestMappings(); ++i) {
      // Important: inclusive range.
      std::uniform_int_distribution<int> distribution(0, i);
      int j = distribution(generator);
      // int j = distribution(tmp_generator);
      if (j < mapping_parameters_.max_num_best_mappings) {
        best_mapping_indices[j] = i;
      }
    }
    std::sort(best_mapping_indices.begin(), best_mapping_indices.end());
  }

  int best_mapping_index = 0;
  int num_best_mappings_reported = 0;
  const int num_best_mappings_to_report =
      std::min(mapping_parameters_.max_num_best_mappings,
               paired_end_mapping_metadata.GetNumBestMappings());

  ProcessBestMappingsForPairedEndReadOnOneDirection(
      kPositive, kNegative, pair_index, read_batch1, read_batch2, barcode_batch,
      reference, best_mapping_indices, best_mapping_index,
      num_best_mappings_reported, force_mapq, paired_end_mapping_metadata,
      mappings_on_diff_ref_seqs);

  if (num_best_mappings_reported != num_best_mappings_to_report) {
    ProcessBestMappingsForPairedEndReadOnOneDirection(
        kNegative, kPositive, pair_index, read_batch1, read_batch2,
        barcode_batch, reference, best_mapping_indices, best_mapping_index,
        num_best_mappings_reported, force_mapq, paired_end_mapping_metadata,
        mappings_on_diff_ref_seqs);
  }

  if (mapping_parameters_.split_alignment &&
      num_best_mappings_reported != num_best_mappings_to_report) {
    ProcessBestMappingsForPairedEndReadOnOneDirection(
        kPositive, kPositive, pair_index, read_batch1, read_batch2,
        barcode_batch, reference, best_mapping_indices, best_mapping_index,
        num_best_mappings_reported, force_mapq, paired_end_mapping_metadata,
        mappings_on_diff_ref_seqs);
  }

  if (mapping_parameters_.split_alignment &&
      num_best_mappings_reported != num_best_mappings_to_report) {
    ProcessBestMappingsForPairedEndReadOnOneDirection(
        kNegative, kNegative, pair_index, read_batch1, read_batch2,
        barcode_batch, reference, best_mapping_indices, best_mapping_index,
        num_best_mappings_reported, force_mapq, paired_end_mapping_metadata,
        mappings_on_diff_ref_seqs);
  }
}

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::ProcessBestMappingsForSingleEndRead(
    const Strand mapping_strand, uint32_t read_index,
    const SequenceBatch &read_batch, const SequenceBatch &barcode_batch,
    const SequenceBatch &reference, const MappingMetadata &mapping_metadata,
    const std::vector<int> &best_mapping_indices, int &best_mapping_index,
    int &num_best_mappings_reported,
    std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs) {
  const std::vector<DraftMapping> &mappings =
      mapping_strand == kPositive ? mapping_metadata.positive_mappings_
                                  : mapping_metadata.negative_mappings_;
  const std::vector<int> &split_sites =
      mapping_strand == kPositive ? mapping_metadata.positive_split_sites_
                                  : mapping_metadata.negative_split_sites_;

  const char *read = read_batch.GetSequenceAt(read_index);
  const uint32_t read_id = read_batch.GetSequenceIdAt(read_index);
  const char *read_name = read_batch.GetSequenceNameAt(read_index);
  const uint32_t read_length = read_batch.GetSequenceLengthAt(read_index);
  const std::string &negative_read =
      read_batch.GetNegativeSequenceAt(read_index);

  MappingInMemory mapping_in_memory;
  mapping_in_memory.read_id = read_id;
  mapping_in_memory.read_name = read_name;
  mapping_in_memory.is_unique = (mapping_metadata.num_best_mappings_ == 1);

  uint64_t barcode_key = 0;
  if (!mapping_parameters_.is_bulk_data) {
    barcode_key = barcode_batch.GenerateSeedFromSequenceAt(
        read_index, /*start_position=*/0,
        barcode_batch.GetSequenceLengthAt(read_index));
  }
  mapping_in_memory.barcode_key = barcode_key;

  mapping_in_memory.strand = mapping_strand;
  mapping_in_memory.read_sequence =
      mapping_strand == kPositive ? read : negative_read.data();
  mapping_in_memory.read_length = read_length;

  for (uint32_t mi = 0; mi < mappings.size(); ++mi) {
    if (mappings[mi].GetNumErrors() > mapping_metadata.min_num_errors_) {
      continue;
    }

    if (best_mapping_index ==
        best_mapping_indices[num_best_mappings_reported]) {
      mapping_in_memory.rid = mappings[mi].GetReferenceSequenceIndex();

      if (mapping_parameters_.split_alignment) {
        mapping_in_memory.read_split_site = split_sites[mi];
      }

      GetRefStartEndPositionForReadFromMapping(mappings[mi], reference,
                                               mapping_in_memory);

      const uint16_t alignment_length = mapping_in_memory.ref_end_position -
                                        mapping_in_memory.ref_start_position +
                                        1;
      const uint8_t mapq = GetMAPQForSingleEndRead(
          mapping_strand, /*num_errors=*/mappings[mi].GetNumErrors(),
          alignment_length, read_length,
          /*max_num_error_difference=*/mapping_parameters_.error_threshold,
          mapping_metadata);
      mapping_in_memory.mapq = mapq;

      if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_SAM) {
        uint16_t flag = mapping_strand == kPositive ? 0 : BAM_FREVERSE;
        if (num_best_mappings_reported >= 1) {
          flag |= BAM_FSECONDARY;
        }
        mapping_in_memory.SAM_flag = flag;
        mapping_in_memory.qual_sequence =
            read_batch.GetSequenceQualAt(read_index);
      }

      EmplaceBackSingleEndMappingRecord(mapping_in_memory,
                                        mappings_on_diff_ref_seqs);

      num_best_mappings_reported++;
      if (num_best_mappings_reported ==
          std::min(mapping_parameters_.max_num_best_mappings,
                   mapping_metadata.num_best_mappings_)) {
        break;
      }
    }

    best_mapping_index++;
  }
}

template <typename MappingRecord>
void MappingGenerator<MappingRecord>::
    GenerateBestMappingsForPairedEndReadOnOneDirection(
        const Strand first_read_strand, const Strand second_read_strand,
        uint32_t pair_index, const SequenceBatch &read_batch1,
        const SequenceBatch &read_batch2, const SequenceBatch &reference,
        PairedEndMappingMetadata &paired_end_mapping_metadata) {
  uint32_t i1 = 0;
  uint32_t i2 = 0;
  uint32_t min_overlap_length = mapping_parameters_.min_read_length;
  uint32_t read1_length = read_batch1.GetSequenceLengthAt(pair_index);
  uint32_t read2_length = read_batch2.GetSequenceLengthAt(pair_index);

  const std::vector<DraftMapping> &mappings1 =
      first_read_strand == kPositive
          ? paired_end_mapping_metadata.mapping_metadata1_.positive_mappings_
          : paired_end_mapping_metadata.mapping_metadata1_.negative_mappings_;
  const std::vector<DraftMapping> &mappings2 =
      second_read_strand == kPositive
          ? paired_end_mapping_metadata.mapping_metadata2_.positive_mappings_
          : paired_end_mapping_metadata.mapping_metadata2_.negative_mappings_;

  std::vector<std::pair<uint32_t, uint32_t>> &best_mappings =
      paired_end_mapping_metadata.GetBestMappings(first_read_strand,
                                                  second_read_strand);
  int &min_sum_errors = paired_end_mapping_metadata.min_sum_errors_;
  int &num_best_mappings = paired_end_mapping_metadata.num_best_mappings_;
  int &second_min_sum_errors =
      paired_end_mapping_metadata.second_min_sum_errors_;
  int &num_second_best_mappings =
      paired_end_mapping_metadata.num_second_best_mappings_;

#ifdef LI_DEBUG
  for (int i = 0; i < mappings1.size(); ++i)
    printf("mappings1 %d %d:%d\n", i,
           (int)(mappings1[i].GetReferenceSequenceIndex()),
           (int)mappings1[i].GetReferenceSequencePosition());
  for (int i = 0; i < mappings1.size(); ++i)
    printf("mappings2 %d %d:%d\n", i,
           (int)(mappings2[i].GetReferenceSequenceIndex()),
           (int)mappings2[i].GetReferenceSequencePosition());
#endif

  if (mapping_parameters_.split_alignment) {
    if (mappings1.size() == 0 || mappings2.size() == 0) {
      return;
    }
    // For split alignment, selecting the pairs whose both single-end are the
    // best.
    for (i1 = 0; i1 < mappings1.size(); ++i1) {
      if (mappings1[i1].GetNumErrors() !=
          paired_end_mapping_metadata.mapping_metadata1_.min_num_errors_) {
        continue;
      }
      for (i2 = 0; i2 < mappings2.size(); ++i2) {
        if (mappings2[i2].GetNumErrors() !=
            paired_end_mapping_metadata.mapping_metadata2_.min_num_errors_) {
          continue;
        }
        best_mappings.emplace_back(i1, i2);
        min_sum_errors =
            paired_end_mapping_metadata.mapping_metadata1_.min_num_errors_ +
            paired_end_mapping_metadata.mapping_metadata2_.min_num_errors_;
        //*second_min_sum_errors = min_num_errors1 + min_num_errors2 + 1;
        num_best_mappings++;
      }
    }

    return;
  }

  while (i1 < mappings1.size() && i2 < mappings2.size()) {
    if ((first_read_strand == kNegative &&
         mappings1[i1].position > mappings2[i2].position +
                                      mapping_parameters_.max_insert_size -
                                      read1_length) ||
        (first_read_strand == kPositive &&
         mappings1[i1].position >
             mappings2[i2].position + read2_length - min_overlap_length)) {
      ++i2;
    } else if ((first_read_strand == kPositive &&
                mappings2[i2].position >
                    mappings1[i1].position +
                        mapping_parameters_.max_insert_size - read2_length) ||
               (first_read_strand == kNegative &&
                mappings2[i2].position > mappings1[i1].position + read1_length -
                                             min_overlap_length)) {
      ++i1;
    } else {
      // Ok, find a pair, we store current ni2 somewhere and keep looking until
      // we go out of the range, then we go back and then move to next pi1 and
      // keep doing the similar thing.
      uint32_t current_i2 = i2;
      while (
          current_i2 < mappings2.size() &&
          ((first_read_strand == kPositive &&
            mappings2[current_i2].position <=
                mappings1[i1].position + mapping_parameters_.max_insert_size -
                    read2_length) ||
           (first_read_strand == kNegative &&
            mappings2[current_i2].position <=
                mappings1[i1].position + read1_length - min_overlap_length))) {
#ifdef LI_DEBUG
        printf(
            "%s passed: %llu %d %llu %d: %d %d %d\n", __func__,
            mappings1[i1].GetReferenceSequenceIndex(),
            int(mappings1[i1].GetReferenceSequencePosition()),
            mappings2[current_i2].GetReferenceSequenceIndex(),
            int(mappings2[current_i2].GetReferenceSequencePosition()),
            mappings1[i1].GetNumErrors() + mappings2[current_i2].GetNumErrors(),
            mappings1[i1].GetNumErrors(), mappings2[current_i2].GetNumErrors());
#endif

        int current_sum_errors =
            mappings1[i1].GetNumErrors() + mappings2[current_i2].GetNumErrors();
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
    ProcessBestMappingsForPairedEndReadOnOneDirection(
        const Strand first_read_strand, const Strand second_read_strand,
        uint32_t pair_index, const SequenceBatch &read_batch1,
        const SequenceBatch &read_batch2, const SequenceBatch &barcode_batch,
        const SequenceBatch &reference,
        const std::vector<int> &best_mapping_indices, int &best_mapping_index,
        int &num_best_mappings_reported, int force_mapq,
        const PairedEndMappingMetadata &paired_end_mapping_metadata,
        std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs) {
  PairedEndMappingInMemory paired_end_mapping_in_memory;

  paired_end_mapping_in_memory.mapping_in_memory1.strand = first_read_strand;
  paired_end_mapping_in_memory.mapping_in_memory2.strand = second_read_strand;

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

  paired_end_mapping_in_memory.mapping_in_memory1.read_id = read_id;
  paired_end_mapping_in_memory.mapping_in_memory2.read_id = read_id;

  paired_end_mapping_in_memory.mapping_in_memory1.read_name = read1_name;
  paired_end_mapping_in_memory.mapping_in_memory2.read_name = read2_name;

  paired_end_mapping_in_memory.mapping_in_memory1.read_length = read1_length;
  paired_end_mapping_in_memory.mapping_in_memory2.read_length = read2_length;

  const MappingMetadata &mapping_metadata1 =
      paired_end_mapping_metadata.mapping_metadata1_;
  const MappingMetadata &mapping_metadata2 =
      paired_end_mapping_metadata.mapping_metadata2_;

  const std::vector<DraftMapping> &mappings1 =
      first_read_strand == kPositive ? mapping_metadata1.positive_mappings_
                                     : mapping_metadata1.negative_mappings_;
  const std::vector<DraftMapping> &mappings2 =
      second_read_strand == kPositive ? mapping_metadata2.positive_mappings_
                                      : mapping_metadata2.negative_mappings_;

  const std::vector<int> &split_sites1 =
      first_read_strand == kPositive ? mapping_metadata1.positive_split_sites_
                                     : mapping_metadata1.negative_split_sites_;
  const std::vector<int> &split_sites2 =
      second_read_strand == kPositive ? mapping_metadata2.positive_split_sites_
                                      : mapping_metadata2.negative_split_sites_;

  const std::vector<std::pair<uint32_t, uint32_t>> &best_mappings =
      paired_end_mapping_metadata.GetBestMappings(first_read_strand,
                                                  second_read_strand);

  const uint8_t is_unique =
      (paired_end_mapping_metadata.num_best_mappings_ == 1 ||
       mapping_metadata1.num_best_mappings_ == 1 ||
       mapping_metadata2.num_best_mappings_ == 1)
          ? 1
          : 0;
  paired_end_mapping_in_memory.is_unique = is_unique;

  uint64_t barcode_key = 0;
  if (!mapping_parameters_.is_bulk_data) {
    barcode_key = barcode_batch.GenerateSeedFromSequenceAt(
        pair_index, /*start_position=*/0,
        barcode_batch.GetSequenceLengthAt(pair_index));
  }
  paired_end_mapping_in_memory.mapping_in_memory1.barcode_key = barcode_key;
  paired_end_mapping_in_memory.mapping_in_memory2.barcode_key = barcode_key;

  for (uint32_t mi = 0; mi < best_mappings.size(); ++mi) {
    const uint32_t i1 = best_mappings[mi].first;
    const uint32_t i2 = best_mappings[mi].second;
    const int current_sum_errors =
        mappings1[i1].GetNumErrors() + mappings2[i2].GetNumErrors();

    if (current_sum_errors > paired_end_mapping_metadata.min_sum_errors_) {
      continue;
    }

    if (best_mapping_index ==
        best_mapping_indices[num_best_mappings_reported]) {
      const uint32_t rid1 = mappings1[i1].GetReferenceSequenceIndex();
      const uint32_t rid2 = mappings2[i2].GetReferenceSequenceIndex();

      paired_end_mapping_in_memory.mapping_in_memory1.rid = rid1;
      paired_end_mapping_in_memory.mapping_in_memory2.rid = rid2;

      paired_end_mapping_in_memory.mapping_in_memory1.read_sequence =
          first_read_strand == kPositive ? read1 : negative_read1.data();
      paired_end_mapping_in_memory.mapping_in_memory2.read_sequence =
          second_read_strand == kPositive ? read2 : negative_read2.data();

      if (mapping_parameters_.split_alignment) {
        paired_end_mapping_in_memory.mapping_in_memory1.read_split_site =
            split_sites1[i1];
        paired_end_mapping_in_memory.mapping_in_memory2.read_split_site =
            split_sites2[i2];
      }

      GetRefStartEndPositionForReadFromMapping(
          mappings1[i1], reference,
          paired_end_mapping_in_memory.mapping_in_memory1);
      GetRefStartEndPositionForReadFromMapping(
          mappings2[i2], reference,
          paired_end_mapping_in_memory.mapping_in_memory2);

      uint8_t mapq1 = 0;
      uint8_t mapq2 = 0;
      const uint8_t mapq = GetMAPQForPairedEndRead(
          first_read_strand, second_read_strand,
          /*read1_num_errors=*/mappings1[i1].GetNumErrors(),
          /*read2_num_errors=*/mappings2[i2].GetNumErrors(),
          paired_end_mapping_in_memory.mapping_in_memory1.GetFragmentLength(),
          paired_end_mapping_in_memory.mapping_in_memory2.GetFragmentLength(),
          read1_length, read2_length, force_mapq, paired_end_mapping_metadata,
          mapq1, mapq2);
      paired_end_mapping_in_memory.mapq = mapq;
      paired_end_mapping_in_memory.mapping_in_memory1.mapq = mapq;
      paired_end_mapping_in_memory.mapping_in_memory2.mapq = mapq;

      if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_SAM) {
        uint16_t flag1 = 3;
        uint16_t flag2 = 3;
        if (first_read_strand == kNegative) {
          flag1 |= BAM_FREVERSE;
          flag2 |= BAM_FMREVERSE;
        }
        if (second_read_strand == kNegative) {
          flag1 |= BAM_FMREVERSE;
          flag2 |= BAM_FREVERSE;
        }
        flag1 |= BAM_FREAD1;
        flag2 |= BAM_FREAD2;
        if (num_best_mappings_reported >= 1) {
          flag1 |= BAM_FSECONDARY;
          flag2 |= BAM_FSECONDARY;
        }
        paired_end_mapping_in_memory.mapping_in_memory1.SAM_flag = flag1;
        paired_end_mapping_in_memory.mapping_in_memory2.SAM_flag = flag2;
        paired_end_mapping_in_memory.mapping_in_memory1.qual_sequence =
            read_batch1.GetSequenceQualAt(pair_index);
        paired_end_mapping_in_memory.mapping_in_memory2.qual_sequence =
            read_batch2.GetSequenceQualAt(pair_index);
        paired_end_mapping_in_memory.mapping_in_memory1.is_unique = is_unique;
        paired_end_mapping_in_memory.mapping_in_memory2.is_unique = is_unique;
      }

      EmplaceBackPairedEndMappingRecord(paired_end_mapping_in_memory,
                                        mappings_on_diff_ref_seqs);

      num_best_mappings_reported++;
      if (num_best_mappings_reported ==
          std::min(mapping_parameters_.max_num_best_mappings,
                   paired_end_mapping_metadata.num_best_mappings_)) {
        break;
      }
    }

    best_mapping_index++;
  }
}

// The computed ref start and end coordinates are left closed and right closed.
template <typename MappingRecord>
void MappingGenerator<MappingRecord>::GetRefStartEndPositionForReadFromMapping(
    const DraftMapping &mapping, const SequenceBatch &reference,
    MappingInMemory &mapping_in_memory) {
  // For now this mat is only used by ksw to generate mappings in SAM format.
  int8_t mat[25];
  // if (output_mapping_in_SAM_) {
  int i, j, k;
  for (i = k = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j)
      mat[k++] = i == j ? mapping_parameters_.match_score
                        : -mapping_parameters_.mismatch_penalty;
    mat[k++] = 0;  // ambiguous base
  }
  for (j = 0; j < 5; ++j) mat[k++] = 0;
  //}

  const uint32_t rid = mapping.GetReferenceSequenceIndex();
  const uint32_t ref_position = mapping.GetReferenceSequencePosition();

  const int full_read_length = mapping_in_memory.read_length;
  int read_length = mapping_in_memory.read_length;

  const int min_num_errors = mapping.GetNumErrors();

  int split_site =
      mapping_in_memory.strand == kPositive ? 0 : mapping_in_memory.read_length;

  int gap_beginning = 0;
  int actual_num_errors = 0;

  if (mapping_parameters_.split_alignment) {
    split_site = mapping_in_memory.read_split_site & 0xffff;
    // Beginning means the 5' end of the read.
    gap_beginning = (mapping_in_memory.read_split_site >> 16) & 0xff;
    // In split alignment, -num_errors is the number of matches.
    actual_num_errors = (mapping_in_memory.read_split_site >> 24) & 0xff;
    read_length = split_site - gap_beginning;
  }

  uint32_t verification_window_start_position =
      ref_position + 1 >
              (uint32_t)(read_length + mapping_parameters_.error_threshold)
          ? ref_position + 1 - read_length - mapping_parameters_.error_threshold
          : 0;

  if (ref_position + mapping_parameters_.error_threshold >=
      reference.GetSequenceLengthAt(rid)) {
    verification_window_start_position = reference.GetSequenceLengthAt(rid) -
                                         mapping_parameters_.error_threshold -
                                         read_length;
  }

  if (verification_window_start_position < 0) {
    verification_window_start_position = 0;
  }

  if (mapping_parameters_.split_alignment) {
    if (split_site < full_read_length &&
        mapping_parameters_.mapping_output_format == MAPPINGFORMAT_SAM &&
        split_site > 3 * mapping_parameters_.error_threshold) {
      split_site -= 3 * mapping_parameters_.error_threshold;
    }
    read_length = split_site - gap_beginning;
  }

  if (mapping_in_memory.strand == kPositive) {
    if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_SAM) {
      mapping_in_memory.n_cigar = 0;

      int mapping_start_position = 0;
      int mapping_end_position = 0;

      ksw_semi_global3(
          read_length + 2 * mapping_parameters_.error_threshold,
          reference.GetSequenceAt(rid) + verification_window_start_position,
          read_length, mapping_in_memory.read_sequence + gap_beginning, 5, mat,
          mapping_parameters_.gap_open_penalties[0],
          mapping_parameters_.gap_extension_penalties[0],
          mapping_parameters_.gap_open_penalties[1],
          mapping_parameters_.gap_extension_penalties[1],
          mapping_parameters_.error_threshold * 2 + 1,
          &(mapping_in_memory.n_cigar), &(mapping_in_memory.cigar),
          &mapping_start_position, &mapping_end_position);

      if (gap_beginning > 0) {
        int new_ref_start_position = AdjustGapBeginning(
            mapping_in_memory.strand, reference.GetSequenceAt(rid),
            mapping_in_memory.read_sequence, &gap_beginning, read_length - 1,
            verification_window_start_position + mapping_start_position,
            verification_window_start_position + mapping_end_position - 1,
            &(mapping_in_memory.n_cigar), &(mapping_in_memory.cigar));
        mapping_start_position =
            new_ref_start_position - verification_window_start_position;
      }

      GenerateNMAndMDTag(
          reference.GetSequenceAt(rid),
          mapping_in_memory.read_sequence + gap_beginning,
          verification_window_start_position + mapping_start_position,
          mapping_in_memory);

      mapping_in_memory.ref_start_position =
          verification_window_start_position + mapping_start_position;
      mapping_in_memory.ref_end_position =
          verification_window_start_position + mapping_end_position - 1;
    } else {
      int mapping_start_position = 0;
      if (!mapping_parameters_.split_alignment) {
        BandedTraceback(
            mapping_parameters_.error_threshold, min_num_errors,
            reference.GetSequenceAt(rid) + verification_window_start_position,
            mapping_in_memory.read_sequence, read_length,
            &mapping_start_position);
      } else {
        BandedTraceback(
            mapping_parameters_.error_threshold, actual_num_errors,
            reference.GetSequenceAt(rid) + verification_window_start_position,
            mapping_in_memory.read_sequence + gap_beginning, read_length,
            &mapping_start_position);
      }

      if (gap_beginning > 0) {
        int new_ref_start_position = AdjustGapBeginning(
            mapping_in_memory.strand, reference.GetSequenceAt(rid),
            mapping_in_memory.read_sequence, &gap_beginning, read_length - 1,
            verification_window_start_position + mapping_start_position,
            ref_position, &(mapping_in_memory.n_cigar),
            &(mapping_in_memory.cigar));

        mapping_start_position =
            new_ref_start_position - verification_window_start_position;
      }

      mapping_in_memory.ref_start_position =
          verification_window_start_position + mapping_start_position;
      mapping_in_memory.ref_end_position = ref_position;
    }

    return;
  }

  //  reversed read looks like:
  //
  //      veri_start_pos       ref_position
  //  ref   --|-------------------|------------------->
  //  read     <-|--read_length---|--gap_beginning--
  //          split_site
  //

  const int read_start_site = full_read_length - split_site;
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_SAM) {
    mapping_in_memory.n_cigar = 0;

    int mapping_start_position = 0;
    int mapping_end_position = 0;

    ksw_semi_global3(read_length + 2 * mapping_parameters_.error_threshold,
                     reference.GetSequenceAt(rid) +
                         verification_window_start_position + read_start_site,
                     read_length,
                     mapping_in_memory.read_sequence + read_start_site, 5, mat,
                     mapping_parameters_.gap_open_penalties[0],
                     mapping_parameters_.gap_extension_penalties[0],
                     mapping_parameters_.gap_open_penalties[1],
                     mapping_parameters_.gap_extension_penalties[1],
                     mapping_parameters_.error_threshold * 2 + 1,
                     &(mapping_in_memory.n_cigar), &(mapping_in_memory.cigar),
                     &mapping_start_position, &mapping_end_position);

    if (gap_beginning > 0) {
      int new_ref_end_position = AdjustGapBeginning(
          mapping_in_memory.strand, reference.GetSequenceAt(rid),
          mapping_in_memory.read_sequence + read_start_site, &gap_beginning,
          read_length - 1,
          verification_window_start_position + mapping_start_position,
          verification_window_start_position + mapping_end_position - 1,
          &(mapping_in_memory.n_cigar), &(mapping_in_memory.cigar));

      // The returned position is right-closed, so need to plus one to match
      // bed convention
      mapping_end_position = new_ref_end_position + 1 -
                             verification_window_start_position -
                             read_start_site;
      read_length = split_site - gap_beginning;
    }

    GenerateNMAndMDTag(reference.GetSequenceAt(rid),
                       mapping_in_memory.read_sequence + read_start_site,
                       verification_window_start_position + read_start_site +
                           mapping_start_position,
                       mapping_in_memory);

    mapping_in_memory.ref_start_position = verification_window_start_position +
                                           read_start_site +
                                           mapping_start_position;
    mapping_in_memory.ref_end_position = verification_window_start_position +
                                         read_start_site +
                                         mapping_end_position - 1;
  } else {
    int mapping_start_position = mapping_parameters_.error_threshold;
    int mapping_end_position =
        ref_position - verification_window_start_position + 1;
    // int n_cigar = 0;
    // uint32_t *cigar;
    // ksw_semi_global3(read_length + 2 * mapping_parameters_.error_threshold,
    // reference.GetSequenceAt(rid) + verification_window_start_position,
    // read_length, negative_read.data() + split_sites[mi], 5, mat,
    // mapping_parameters_.gap_open_penalties[0],
    // mapping_parameters_.gap_extension_penalties[0],
    // mapping_parameters_.gap_open_penalties[1],
    // mapping_parameters_.gap_extension_penalties[1],
    // mapping_parameters_.error_threshold * 2 + 1, &n_cigar, &cigar,
    // &mapping_start_position, &mapping_end_position); mapq =
    // GetMAPQForSingleEndRead(mapping_parameters_.error_threshold, 0, 0,
    // mapping_end_position - mapping_start_position + 1, min_num_errors,
    // num_best_mappings, second_min_num_errors, num_second_best_mappings);
    // uint32_t fragment_start_position = verification_window_start_position +
    // mapping_start_position; uint16_t fragment_length = mapping_end_position
    // - mapping_start_position + 1;
    if (!mapping_parameters_.split_alignment) {
      BandedTraceback(
          mapping_parameters_.error_threshold, min_num_errors,
          reference.GetSequenceAt(rid) + verification_window_start_position,
          mapping_in_memory.read_sequence + read_start_site, read_length,
          &mapping_start_position);
    } else {
      // BandedTracebackToEnd(mapping_parameters_.error_threshold,actual_num_errors,
      // reference.GetSequenceAt(rid)
      // + verification_window_start_position, read + read_start_site,
      // read_length, &mapping_end_position);
      BandedAlignPatternToText(
          mapping_parameters_.error_threshold,
          reference.GetSequenceAt(rid) + verification_window_start_position,
          mapping_in_memory.read_sequence + read_start_site, read_length,
          &mapping_end_position);
      // seems banded align's mapping end position is included?
      mapping_end_position += 1;
    }

    if (gap_beginning > 0) {
      int new_ref_end_position = AdjustGapBeginning(
          mapping_in_memory.strand, reference.GetSequenceAt(rid),
          mapping_in_memory.read_sequence + read_start_site, &gap_beginning,
          read_length - 1,
          verification_window_start_position + mapping_start_position,
          verification_window_start_position + mapping_end_position,
          &(mapping_in_memory.n_cigar), &(mapping_in_memory.cigar));

      // The returned position is right-closed, so need to plus one to match
      // bed convention.
      mapping_end_position =
          new_ref_end_position - verification_window_start_position + 1;
      read_length = split_site - gap_beginning;
    }

    mapping_in_memory.ref_start_position =
        verification_window_start_position + mapping_start_position;
    mapping_in_memory.ref_end_position =
        verification_window_start_position + mapping_end_position - 1;
  }
}

template <typename MappingRecord>
uint8_t MappingGenerator<MappingRecord>::GetMAPQForSingleEndRead(
    const Strand strand, int num_errors, uint16_t alignment_length,
    int read_length, int max_num_error_difference,
    const MappingMetadata &mapping_metadata) {
  int mapq_coef_length = 50;
  int mapq_coef_fraction = log(mapq_coef_length);

  if (!mapping_parameters_.split_alignment) {
    alignment_length =
        alignment_length > read_length ? alignment_length : read_length;
  }

  double alignment_identity = 1 - (double)num_errors / alignment_length;

  if (mapping_parameters_.split_alignment) {
    alignment_identity = (double)(-num_errors) / alignment_length;
    if (alignment_identity > 1) alignment_identity = 1;
  }

  int mapq = 0;
  int second_min_num_errors = mapping_metadata.second_min_num_errors_;

  if (mapping_metadata.num_best_mappings_ > 1) {
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
    if (second_min_num_errors > num_errors + max_num_error_difference) {
      second_min_num_errors = num_errors + max_num_error_difference;
    }

    double tmp = alignment_length < mapq_coef_length
                     ? 1.0
                     : mapq_coef_fraction / log(alignment_length);
    tmp *= alignment_identity * alignment_identity;
    mapq = 5 * 6.02 * (second_min_num_errors - num_errors) * tmp * tmp + 0.499;
  }

  if (mapping_metadata.num_second_best_mappings_ > 0) {
    mapq -= (int)(4.343 * log(mapping_metadata.num_second_best_mappings_ + 1) +
                  0.499);
  }

  if (mapq > 60) {
    mapq = 60;
  }
  if (mapq < 0) {
    mapq = 0;
  }

  if (mapping_metadata.repetitive_seed_length_ > 0) {
    double frac_rep =
        (mapping_metadata.repetitive_seed_length_) / (double)read_length;
    if (mapping_metadata.repetitive_seed_length_ >= (uint32_t)read_length) {
      frac_rep = 0.999;
    }
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

  if (mapping_parameters_.split_alignment &&
      alignment_length < read_length - mapping_parameters_.error_threshold &&
      second_min_num_errors != num_errors) {
    if (mapping_metadata.repetitive_seed_length_ >= alignment_length &&
        mapping_metadata.repetitive_seed_length_ < (uint32_t)read_length &&
        alignment_length < read_length / 3) {
      mapq = 0;
    }
    const int diff = second_min_num_errors - num_errors;
    const uint32_t num_candidates =
        strand == kPositive ? mapping_metadata.positive_candidates_.size()
                            : mapping_metadata.negative_candidates_.size();
    if (second_min_num_errors - num_errors <=
            mapping_parameters_.error_threshold * 3 / 4 &&
        num_candidates >= 5) {
      mapq -= (num_candidates / 5 / diff);
    }
    if (mapq < 0) {
      mapq = 0;
    }
    if (mapping_metadata.num_second_best_mappings_ > 0 &&
        second_min_num_errors - num_errors <=
            mapping_parameters_.error_threshold * 3 / 4) {
      mapq /= (mapping_metadata.num_second_best_mappings_ / diff + 1);
    }
  }

  return (uint8_t)mapq;
}

#define raw_mapq(diff, a) ((int)(5 * 6.02 * (diff) / (a) + .499))

template <typename MappingRecord>
uint8_t MappingGenerator<MappingRecord>::GetMAPQForPairedEndRead(
    const Strand first_read_strand, const Strand second_read_strand,
    int read1_num_errors, int read2_num_errors, uint16_t read1_alignment_length,
    uint16_t read2_alignment_length, int read1_length, int read2_length,
    int force_mapq, const PairedEndMappingMetadata &paired_end_mapping_metadata,
    uint8_t &mapq1, uint8_t &mapq2) {
  const MappingMetadata &mapping_metadata1 =
      paired_end_mapping_metadata.mapping_metadata1_;
  const MappingMetadata &mapping_metadata2 =
      paired_end_mapping_metadata.mapping_metadata2_;

#ifdef CHROMAP_DEBUG
  std::cerr
      << " rl1:"
      << (int)paired_end_mapping_metadata.mapping_metadata1_
             .repetitive_seed_length_
      << " rl2:"
      << (int)paired_end_mapping_metadata.mapping_metadata2_
             .repetitive_seed_length_
      << " pal:" << (int)read1_alignment_length
      << " nal:" << (int)read2_alignment_length
      << " me:" << paired_end_mapping_metadata.min_sum_errors_
      << " #bm:" << paired_end_mapping_metadata.num_best_mappings_
      << " sme:" << paired_end_mapping_metadata.second_min_sum_errors_
      << " #sbm:" << paired_end_mapping_metadata.num_second_best_mappings_
      << " ne1:" << read1_num_errors << " ne2:" << read2_num_errors << " me1:"
      << paired_end_mapping_metadata.mapping_metadata1_.min_num_errors_
      << " me2:"
      << paired_end_mapping_metadata.mapping_metadata2_.min_num_errors_
      << " #bm1:"
      << paired_end_mapping_metadata.mapping_metadata1_.num_best_mappings_
      << " #bm2:"
      << paired_end_mapping_metadata.mapping_metadata2_.num_best_mappings_
      << " sme1:"
      << paired_end_mapping_metadata.mapping_metadata1_.second_min_num_errors_
      << " sme2:"
      << paired_end_mapping_metadata.mapping_metadata2_.second_min_num_errors_
      << " #sbm1:"
      << paired_end_mapping_metadata.mapping_metadata1_
             .num_second_best_mappings_
      << " #sbm2:"
      << paired_end_mapping_metadata.mapping_metadata2_
             .num_second_best_mappings_
      << "\n";
#endif

  uint8_t mapq_pe = 0;
  int min_num_unpaired_sum_errors =
      mapping_metadata1.min_num_errors_ + mapping_metadata2.min_num_errors_ + 3;

  if (paired_end_mapping_metadata.num_best_mappings_ <= 1) {
    int adjusted_second_min_sum_errors =
        paired_end_mapping_metadata.second_min_sum_errors_ <
                min_num_unpaired_sum_errors
            ? paired_end_mapping_metadata.second_min_sum_errors_
            : min_num_unpaired_sum_errors;

    mapq_pe = raw_mapq(adjusted_second_min_sum_errors -
                           paired_end_mapping_metadata.min_sum_errors_,
                       1);

#ifdef CHROMAP_DEBUG
    std::cerr << "mapqpe: " << (int)mapq_pe << "\n";
#endif

    if (paired_end_mapping_metadata.num_second_best_mappings_ > 0) {
      mapq_pe -=
          (int)(4.343 *
                    log(paired_end_mapping_metadata.num_second_best_mappings_ +
                        1) +
                0.499);
    }

    if (mapq_pe > 60) {
      mapq_pe = 60;
    }
    if (mapq_pe < 0) {
      mapq_pe = 0;
    }

#ifdef CHROMAP_DEBUG
    std::cerr << "mapqpe: " << (int)mapq_pe << "\n";
#endif

    int repetitive_seed_length = mapping_metadata1.repetitive_seed_length_ +
                                 mapping_metadata2.repetitive_seed_length_;

    if (repetitive_seed_length > 0) {
      double total_read_length = read1_length + read2_length;
      double frac_rep = (double)repetitive_seed_length / total_read_length;
      if (repetitive_seed_length >= total_read_length) {
        frac_rep = 0.999;
      }

      double alignment_identity1 =
          1 - (double)read1_num_errors / (read1_length > read1_alignment_length
                                              ? read1_length
                                              : read1_alignment_length);

      double alignment_identity2 =
          1 - (double)read2_num_errors / (read2_length > read2_alignment_length
                                              ? read2_length
                                              : read2_alignment_length);

      double alignment_identity = alignment_identity1 < alignment_identity2
                                      ? alignment_identity1
                                      : alignment_identity2;

      if (alignment_identity <= 0.95) {
        mapq_pe = mapq_pe * (1 - sqrt(frac_rep)) + 0.499;
      } else if (alignment_identity <= 0.97) {
        mapq_pe = mapq_pe * (1 - frac_rep) + 0.499;
      } else if (alignment_identity >= 0.999) {
        mapq_pe =
            mapq_pe * (1 - frac_rep * frac_rep * frac_rep * frac_rep) + 0.499;
      } else {
        mapq_pe = mapq_pe * (1 - frac_rep * frac_rep) + 0.499;
      }
    }
  }

  mapq1 = GetMAPQForSingleEndRead(
      first_read_strand, read1_num_errors, read1_alignment_length, read1_length,
      /*max_num_error_difference=*/2, mapping_metadata1);

  mapq2 = GetMAPQForSingleEndRead(second_read_strand, read2_num_errors,
                                  read2_alignment_length, read2_length,
                                  /*max_num_error_difference=*/2,
                                  mapping_metadata2);

#ifdef CHROMAP_DEBUG
  std::cerr << " 1:" << (int)mapq1 << " 2:" << (int)mapq2
            << " mapq_pe:" << (int)mapq_pe << "\n";
#endif

  if (!mapping_parameters_.split_alignment) {
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

#ifdef CHROMAP_DEBUG
  std::cerr << " 1:" << (int)mapq1 << " 2:" << (int)mapq2 << "\n\n";
#endif

  uint8_t mapq = mapq1 < mapq2 ? mapq1 : mapq2;

  if (mapq < 60 && force_mapq >= 0 && force_mapq < mapq) {
    mapq = force_mapq;
  }

  return (uint8_t)mapq;
}

}  // namespace chromap

#endif  // MAPPING_GENERATOR_H_
