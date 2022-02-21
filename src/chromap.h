#ifndef CHROMAP_H_
#define CHROMAP_H_

#include <memory>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#include "feature_barcode_matrix.h"
#include "index.h"
#include "index_parameters.h"
#include "khash.h"
#include "mapping_metadata.h"
#include "mapping_parameters.h"
#include "mapping_processor.h"
#include "mapping_writer.h"
#include "paired_end_mapping_metadata.h"
#include "sequence_batch.h"
#include "temp_mapping.h"
#include "utils.h"

#define CHROMAP_VERSION "0.1.6-r310"

namespace chromap {

class ChromapDriver {
 public:
  ChromapDriver() {}
  ~ChromapDriver() {}
  void ParseArgsAndRun(int argc, char *argv[]);
};

template <typename MappingRecord>
class Chromap {
 public:
  // For index construction
  Chromap(const IndexParameters &index_parameters)
      : index_parameters_(index_parameters) {
    barcode_lookup_table_ = NULL;
    barcode_whitelist_lookup_table_ = NULL;
  }

  // For mapping
  Chromap(const MappingParameters &mapping_parameters)
      : mapping_parameters_(mapping_parameters) {
    barcode_lookup_table_ = kh_init(k64_seq);
    barcode_whitelist_lookup_table_ = kh_init(k64_seq);

    ParseReadFormat(mapping_parameters.read_format);
  }

  ~Chromap() {
    if (barcode_whitelist_lookup_table_ != NULL) {
      kh_destroy(k64_seq, barcode_whitelist_lookup_table_);
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

  void ConstructIndex();

  void MapSingleEndReads();

  void MapPairedEndReads();

 private:
  uint32_t SampleInputBarcodesAndExamineLength();

  uint32_t LoadSingleEndReadsWithBarcodes(SequenceBatch &read_batch,
                                          SequenceBatch &barcode_batch);

  uint32_t LoadPairedEndReadsWithBarcodes(SequenceBatch &read_batch1,
                                          SequenceBatch &read_batch2,
                                          SequenceBatch &barcode_batch);

  void TrimAdapterForPairedEndRead(uint32_t pair_index,
                                   SequenceBatch &read_batch1,
                                   SequenceBatch &read_batch2);

  bool PairedEndReadWithBarcodeIsDuplicate(uint32_t pair_index,
                                           const SequenceBatch &barcode_batch,
                                           const SequenceBatch &read_batch1,
                                           const SequenceBatch &read_batch2);

  void LoadBarcodeWhitelist();

  void ComputeBarcodeAbundance(uint64_t max_num_sample_barcodes);

  void UpdateBarcodeAbundance(uint32_t num_loaded_barcodes,
                              const SequenceBatch &barcode_batch);

  bool CorrectBarcodeAt(uint32_t barcode_index, SequenceBatch &barcode_batch,
                        uint64_t &num_barcode_in_whitelist,
                        uint64_t &num_corrected_barcode);

  void OutputBarcodeStatistics();

  void OutputMappingStatistics();

  void OutputMappingStatistics(uint32_t num_reference_sequences,
                               const std::vector<std::vector<MappingRecord> >
                                   &mappings_on_diff_ref_seqs);

  void GenerateCustomizedRidRank(const std::string &rid_order_path,
                                 uint32_t num_reference_sequences,
                                 const SequenceBatch &reference,
                                 std::vector<int> &rid_rank);

  void RerankCandidatesRid(std::vector<Candidate> &candidates);

  void ParseReadFormat(const std::string &read_format);

  // Parameters
  const IndexParameters index_parameters_;
  const MappingParameters mapping_parameters_;
  // Default batch size, # reads for single-end reads, # read pairs for
  // paired-end reads.
  const uint32_t read_batch_size_ = 500000;
  // 0-start, 1-end (includsive), 2-strand(-1:minus, 1:plus)
  int barcode_format_[3];
  int read1_format_[3];
  int read2_format_[3];
  std::vector<int> custom_rid_rank_;
  std::vector<int> pairs_custom_rid_rank_;
  khash_t(k64_seq) * barcode_whitelist_lookup_table_;
  // For identical read dedupe
  khash_t(k64_seq) * barcode_lookup_table_;
  std::vector<khash_t(k128) *> read_lookup_tables_;
  // For mapping
  const int min_unique_mapping_mapq_ = 4;
  // For mapping stats.
  uint64_t num_candidates_ = 0;
  uint64_t num_mappings_ = 0;
  uint64_t num_mapped_reads_ = 0;
  uint64_t num_uniquely_mapped_reads_ = 0;
  uint64_t num_reads_ = 0;
  uint64_t num_duplicated_reads_ = 0;  // # identical reads
  // For barcode stats
  const uint64_t initial_num_sample_barcodes_ = 20000000;
  uint64_t num_sample_barcodes_ = 0;
  uint64_t num_barcode_in_whitelist_ = 0;
  uint64_t num_corrected_barcode_ = 0;
  uint32_t barcode_length_ = 0;
  bool skip_barcode_check_ = false;
};

}  // namespace chromap

#endif  // CHROMAP_H_
