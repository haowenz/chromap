#ifndef OUTPUTTOOLS_H_
#define OUTPUTTOOLS_H_

#include <assert.h>

#include <cinttypes>
#include <cstring>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "barcode_translator.h"
#include "bed_mapping.h"
#include "mapping.h"
#include "mapping_parameters.h"
#include "paf_mapping.h"
#include "pairs_mapping.h"
#include "sam_mapping.h"
#include "sequence_batch.h"

namespace chromap {

template <typename MappingRecord>
bool ReadIdLess(const std::pair<uint32_t, MappingRecord> &a,
                const std::pair<uint32_t, MappingRecord> &b) {
  return a.second.read_id_ < b.second.read_id_;
}

template <typename MappingRecord>
class OutputTools {
 public:
  OutputTools() {}
  ~OutputTools() {}

  // Output the mappings in a temp file.
  inline void OutputTempMapping(
      const std::string &temp_mapping_output_file_path,
      uint32_t num_reference_sequences,
      const std::vector<std::vector<MappingRecord> > &mappings) {
    FILE *temp_mapping_output_file =
        fopen(temp_mapping_output_file_path.c_str(), "wb");
    assert(temp_mapping_output_file != NULL);
    for (size_t ri = 0; ri < num_reference_sequences; ++ri) {
      // make sure mappings[ri] exists even if its size is 0
      size_t num_mappings = mappings[ri].size();
      fwrite(&num_mappings, sizeof(size_t), 1, temp_mapping_output_file);
      if (mappings[ri].size() > 0) {
        fwrite(mappings[ri].data(), sizeof(MappingRecord), mappings[ri].size(),
               temp_mapping_output_file);
      }
    }
    fclose(temp_mapping_output_file);
  }

  inline void InitializeMappingOutput(
      uint32_t cell_barcode_length, const std::string &mapping_output_file_path,
      MappingOutputFormat &format) {
    cell_barcode_length_ = cell_barcode_length;
    mapping_output_file_path_ = mapping_output_file_path;
    mapping_output_file_ = fopen(mapping_output_file_path_.c_str(), "w");
    assert(mapping_output_file_ != NULL);
    mapping_output_format_ = format;
  }

  inline void SetBarcodeLength(uint32_t cell_barcode_length) {
    cell_barcode_length_ = cell_barcode_length;
  }

  inline void FinalizeMappingOutput() { fclose(mapping_output_file_); }

  inline void AppendMappingOutput(const std::string &line) {
    fprintf(mapping_output_file_, "%s", line.data());
  }

  void OutputHeader(uint32_t num_reference_sequences,
                    const SequenceBatch &reference);

  void AppendMapping(uint32_t rid, const SequenceBatch &reference,
                     const MappingRecord &mapping);

  inline uint32_t GetNumMappings() const { return num_mappings_; }

  inline std::string Seed2Sequence(uint64_t seed, uint32_t seed_length) const {
    std::string sequence;
    sequence.reserve(seed_length);
    uint64_t mask = 3;
    for (uint32_t i = 0; i < seed_length; ++i) {
      sequence.push_back(SequenceBatch::Uint8ToChar(
          (seed >> ((seed_length - 1 - i) * 2)) & mask));
    }
    return sequence;
  }

  inline void SetPairsCustomRidRank(const std::vector<int> &custom_rid_rank) {
    custom_rid_rank_ = custom_rid_rank;
  }

  inline void SetBarcodeTranslateTable(std::string &file) {
    barcode_translator_.SetTranslateTable(file);
  }

  std::vector<int> custom_rid_rank_;  // for pairs
 protected:
  std::string mapping_output_file_path_;
  FILE *mapping_output_file_;
  // TODO(Haowen): use this variable to decide output in BED or TagAlign. It
  // should be removed later.
  MappingOutputFormat mapping_output_format_ = MAPPINGFORMAT_BED;
  uint32_t num_mappings_;
  uint32_t cell_barcode_length_ = 16;
  BarcodeTranslator barcode_translator_;
};

// Specialization for BED format.
template <>
void OutputTools<MappingWithBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference);

template <>
void OutputTools<MappingWithBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const MappingWithBarcode &mapping);

template <>
void OutputTools<MappingWithoutBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference);

template <>
void OutputTools<MappingWithoutBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const MappingWithoutBarcode &mapping);

// Specialization for BEDPE format.
template <>
void OutputTools<PairedEndMappingWithoutBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference);

template <>
void OutputTools<PairedEndMappingWithoutBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedEndMappingWithoutBarcode &mapping);

template <>
void OutputTools<PairedEndMappingWithBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference);

template <>
void OutputTools<PairedEndMappingWithBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedEndMappingWithBarcode &mapping);

// Specialization for PAF format.
template <>
void OutputTools<PAFMapping>::OutputHeader(uint32_t num_reference_sequences,
                                           const SequenceBatch &reference);

template <>
void OutputTools<PAFMapping>::AppendMapping(uint32_t rid,
                                            const SequenceBatch &reference,
                                            const PAFMapping &mapping);

template <>
void OutputTools<PAFMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PAFMapping> > &mappings);

// Specialization for PairedPAF format.
template <>
void OutputTools<PairedPAFMapping>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference);

template <>
void OutputTools<PairedPAFMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PairedPAFMapping> > &mappings);

template <>
void OutputTools<PairedPAFMapping>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedPAFMapping &mapping);

// Specialization for SAM format.
template <>
void OutputTools<SAMMapping>::OutputHeader(uint32_t num_reference_sequences,
                                           const SequenceBatch &reference);

template <>
void OutputTools<SAMMapping>::AppendMapping(uint32_t rid,
                                            const SequenceBatch &reference,
                                            const SAMMapping &mapping);

template <>
void OutputTools<SAMMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<SAMMapping> > &mappings);

// Specialization for pairs format.
template <>
void OutputTools<PairsMapping>::OutputHeader(uint32_t num_reference_sequences,
                                             const SequenceBatch &reference);

template <>
void OutputTools<PairsMapping>::AppendMapping(uint32_t rid,
                                              const SequenceBatch &reference,
                                              const PairsMapping &mapping);

template <>
void OutputTools<PairsMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PairsMapping> > &mappings);

}  // namespace chromap

#endif  // OUTPUTTOOLS_H_
