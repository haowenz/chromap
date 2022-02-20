#ifndef MAPPING_WRITER_H_
#define MAPPING_WRITER_H_

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
class MappingWriter {
 public:
  MappingWriter() = delete;

  MappingWriter(const std::string &mapping_output_file_path,
                MappingOutputFormat &mapping_output_format,
                uint32_t cell_barcode_length,
                const std::string *barcode_translate_table_file_path,
                const std::vector<int> *custom_rid_rank)
      : mapping_output_file_path_(mapping_output_file_path),
        mapping_output_format_(mapping_output_format),
        cell_barcode_length_(cell_barcode_length),
        custom_rid_rank_(custom_rid_rank != nullptr ? *custom_rid_rank
                                                    : std::vector<int>()) {
    if (barcode_translate_table_file_path != nullptr) {
      barcode_translator_.SetTranslateTable(*barcode_translate_table_file_path);
    }

    mapping_output_file_ = fopen(mapping_output_file_path_.c_str(), "w");
    assert(mapping_output_file_ != nullptr);
  }

  ~MappingWriter() { fclose(mapping_output_file_); }

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

  void OutputHeader(uint32_t num_reference_sequences,
                    const SequenceBatch &reference);

  void AppendMapping(uint32_t rid, const SequenceBatch &reference,
                     const MappingRecord &mapping);

 protected:
  inline void AppendMappingOutput(const std::string &line) {
    fprintf(mapping_output_file_, "%s", line.data());
  }

  const std::string mapping_output_file_path_;
  // TODO(Haowen): use this variable to decide output in BED or TagAlign. It
  // should be removed later.
  const MappingOutputFormat mapping_output_format_;
  const uint32_t cell_barcode_length_;
  FILE *mapping_output_file_ = nullptr;
  BarcodeTranslator barcode_translator_;
  // for pairs
  const std::vector<int> custom_rid_rank_;
};

// Specialization for BED format.
template <>
void MappingWriter<MappingWithBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference);

template <>
void MappingWriter<MappingWithBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const MappingWithBarcode &mapping);

template <>
void MappingWriter<MappingWithoutBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference);

template <>
void MappingWriter<MappingWithoutBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const MappingWithoutBarcode &mapping);

// Specialization for BEDPE format.
template <>
void MappingWriter<PairedEndMappingWithoutBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference);

template <>
void MappingWriter<PairedEndMappingWithoutBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedEndMappingWithoutBarcode &mapping);

template <>
void MappingWriter<PairedEndMappingWithBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference);

template <>
void MappingWriter<PairedEndMappingWithBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedEndMappingWithBarcode &mapping);

// Specialization for PAF format.
template <>
void MappingWriter<PAFMapping>::OutputHeader(uint32_t num_reference_sequences,
                                             const SequenceBatch &reference);

template <>
void MappingWriter<PAFMapping>::AppendMapping(uint32_t rid,
                                              const SequenceBatch &reference,
                                              const PAFMapping &mapping);

template <>
void MappingWriter<PAFMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PAFMapping> > &mappings);

// Specialization for PairedPAF format.
template <>
void MappingWriter<PairedPAFMapping>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference);

template <>
void MappingWriter<PairedPAFMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PairedPAFMapping> > &mappings);

template <>
void MappingWriter<PairedPAFMapping>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedPAFMapping &mapping);

// Specialization for SAM format.
template <>
void MappingWriter<SAMMapping>::OutputHeader(uint32_t num_reference_sequences,
                                             const SequenceBatch &reference);

template <>
void MappingWriter<SAMMapping>::AppendMapping(uint32_t rid,
                                              const SequenceBatch &reference,
                                              const SAMMapping &mapping);

template <>
void MappingWriter<SAMMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<SAMMapping> > &mappings);

// Specialization for pairs format.
template <>
void MappingWriter<PairsMapping>::OutputHeader(uint32_t num_reference_sequences,
                                               const SequenceBatch &reference);

template <>
void MappingWriter<PairsMapping>::AppendMapping(uint32_t rid,
                                                const SequenceBatch &reference,
                                                const PairsMapping &mapping);

template <>
void MappingWriter<PairsMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PairsMapping> > &mappings);

}  // namespace chromap

#endif  // MAPPING_WRITER_H_
