#ifndef OUTPUTTOOLS_H_
#define OUTPUTTOOLS_H_

#include <assert.h>

#include <cinttypes>
#include <cstring>
#include <functional>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "bed_mapping.h"
#include "mapping.h"
#include "paf_mapping.h"
#include "pairs_mapping.h"
#include "sam_mapping.h"
#include "sequence_batch.h"
#include "khash.h"

namespace chromap {

KHASH_INIT(k64_str, uint64_t, std::string, 1, kh_int64_hash_func, kh_int64_hash_equal);

enum MappingOutputFormat {
  MAPPINGFORMAT_UNKNOWN,
  MAPPINGFORMAT_BED,
  MAPPINGFORMAT_TAGALIGN,
  MAPPINGFORMAT_PAF,
  MAPPINGFORMAT_SAM,
  MAPPINGFORMAT_PAIRS
};

template <typename MappingRecord>
bool ReadIdLess(const std::pair<uint32_t, MappingRecord> &a,
                const std::pair<uint32_t, MappingRecord> &b) {
  return a.second.read_id_ < b.second.read_id_;
}

// The class for handling barcode convertion.
class BarcodeTranslator
{
private:
  khash_t(k64_str) *barcode_translate_table;
  int from_bc_length;
  
  std::string Seed2Sequence(uint64_t seed, uint32_t seed_length) const {
    std::string sequence;
    sequence.reserve(seed_length);
    uint64_t mask = 3;
    for (uint32_t i = 0; i < seed_length; ++i) {
      sequence.push_back(SequenceBatch::Uint8ToChar(
          (seed >> ((seed_length - 1 - i) * 2)) & mask));
    }
    return sequence;
  }
  
  void ProcessTranslateFileLine(std::string &line) {
    int i;
    int len = line.length();
    std::string to;
    for (i = 0; i < len; ++i) {
      if (line[i] == ',' || line[i] == '\t')
        break;
    }

    to = line.substr(0, i);
    //from = line.substr(i + 1, len - i - 1);
    from_bc_length = len - i - 1;
    uint64_t from_seed = SequenceBatch::GenerateSeedFromSequence(line.c_str(), len, i + 1, from_bc_length);
    
    int khash_return_code;
    khiter_t barcode_translate_table_iter = kh_put(k64_str, barcode_translate_table, from_seed, &khash_return_code);
    kh_value(barcode_translate_table, barcode_translate_table_iter) = to;
  }

public:
  BarcodeTranslator() {
    barcode_translate_table = NULL;
    from_bc_length = -1;
  }
  
  ~BarcodeTranslator() {
    if (barcode_translate_table != NULL) {
      kh_destroy(k64_str, barcode_translate_table);
    }
  }

  void SetTranslateTable(const std::string &file) {
    barcode_translate_table = kh_init(k64_str);
    std::ifstream file_stream(file);
    std::string file_line;

    while (getline(file_stream, file_line)) {
      ProcessTranslateFileLine(file_line); 
    }
  }
  
  std::string Translate(uint64_t bc, uint32_t bc_length) {
    if (barcode_translate_table == NULL) {
      return Seed2Sequence(bc, bc_length);
    } 

    std::string ret;  
    uint64_t i;
    for (i = 0; i < bc_length / from_bc_length; ++i) {
      uint64_t seed = (bc << (i * from_bc_length)) >> (bc_length - from_bc_length);
      khiter_t barcode_translate_table_iter = kh_get(k64_str, barcode_translate_table, seed);
      if (barcode_translate_table_iter == kh_end(barcode_translate_table)) {
        std::cerr << "Barcode does not exist in the translation table.\n" << std::endl;
        exit(-1);
      }
      const std::string &bc_to = kh_value(barcode_translate_table, barcode_translate_table_iter);
      if (i == 0) {
        ret = bc_to;
      } else {
        ret += "-" + bc_to;
      }
    }
    return ret;
  }
};

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

  // Below are functions to output feature matrix.
  inline void InitializeMatrixOutput(const std::string &matrix_output_prefix) {
    matrix_output_prefix_ = matrix_output_prefix;
    matrix_output_file_ =
        fopen((matrix_output_prefix_ + "_matrix.mtx").c_str(), "w");
    assert(matrix_output_file_ != NULL);
    peak_output_file_ =
        fopen((matrix_output_prefix_ + "_peaks.bed").c_str(), "w");
    assert(peak_output_file_ != NULL);
    barcode_output_file_ =
        fopen((matrix_output_prefix_ + "_barcode.tsv").c_str(), "w");
    assert(barcode_output_file_ != NULL);
  }

  void OutputPeaks(uint32_t bin_size, uint32_t num_sequences,
                   const SequenceBatch &reference) {
    for (uint32_t rid = 0; rid < num_sequences; ++rid) {
      uint32_t sequence_length = reference.GetSequenceLengthAt(rid);
      const char *sequence_name = reference.GetSequenceNameAt(rid);
      for (uint32_t position = 0; position < sequence_length;
           position += bin_size) {
        fprintf(peak_output_file_, "%s\t%u\t%u\n", sequence_name, position + 1,
                position + bin_size);
      }
    }
  }

  void OutputPeaks(uint32_t peak_start_position, uint16_t peak_length,
                   uint32_t rid, const SequenceBatch &reference) {
    const char *sequence_name = reference.GetSequenceNameAt(rid);
    fprintf(peak_output_file_, "%s\t%u\t%u\n", sequence_name,
            peak_start_position + 1, peak_start_position + peak_length);
  }

  void AppendBarcodeOutput(uint64_t barcode_key) {
    fprintf(barcode_output_file_, "%s-1\n",
            barcode_translator.Translate(barcode_key, cell_barcode_length_).data());
            //Seed2Sequence(barcode_key, cell_barcode_length_).data());
  }

  void WriteMatrixOutputHead(uint64_t num_peaks, uint64_t num_barcodes,
                             uint64_t num_lines) {
    fprintf(matrix_output_file_, "%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\n",
            num_peaks, num_barcodes, num_lines);
  }

  void AppendMatrixOutput(uint32_t peak_index, uint32_t barcode_index,
                          uint32_t num_mappings) {
    fprintf(matrix_output_file_, "%u\t%u\t%u\n", peak_index, barcode_index,
            num_mappings);
  }

  inline void FinalizeMatrixOutput() {
    fclose(matrix_output_file_);
    fclose(peak_output_file_);
    fclose(barcode_output_file_);
  }

  inline void SetPairsCustomRidRank(const std::vector<int> &custom_rid_rank) {
    custom_rid_rank_ = custom_rid_rank;
  }

  inline void SetBarcodeTranslateTable(std::string &file) {
    barcode_translator.SetTranslateTable(file);
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
  std::string matrix_output_prefix_;
  FILE *peak_output_file_;
  FILE *barcode_output_file_;
  FILE *matrix_output_file_;
  BarcodeTranslator barcode_translator;
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
