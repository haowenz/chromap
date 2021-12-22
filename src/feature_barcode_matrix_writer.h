#ifndef FEATUREBARCODEMATRIXWRITER_H_
#define FEATUREBARCODEMATRIXWRITER_H_

#include <assert.h>

#include <cinttypes>
#include <cstring>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "barcode_translator.h"
#include "sequence_batch.h"

namespace chromap {

// The code here is not working properly since the barcode length is not set.
// But this feature is not used in the realse for now so this is fine.
class FeatureBarcodeMatrixWriter {
 public:
  FeatureBarcodeMatrixWriter() {}
  ~FeatureBarcodeMatrixWriter() {}

  inline void InitializeMatrixOutput(const std::string &matrix_output_prefix) {
    matrix_output_prefix_ = matrix_output_prefix;
    matrix_output_file_ =
        fopen((matrix_output_prefix_ + "_matrix.mtx").c_str(), "w");
    assert(matrix_output_file_ != nullptr);
    peak_output_file_ =
        fopen((matrix_output_prefix_ + "_peaks.bed").c_str(), "w");
    assert(peak_output_file_ != nullptr);
    barcode_output_file_ =
        fopen((matrix_output_prefix_ + "_barcode.tsv").c_str(), "w");
    assert(barcode_output_file_ != nullptr);
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
            barcode_translator_.Translate(barcode_key, cell_barcode_length_)
                .data());
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

  inline void SetBarcodeTranslateTable(const std::string &file) {
    barcode_translator_.SetTranslateTable(file);
  }

  inline void SetBarcodeLength(uint32_t cell_barcode_length) {
    cell_barcode_length_ = cell_barcode_length;
  }

 protected:
  uint32_t cell_barcode_length_ = 16;
  std::string matrix_output_prefix_;
  FILE *peak_output_file_ = nullptr;
  FILE *barcode_output_file_ = nullptr;
  FILE *matrix_output_file_ = nullptr;
  BarcodeTranslator barcode_translator_;
};

}  // namespace chromap

#endif  // FEATUREBARCODEMATRIXWRITER_H_
