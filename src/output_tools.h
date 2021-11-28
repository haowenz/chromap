#ifndef OUTPUTTOOLS_H_
#define OUTPUTTOOLS_H_

#include <assert.h>

#include <cinttypes>
#include <cstring>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "bed_mapping.h"
#include "mapping.h"
#include "paf_mapping.h"
#include "pairs_mapping.h"
#include "sam_mapping.h"
#include "sequence_batch.h"

namespace chromap {

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

template <typename MappingRecord>
struct TempMappingFileHandle {
  std::string file_path;
  FILE *file;
  bool all_loaded;
  uint32_t num_mappings;
  uint32_t block_size;
  uint32_t current_rid;
  uint32_t current_mapping_index;
  uint32_t num_mappings_on_current_rid;
  uint32_t num_loaded_mappings_on_current_rid;
  std::vector<MappingRecord>
      mappings;  // this vector only keep mappings on the same ref seq
  inline void InitializeTempMappingLoading(uint32_t num_reference_sequences) {
    file = fopen(file_path.c_str(), "rb");
    assert(file != NULL);
    all_loaded = false;
    current_rid = 0;
    fread(&num_mappings_on_current_rid, sizeof(size_t), 1, file);
    mappings.resize(block_size);
    num_loaded_mappings_on_current_rid = 0;
    // std::cerr << "Block size: " << block_size << ", initialize temp file " <<
    // file_path << "\n";
  }
  inline void FinalizeTempMappingLoading() { fclose(file); }
  inline void LoadTempMappingBlock(uint32_t num_reference_sequences) {
    num_mappings = 0;
    while (num_mappings == 0) {
      // Only keep mappings on one ref seq, which means # mappings in buffer can
      // be less than block size Two cases: current ref seq has remainings or
      // not
      if (num_loaded_mappings_on_current_rid < num_mappings_on_current_rid) {
        // Check if # remains larger than block size
        uint32_t num_mappings_to_load_on_current_rid =
            num_mappings_on_current_rid - num_loaded_mappings_on_current_rid;
        if (num_mappings_to_load_on_current_rid > block_size) {
          num_mappings_to_load_on_current_rid = block_size;
        }
        // std::cerr << num_mappings_to_load_on_current_rid << " " <<
        // num_loaded_mappings_on_current_rid << " " <<
        // num_mappings_on_current_rid << "\n"; std::cerr << mappings.size() <<
        // "\n";
        fread(mappings.data(), sizeof(MappingRecord),
              num_mappings_to_load_on_current_rid, file);
        // std::cerr << "Load mappings\n";
        num_loaded_mappings_on_current_rid +=
            num_mappings_to_load_on_current_rid;
        num_mappings = num_mappings_to_load_on_current_rid;
      } else {
        // Move to next rid
        ++current_rid;
        if (current_rid < num_reference_sequences) {
          // std::cerr << "Load size\n";
          fread(&num_mappings_on_current_rid, sizeof(size_t), 1, file);
          // std::cerr << "Load size " << num_mappings_on_current_rid << "\n";
          num_loaded_mappings_on_current_rid = 0;
        } else {
          all_loaded = true;
          break;
        }
      }
    }
    current_mapping_index = 0;
  }
  inline void Next(uint32_t num_reference_sequences) {
    ++current_mapping_index;
    if (current_mapping_index >= num_mappings) {
      LoadTempMappingBlock(num_reference_sequences);
    }
  }
};

template <>
inline void TempMappingFileHandle<PAFMapping>::LoadTempMappingBlock(
    uint32_t num_reference_sequences) {
  num_mappings = 0;
  while (num_mappings == 0) {
    // Only keep mappings on one ref seq, which means # mappings in buffer can
    // be less than block size Two cases: current ref seq has remainings or not
    if (num_loaded_mappings_on_current_rid < num_mappings_on_current_rid) {
      // Check if # remains larger than block size
      uint32_t num_mappings_to_load_on_current_rid =
          num_mappings_on_current_rid - num_loaded_mappings_on_current_rid;
      if (num_mappings_to_load_on_current_rid > block_size) {
        num_mappings_to_load_on_current_rid = block_size;
      }
      // std::cerr << num_mappings_to_load_on_current_rid << " " <<
      // num_loaded_mappings_on_current_rid << " " <<
      // num_mappings_on_current_rid
      // << "\n"; std::cerr << mappings.size() << "\n";
      for (size_t mi = 0; mi < num_mappings_to_load_on_current_rid; ++mi) {
        mappings[mi].LoadFromFile(file);
      }
      // fread(mappings.data(), sizeof(MappingRecord),
      // num_mappings_to_load_on_current_rid, file); std::cerr << "Load
      // mappings\n";
      num_loaded_mappings_on_current_rid += num_mappings_to_load_on_current_rid;
      num_mappings = num_mappings_to_load_on_current_rid;
    } else {
      // Move to next rid
      ++current_rid;
      if (current_rid < num_reference_sequences) {
        // std::cerr << "Load size\n";
        fread(&num_mappings_on_current_rid, sizeof(size_t), 1, file);
        // std::cerr << "Load size " << num_mappings_on_current_rid << "\n";
        num_loaded_mappings_on_current_rid = 0;
      } else {
        all_loaded = true;
        break;
      }
    }
  }
  current_mapping_index = 0;
}

template <>
inline void TempMappingFileHandle<PairedPAFMapping>::LoadTempMappingBlock(
    uint32_t num_reference_sequences) {
  num_mappings = 0;
  while (num_mappings == 0) {
    // Only keep mappings on one ref seq, which means # mappings in buffer can
    // be less than block size Two cases: current ref seq has remainings or not
    if (num_loaded_mappings_on_current_rid < num_mappings_on_current_rid) {
      // Check if # remains larger than block size
      uint32_t num_mappings_to_load_on_current_rid =
          num_mappings_on_current_rid - num_loaded_mappings_on_current_rid;
      if (num_mappings_to_load_on_current_rid > block_size) {
        num_mappings_to_load_on_current_rid = block_size;
      }
      // std::cerr << num_mappings_to_load_on_current_rid << " " <<
      // num_loaded_mappings_on_current_rid << " " <<
      // num_mappings_on_current_rid
      // << "\n"; std::cerr << mappings.size() << "\n";
      for (size_t mi = 0; mi < num_mappings_to_load_on_current_rid; ++mi) {
        mappings[mi].LoadFromFile(file);
      }
      // fread(mappings.data(), sizeof(MappingRecord),
      // num_mappings_to_load_on_current_rid, file); std::cerr << "Load
      // mappings\n";
      num_loaded_mappings_on_current_rid += num_mappings_to_load_on_current_rid;
      num_mappings = num_mappings_to_load_on_current_rid;
    } else {
      // Move to next rid
      ++current_rid;
      if (current_rid < num_reference_sequences) {
        // std::cerr << "Load size\n";
        fread(&num_mappings_on_current_rid, sizeof(size_t), 1, file);
        // std::cerr << "Load size " << num_mappings_on_current_rid << "\n";
        num_loaded_mappings_on_current_rid = 0;
      } else {
        all_loaded = true;
        break;
      }
    }
  }
  current_mapping_index = 0;
}

template <>
inline void TempMappingFileHandle<SAMMapping>::LoadTempMappingBlock(
    uint32_t num_reference_sequences) {
  num_mappings = 0;
  while (num_mappings == 0) {
    // Only keep mappings on one ref seq, which means # mappings in buffer can
    // be less than block size Two cases: current ref seq has remainings or not
    if (num_loaded_mappings_on_current_rid < num_mappings_on_current_rid) {
      // Check if # remains larger than block size
      uint32_t num_mappings_to_load_on_current_rid =
          num_mappings_on_current_rid - num_loaded_mappings_on_current_rid;
      if (num_mappings_to_load_on_current_rid > block_size) {
        num_mappings_to_load_on_current_rid = block_size;
      }
      // std::cerr << num_mappings_to_load_on_current_rid << " " <<
      // num_loaded_mappings_on_current_rid << " " <<
      // num_mappings_on_current_rid
      // << "\n"; std::cerr << mappings.size() << "\n";
      for (size_t mi = 0; mi < num_mappings_to_load_on_current_rid; ++mi) {
        mappings[mi].LoadFromFile(file);
      }
      // fread(mappings.data(), sizeof(MappingRecord),
      // num_mappings_to_load_on_current_rid, file); std::cerr << "Load
      // mappings\n";
      num_loaded_mappings_on_current_rid += num_mappings_to_load_on_current_rid;
      num_mappings = num_mappings_to_load_on_current_rid;
    } else {
      // Move to next rid
      ++current_rid;
      if (current_rid < num_reference_sequences) {
        // std::cerr << "Load size\n";
        fread(&num_mappings_on_current_rid, sizeof(size_t), 1, file);
        // std::cerr << "Load size " << num_mappings_on_current_rid << "\n";
        num_loaded_mappings_on_current_rid = 0;
      } else {
        all_loaded = true;
        break;
      }
    }
  }
  current_mapping_index = 0;
}

template <>
inline void TempMappingFileHandle<PairsMapping>::LoadTempMappingBlock(
    uint32_t num_reference_sequences) {
  num_mappings = 0;
  while (num_mappings == 0) {
    // Only keep mappings on one ref seq, which means # mappings in buffer can
    // be less than block size Two cases: current ref seq has remainings or not
    if (num_loaded_mappings_on_current_rid < num_mappings_on_current_rid) {
      // Check if # remains larger than block size
      uint32_t num_mappings_to_load_on_current_rid =
          num_mappings_on_current_rid - num_loaded_mappings_on_current_rid;
      if (num_mappings_to_load_on_current_rid > block_size) {
        num_mappings_to_load_on_current_rid = block_size;
      }
      // std::cerr << num_mappings_to_load_on_current_rid << " " <<
      // num_loaded_mappings_on_current_rid << " " <<
      // num_mappings_on_current_rid
      // << "\n"; std::cerr << mappings.size() << "\n";
      for (size_t mi = 0; mi < num_mappings_to_load_on_current_rid; ++mi) {
        mappings[mi].LoadFromFile(file);
      }
      // fread(mappings.data(), sizeof(MappingRecord),
      // num_mappings_to_load_on_current_rid, file); std::cerr << "Load
      // mappings\n";
      num_loaded_mappings_on_current_rid += num_mappings_to_load_on_current_rid;
      num_mappings = num_mappings_to_load_on_current_rid;
    } else {
      // Move to next rid
      ++current_rid;
      if (current_rid < num_reference_sequences) {
        // std::cerr << "Load size\n";
        fread(&num_mappings_on_current_rid, sizeof(size_t), 1, file);
        // std::cerr << "Load size " << num_mappings_on_current_rid << "\n";
        num_loaded_mappings_on_current_rid = 0;
      } else {
        all_loaded = true;
        break;
      }
    }
  }
  current_mapping_index = 0;
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
            Seed2Sequence(barcode_key, cell_barcode_length_).data());
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
