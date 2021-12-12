#include "output_tools.h"

namespace chromap {

// Specialization for BED format.
template <>
void OutputTools<MappingWithBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void OutputTools<MappingWithBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const MappingWithBarcode &mapping) {
  if (mapping_output_format_ == MAPPINGFORMAT_BED) {
    std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t mapping_end_position = mapping.GetEndPosition();
    this->AppendMappingOutput(
        std::string(reference_sequence_name) + "\t" +
        std::to_string(mapping.GetStartPosition()) + "\t" +
        std::to_string(mapping_end_position) + "\t" +
        barcode_translator.Translate(mapping.cell_barcode_, cell_barcode_length_) + "\n");
  } else {
    std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t mapping_end_position = mapping.GetEndPosition();
    this->AppendMappingOutput(std::string(reference_sequence_name) + "\t" +
                              std::to_string(mapping.GetStartPosition()) +
                              "\t" + std::to_string(mapping_end_position) +
                              "\tN\t" + std::to_string(mapping.mapq_) + "\t" +
                              strand + "\n");
  }
}

template <>
void OutputTools<MappingWithoutBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void OutputTools<MappingWithoutBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const MappingWithoutBarcode &mapping) {
  if (mapping_output_format_ == MAPPINGFORMAT_BED) {
    std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t mapping_end_position = mapping.GetEndPosition();
    this->AppendMappingOutput(std::string(reference_sequence_name) + "\t" +
                              std::to_string(mapping.GetStartPosition()) +
                              "\t" + std::to_string(mapping_end_position) +
                              "\tN\t" + std::to_string(mapping.mapq_) + "\t" +
                              strand + "\n");
  } else {
    std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t mapping_end_position = mapping.GetEndPosition();
    this->AppendMappingOutput(std::string(reference_sequence_name) + "\t" +
                              std::to_string(mapping.GetStartPosition()) +
                              "\t" + std::to_string(mapping_end_position) +
                              "\tN\t" + std::to_string(mapping.mapq_) + "\t" +
                              strand + "\n");
  }
}

// Specialization for BEDPE format.
template <>
void OutputTools<PairedEndMappingWithoutBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void OutputTools<PairedEndMappingWithoutBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedEndMappingWithoutBarcode &mapping) {
  if (mapping_output_format_ == MAPPINGFORMAT_BED) {
    std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t mapping_end_position = mapping.GetEndPosition();
    this->AppendMappingOutput(std::string(reference_sequence_name) + "\t" +
                              std::to_string(mapping.GetStartPosition()) +
                              "\t" + std::to_string(mapping_end_position) +
                              "\tN\t" + std::to_string(mapping.mapq_) + "\t" +
                              strand + "\n");
  } else {
    bool positive_strand = mapping.IsPositiveStrand();
    uint32_t positive_read_end =
        mapping.fragment_start_position_ + mapping.positive_alignment_length_;
    uint32_t negative_read_end =
        mapping.fragment_start_position_ + mapping.fragment_length_;
    uint32_t negative_read_start =
        negative_read_end - mapping.negative_alignment_length_;
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    if (positive_strand) {
      this->AppendMappingOutput(
          std::string(reference_sequence_name) + "\t" +
          std::to_string(mapping.fragment_start_position_) + "\t" +
          std::to_string(positive_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t+\n" +
          std::string(reference_sequence_name) + "\t" +
          std::to_string(negative_read_start) + "\t" +
          std::to_string(negative_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t-\n");
    } else {
      this->AppendMappingOutput(
          std::string(reference_sequence_name) + "\t" +
          std::to_string(negative_read_start) + "\t" +
          std::to_string(negative_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t-\n" +
          std::string(reference_sequence_name) + "\t" +
          std::to_string(mapping.fragment_start_position_) + "\t" +
          std::to_string(positive_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t+\n");
    }
  }
}

template <>
void OutputTools<PairedEndMappingWithBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void OutputTools<PairedEndMappingWithBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedEndMappingWithBarcode &mapping) {
  if (mapping_output_format_ == MAPPINGFORMAT_BED) {
    std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t mapping_end_position = mapping.GetEndPosition();
    this->AppendMappingOutput(
        std::string(reference_sequence_name) + "\t" +
        std::to_string(mapping.GetStartPosition()) + "\t" +
        std::to_string(mapping_end_position) + "\t" +
        barcode_translator.Translate(mapping.cell_barcode_, cell_barcode_length_) + "\t" +
        std::to_string(mapping.num_dups_) + "\n");
  } else {
    bool positive_strand = mapping.IsPositiveStrand();
    uint32_t positive_read_end =
        mapping.fragment_start_position_ + mapping.positive_alignment_length_;
    uint32_t negative_read_end =
        mapping.fragment_start_position_ + mapping.fragment_length_;
    uint32_t negative_read_start =
        negative_read_end - mapping.negative_alignment_length_;
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    if (positive_strand) {
      this->AppendMappingOutput(
          std::string(reference_sequence_name) + "\t" +
          std::to_string(mapping.fragment_start_position_) + "\t" +
          std::to_string(positive_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t+\n" +
          std::string(reference_sequence_name) + "\t" +
          std::to_string(negative_read_start) + "\t" +
          std::to_string(negative_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t-\n");
    } else {
      this->AppendMappingOutput(
          std::string(reference_sequence_name) + "\t" +
          std::to_string(negative_read_start) + "\t" +
          std::to_string(negative_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t-\n" +
          std::string(reference_sequence_name) + "\t" +
          std::to_string(mapping.fragment_start_position_) + "\t" +
          std::to_string(positive_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t+\n");
    }
  }
}

// Specialization for PAF format.
template <>
void OutputTools<PAFMapping>::OutputHeader(uint32_t num_reference_sequences,
                                           const SequenceBatch &reference) {}

template <>
void OutputTools<PAFMapping>::AppendMapping(uint32_t rid,
                                            const SequenceBatch &reference,
                                            const PAFMapping &mapping) {
  const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
  uint32_t reference_sequence_length = reference.GetSequenceLengthAt(rid);
  std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
  uint32_t mapping_end_position =
      mapping.fragment_start_position_ + mapping.fragment_length_;
  this->AppendMappingOutput(
      mapping.read_name_ + "\t" + std::to_string(mapping.read_length_) + "\t" +
      std::to_string(0) + "\t" + std::to_string(mapping.read_length_) + "\t" +
      strand + "\t" + std::string(reference_sequence_name) + "\t" +
      std::to_string(reference_sequence_length) + "\t" +
      std::to_string(mapping.fragment_start_position_) + "\t" +
      std::to_string(mapping_end_position) + "\t" +
      std::to_string(mapping.read_length_) + "\t" +
      std::to_string(mapping.fragment_length_) + "\t" +
      std::to_string(mapping.mapq_) + "\n");
}

template <>
void OutputTools<PAFMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PAFMapping> > &mappings) {
  FILE *temp_mapping_output_file =
      fopen(temp_mapping_output_file_path.c_str(), "wb");
  assert(temp_mapping_output_file != NULL);
  for (size_t ri = 0; ri < num_reference_sequences; ++ri) {
    // make sure mappings[ri] exists even if its size is 0
    size_t num_mappings = mappings[ri].size();
    fwrite(&num_mappings, sizeof(size_t), 1, temp_mapping_output_file);
    if (mappings[ri].size() > 0) {
      for (size_t mi = 0; mi < num_mappings; ++mi) {
        mappings[ri][mi].WriteToFile(temp_mapping_output_file);
      }
      // fwrite(mappings[ri].data(), sizeof(MappingRecord), mappings[ri].size(),
      // temp_mapping_output_file);
    }
  }
  fclose(temp_mapping_output_file);
}

// Specialization for PairedPAF format.
template <>
void OutputTools<PairedPAFMapping>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void OutputTools<PairedPAFMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PairedPAFMapping> > &mappings) {
  FILE *temp_mapping_output_file =
      fopen(temp_mapping_output_file_path.c_str(), "wb");
  assert(temp_mapping_output_file != NULL);
  for (size_t ri = 0; ri < num_reference_sequences; ++ri) {
    // make sure mappings[ri] exists even if its size is 0
    size_t num_mappings = mappings[ri].size();
    fwrite(&num_mappings, sizeof(size_t), 1, temp_mapping_output_file);
    if (mappings[ri].size() > 0) {
      for (size_t mi = 0; mi < num_mappings; ++mi) {
        mappings[ri][mi].WriteToFile(temp_mapping_output_file);
      }
      // fwrite(mappings[ri].data(), sizeof(MappingRecord), mappings[ri].size(),
      // temp_mapping_output_file);
    }
  }
  fclose(temp_mapping_output_file);
}

template <>
void OutputTools<PairedPAFMapping>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedPAFMapping &mapping) {
  bool positive_strand = mapping.IsPositiveStrand();
  uint32_t positive_read_end =
      mapping.fragment_start_position_ + mapping.positive_alignment_length_;
  uint32_t negative_read_end =
      mapping.fragment_start_position_ + mapping.fragment_length_;
  uint32_t negative_read_start =
      negative_read_end - mapping.negative_alignment_length_;
  const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
  uint32_t reference_sequence_length = reference.GetSequenceLengthAt(rid);
  if (positive_strand) {
    this->AppendMappingOutput(
        mapping.read1_name_ + "\t" + std::to_string(mapping.read1_length_) +
        "\t" + std::to_string(0) + "\t" +
        std::to_string(mapping.read1_length_) + "\t" + "+" + "\t" +
        std::string(reference_sequence_name) + "\t" +
        std::to_string(reference_sequence_length) + "\t" +
        std::to_string(mapping.fragment_start_position_) + "\t" +
        std::to_string(positive_read_end) + "\t" +
        std::to_string(mapping.read1_length_) + "\t" +
        std::to_string(mapping.positive_alignment_length_) + "\t" +
        std::to_string(mapping.mapq1_) + "\n");
    this->AppendMappingOutput(
        mapping.read2_name_ + "\t" + std::to_string(mapping.read2_length_) +
        "\t" + std::to_string(0) + "\t" +
        std::to_string(mapping.read2_length_) + "\t" + "-" + "\t" +
        std::string(reference_sequence_name) + "\t" +
        std::to_string(reference_sequence_length) + "\t" +
        std::to_string(negative_read_start) + "\t" +
        std::to_string(negative_read_end) + "\t" +
        std::to_string(mapping.read2_length_) + "\t" +
        std::to_string(mapping.negative_alignment_length_) + "\t" +
        std::to_string(mapping.mapq2_) + "\n");
  } else {
    this->AppendMappingOutput(
        mapping.read1_name_ + "\t" + std::to_string(mapping.read1_length_) +
        "\t" + std::to_string(0) + "\t" +
        std::to_string(mapping.read1_length_) + "\t" + "-" + "\t" +
        std::string(reference_sequence_name) + "\t" +
        std::to_string(reference_sequence_length) + "\t" +
        std::to_string(negative_read_start) + "\t" +
        std::to_string(negative_read_end) + "\t" +
        std::to_string(mapping.read1_length_) + "\t" +
        std::to_string(mapping.negative_alignment_length_) + "\t" +
        std::to_string(mapping.mapq1_) + "\n");
    this->AppendMappingOutput(
        mapping.read2_name_ + "\t" + std::to_string(mapping.read2_length_) +
        "\t" + std::to_string(0) + "\t" +
        std::to_string(mapping.read2_length_) + "\t" + "+" + "\t" +
        std::string(reference_sequence_name) + "\t" +
        std::to_string(reference_sequence_length) + "\t" +
        std::to_string(mapping.fragment_start_position_) + "\t" +
        std::to_string(positive_read_end) + "\t" +
        std::to_string(mapping.read2_length_) + "\t" +
        std::to_string(mapping.positive_alignment_length_) + "\t" +
        std::to_string(mapping.mapq2_) + "\n");
  }
}

// Specialization for SAM format.
template <>
void OutputTools<SAMMapping>::OutputHeader(uint32_t num_reference_sequences,
                                           const SequenceBatch &reference) {
  for (uint32_t rid = 0; rid < num_reference_sequences; ++rid) {
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t reference_sequence_length = reference.GetSequenceLengthAt(rid);
    this->AppendMappingOutput(
        "@SQ\tSN:" + std::string(reference_sequence_name) +
        "\tLN:" + std::to_string(reference_sequence_length) + "\n");
  }
}

template <>
void OutputTools<SAMMapping>::AppendMapping(uint32_t rid,
                                            const SequenceBatch &reference,
                                            const SAMMapping &mapping) {
  // const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
  // uint32_t reference_sequence_length = reference.GetSequenceLengthAt(rid);
  // std::string strand = (mapping.direction & 1) == 1 ? "+" : "-";
  // uint32_t mapping_end_position = mapping.fragment_start_position +
  // mapping.fragment_length;
  const char *reference_sequence_name =
      (mapping.flag_ & BAM_FUNMAP) > 0 ? "*" : reference.GetSequenceNameAt(rid);
  const uint32_t mapping_start_position = mapping.GetStartPosition();
  this->AppendMappingOutput(
      mapping.read_name_ + "\t" + std::to_string(mapping.flag_) + "\t" +
      std::string(reference_sequence_name) + "\t" +
      std::to_string(mapping_start_position) + "\t" +
      std::to_string(mapping.mapq_) + "\t" + mapping.GenerateCigarString() +
      "\t*\t" + std::to_string(0) + "\t" + std::to_string(0) + "\t" +
      mapping.sequence_ + "\t" + mapping.sequence_qual_ + "\t" +
      mapping.GenerateIntTagString("NM", mapping.NM_) +
      "\tMD:Z:" + mapping.MD_);
  if (cell_barcode_length_ > 0) {
    this->AppendMappingOutput(
        "\tCB:Z:" + barcode_translator.Translate(mapping.cell_barcode_, cell_barcode_length_));
  }
  this->AppendMappingOutput("\n");
}

template <>
void OutputTools<SAMMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<SAMMapping> > &mappings) {
  FILE *temp_mapping_output_file =
      fopen(temp_mapping_output_file_path.c_str(), "wb");
  assert(temp_mapping_output_file != NULL);
  for (size_t ri = 0; ri < num_reference_sequences; ++ri) {
    // make sure mappings[ri] exists even if its size is 0
    size_t num_mappings = mappings[ri].size();
    fwrite(&num_mappings, sizeof(size_t), 1, temp_mapping_output_file);
    if (mappings[ri].size() > 0) {
      for (size_t mi = 0; mi < num_mappings; ++mi) {
        mappings[ri][mi].WriteToFile(temp_mapping_output_file);
      }
      // fwrite(mappings[ri].data(), sizeof(MappingRecord), mappings[ri].size(),
      // temp_mapping_output_file);
    }
  }
  fclose(temp_mapping_output_file);
}

// Specialization for pairs format.
template <>
void OutputTools<PairsMapping>::OutputHeader(uint32_t num_reference_sequences,
                                             const SequenceBatch &reference) {
  std::vector<uint32_t> rid_order;
  rid_order.resize(num_reference_sequences);
  uint32_t i;
  for (i = 0; i < num_reference_sequences; ++i) {
    rid_order[this->custom_rid_rank_[i]] = i;
  }
  this->AppendMappingOutput("## pairs format v1.0.0\n#shape: upper triangle\n");
  for (i = 0; i < num_reference_sequences; ++i) {
    uint32_t rid = rid_order[i];
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t reference_sequence_length = reference.GetSequenceLengthAt(rid);
    this->AppendMappingOutput(
        "#chromsize: " + std::string(reference_sequence_name) + " " +
        std::to_string(reference_sequence_length) + "\n");
  }
  this->AppendMappingOutput(
      "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n");
}

template <>
void OutputTools<PairsMapping>::AppendMapping(uint32_t rid,
                                              const SequenceBatch &reference,
                                              const PairsMapping &mapping) {
  const char *reference_sequence_name1 =
      reference.GetSequenceNameAt(mapping.rid1_);
  const char *reference_sequence_name2 =
      reference.GetSequenceNameAt(mapping.rid2_);
  this->AppendMappingOutput(mapping.read_name_ + "\t" +
                            std::string(reference_sequence_name1) + "\t" +
                            std::to_string(mapping.GetPosition(1)) + "\t" +
                            std::string(reference_sequence_name2) + "\t" +
                            std::to_string(mapping.GetPosition(2)) + "\t" +
                            std::string(1, mapping.GetDirection(1)) + "\t" +
                            std::string(1, mapping.GetDirection(2)) + "\tUU\n");
}

template <>
void OutputTools<PairsMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PairsMapping> > &mappings) {
  FILE *temp_mapping_output_file =
      fopen(temp_mapping_output_file_path.c_str(), "wb");
  assert(temp_mapping_output_file != NULL);
  for (size_t ri = 0; ri < num_reference_sequences; ++ri) {
    // make sure mappings[ri] exists even if its size is 0
    size_t num_mappings = mappings[ri].size();
    fwrite(&num_mappings, sizeof(size_t), 1, temp_mapping_output_file);
    if (mappings[ri].size() > 0) {
      for (size_t mi = 0; mi < num_mappings; ++mi) {
        mappings[ri][mi].WriteToFile(temp_mapping_output_file);
      }
      // fwrite(mappings[ri].data(), sizeof(MappingRecord), mappings[ri].size(),
      // temp_mapping_output_file);
    }
  }
  fclose(temp_mapping_output_file);
}

}  // namespace chromap
