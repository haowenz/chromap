#include "mapping_writer.h"

namespace chromap {

// Specialization for BED format.
template <>
void MappingWriter<MappingWithBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void MappingWriter<MappingWithBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const MappingWithBarcode &mapping) {
  AppendCStringToLine(reference.GetSequenceNameAt(rid));
  AppendCharToLine('\t');
  AppendUnsignedToLine(mapping.GetStartPosition());
  AppendCharToLine('\t');
  AppendUnsignedToLine(mapping.GetEndPosition());
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BED) {
    AppendCharToLine('\t');
    AppendBarcodeToLine(mapping.cell_barcode_);
    AppendCharToLine('\t');
    AppendUnsignedToLine(mapping.num_dups_);
  } else {
    AppendTextToLine("\tN\t", 3);
    AppendUnsignedToLine(mapping.mapq_);
    AppendCharToLine('\t');
    AppendCharToLine(mapping.IsPositiveStrand() ? '+' : '-');
  }
  AppendCharToLine('\n');
  WriteLine();
}

template <>
void MappingWriter<MappingWithoutBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void MappingWriter<MappingWithoutBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const MappingWithoutBarcode &mapping) {
  // BED and TagAlign lines are the same.
  AppendCStringToLine(reference.GetSequenceNameAt(rid));
  AppendCharToLine('\t');
  AppendUnsignedToLine(mapping.GetStartPosition());
  AppendCharToLine('\t');
  AppendUnsignedToLine(mapping.GetEndPosition());
  AppendTextToLine("\tN\t", 3);
  AppendUnsignedToLine(mapping.mapq_);
  AppendCharToLine('\t');
  AppendCharToLine(mapping.IsPositiveStrand() ? '+' : '-');
  AppendCharToLine('\t');
  AppendUnsignedToLine(mapping.num_dups_);
  AppendCharToLine('\n');
  WriteLine();
}

// Specialization for BEDPE format.
template <>
void MappingWriter<PairedEndMappingWithoutBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void MappingWriter<PairedEndMappingWithoutBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedEndMappingWithoutBarcode &mapping) {
  const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BED) {
    AppendCStringToLine(reference_sequence_name);
    AppendCharToLine('\t');
    AppendUnsignedToLine(mapping.GetStartPosition());
    AppendCharToLine('\t');
    AppendUnsignedToLine(mapping.GetEndPosition());
    AppendTextToLine("\tN\t", 3);
    AppendUnsignedToLine(mapping.mapq_);
    AppendCharToLine('\t');
    AppendCharToLine(mapping.IsPositiveStrand() ? '+' : '-');
    AppendCharToLine('\t');
    AppendUnsignedToLine(mapping.num_dups_);
    AppendCharToLine('\n');
  } else {
    bool positive_strand = mapping.IsPositiveStrand();
    uint32_t positive_read_end =
        mapping.fragment_start_position_ + mapping.positive_alignment_length_;
    uint32_t negative_read_end =
        mapping.fragment_start_position_ + mapping.fragment_length_;
    uint32_t negative_read_start =
        negative_read_end - mapping.negative_alignment_length_;
    for (int mate = 0; mate < 2; ++mate) {
      // The positive-strand read first when the fragment is on the positive
      // strand.
      const bool positive_read = (mate == 0) == positive_strand;
      AppendCStringToLine(reference_sequence_name);
      AppendCharToLine('\t');
      AppendUnsignedToLine(positive_read ? mapping.fragment_start_position_
                                         : negative_read_start);
      AppendCharToLine('\t');
      AppendUnsignedToLine(positive_read ? positive_read_end
                                         : negative_read_end);
      AppendTextToLine("\tN\t", 3);
      AppendUnsignedToLine(mapping.mapq_);
      AppendCharToLine('\t');
      AppendCharToLine(positive_read ? '+' : '-');
      if (mate == 1) {
        AppendCharToLine('\t');
        AppendUnsignedToLine(mapping.num_dups_);
      }
      AppendCharToLine('\n');
    }
  }
  WriteLine();
}

template <>
void MappingWriter<PairedEndMappingWithBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void MappingWriter<PairedEndMappingWithBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedEndMappingWithBarcode &mapping) {
  const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BED) {
    AppendCStringToLine(reference_sequence_name);
    AppendCharToLine('\t');
    AppendUnsignedToLine(mapping.GetStartPosition());
    AppendCharToLine('\t');
    AppendUnsignedToLine(mapping.GetEndPosition());
    AppendCharToLine('\t');
    AppendBarcodeToLine(mapping.cell_barcode_);
    AppendCharToLine('\t');
    AppendUnsignedToLine(mapping.num_dups_);
    AppendCharToLine('\n');
  } else {
    bool positive_strand = mapping.IsPositiveStrand();
    uint32_t positive_read_end =
        mapping.fragment_start_position_ + mapping.positive_alignment_length_;
    uint32_t negative_read_end =
        mapping.fragment_start_position_ + mapping.fragment_length_;
    uint32_t negative_read_start =
        negative_read_end - mapping.negative_alignment_length_;
    for (int mate = 0; mate < 2; ++mate) {
      // The positive-strand read first when the fragment is on the positive
      // strand.
      const bool positive_read = (mate == 0) == positive_strand;
      AppendCStringToLine(reference_sequence_name);
      AppendCharToLine('\t');
      AppendUnsignedToLine(positive_read ? mapping.fragment_start_position_
                                         : negative_read_start);
      AppendCharToLine('\t');
      AppendUnsignedToLine(positive_read ? positive_read_end
                                         : negative_read_end);
      AppendTextToLine("\tN\t", 3);
      AppendUnsignedToLine(mapping.mapq_);
      AppendCharToLine('\t');
      AppendCharToLine(positive_read ? '+' : '-');
      AppendCharToLine('\n');
    }
  }
  WriteLine();
}

// Specialization for PAF format.
template <>
void MappingWriter<PAFMapping>::OutputHeader(uint32_t num_reference_sequences,
                                             const SequenceBatch &reference) {}

template <>
void MappingWriter<PAFMapping>::AppendMapping(uint32_t rid,
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
void MappingWriter<PAFMapping>::OutputTempMapping(
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
void MappingWriter<PairedPAFMapping>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void MappingWriter<PairedPAFMapping>::OutputTempMapping(
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
void MappingWriter<PairedPAFMapping>::AppendMapping(
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
void MappingWriter<SAMMapping>::OutputHeader(uint32_t num_reference_sequences,
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
void MappingWriter<SAMMapping>::AppendMapping(uint32_t rid,
                                              const SequenceBatch &reference,
                                              const SAMMapping &mapping) {
  // const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
  // uint32_t reference_sequence_length = reference.GetSequenceLengthAt(rid);
  // std::string strand = (mapping.direction & 1) == 1 ? "+" : "-";
  // uint32_t mapping_end_position = mapping.fragment_start_position +
  // mapping.fragment_length;
  const char *reference_sequence_name =
      (mapping.flag_ & BAM_FUNMAP) > 0 ? "*" : reference.GetSequenceNameAt(rid);
  const char *mate_ref_sequence_name =
      mapping.mrid_ < 0 ? "*" : 
      ((uint32_t)mapping.mrid_ == rid ? "=" : reference.GetSequenceNameAt(mapping.mrid_));
  const uint32_t mapping_start_position = mapping.GetStartPosition();
  const uint32_t mate_mapping_start_position = mapping.mrid_ < 0 ? 0 : (mapping.mpos_ + 1);
  this->AppendMappingOutput(
      mapping.read_name_ + "\t" + std::to_string(mapping.flag_) + "\t" +
      std::string(reference_sequence_name) + "\t" +
      std::to_string(mapping_start_position) + "\t" +
      std::to_string(mapping.mapq_) + "\t" + mapping.GenerateCigarString() +
      "\t" + std::string(mate_ref_sequence_name) + "\t" + 
      std::to_string(mate_mapping_start_position) + "\t" + 
      std::to_string(mapping.tlen_) + "\t" +
      mapping.sequence_ + "\t" + mapping.sequence_qual_ + "\t" +
      mapping.GenerateIntTagString("NM", mapping.NM_) +
      "\tMD:Z:" + mapping.MD_);
  if (cell_barcode_length_ > 0) {
    this->AppendMappingOutput("\tCB:Z:" +
                              barcode_translator_.Translate(
                                  mapping.cell_barcode_, cell_barcode_length_));
  }
  this->AppendMappingOutput("\n");
}

template <>
void MappingWriter<SAMMapping>::OutputTempMapping(
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
void MappingWriter<PairsMapping>::OutputHeader(uint32_t num_reference_sequences,
                                               const SequenceBatch &reference) {
  std::vector<uint32_t> rid_order;
  rid_order.resize(num_reference_sequences);
  uint32_t i;
  for (i = 0; i < num_reference_sequences; ++i) {
    rid_order[pairs_custom_rid_rank_[i]] = i;
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
      "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2\n");
}

template <>
void MappingWriter<PairsMapping>::AppendMapping(uint32_t rid,
                                                const SequenceBatch &reference,
                                                const PairsMapping &mapping) {
  AppendStringToLine(mapping.read_name_);
  AppendCharToLine('\t');
  AppendCStringToLine(reference.GetSequenceNameAt(mapping.rid1_));
  AppendCharToLine('\t');
  AppendSignedToLine(mapping.GetPosition(1));
  AppendCharToLine('\t');
  AppendCStringToLine(reference.GetSequenceNameAt(mapping.rid2_));
  AppendCharToLine('\t');
  AppendSignedToLine(mapping.GetPosition(2));
  AppendCharToLine('\t');
  AppendCharToLine(mapping.GetStrand(1));
  AppendCharToLine('\t');
  AppendCharToLine(mapping.GetStrand(2));
  AppendTextToLine("\tUU\t", 4);
  AppendUnsignedToLine(mapping.mapq_);  // mapq1
  AppendCharToLine('\t');
  AppendUnsignedToLine(mapping.mapq_);  // mapq2
  AppendCharToLine('\n');
  WriteLine();
}

template <>
void MappingWriter<PairsMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PairsMapping> > &mappings) {
  FILE *temp_mapping_output_file =
      fopen(temp_mapping_output_file_path.c_str(), "wb");
  assert(temp_mapping_output_file != NULL);
  // The bytes of a size_t count per ref seq and a WriteToFile call per
  // mapping, collected in a buffer written with one fwrite whenever the next
  // item would not fit. A record is at most 65,575 bytes (its name length is a
  // uint16_t), so every item fits in an empty buffer.
  std::vector<char> buffer((size_t)1 << 20);
  size_t num_buffered_bytes = 0;
  for (size_t ri = 0; ri < num_reference_sequences; ++ri) {
    // make sure mappings[ri] exists even if its size is 0
    size_t num_mappings = mappings[ri].size();
    if (num_buffered_bytes + sizeof(size_t) > buffer.size()) {
      fwrite(buffer.data(), sizeof(char), num_buffered_bytes,
             temp_mapping_output_file);
      num_buffered_bytes = 0;
    }
    memcpy(buffer.data() + num_buffered_bytes, &num_mappings, sizeof(size_t));
    num_buffered_bytes += sizeof(size_t);
    for (size_t mi = 0; mi < num_mappings; ++mi) {
      const size_t num_record_bytes = mappings[ri][mi].GetSerializedSize();
      if (num_buffered_bytes + num_record_bytes > buffer.size()) {
        fwrite(buffer.data(), sizeof(char), num_buffered_bytes,
               temp_mapping_output_file);
        num_buffered_bytes = 0;
      }
      mappings[ri][mi].SerializeTo(buffer.data() + num_buffered_bytes);
      num_buffered_bytes += num_record_bytes;
    }
  }
  fwrite(buffer.data(), sizeof(char), num_buffered_bytes,
         temp_mapping_output_file);
  fclose(temp_mapping_output_file);
}

}  // namespace chromap
