#ifndef PAFMAPPING_H_
#define PAFMAPPING_H_

#include <string>

#include "mapping.h"

namespace chromap {

// When direction = 1, strand is positive
class PAFMapping : public Mapping {
 public:
  uint32_t read_id_;
  std::string read_name_;
  uint16_t read_length_;
  uint32_t fragment_start_position_;
  uint16_t fragment_length_;
  uint8_t mapq_ : 6, direction_ : 1, is_unique_ : 1;
  uint8_t num_dups_;
  PAFMapping() : num_dups_(0) {}
  PAFMapping(uint32_t read_id, const std::string &read_name,
             uint16_t read_length, uint32_t fragment_start_position,
             uint16_t fragment_length, uint8_t mapq, uint8_t direction,
             uint8_t is_unique, uint8_t num_dups)
      : read_id_(read_id),
        read_name_(read_name),
        read_length_(read_length),
        fragment_start_position_(fragment_start_position),
        fragment_length_(fragment_length),
        mapq_(mapq),
        direction_(direction),
        is_unique_(is_unique),
        num_dups_(num_dups){};
  bool operator<(const PAFMapping &m) const {
    return std::tie(fragment_start_position_, fragment_length_, mapq_,
                    direction_, is_unique_, read_id_, read_length_) <
           std::tie(m.fragment_start_position_, m.fragment_length_, m.mapq_,
                    m.direction_, m.is_unique_, m.read_id_, m.read_length_);
  }
  bool operator==(const PAFMapping &m) const {
    return std::tie(fragment_start_position_) ==
           std::tie(m.fragment_start_position_);
  }
  bool IsSamePosition(const PAFMapping &m) const {
    return std::tie(fragment_start_position_) ==
           std::tie(m.fragment_start_position_);
  }
  uint64_t GetBarcode() const { return 0; }
  void Tn5Shift() {
    if (direction_ == 1) {
      fragment_start_position_ += 4;
    } else {
      fragment_length_ -= 5;
    }
  }
  bool IsPositiveStrand() const { return direction_ > 0 ? true : false; }
  uint32_t GetStartPosition() const {  // inclusive
    return fragment_start_position_;
  }
  uint32_t GetEndPosition() const {  // exclusive
    return fragment_start_position_ + fragment_length_;
  }
  uint16_t GetByteSize() const {
    return 2 * sizeof(uint32_t) + 2 * sizeof(uint16_t) + 2 * sizeof(uint8_t) +
           read_name_.length() * sizeof(char);
  }
  size_t WriteToFile(FILE *temp_mapping_output_file) const {
    size_t num_written_bytes = 0;
    num_written_bytes +=
        fwrite(&read_id_, sizeof(uint32_t), 1, temp_mapping_output_file);
    uint16_t read_name_length = read_name_.length();
    num_written_bytes += fwrite(&read_name_length, sizeof(uint16_t), 1,
                                temp_mapping_output_file);
    num_written_bytes += fwrite(read_name_.data(), sizeof(char),
                                read_name_length, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&read_length_, sizeof(uint16_t), 1, temp_mapping_output_file);
    num_written_bytes += fwrite(&fragment_start_position_, sizeof(uint32_t), 1,
                                temp_mapping_output_file);
    num_written_bytes += fwrite(&fragment_length_, sizeof(uint16_t), 1,
                                temp_mapping_output_file);
    uint8_t mapq_direction_is_unique =
        (mapq_ << 2) | (direction_ << 1) | is_unique_;
    num_written_bytes += fwrite(&mapq_direction_is_unique, sizeof(uint8_t), 1,
                                temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&num_dups_, sizeof(uint8_t), 1, temp_mapping_output_file);
    return num_written_bytes;
  }
  size_t LoadFromFile(FILE *temp_mapping_output_file) {
    size_t num_read_bytes = 0;
    num_read_bytes +=
        fread(&read_id_, sizeof(uint32_t), 1, temp_mapping_output_file);
    uint16_t read_name_length = 0;
    num_read_bytes +=
        fread(&read_name_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    read_name_ = std::string(read_name_length, '\0');
    num_read_bytes += fread(&(read_name_[0]), sizeof(char), read_name_length,
                            temp_mapping_output_file);
    num_read_bytes +=
        fread(&read_length_, sizeof(uint16_t), 1, temp_mapping_output_file);
    num_read_bytes += fread(&fragment_start_position_, sizeof(uint32_t), 1,
                            temp_mapping_output_file);
    num_read_bytes +=
        fread(&fragment_length_, sizeof(uint16_t), 1, temp_mapping_output_file);
    uint8_t mapq_direction_is_unique = 0;
    num_read_bytes += fread(&mapq_direction_is_unique, sizeof(uint8_t), 1,
                            temp_mapping_output_file);
    mapq_ = (mapq_direction_is_unique >> 2);
    direction_ = (mapq_direction_is_unique >> 1) & 1;
    is_unique_ = mapq_direction_is_unique & 1;
    num_read_bytes +=
        fread(&num_dups_, sizeof(uint8_t), 1, temp_mapping_output_file);
    return num_read_bytes;
  }
};

class PairedPAFMapping : public Mapping {
 public:
  uint32_t read_id_;
  std::string read1_name_;
  std::string read2_name_;
  uint16_t read1_length_;
  uint16_t read2_length_;
  uint32_t fragment_start_position_;
  uint16_t fragment_length_;
  uint16_t positive_alignment_length_;
  uint16_t negative_alignment_length_;
  uint8_t mapq_;
  uint16_t mapq1_ : 6, mapq2_ : 6, direction_ : 1, is_unique_ : 1,
      reserved_ : 2;
  uint8_t num_dups_;
  // uint8_t mapq; // least significant bit saves the direction of mapping
  PairedPAFMapping() : num_dups_(0) {}
  PairedPAFMapping(uint32_t read_id, std::string read1_name,
                   std::string read2_name, uint16_t read1_length,
                   uint16_t read2_length, uint32_t fragment_start_position,
                   uint16_t fragment_length, uint16_t positive_alignment_length,
                   uint16_t negative_alignment_length, uint8_t mapq,
                   uint16_t mapq1, uint16_t mapq2, uint16_t direction,
                   uint16_t is_unique, uint8_t num_dups)
      : read_id_(read_id),
        read1_name_(read1_name),
        read2_name_(read2_name),
        fragment_start_position_(fragment_start_position),
        fragment_length_(fragment_length),
        positive_alignment_length_(positive_alignment_length),
        negative_alignment_length_(negative_alignment_length),
        mapq_(mapq),
        mapq1_(mapq1),
        mapq2_(mapq2),
        direction_(direction),
        is_unique_(is_unique),
        num_dups_(num_dups) {}
  bool operator<(const PairedPAFMapping &m) const {
    return std::tie(fragment_start_position_, fragment_length_, mapq1_, mapq2_,
                    direction_, is_unique_, read_id_,
                    positive_alignment_length_, negative_alignment_length_) <
           std::tie(m.fragment_start_position_, m.fragment_length_, m.mapq1_,
                    m.mapq2_, m.direction_, m.is_unique_, m.read_id_,
                    m.positive_alignment_length_, m.negative_alignment_length_);
  }
  bool operator==(const PairedPAFMapping &m) const {
    return std::tie(fragment_start_position_, fragment_length_) ==
           std::tie(m.fragment_start_position_, m.fragment_length_);
  }
  bool IsSamePosition(const PairedPAFMapping &m) const {
    return std::tie(fragment_start_position_, fragment_length_) ==
           std::tie(m.fragment_start_position_, m.fragment_length_);
  }
  uint64_t GetBarcode() const { return 0; }
  void Tn5Shift() {
    fragment_start_position_ += 4;
    positive_alignment_length_ -= 4;
    fragment_length_ -= 9;
    negative_alignment_length_ -= 5;
  }
  bool IsPositiveStrand() const { return direction_ > 0 ? true : false; }
  uint32_t GetStartPosition() const {  // inclusive
    return fragment_start_position_;
  }
  uint32_t GetEndPosition() const {  // exclusive
    return fragment_start_position_ + fragment_length_;
  }
  uint16_t GetByteSize() const {
    return 2 * sizeof(uint32_t) + 6 * sizeof(uint16_t) + 2 * sizeof(uint8_t) +
           (read1_name_.length() + read2_name_.length()) * sizeof(char);
  }
  size_t WriteToFile(FILE *temp_mapping_output_file) const {
    size_t num_written_bytes = 0;
    num_written_bytes +=
        fwrite(&read_id_, sizeof(uint32_t), 1, temp_mapping_output_file);
    uint16_t read1_name_length = read1_name_.length();
    num_written_bytes += fwrite(&read1_name_length, sizeof(uint16_t), 1,
                                temp_mapping_output_file);
    num_written_bytes += fwrite(read1_name_.data(), sizeof(char),
                                read1_name_length, temp_mapping_output_file);
    uint16_t read2_name_length = read2_name_.length();
    num_written_bytes += fwrite(&read2_name_length, sizeof(uint16_t), 1,
                                temp_mapping_output_file);
    num_written_bytes += fwrite(read2_name_.data(), sizeof(char),
                                read2_name_length, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&read1_length_, sizeof(uint16_t), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&read2_length_, sizeof(uint16_t), 1, temp_mapping_output_file);
    num_written_bytes += fwrite(&fragment_start_position_, sizeof(uint32_t), 1,
                                temp_mapping_output_file);
    num_written_bytes += fwrite(&fragment_length_, sizeof(uint16_t), 1,
                                temp_mapping_output_file);
    num_written_bytes += fwrite(&positive_alignment_length_, sizeof(uint16_t),
                                1, temp_mapping_output_file);
    num_written_bytes += fwrite(&negative_alignment_length_, sizeof(uint16_t),
                                1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&mapq_, sizeof(uint8_t), 1, temp_mapping_output_file);
    uint16_t mapq1_mapq2_direction_is_unique =
        (mapq1_ << 10) | (mapq2_ << 4) | (direction_ << 3) | (is_unique_ << 2);
    num_written_bytes += fwrite(&mapq1_mapq2_direction_is_unique,
                                sizeof(uint16_t), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&num_dups_, sizeof(uint8_t), 1, temp_mapping_output_file);
    return num_written_bytes;
  }
  size_t LoadFromFile(FILE *temp_mapping_output_file) {
    size_t num_read_bytes = 0;
    num_read_bytes +=
        fread(&read_id_, sizeof(uint32_t), 1, temp_mapping_output_file);
    uint16_t read1_name_length = 0;
    num_read_bytes += fread(&read1_name_length, sizeof(uint16_t), 1,
                            temp_mapping_output_file);
    read1_name_ = std::string(read1_name_length, '\0');
    num_read_bytes += fread(&(read1_name_[0]), sizeof(char), read1_name_length,
                            temp_mapping_output_file);
    uint16_t read2_name_length = 0;
    num_read_bytes += fread(&read2_name_length, sizeof(uint16_t), 1,
                            temp_mapping_output_file);
    read2_name_ = std::string(read2_name_length, '\0');
    num_read_bytes += fread(&(read2_name_[0]), sizeof(char), read2_name_length,
                            temp_mapping_output_file);
    num_read_bytes +=
        fread(&read1_length_, sizeof(uint16_t), 1, temp_mapping_output_file);
    num_read_bytes +=
        fread(&read2_length_, sizeof(uint16_t), 1, temp_mapping_output_file);
    num_read_bytes += fread(&fragment_start_position_, sizeof(uint32_t), 1,
                            temp_mapping_output_file);
    num_read_bytes +=
        fread(&fragment_length_, sizeof(uint16_t), 1, temp_mapping_output_file);
    num_read_bytes += fread(&positive_alignment_length_, sizeof(uint16_t), 1,
                            temp_mapping_output_file);
    num_read_bytes += fread(&negative_alignment_length_, sizeof(uint16_t), 1,
                            temp_mapping_output_file);
    num_read_bytes +=
        fread(&mapq_, sizeof(uint8_t), 1, temp_mapping_output_file);
    uint16_t mapq1_mapq2_direction_is_unique = 0;
    num_read_bytes += fread(&mapq1_mapq2_direction_is_unique, sizeof(uint16_t),
                            1, temp_mapping_output_file);
    mapq1_ = (mapq1_mapq2_direction_is_unique >> 10);
    mapq2_ = ((mapq1_mapq2_direction_is_unique << 6) >> 10);
    direction_ = (mapq1_mapq2_direction_is_unique >> 3) & 1;
    is_unique_ = (mapq1_mapq2_direction_is_unique >> 2) & 1;
    num_read_bytes +=
        fread(&num_dups_, sizeof(uint8_t), 1, temp_mapping_output_file);
    return num_read_bytes;
  }
};

}  // namespace chromap

#endif  // PAFMAPPING_H_
