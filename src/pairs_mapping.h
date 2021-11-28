#ifndef PAIRSMAPPING_H_
#define PAIRSMAPPING_H_

#include <string>

#include "mapping.h"

namespace chromap {

// Format for pairtools for HiC data.
class PairsMapping : public Mapping {
 public:
  uint32_t read_id_;
  std::string read_name_;
  uint64_t cell_barcode_;
  int rid1_;
  int rid2_;
  uint32_t pos1_;
  uint32_t pos2_;
  int direction1_;  // 1-positive. 0-negative
  int direction2_;
  uint16_t mapq_ : 8, is_unique_ : 1, num_dups_ : 7;

  PairsMapping() : num_dups_(0) {}
  PairsMapping(uint32_t read_id, std::string read_name, uint64_t cell_barcode,
               int rid1, int rid2, uint32_t pos1, uint32_t pos2, int direction1,
               int direction2, uint8_t mapq, uint8_t is_unique,
               uint8_t num_dups)
      : read_id_(read_id),
        read_name_(read_name),
        cell_barcode_(cell_barcode),
        rid1_(rid1),
        rid2_(rid2),
        pos1_(pos1),
        pos2_(pos2),
        direction1_(direction1),
        direction2_(direction2),
        mapq_(mapq),
        is_unique_(is_unique),
        num_dups_(num_dups) {}
  bool operator<(const PairsMapping &m) const {
    return std::tie(rid1_, rid2_, pos1_, pos2_, mapq_, read_id_) <
           std::tie(m.rid1_, m.rid2_, m.pos1_, m.pos2_, m.mapq_, m.read_id_);
  }
  bool operator==(const PairsMapping &m) const {
    return std::tie(rid1_, pos1_, rid2_, pos2_) ==
           std::tie(m.rid1_, m.pos1_, m.rid2_, m.pos2_);
    // return std::tie(pos1, pos2, rid1, rid2, is_rev1, is_rev2) ==
    // std::tie(m.pos1, m.pos2, m.rid1, m.rid2, m.is_rev1, m.is_rev2);
  }
  bool IsSamePosition(const PairsMapping &m) const {
    return std::tie(rid1_, pos1_, rid2_, pos2_) ==
           std::tie(m.rid1_, m.pos1_, m.rid2_, m.pos2_);
  }
  uint64_t GetBarcode() const { return 0; }
  void Tn5Shift() {
    // We don't support Tn5 shift in SAM format because it has other fields that
    // depend mapping position.
  }

  int GetPosition(int idx) const {
    if (idx == 2) {
      return pos2_ + 1;
    }
    return pos1_ + 1;
  }

  char GetDirection(int idx) const {
    int d = direction1_;
    if (idx == 2) {
      d = direction2_;
    }
    return d > 0 ? '+' : '-';
  }

  bool IsPositiveStrand() const { return direction1_ > 0 ? true : false; }
  uint32_t GetStartPosition() const {  // inclusive
    return pos1_;
  }
  uint32_t GetEndPosition() const {  // exclusive
    return pos2_;
  }
  uint16_t GetByteSize() const {
    return 5 * sizeof(uint32_t) + 1 * sizeof(uint16_t) + 4 * sizeof(int) +
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
        fwrite(&cell_barcode_, sizeof(uint64_t), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&rid1_, sizeof(int), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&rid2_, sizeof(int), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&pos1_, sizeof(uint32_t), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&pos2_, sizeof(uint32_t), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&direction1_, sizeof(int), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&direction2_, sizeof(int), 1, temp_mapping_output_file);
    uint16_t mapq_unique_dups = (mapq_ << 8) | (is_unique_ << 7) | num_dups_;
    num_written_bytes += fwrite(&mapq_unique_dups, sizeof(uint16_t), 1,
                                temp_mapping_output_file);
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
        fread(&cell_barcode_, sizeof(uint64_t), 1, temp_mapping_output_file);
    num_read_bytes += fread(&rid1_, sizeof(int), 1, temp_mapping_output_file);
    num_read_bytes += fread(&rid2_, sizeof(int), 1, temp_mapping_output_file);
    num_read_bytes +=
        fread(&pos1_, sizeof(uint32_t), 1, temp_mapping_output_file);
    num_read_bytes +=
        fread(&pos2_, sizeof(uint32_t), 1, temp_mapping_output_file);
    num_read_bytes +=
        fread(&direction1_, sizeof(int), 1, temp_mapping_output_file);
    num_read_bytes +=
        fread(&direction2_, sizeof(int), 1, temp_mapping_output_file);
    uint16_t mapq_unique_dups = 0;
    num_read_bytes +=
        fread(&mapq_unique_dups, sizeof(uint16_t), 1, temp_mapping_output_file);
    mapq_ = (mapq_unique_dups >> 8);
    is_unique_ = (mapq_unique_dups >> 7) & 1;
    num_dups_ = ((mapq_unique_dups << 9) >> 9);
    return num_read_bytes;
  }
};

}  // namespace chromap

#endif  // PAIRSMAPPING_H_
