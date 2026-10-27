#ifndef PAIRSMAPPING_H_
#define PAIRSMAPPING_H_

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

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
  int strand1_;  // 1-positive. 0-negative
  int strand2_;
  uint16_t mapq_ : 8, is_unique_ : 1, num_dups_ : 7;

  PairsMapping() : num_dups_(0) {}
  PairsMapping(uint32_t read_id, std::string read_name, uint64_t cell_barcode,
               int rid1, int rid2, uint32_t pos1, uint32_t pos2, int strand1,
               int strand2, uint8_t mapq, uint8_t is_unique, uint8_t num_dups)
      : read_id_(read_id),
        read_name_(read_name),
        cell_barcode_(cell_barcode),
        rid1_(rid1),
        rid2_(rid2),
        pos1_(pos1),
        pos2_(pos2),
        strand1_(strand1),
        strand2_(strand2),
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

  char GetStrand(int idx) const {
    int d = strand1_;
    if (idx == 2) {
      d = strand2_;
    }
    return d > 0 ? '+' : '-';
  }

  bool IsPositiveStrand() const { return strand1_ > 0 ? true : false; }
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
  // A record's bytes in a temp file, in this order: read_id_ (uint32_t), the
  // read name's length (uint16_t) and characters, cell_barcode_ (uint64_t),
  // rid1_, rid2_ (int), pos1_, pos2_ (uint32_t), strand1_, strand2_ (int), and
  // mapq_, is_unique_ and num_dups_ packed in a uint16_t. The head is the
  // fields before the name, the tail those after it.
  static constexpr size_t kSerializedHeadSize =
      sizeof(uint32_t) + sizeof(uint16_t);
  static constexpr size_t kSerializedTailSize =
      sizeof(uint64_t) + 4 * sizeof(int) + 2 * sizeof(uint32_t) +
      sizeof(uint16_t);

  // The length of the record whose head starts at head.
  static size_t GetSerializedSize(const char *head) {
    uint16_t read_name_length = 0;
    memcpy(&read_name_length, head + sizeof(uint32_t), sizeof(uint16_t));
    return kSerializedHeadSize + read_name_length + kSerializedTailSize;
  }

  size_t GetSerializedSize() const {
    const uint16_t read_name_length = read_name_.length();
    return kSerializedHeadSize + read_name_length + kSerializedTailSize;
  }

  // Writes the record's GetSerializedSize() bytes at out and returns the end.
  char *SerializeTo(char *out) const {
    const uint16_t read_name_length = read_name_.length();
    memcpy(out, &read_id_, sizeof(uint32_t));
    out += sizeof(uint32_t);
    memcpy(out, &read_name_length, sizeof(uint16_t));
    out += sizeof(uint16_t);
    memcpy(out, read_name_.data(), read_name_length);
    out += read_name_length;
    memcpy(out, &cell_barcode_, sizeof(uint64_t));
    out += sizeof(uint64_t);
    memcpy(out, &rid1_, sizeof(int));
    out += sizeof(int);
    memcpy(out, &rid2_, sizeof(int));
    out += sizeof(int);
    memcpy(out, &pos1_, sizeof(uint32_t));
    out += sizeof(uint32_t);
    memcpy(out, &pos2_, sizeof(uint32_t));
    out += sizeof(uint32_t);
    memcpy(out, &strand1_, sizeof(int));
    out += sizeof(int);
    memcpy(out, &strand2_, sizeof(int));
    out += sizeof(int);
    const uint16_t mapq_unique_dups =
        (mapq_ << 8) | (is_unique_ << 7) | num_dups_;
    memcpy(out, &mapq_unique_dups, sizeof(uint16_t));
    out += sizeof(uint16_t);
    return out;
  }

  // Reads a record written by SerializeTo at in and returns its end.
  const char *DeserializeFrom(const char *in) {
    memcpy(&read_id_, in, sizeof(uint32_t));
    in += sizeof(uint32_t);
    uint16_t read_name_length = 0;
    memcpy(&read_name_length, in, sizeof(uint16_t));
    in += sizeof(uint16_t);
    read_name_.assign(in, read_name_length);
    in += read_name_length;
    memcpy(&cell_barcode_, in, sizeof(uint64_t));
    in += sizeof(uint64_t);
    memcpy(&rid1_, in, sizeof(int));
    in += sizeof(int);
    memcpy(&rid2_, in, sizeof(int));
    in += sizeof(int);
    memcpy(&pos1_, in, sizeof(uint32_t));
    in += sizeof(uint32_t);
    memcpy(&pos2_, in, sizeof(uint32_t));
    in += sizeof(uint32_t);
    memcpy(&strand1_, in, sizeof(int));
    in += sizeof(int);
    memcpy(&strand2_, in, sizeof(int));
    in += sizeof(int);
    uint16_t mapq_unique_dups = 0;
    memcpy(&mapq_unique_dups, in, sizeof(uint16_t));
    in += sizeof(uint16_t);
    mapq_ = (mapq_unique_dups >> 8);
    is_unique_ = (mapq_unique_dups >> 7) & 1;
    num_dups_ = ((mapq_unique_dups << 9) >> 9);
    return in;
  }

  // The count a field-by-field fwrite or fread of the record would return
  // when its first num_bytes bytes were transferred: one per whole fixed-size
  // field and one per name character.
  static size_t CountSerializedItems(size_t num_bytes,
                                     uint16_t read_name_length) {
    const size_t field_sizes[] = {sizeof(uint32_t), sizeof(uint16_t), 0,
                                  sizeof(uint64_t), sizeof(int),
                                  sizeof(int),      sizeof(uint32_t),
                                  sizeof(uint32_t), sizeof(int),
                                  sizeof(int),      sizeof(uint16_t)};
    size_t num_items = 0;
    size_t field_end = 0;
    for (size_t fi = 0; fi < sizeof(field_sizes) / sizeof(size_t); ++fi) {
      if (fi == 2) {
        const size_t num_name_bytes =
            num_bytes - field_end < read_name_length ? num_bytes - field_end
                                                     : read_name_length;
        num_items += num_name_bytes;
        field_end += read_name_length;
        if (field_end > num_bytes) {
          break;
        }
        continue;
      }
      field_end += field_sizes[fi];
      if (field_end > num_bytes) {
        break;
      }
      ++num_items;
    }
    return num_items;
  }

  // One fwrite of the record's bytes. Returns the count the field-by-field
  // fwrite calls of the same bytes returned.
  size_t WriteToFile(FILE *temp_mapping_output_file) const {
    const size_t kStackBufferSize = 512;
    char stack_buffer[kStackBufferSize];
    std::vector<char> heap_buffer;
    const size_t num_bytes = GetSerializedSize();
    char *buffer = stack_buffer;
    if (num_bytes > kStackBufferSize) {
      heap_buffer.resize(num_bytes);
      buffer = heap_buffer.data();
    }
    SerializeTo(buffer);
    const size_t num_written_bytes =
        fwrite(buffer, sizeof(char), num_bytes, temp_mapping_output_file);
    return CountSerializedItems(num_written_bytes, read_name_.length());
  }

  // Reads exactly one record with two freads, the head and then the rest.
  // Bytes missing at the end of the file read as zero. Returns the count the
  // field-by-field fread calls of the same bytes returned.
  size_t LoadFromFile(FILE *temp_mapping_output_file) {
    const size_t kStackBufferSize = 512;
    char stack_buffer[kStackBufferSize];
    std::vector<char> heap_buffer;
    char *buffer = stack_buffer;
    memset(buffer, 0, kSerializedHeadSize);
    size_t num_read_bytes = fread(buffer, sizeof(char), kSerializedHeadSize,
                                  temp_mapping_output_file);
    const size_t num_bytes = GetSerializedSize(buffer);
    if (num_bytes > kStackBufferSize) {
      heap_buffer.resize(num_bytes);
      memcpy(heap_buffer.data(), buffer, kSerializedHeadSize);
      buffer = heap_buffer.data();
    }
    memset(buffer + kSerializedHeadSize, 0, num_bytes - kSerializedHeadSize);
    if (num_read_bytes == kSerializedHeadSize) {
      num_read_bytes += fread(buffer + kSerializedHeadSize, sizeof(char),
                              num_bytes - kSerializedHeadSize,
                              temp_mapping_output_file);
    }
    DeserializeFrom(buffer);
    return CountSerializedItems(
        num_read_bytes,
        (uint16_t)(num_bytes - kSerializedHeadSize - kSerializedTailSize));
  }
};

}  // namespace chromap

#endif  // PAIRSMAPPING_H_
