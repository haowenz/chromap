#ifndef BEDMAPPING_H_
#define BEDMAPPING_H_

#include <string>

#include "mapping.h"

namespace chromap {

class MappingWithBarcode : public Mapping {
 public:
  uint32_t read_id_;
  uint64_t cell_barcode_;
  uint32_t fragment_start_position_;
  uint16_t fragment_length_;
  uint8_t mapq_ : 6, direction_ : 1, is_unique_ : 1;
  uint8_t num_dups_;
  // uint8_t mapq;
  MappingWithBarcode() : num_dups_(0) {}
  MappingWithBarcode(uint32_t read_id, uint64_t cell_barcode,
                     uint32_t fragment_start_position, uint16_t fragment_length,
                     uint8_t mapq, uint8_t direction, uint8_t is_unique,
                     uint8_t num_dups)
      : read_id_(read_id),
        cell_barcode_(cell_barcode),
        fragment_start_position_(fragment_start_position),
        fragment_length_(fragment_length),
        mapq_(mapq),
        direction_(direction),
        is_unique_(is_unique),
        num_dups_(num_dups) {}
  bool operator<(const MappingWithBarcode &m) const {
    return std::tie(fragment_start_position_, fragment_length_, cell_barcode_,
                    mapq_, direction_, is_unique_, read_id_) <
           std::tie(m.fragment_start_position_, m.fragment_length_,
                    m.cell_barcode_, m.mapq_, m.direction_, m.is_unique_,
                    m.read_id_);
  }
  bool operator==(const MappingWithBarcode &m) const {
    return std::tie(cell_barcode_, fragment_start_position_) ==
           std::tie(m.cell_barcode_, m.fragment_start_position_);
  }
  bool IsSamePosition(const MappingWithBarcode &m) const {
    return std::tie(fragment_start_position_) ==
           std::tie(m.fragment_start_position_);
  }
  uint64_t GetBarcode() const { return cell_barcode_; }
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
};

class MappingWithoutBarcode : public Mapping {
 public:
  uint32_t read_id_;
  uint32_t fragment_start_position_;
  uint16_t fragment_length_;
  // uint8_t mapq;
  uint8_t mapq_ : 6, direction_ : 1, is_unique_ : 1;
  uint16_t num_dups_; // Need higher limit in bulk setting

  MappingWithoutBarcode() : num_dups_(0) {}
  MappingWithoutBarcode(uint32_t read_id, uint32_t fragment_start_position,
                        uint16_t fragment_length, uint16_t mapq,
                        uint8_t direction, uint8_t is_unique, uint8_t num_dups)
      : read_id_(read_id),
        fragment_start_position_(fragment_start_position),
        fragment_length_(fragment_length),
        mapq_(mapq),
        direction_(direction),
        is_unique_(is_unique),
        num_dups_(num_dups) {}

  bool operator<(const MappingWithoutBarcode &m) const {
    return std::tie(fragment_start_position_, fragment_length_, mapq_,
                    direction_, is_unique_, read_id_) <
           std::tie(m.fragment_start_position_, m.fragment_length_, m.mapq_,
                    m.direction_, m.is_unique_, m.read_id_);
  }
  bool operator==(const MappingWithoutBarcode &m) const {
    return std::tie(fragment_start_position_) ==
           std::tie(m.fragment_start_position_);
  }
  bool IsSamePosition(const MappingWithoutBarcode &m) const {
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
};

class PairedEndMappingWithBarcode : public Mapping {
 public:
  uint32_t read_id_;
  uint64_t cell_barcode_;
  uint32_t fragment_start_position_;
  uint16_t fragment_length_;
  uint8_t mapq_ : 6, direction_ : 1, is_unique_ : 1;
  uint8_t num_dups_;
  // uint8_t mapq;
  uint16_t positive_alignment_length_;
  uint16_t negative_alignment_length_;
  PairedEndMappingWithBarcode() : num_dups_(0) {}
  PairedEndMappingWithBarcode(uint32_t read_id, uint64_t cell_barcode,
                              uint32_t fragment_start_position,
                              uint16_t fragment_length, uint8_t mapq,
                              uint8_t direction, uint8_t is_unique,
                              uint8_t num_dups,
                              uint16_t positive_alignment_length,
                              uint16_t negative_alignment_length)
      : read_id_(read_id),
        cell_barcode_(cell_barcode),
        fragment_start_position_(fragment_start_position),
        fragment_length_(fragment_length),
        mapq_(mapq),
        direction_(direction),
        is_unique_(is_unique),
        num_dups_(num_dups),
        positive_alignment_length_(positive_alignment_length),
        negative_alignment_length_(negative_alignment_length) {}
  bool operator<(const PairedEndMappingWithBarcode &m) const {
    return std::tie(fragment_start_position_, fragment_length_, cell_barcode_,
                    mapq_, direction_, is_unique_, read_id_,
                    positive_alignment_length_, negative_alignment_length_) <
           std::tie(m.fragment_start_position_, m.fragment_length_,
                    m.cell_barcode_, m.mapq_, m.direction_, m.is_unique_,
                    m.read_id_, m.positive_alignment_length_,
                    m.negative_alignment_length_);
  }
  bool operator==(const PairedEndMappingWithBarcode &m) const {
    return std::tie(cell_barcode_, fragment_start_position_,
                    fragment_length_) == std::tie(m.cell_barcode_,
                                                  m.fragment_start_position_,
                                                  m.fragment_length_);
  }
  bool IsSamePosition(const PairedEndMappingWithBarcode &m) const {
    return std::tie(fragment_start_position_, fragment_length_) ==
           std::tie(m.fragment_start_position_, m.fragment_length_);
  }
  uint64_t GetBarcode() const { return cell_barcode_; }
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
};

class PairedEndMappingWithoutBarcode : public Mapping {
 public:
  uint32_t read_id_;
  uint32_t fragment_start_position_;
  uint16_t fragment_length_;
  uint8_t mapq_ : 6, direction_ : 1, is_unique_ : 1;
  uint8_t num_dups_;
  // uint8_t mapq;
  uint16_t positive_alignment_length_;
  uint16_t negative_alignment_length_;
  PairedEndMappingWithoutBarcode() : num_dups_(0) {}
  PairedEndMappingWithoutBarcode(uint32_t read_id,
                                 uint32_t fragment_start_position,
                                 uint16_t fragment_length, uint8_t mapq,
                                 uint8_t direction, uint8_t is_unique,
                                 uint16_t num_dups,
                                 uint16_t positive_alignment_length,
                                 uint16_t negative_alignment_length)
      : read_id_(read_id),
        fragment_start_position_(fragment_start_position),
        fragment_length_(fragment_length),
        mapq_(mapq),
        direction_(direction),
        is_unique_(is_unique),
        num_dups_(num_dups),
        positive_alignment_length_(positive_alignment_length),
        negative_alignment_length_(negative_alignment_length) {}

  bool operator<(const PairedEndMappingWithoutBarcode &m) const {
    return std::tie(fragment_start_position_, fragment_length_, mapq_,
                    direction_, is_unique_, read_id_,
                    positive_alignment_length_, negative_alignment_length_) <
           std::tie(m.fragment_start_position_, m.fragment_length_, m.mapq_,
                    m.direction_, m.is_unique_, m.read_id_,
                    m.positive_alignment_length_, m.negative_alignment_length_);
  }
  bool operator==(const PairedEndMappingWithoutBarcode &m) const {
    return std::tie(fragment_start_position_, fragment_length_) ==
           std::tie(m.fragment_start_position_, m.fragment_length_);
  }
  bool IsSamePosition(const PairedEndMappingWithoutBarcode &m) const {
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
};

}  // namespace chromap

#endif  // BEDMAPPING_H_
