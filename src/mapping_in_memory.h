#ifndef MAPPING_IN_MEMORY_H_
#define MAPPING_IN_MEMORY_H_

#include <stdint.h>

#include <string>

#include "utils.h"

namespace chromap {

// Regardless of mapping format, this struct can temporarily hold a mapping in
// memory for easily passing it into the several functions before pushing it
// to the result mapping vector. It never owns the read or the read qual. It
// owns the cigar before push the mapping to the vector. (For now, the cigar
// memory is released once a SAMMapping is created.) Since this struct is large,
// we should never create a huge vector of this struct.
struct MappingInMemory {
  uint32_t read_id = 0;
  int read_split_site = 0;
  int read_length = 0;

  uint32_t rid = 0;
  uint32_t ref_start_position = 0;
  uint32_t ref_end_position = 0;

  uint64_t barcode_key = 0;

  Strand strand = kPositive;
  bool is_unique = true;
  uint8_t mapq = 0;

  // It does NOT own read or read qual.
  const char *read_name = nullptr;
  const char *read_sequence = nullptr;
  const char *qual_sequence = nullptr;

  // SAM fields or tags.
  uint16_t SAM_flag = 0;

  uint32_t *cigar = nullptr;
  int n_cigar = 0;

  int NM = 0;

  std::string MD_tag;

  inline uint8_t GetStrand() const { return (strand == kPositive ? 1 : 0); }

  inline uint32_t GetFragmentStartPosition() const {
    return ref_start_position;
  }

  // TODO(Haowen): change this to alignment length.
  inline uint16_t GetFragmentLength() const {
    return ref_end_position - ref_start_position + 1;
  }

  inline uint16_t GetAlignmentLength() const {
    return ref_end_position - ref_start_position + 1;
  }
};

struct PairedEndMappingInMemory {
  MappingInMemory mapping_in_memory1;
  MappingInMemory mapping_in_memory2;
  uint8_t is_unique;
  uint8_t mapq;

  inline uint8_t GetStrand() const {
    return (mapping_in_memory1.strand == kPositive ? 1 : 0);
  }

  inline uint32_t GetReadId() const { return mapping_in_memory1.read_id; }

  inline uint64_t GetBarcode() const { return mapping_in_memory1.barcode_key; }

  inline uint32_t GetFragmentStartPosition() const {
    if (mapping_in_memory1.strand == kPositive) {
      return mapping_in_memory1.GetFragmentStartPosition();
    }

    return mapping_in_memory2.GetFragmentStartPosition();
  }

  inline int GetFragmentLength() const {
    if (mapping_in_memory1.strand == kPositive) {
      return mapping_in_memory2.ref_end_position -
             mapping_in_memory1.ref_start_position + 1;
    }
    return mapping_in_memory1.ref_end_position -
           mapping_in_memory2.ref_start_position + 1;
  }

  inline uint32_t GetPositiveAlignmentLength() const {
    if (mapping_in_memory1.strand == kPositive) {
      return mapping_in_memory1.GetAlignmentLength();
    }
    return mapping_in_memory2.GetAlignmentLength();
  }

  inline uint32_t GetNegativeAlignmentLength() const {
    if (mapping_in_memory1.strand == kNegative) {
      return mapping_in_memory1.GetAlignmentLength();
    }
    return mapping_in_memory2.GetAlignmentLength();
  }
};

}  // namespace chromap

#endif  // MAPPING_IN_MEMORY_H_
