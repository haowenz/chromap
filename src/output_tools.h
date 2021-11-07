#ifndef OUTPUTTOOLS_H_
#define OUTPUTTOOLS_H_

#include <assert.h>
#include <inttypes.h>

#include <cstring>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

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

/****************************
 **** CIGAR related macros ***
 *****************************/

#define BAM_CMATCH 0
#define BAM_CINS 1
#define BAM_CDEL 2
#define BAM_CREF_SKIP 3
#define BAM_CSOFT_CLIP 4
#define BAM_CHARD_CLIP 5
#define BAM_CPAD 6
#define BAM_CEQUAL 7
#define BAM_CDIFF 8
#define BAM_CBACK 9

#define BAM_CIGAR_STR "MIDNSHP=XB"
#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK 0xf
#define BAM_CIGAR_TYPE 0x3C1A7

/*! @abstract Table for converting a CIGAR operator character to BAM_CMATCH etc.
 * Result is operator code or -1. Be sure to cast the index if it is a plain
 *char: int op = bam_cigar_table[(unsigned char) ch];
 **/
// extern const int8_t bam_cigar_table[256];
const int8_t bam_cigar_table[256] = {
    // 0 .. 47
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,

    // 48 .. 63  (including =)
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, BAM_CEQUAL, -1, -1,

    // 64 .. 79  (including MIDNHB)
    -1, -1, BAM_CBACK, -1, BAM_CDEL, -1, -1, -1, BAM_CHARD_CLIP, BAM_CINS, -1,
    -1, -1, BAM_CMATCH, BAM_CREF_SKIP, -1,

    // 80 .. 95  (including SPX)
    BAM_CPAD, -1, -1, BAM_CSOFT_CLIP, -1, -1, -1, -1, BAM_CDIFF, -1, -1, -1, -1,
    -1, -1, -1,

    // 96 .. 127
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,

    // 128 .. 255
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
#define bam_cigar_op(c) ((c)&BAM_CIGAR_MASK)
#define bam_cigar_oplen(c) ((c) >> BAM_CIGAR_SHIFT)
// Note that BAM_CIGAR_STR is padded to length 16 bytes below so that
// the array look-up will not fall off the end.  '?' is chosen as the
// padding character so it's easy to spot if one is emitted, and will
// result in a parsing failure (in sam_parse1(), at least) if read.
#define bam_cigar_opchr(c) (BAM_CIGAR_STR "??????"[bam_cigar_op(c)])
#define bam_cigar_gen(l, o) ((l) << BAM_CIGAR_SHIFT | (o))

/* bam_cigar_type returns a bit flag with:
 *   bit 1 set if the cigar operation consumes the query
 *   bit 2 set if the cigar operation consumes the reference
 *
 * For reference, the unobfuscated truth table for this function is:
 * BAM_CIGAR_TYPE  QUERY  REFERENCE
 * --------------------------------
 * BAM_CMATCH      1      1
 * BAM_CINS        1      0
 * BAM_CDEL        0      1
 * BAM_CREF_SKIP   0      1
 * BAM_CSOFT_CLIP  1      0
 * BAM_CHARD_CLIP  0      0
 * BAM_CPAD        0      0
 * BAM_CEQUAL      1      1
 * BAM_CDIFF       1      1
 * BAM_CBACK       0      0
 * --------------------------------
 */
#define bam_cigar_type(o) (BAM_CIGAR_TYPE >> ((o) << 1) & 3)
// bit 1: consume query; bit 2: consume reference

/*! @abstract the read is paired in sequencing, no matter whether it is mapped
 * in a pair */
#define BAM_FPAIRED 1
/*! @abstract the read is mapped in a proper pair */
#define BAM_FPROPER_PAIR 2
/*! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
#define BAM_FUNMAP 4
/*! @abstract the mate is unmapped */
#define BAM_FMUNMAP 8
/*! @abstract the read is mapped to the reverse strand */
#define BAM_FREVERSE 16
/*! @abstract the mate is mapped to the reverse strand */
#define BAM_FMREVERSE 32
/*! @abstract this is read1 */
#define BAM_FREAD1 64
/*! @abstract this is read2 */
#define BAM_FREAD2 128
/*! @abstract not primary alignment */
#define BAM_FSECONDARY 256
/*! @abstract QC failure */
#define BAM_FQCFAIL 512
/*! @abstract optical or PCR duplicate */
#define BAM_FDUP 1024
/*! @abstract supplementary alignment */
#define BAM_FSUPPLEMENTARY 2048

template <typename MappingRecord>
bool ReadIdLess(const std::pair<uint32_t, MappingRecord> &a,
                const std::pair<uint32_t, MappingRecord> &b) {
  return a.second.read_id < b.second.read_id;
}

// When direction = 1, strand is positive
struct PAFMapping {
  uint32_t read_id;
  std::string read_name;
  uint16_t read_length;
  uint32_t fragment_start_position;
  uint16_t fragment_length;
  uint8_t mapq : 6, direction : 1, is_unique : 1;
  uint8_t num_dups;
  bool operator<(const PAFMapping &m) const {
    return std::tie(fragment_start_position, fragment_length, mapq, direction,
                    is_unique, read_id, read_length) <
           std::tie(m.fragment_start_position, m.fragment_length, m.mapq,
                    m.direction, m.is_unique, m.read_id, m.read_length);
  }
  bool operator==(const PAFMapping &m) const {
    return std::tie(fragment_start_position) ==
           std::tie(m.fragment_start_position);
  }
  bool HasSamePosition(const PAFMapping &m) const {
    return std::tie(fragment_start_position) ==
           std::tie(m.fragment_start_position);
  }
  uint32_t GetBarcode() const { return 0; }
  void Tn5Shift() {
    if (direction == 1) {
      fragment_start_position += 4;
    } else {
      fragment_length -= 5;
    }
  }
  bool IsPositive() const { return direction > 0 ? true : false; }
  uint32_t GetStartPosition() const {  // inclusive
    return fragment_start_position;
  }
  uint32_t GetEndPosition() const {  // exclusive
    return fragment_start_position + fragment_length;
  }
  uint16_t GetStructSize() const {
    return 2 * sizeof(uint32_t) + 2 * sizeof(uint16_t) + 2 * sizeof(uint8_t) +
           read_name.length() * sizeof(char);
  }
  size_t WriteToFile(FILE *temp_mapping_output_file) const {
    size_t num_written_bytes = 0;
    num_written_bytes +=
        fwrite(&read_id, sizeof(uint32_t), 1, temp_mapping_output_file);
    uint16_t read_name_length = read_name.length();
    num_written_bytes += fwrite(&read_name_length, sizeof(uint16_t), 1,
                                temp_mapping_output_file);
    num_written_bytes += fwrite(read_name.data(), sizeof(char),
                                read_name_length, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&read_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    num_written_bytes += fwrite(&fragment_start_position, sizeof(uint32_t), 1,
                                temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&fragment_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    uint8_t mapq_direction_is_unique =
        (mapq << 2) | (direction << 1) | is_unique;
    num_written_bytes += fwrite(&mapq_direction_is_unique, sizeof(uint8_t), 1,
                                temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&num_dups, sizeof(uint8_t), 1, temp_mapping_output_file);
    return num_written_bytes;
  }
  size_t LoadFromFile(FILE *temp_mapping_output_file) {
    size_t num_read_bytes = 0;
    num_read_bytes +=
        fread(&read_id, sizeof(uint32_t), 1, temp_mapping_output_file);
    uint16_t read_name_length = 0;
    num_read_bytes +=
        fread(&read_name_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    read_name = std::string(read_name_length, '\0');
    num_read_bytes += fread(&(read_name[0]), sizeof(char), read_name_length,
                            temp_mapping_output_file);
    num_read_bytes +=
        fread(&read_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    num_read_bytes += fread(&fragment_start_position, sizeof(uint32_t), 1,
                            temp_mapping_output_file);
    num_read_bytes +=
        fread(&fragment_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    uint8_t mapq_direction_is_unique = 0;
    num_read_bytes += fread(&mapq_direction_is_unique, sizeof(uint8_t), 1,
                            temp_mapping_output_file);
    mapq = (mapq_direction_is_unique >> 2);
    direction = (mapq_direction_is_unique >> 1) & 1;
    is_unique = mapq_direction_is_unique & 1;
    num_read_bytes +=
        fread(&num_dups, sizeof(uint8_t), 1, temp_mapping_output_file);
    return num_read_bytes;
  }
};

struct PairedPAFMapping {
  uint32_t read_id;
  std::string read1_name;
  std::string read2_name;
  uint16_t read1_length;
  uint16_t read2_length;
  uint32_t fragment_start_position;
  uint16_t fragment_length;
  uint16_t positive_alignment_length;
  uint16_t negative_alignment_length;
  uint8_t mapq;
  uint16_t mapq1 : 6, mapq2 : 6, direction : 1, is_unique : 1, reserved : 2;
  uint8_t num_dups;
  // uint8_t mapq; // least significant bit saves the direction of mapping
  bool operator<(const PairedPAFMapping &m) const {
    return std::tie(fragment_start_position, fragment_length, mapq1, mapq2,
                    direction, is_unique, read_id, positive_alignment_length,
                    negative_alignment_length) <
           std::tie(m.fragment_start_position, m.fragment_length, m.mapq1,
                    m.mapq2, m.direction, m.is_unique, m.read_id,
                    m.positive_alignment_length, m.negative_alignment_length);
  }
  bool operator==(const PairedPAFMapping &m) const {
    return std::tie(fragment_start_position, fragment_length) ==
           std::tie(m.fragment_start_position, m.fragment_length);
  }
  bool HasSamePosition(const PairedPAFMapping &m) const {
    return std::tie(fragment_start_position, fragment_length) ==
           std::tie(m.fragment_start_position, m.fragment_length);
  }
  uint32_t GetBarcode() const { return 0; }
  void Tn5Shift() {
    fragment_start_position += 4;
    positive_alignment_length -= 4;
    fragment_length -= 9;
    negative_alignment_length -= 5;
  }
  bool IsPositive() const { return direction > 0 ? true : false; }
  uint32_t GetStartPosition() const {  // inclusive
    return fragment_start_position;
  }
  uint32_t GetEndPosition() const {  // exclusive
    return fragment_start_position + fragment_length;
  }
  uint16_t GetStructSize() const {
    return 2 * sizeof(uint32_t) + 6 * sizeof(uint16_t) + 2 * sizeof(uint8_t) +
           (read1_name.length() + read2_name.length()) * sizeof(char);
  }
  size_t WriteToFile(FILE *temp_mapping_output_file) const {
    size_t num_written_bytes = 0;
    num_written_bytes +=
        fwrite(&read_id, sizeof(uint32_t), 1, temp_mapping_output_file);
    uint16_t read1_name_length = read1_name.length();
    num_written_bytes += fwrite(&read1_name_length, sizeof(uint16_t), 1,
                                temp_mapping_output_file);
    num_written_bytes += fwrite(read1_name.data(), sizeof(char),
                                read1_name_length, temp_mapping_output_file);
    uint16_t read2_name_length = read2_name.length();
    num_written_bytes += fwrite(&read2_name_length, sizeof(uint16_t), 1,
                                temp_mapping_output_file);
    num_written_bytes += fwrite(read2_name.data(), sizeof(char),
                                read2_name_length, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&read1_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&read2_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    num_written_bytes += fwrite(&fragment_start_position, sizeof(uint32_t), 1,
                                temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&fragment_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    num_written_bytes += fwrite(&positive_alignment_length, sizeof(uint16_t), 1,
                                temp_mapping_output_file);
    num_written_bytes += fwrite(&negative_alignment_length, sizeof(uint16_t), 1,
                                temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&mapq, sizeof(uint8_t), 1, temp_mapping_output_file);
    uint16_t mapq1_mapq2_direction_is_unique =
        (mapq1 << 10) | (mapq2 << 4) | (direction << 3) | (is_unique << 2);
    num_written_bytes += fwrite(&mapq1_mapq2_direction_is_unique,
                                sizeof(uint16_t), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&num_dups, sizeof(uint8_t), 1, temp_mapping_output_file);
    return num_written_bytes;
  }
  size_t LoadFromFile(FILE *temp_mapping_output_file) {
    size_t num_read_bytes = 0;
    num_read_bytes +=
        fread(&read_id, sizeof(uint32_t), 1, temp_mapping_output_file);
    uint16_t read1_name_length = 0;
    num_read_bytes += fread(&read1_name_length, sizeof(uint16_t), 1,
                            temp_mapping_output_file);
    read1_name = std::string(read1_name_length, '\0');
    num_read_bytes += fread(&(read1_name[0]), sizeof(char), read1_name_length,
                            temp_mapping_output_file);
    uint16_t read2_name_length = 0;
    num_read_bytes += fread(&read2_name_length, sizeof(uint16_t), 1,
                            temp_mapping_output_file);
    read2_name = std::string(read2_name_length, '\0');
    num_read_bytes += fread(&(read2_name[0]), sizeof(char), read2_name_length,
                            temp_mapping_output_file);
    num_read_bytes +=
        fread(&read1_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    num_read_bytes +=
        fread(&read2_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    num_read_bytes += fread(&fragment_start_position, sizeof(uint32_t), 1,
                            temp_mapping_output_file);
    num_read_bytes +=
        fread(&fragment_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    num_read_bytes += fread(&positive_alignment_length, sizeof(uint16_t), 1,
                            temp_mapping_output_file);
    num_read_bytes += fread(&negative_alignment_length, sizeof(uint16_t), 1,
                            temp_mapping_output_file);
    num_read_bytes +=
        fread(&mapq, sizeof(uint8_t), 1, temp_mapping_output_file);
    uint16_t mapq1_mapq2_direction_is_unique = 0;
    num_read_bytes += fread(&mapq1_mapq2_direction_is_unique, sizeof(uint16_t),
                            1, temp_mapping_output_file);
    mapq1 = (mapq1_mapq2_direction_is_unique >> 10);
    mapq2 = ((mapq1_mapq2_direction_is_unique << 6) >> 10);
    direction = (mapq1_mapq2_direction_is_unique >> 3) & 1;
    is_unique = (mapq1_mapq2_direction_is_unique >> 2) & 1;
    num_read_bytes +=
        fread(&num_dups, sizeof(uint8_t), 1, temp_mapping_output_file);
    return num_read_bytes;
  }
};

struct SAMMapping {
  uint32_t read_id;
  std::string read_name;
  uint64_t cell_barcode;
  // uint16_t read_length;
  // uint32_t fragment_start_position;
  // uint16_t fragment_length;
  // uint8_t direction : 1, is_unique : 1;
  uint8_t num_dups;
  // uint16_t positive_alignment_length;
  // uint16_t negative_alignment_length;

  int64_t pos;  // forward strand 5'-end mapping position (inclusive)
  int rid;      // reference sequence index in bntseq_t; <0 for unmapped
  int flag;     // extra flag
  uint32_t is_rev : 1, is_alt : 1, is_unique : 1, mapq : 7,
      NM : 22;  // is_rev: whether on the reverse strand; mapq: mapping quality;
                // NM: edit distance
  int n_cigar;  // number of CIGAR operations
  uint32_t *cigar;  // CIGAR in the BAM encoding: opLen<<4|op; op to integer
                    // mapping: MIDSH=>01234
  std::string MD;
  std::string sequence;
  std::string sequence_qual;
  // char *XA;        // alternative mappings
  // int score, sub, alt_sc;
  bool operator<(const SAMMapping &m) const {
    int read1_flag = flag & BAM_FREAD1;
    int m_read1_flag = m.flag & BAM_FREAD1;
    return std::tie(rid, pos, cell_barcode, mapq, read_id, read1_flag) <
           std::tie(m.rid, m.pos, m.cell_barcode, m.mapq, m.read_id,
                    m_read1_flag);
  }
  bool operator==(const SAMMapping &m) const {
    return std::tie(pos, rid, cell_barcode, is_rev) ==
           std::tie(m.pos, m.rid, m.cell_barcode, m.is_rev);
  }
  bool HasSamePosition(const SAMMapping &m) const {
    return std::tie(pos, rid, is_rev) == std::tie(m.pos, m.rid, m.is_rev);
  }
  uint64_t GetBarcode() const { return cell_barcode; }
  void Tn5Shift() {
    // We don't support Tn5 shift in SAM format because it has other fields that
    // depend mapping position.
  }
  bool IsPositive() const { return is_rev > 0 ? true : false; }
  // For now for convenience, we assume cigar should not be accessed after
  // generating the cigar string for output
  std::string GenerateCigarString() const {
    if (n_cigar == 0) {
      return "*";
    }
    std::string cigar_string = "";
    for (int ci = 0; ci < n_cigar; ++ci) {
      uint32_t op = bam_cigar_op(cigar[ci]);
      uint32_t op_length = bam_cigar_oplen(cigar[ci]);
      // std::cerr << op << " " << op_length << "\n";
      cigar_string.append(std::to_string(op_length));
      // cigar_string.append(std::to_string((BAM_CIGAR_STR[op])));
      cigar_string.push_back((BAM_CIGAR_STR[op]));
    }
    delete[] cigar;
    return cigar_string;
  }
  std::string GenerateIntTagString(const std::string &tag, int value) const {
    std::string tag_string = tag;
    tag_string.append(":i:" + std::to_string(value));
    return tag_string;
  }
  uint32_t GetAlignmentLength() const {
    uint32_t alignment_length = 0;
    for (int ci = 0; ci < n_cigar; ++ci) {
      uint32_t op = bam_cigar_op(cigar[ci]);
      uint32_t op_length = bam_cigar_oplen(cigar[ci]);
      if ((bam_cigar_type(op) & 0x2) > 0) {
        alignment_length += op_length;
      }
    }
    return alignment_length;
  }
  uint32_t GetStartPosition() const {  // inclusive
    return pos + 1;
    /*if (IsPositive()) {
      return pos + 1;
    } else {
      return pos + 1 - GetAlignmentLength() + 1;
    }*/
  }
  uint32_t GetEndPosition() const {  // exclusive
    return pos + GetAlignmentLength();
    /*if (IsPositive()) {
      return pos + GetAlignmentLength();
    } else {
      return pos + 1;
    }*/
  }
  uint16_t GetStructSize() const {
    return 2 * sizeof(uint32_t) + 2 * sizeof(uint16_t) + 2 * sizeof(uint8_t) +
           (read_name.length() + MD.length()) * sizeof(char) +
           n_cigar * sizeof(uint32_t);
  }
  size_t WriteToFile(FILE *temp_mapping_output_file) const {
    size_t num_written_bytes = 0;
    num_written_bytes +=
        fwrite(&read_id, sizeof(uint32_t), 1, temp_mapping_output_file);
    uint16_t read_name_length = read_name.length();
    num_written_bytes += fwrite(&read_name_length, sizeof(uint16_t), 1,
                                temp_mapping_output_file);
    num_written_bytes += fwrite(read_name.data(), sizeof(char),
                                read_name_length, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&cell_barcode, sizeof(uint64_t), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&num_dups, sizeof(uint8_t), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&pos, sizeof(int64_t), 1, temp_mapping_output_file);
    num_written_bytes += fwrite(&rid, sizeof(int), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&flag, sizeof(int), 1, temp_mapping_output_file);
    uint32_t rev_alt_unique_mapq_NM =
        (is_rev << 31) | (is_alt << 30) | (is_unique << 29) | (mapq << 22) | NM;
    num_written_bytes += fwrite(&rev_alt_unique_mapq_NM, sizeof(uint32_t), 1,
                                temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&n_cigar, sizeof(int), 1, temp_mapping_output_file);
    if (n_cigar > 0) {
      num_written_bytes +=
          fwrite(cigar, sizeof(uint32_t), n_cigar, temp_mapping_output_file);
      delete[] cigar;
    }
    uint16_t MD_length = MD.length();
    num_written_bytes +=
        fwrite(&MD_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    if (MD_length > 0) {
      num_written_bytes +=
          fwrite(MD.data(), sizeof(char), MD_length, temp_mapping_output_file);
    }
    uint16_t sequence_length = sequence.length();
    num_written_bytes +=
        fwrite(&sequence_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    num_written_bytes += fwrite(sequence.data(), sizeof(char), sequence_length,
                                temp_mapping_output_file);
    num_written_bytes += fwrite(sequence_qual.data(), sizeof(char),
                                sequence_length, temp_mapping_output_file);
    return num_written_bytes;
  }
  size_t LoadFromFile(FILE *temp_mapping_output_file) {
    int num_read_bytes = 0;
    num_read_bytes +=
        fread(&read_id, sizeof(uint32_t), 1, temp_mapping_output_file);
    uint16_t read_name_length = 0;
    num_read_bytes +=
        fread(&read_name_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    read_name = std::string(read_name_length, '\0');
    num_read_bytes += fread(&(read_name[0]), sizeof(char), read_name_length,
                            temp_mapping_output_file);
    num_read_bytes +=
        fread(&cell_barcode, sizeof(uint64_t), 1, temp_mapping_output_file);
    num_read_bytes +=
        fread(&num_dups, sizeof(uint8_t), 1, temp_mapping_output_file);
    num_read_bytes += fread(&pos, sizeof(int64_t), 1, temp_mapping_output_file);
    num_read_bytes += fread(&rid, sizeof(int), 1, temp_mapping_output_file);
    num_read_bytes += fread(&flag, sizeof(int), 1, temp_mapping_output_file);
    uint32_t rev_alt_unique_mapq_NM = 0;
    num_read_bytes += fread(&rev_alt_unique_mapq_NM, sizeof(uint32_t), 1,
                            temp_mapping_output_file);
    is_rev = (rev_alt_unique_mapq_NM >> 31);
    is_alt = (rev_alt_unique_mapq_NM >> 30) & 1;
    is_unique = (rev_alt_unique_mapq_NM >> 29) & 1;
    mapq = ((rev_alt_unique_mapq_NM << 3) >> 25);
    NM = ((rev_alt_unique_mapq_NM << 10) >> 10);
    num_read_bytes += fread(&n_cigar, sizeof(int), 1, temp_mapping_output_file);
    if (n_cigar > 0) {
      cigar = new uint32_t[n_cigar];
      num_read_bytes +=
          fread(cigar, sizeof(uint32_t), n_cigar, temp_mapping_output_file);
    }
    uint16_t MD_length = 0;
    num_read_bytes +=
        fread(&MD_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    if (MD_length > 0) {
      MD = std::string(MD_length, '\0');
      num_read_bytes +=
          fread(&(MD[0]), sizeof(char), MD_length, temp_mapping_output_file);
    }
    uint16_t sequence_length = 0;
    num_read_bytes +=
        fread(&sequence_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    if (sequence_length > 0) {
      sequence = std::string(sequence_length, '\0');
      sequence_qual = std::string(sequence_length, '\0');
      num_read_bytes += fread(&(sequence[0]), sizeof(char), sequence_length,
                              temp_mapping_output_file);
      num_read_bytes += fread(&(sequence_qual[0]), sizeof(char),
                              sequence_length, temp_mapping_output_file);
    }
    return num_read_bytes;
  }
};

// TODO(Haowen) : Add PairedSAMMapping.

// Format for pairtools for HiC data.
struct PairsMapping {
  uint32_t read_id;
  std::string read_name;
  uint64_t cell_barcode;
  int rid1;
  int rid2;
  uint32_t pos1;
  uint32_t pos2;
  int direction1;  // 1-positive. 0-negative
  int direction2;
  uint16_t mapq : 8, is_unique : 1, num_dups : 7;

  bool operator<(const PairsMapping &m) const {
    return std::tie(rid1, rid2, pos1, pos2, mapq, read_id) <
           std::tie(m.rid1, m.rid2, m.pos1, m.pos2, m.mapq, m.read_id);
  }
  bool operator==(const PairsMapping &m) const {
    return std::tie(rid1, pos1, rid2, pos2) ==
           std::tie(m.rid1, m.pos1, m.rid2, m.pos2);
    // return std::tie(pos1, pos2, rid1, rid2, is_rev1, is_rev2) ==
    // std::tie(m.pos1, m.pos2, m.rid1, m.rid2, m.is_rev1, m.is_rev2);
  }
  bool HasSamePosition(const PairsMapping &m) const {
    return std::tie(rid1, pos1, rid2, pos2) ==
           std::tie(m.rid1, m.pos1, m.rid2, m.pos2);
  }
  uint32_t GetBarcode() const { return 0; }
  void Tn5Shift() {
    // We don't support Tn5 shift in SAM format because it has other fields that
    // depend mapping position.
  }

  int GetPosition(int idx) const {
    if (idx == 2) {
      return pos2 + 1;
    }
    return pos1 + 1;
  }

  char GetDirection(int idx) const {
    int d = direction1;
    if (idx == 2) {
      d = direction2;
    }
    return d > 0 ? '+' : '-';
  }

  bool IsPositive() const { return direction1 > 0 ? true : false; }
  uint32_t GetStartPosition() const {  // inclusive
    return pos1;
  }
  uint32_t GetEndPosition() const {  // exclusive
    return pos2;
  }
  uint16_t GetStructSize() const {
    return 5 * sizeof(uint32_t) + 1 * sizeof(uint16_t) + 4 * sizeof(int) +
           read_name.length() * sizeof(char);
  }
  size_t WriteToFile(FILE *temp_mapping_output_file) const {
    size_t num_written_bytes = 0;
    num_written_bytes +=
        fwrite(&read_id, sizeof(uint32_t), 1, temp_mapping_output_file);
    uint16_t read_name_length = read_name.length();
    num_written_bytes += fwrite(&read_name_length, sizeof(uint16_t), 1,
                                temp_mapping_output_file);
    num_written_bytes += fwrite(read_name.data(), sizeof(char),
                                read_name_length, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&cell_barcode, sizeof(uint64_t), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&rid1, sizeof(int), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&rid2, sizeof(int), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&pos1, sizeof(uint32_t), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&pos2, sizeof(uint32_t), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&direction1, sizeof(int), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&direction2, sizeof(int), 1, temp_mapping_output_file);
    uint16_t mapq_unique_dups = (mapq << 8) | (is_unique << 7) | num_dups;
    num_written_bytes += fwrite(&mapq_unique_dups, sizeof(uint16_t), 1,
                                temp_mapping_output_file);
    return num_written_bytes;
  }
  size_t LoadFromFile(FILE *temp_mapping_output_file) {
    size_t num_read_bytes = 0;
    num_read_bytes +=
        fread(&read_id, sizeof(uint32_t), 1, temp_mapping_output_file);
    uint16_t read_name_length = 0;
    num_read_bytes +=
        fread(&read_name_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    read_name = std::string(read_name_length, '\0');
    num_read_bytes += fread(&(read_name[0]), sizeof(char), read_name_length,
                            temp_mapping_output_file);
    num_read_bytes +=
        fread(&cell_barcode, sizeof(uint64_t), 1, temp_mapping_output_file);
    num_read_bytes += fread(&rid1, sizeof(int), 1, temp_mapping_output_file);
    num_read_bytes += fread(&rid2, sizeof(int), 1, temp_mapping_output_file);
    num_read_bytes +=
        fread(&pos1, sizeof(uint32_t), 1, temp_mapping_output_file);
    num_read_bytes +=
        fread(&pos2, sizeof(uint32_t), 1, temp_mapping_output_file);
    num_read_bytes +=
        fread(&direction1, sizeof(int), 1, temp_mapping_output_file);
    num_read_bytes +=
        fread(&direction2, sizeof(int), 1, temp_mapping_output_file);
    uint16_t mapq_unique_dups = 0;
    num_read_bytes +=
        fread(&mapq_unique_dups, sizeof(uint16_t), 1, temp_mapping_output_file);
    mapq = (mapq_unique_dups >> 8);
    is_unique = (mapq_unique_dups >> 7) & 1;
    num_dups = ((mapq_unique_dups << 9) >> 9);
    return num_read_bytes;
  }
};

struct MappingWithBarcode {
  uint32_t read_id;
  uint64_t cell_barcode;
  uint32_t fragment_start_position;
  uint16_t fragment_length;
  uint8_t mapq : 6, direction : 1, is_unique : 1;
  uint8_t num_dups;
  // uint8_t mapq;
  bool operator<(const MappingWithBarcode &m) const {
    return std::tie(fragment_start_position, fragment_length, cell_barcode,
                    mapq, direction, is_unique, read_id) <
           std::tie(m.fragment_start_position, m.fragment_length,
                    m.cell_barcode, m.mapq, m.direction, m.is_unique,
                    m.read_id);
  }
  bool operator==(const MappingWithBarcode &m) const {
    return std::tie(cell_barcode, fragment_start_position) ==
           std::tie(m.cell_barcode, m.fragment_start_position);
  }
  bool HasSamePosition(const MappingWithBarcode &m) const {
    return std::tie(fragment_start_position) ==
           std::tie(m.fragment_start_position);
  }
  uint64_t GetBarcode() const { return cell_barcode; }
  void Tn5Shift() {
    if (direction == 1) {
      fragment_start_position += 4;
    } else {
      fragment_length -= 5;
    }
  }
  bool IsPositive() const { return direction > 0 ? true : false; }
  uint32_t GetStartPosition() const {  // inclusive
    return fragment_start_position;
  }
  uint32_t GetEndPosition() const {  // exclusive
    return fragment_start_position + fragment_length;
  }
};

struct MappingWithoutBarcode {
  uint32_t read_id;
  uint32_t fragment_start_position;
  uint16_t fragment_length;
  // uint8_t mapq;
  uint8_t mapq : 6, direction : 1, is_unique : 1;
  uint8_t num_dups;
  bool operator<(const MappingWithoutBarcode &m) const {
    return std::tie(fragment_start_position, fragment_length, mapq, direction,
                    is_unique, read_id) <
           std::tie(m.fragment_start_position, m.fragment_length, m.mapq,
                    m.direction, m.is_unique, m.read_id);
  }
  bool operator==(const MappingWithoutBarcode &m) const {
    return std::tie(fragment_start_position) ==
           std::tie(m.fragment_start_position);
  }
  bool HasSamePosition(const MappingWithBarcode &m) const {
    return std::tie(fragment_start_position) ==
           std::tie(m.fragment_start_position);
  }
  uint32_t GetBarcode() const { return 0; }
  void Tn5Shift() {
    if (direction == 1) {
      fragment_start_position += 4;
    } else {
      fragment_length -= 5;
    }
  }
  bool IsPositive() const { return direction > 0 ? true : false; }
  uint32_t GetStartPosition() const {  // inclusive
    return fragment_start_position;
  }
  uint32_t GetEndPosition() const {  // exclusive
    return fragment_start_position + fragment_length;
  }
};

struct PairedEndMappingWithBarcode {
  uint32_t read_id;
  uint64_t cell_barcode;
  uint32_t fragment_start_position;
  uint16_t fragment_length;
  uint8_t mapq : 6, direction : 1, is_unique : 1;
  uint8_t num_dups;
  // uint8_t mapq;
  uint16_t positive_alignment_length;
  uint16_t negative_alignment_length;
  bool operator<(const PairedEndMappingWithBarcode &m) const {
    return std::tie(fragment_start_position, fragment_length, cell_barcode,
                    mapq, direction, is_unique, read_id,
                    positive_alignment_length, negative_alignment_length) <
           std::tie(m.fragment_start_position, m.fragment_length,
                    m.cell_barcode, m.mapq, m.direction, m.is_unique, m.read_id,
                    m.positive_alignment_length, m.negative_alignment_length);
  }
  bool operator==(const PairedEndMappingWithBarcode &m) const {
    return std::tie(cell_barcode, fragment_start_position, fragment_length) ==
           std::tie(m.cell_barcode, m.fragment_start_position,
                    m.fragment_length);
  }
  bool HasSamePosition(const PairedEndMappingWithBarcode &m) const {
    return std::tie(fragment_start_position, fragment_length) ==
           std::tie(m.fragment_start_position, m.fragment_length);
  }
  uint64_t GetBarcode() const { return cell_barcode; }
  void Tn5Shift() {
    fragment_start_position += 4;
    positive_alignment_length -= 4;
    fragment_length -= 9;
    negative_alignment_length -= 5;
  }
  bool IsPositive() const { return direction > 0 ? true : false; }
  uint32_t GetStartPosition() const {  // inclusive
    return fragment_start_position;
  }
  uint32_t GetEndPosition() const {  // exclusive
    return fragment_start_position + fragment_length;
  }
};

struct PairedEndMappingWithoutBarcode {
  uint32_t read_id;
  uint32_t fragment_start_position;
  uint16_t fragment_length;
  uint8_t mapq : 6, direction : 1, is_unique : 1;
  uint8_t num_dups;
  // uint8_t mapq;
  uint16_t positive_alignment_length;
  uint16_t negative_alignment_length;
  bool operator<(const PairedEndMappingWithoutBarcode &m) const {
    return std::tie(fragment_start_position, fragment_length, mapq, direction,
                    is_unique, read_id, positive_alignment_length,
                    negative_alignment_length) <
           std::tie(m.fragment_start_position, m.fragment_length, m.mapq,
                    m.direction, m.is_unique, m.read_id,
                    m.positive_alignment_length, m.negative_alignment_length);
  }
  bool operator==(const PairedEndMappingWithoutBarcode &m) const {
    return std::tie(fragment_start_position, fragment_length) ==
           std::tie(m.fragment_start_position, m.fragment_length);
  }
  bool HasSamePosition(const PairedEndMappingWithoutBarcode &m) const {
    return std::tie(fragment_start_position, fragment_length) ==
           std::tie(m.fragment_start_position, m.fragment_length);
  }
  uint32_t GetBarcode() const { return 0; }
  void Tn5Shift() {
    fragment_start_position += 4;
    positive_alignment_length -= 4;
    fragment_length -= 9;
    negative_alignment_length -= 5;
  }
  bool IsPositive() const { return direction > 0 ? true : false; }
  uint32_t GetStartPosition() const {  // inclusive
    return fragment_start_position;
  }
  uint32_t GetEndPosition() const {  // exclusive
    return fragment_start_position + fragment_length;
  }
};

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
