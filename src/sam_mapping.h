#ifndef SAMMAPPING_H_
#define SAMMAPPING_H_

#include <string>
#include <tuple>
#include <vector>

#include "mapping.h"

namespace chromap {

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

class SAMMapping : public Mapping {
 public:
  uint32_t read_id_;
  std::string read_name_;
  uint64_t cell_barcode_;
  // uint16_t read_length;
  // uint32_t fragment_start_position;
  // uint16_t fragment_length;
  // uint8_t direction : 1, is_unique : 1;
  uint8_t num_dups_;
  // uint16_t positive_alignment_length;
  // uint16_t negative_alignment_length;

  int64_t pos_;  // forward strand 5'-end mapping position (inclusive)
  int rid_;      // reference sequence index in bntseq_t; <0 for unmapped
  int64_t mpos_;  // forward strand 5'-end mapping position for mate (inclusive)
  int mrid_;      // reference sequence index in bntseq_t; <0 for unmapped
  int tlen_;      // template length
  int flag_;     // extra flag
  uint32_t is_rev_ : 1, is_alt_ : 1, is_unique_ : 1, mapq_ : 7,
      NM_ : 22;      // is_rev: whether on the reverse strand; mapq: mapping
                     // quality; NM: edit distance
  int n_cigar_ = 0;  // number of CIGAR operations
  std::vector<uint32_t> cigar_;  // CIGAR in the BAM encoding: opLen<<4|op; op
                                 // to integer mapping: MIDSH=>01234
  std::string MD_;
  std::string sequence_;
  std::string sequence_qual_;
  // char *XA;        // alternative mappings
  // int score, sub, alt_sc;

  SAMMapping() {}

  SAMMapping(uint32_t read_id, const std::string &read_name,
             uint64_t cell_barcode, uint8_t num_dups, int64_t pos, int rid,
             int64_t mpos, int mrid, int tlen, 
             int flag, uint8_t is_rev, uint8_t is_alt, uint8_t is_unique,
             uint8_t mapq, uint32_t NM, int n_cigar, uint32_t *cigar,
             const std::string &MD_tag, const std::string &sequence,
             const std::string &sequence_qual)
      : read_id_(read_id),
        read_name_(read_name),
        cell_barcode_(cell_barcode),
        num_dups_(num_dups),
        pos_(pos),
        rid_(rid),
        mpos_(mpos),
        mrid_(mrid),
        tlen_(tlen),
        flag_(flag),
        is_rev_(is_rev),
        is_alt_(is_alt),
        is_unique_(is_unique),
        mapq_(mapq),
        NM_(NM),
        n_cigar_(n_cigar),
        MD_(MD_tag) {
    cigar_ = std::vector<uint32_t>(cigar, cigar + n_cigar);
    free(cigar);

    if (!IsPositiveStrand()) {
      for (uint32_t i = 0; i < sequence_qual.length(); ++i) {
        sequence_qual_.push_back(sequence_qual[sequence_qual.length() - 1 - i]);
      }
    } else {
      sequence_qual_ = sequence_qual;
    }

    uint32_t sequence_length_deduced_from_cigar = GetSequenceLength();
    if (sequence_length_deduced_from_cigar != sequence.length()) {
      sequence_ = sequence.substr(0, sequence_length_deduced_from_cigar);
      sequence_qual_ =
          sequence_qual_.substr(0, sequence_length_deduced_from_cigar);
    } else {
      sequence_ = sequence;
    }
  }

  bool operator<(const SAMMapping &m) const {
    int read1_flag = flag_ & BAM_FREAD1;
    int m_read1_flag = m.flag_ & BAM_FREAD1;
    return std::tie(rid_, pos_, cell_barcode_, mapq_, mrid_, mpos_, read_id_, read1_flag) <
           std::tie(m.rid_, m.pos_, m.cell_barcode_, m.mapq_, m.mrid_, m.mpos_, m.read_id_,
                    m_read1_flag);
  }
  bool operator==(const SAMMapping &m) const {
    return std::tie(pos_, rid_, cell_barcode_, is_rev_, mrid_, mpos_) ==
           std::tie(m.pos_, m.rid_, m.cell_barcode_, m.is_rev_, mrid_, mpos_);
  }
  bool IsSamePosition(const SAMMapping &m) const {
    return std::tie(pos_, rid_, is_rev_) == std::tie(m.pos_, m.rid_, m.is_rev_);
  }
  uint64_t GetBarcode() const { return cell_barcode_; }
  void Tn5Shift() {
    // We don't support Tn5 shift in SAM format because it has other fields that
    // depend mapping position.
  }
  // TODO(Haowen): I have to change the variable names or this function to make
  // the meaning consistent.
  bool IsPositiveStrand() const { return is_rev_ > 0 ? true : false; }
  // For now for convenience, we assume cigar should not be accessed after
  // generating the cigar string for output
  std::string GenerateCigarString() const {
    if (n_cigar_ == 0) {
      return "*";
    }
    std::string cigar_string = "";
    for (int ci = 0; ci < n_cigar_; ++ci) {
      uint32_t op = bam_cigar_op(cigar_[ci]);
      uint32_t op_length = bam_cigar_oplen(cigar_[ci]);
      // std::cerr << op << " " << op_length << "\n";
      cigar_string.append(std::to_string(op_length));
      // cigar_string.append(std::to_string((BAM_CIGAR_STR[op])));
      cigar_string.push_back((BAM_CIGAR_STR[op]));
    }
    return cigar_string;
  }
  std::string GenerateIntTagString(const std::string &tag, int value) const {
    std::string tag_string = tag;
    tag_string.append(":i:" + std::to_string(value));
    return tag_string;
  }
  uint32_t GetAlignmentLength() const {
    uint32_t alignment_length = 0;
    for (int ci = 0; ci < n_cigar_; ++ci) {
      uint32_t op = bam_cigar_op(cigar_[ci]);
      uint32_t op_length = bam_cigar_oplen(cigar_[ci]);
      if ((bam_cigar_type(op) & 0x2) > 0) {
        alignment_length += op_length;
      }
    }
    return alignment_length;
  }

  uint32_t GetSequenceLength() const {
    uint32_t sequence_length = 0;
    for (int ci = 0; ci < n_cigar_; ++ci) {
      uint32_t op = bam_cigar_op(cigar_[ci]);
      uint32_t op_length = bam_cigar_oplen(cigar_[ci]);
      if ((bam_cigar_type(op) & 0x1) > 0) {
        sequence_length += op_length;
      }
    }
    return sequence_length;
  }

  uint32_t GetStartPosition() const {  // inclusive
    return pos_ + 1;
    /*if (IsPositiveStrand()) {
      return pos + 1;
    } else {
      return pos + 1 - GetAlignmentLength() + 1;
    }*/
  }
  uint32_t GetEndPosition() const {  // exclusive
    return pos_ + GetAlignmentLength();
    /*if (IsPositiveStrand()) {
      return pos + GetAlignmentLength();
    } else {
      return pos + 1;
    }*/
  }
  uint16_t GetByteSize() const {
    return 2 * sizeof(uint32_t) + 2 * sizeof(uint16_t) + 2 * sizeof(uint8_t) +
           (read_name_.length() + MD_.length()) * sizeof(char) +
           n_cigar_ * sizeof(uint32_t);
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
        fwrite(&num_dups_, sizeof(uint8_t), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&pos_, sizeof(int64_t), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&rid_, sizeof(int), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&flag_, sizeof(int), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&mpos_, sizeof(int64_t), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&mrid_, sizeof(int), 1, temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&tlen_, sizeof(int), 1, temp_mapping_output_file);
    uint32_t rev_alt_unique_mapq_NM = (is_rev_ << 31) | (is_alt_ << 30) |
                                      (is_unique_ << 29) | (mapq_ << 22) | NM_;
    num_written_bytes += fwrite(&rev_alt_unique_mapq_NM, sizeof(uint32_t), 1,
                                temp_mapping_output_file);
    num_written_bytes +=
        fwrite(&n_cigar_, sizeof(int), 1, temp_mapping_output_file);
    if (n_cigar_ > 0) {
      num_written_bytes += fwrite(cigar_.data(), sizeof(uint32_t), n_cigar_,
                                  temp_mapping_output_file);
    }
    uint16_t MD_length = MD_.length();
    num_written_bytes +=
        fwrite(&MD_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    if (MD_length > 0) {
      num_written_bytes +=
          fwrite(MD_.data(), sizeof(char), MD_length, temp_mapping_output_file);
    }
    uint16_t sequence_length = sequence_.length();
    num_written_bytes +=
        fwrite(&sequence_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    num_written_bytes += fwrite(sequence_.data(), sizeof(char), sequence_length,
                                temp_mapping_output_file);
    num_written_bytes += fwrite(sequence_qual_.data(), sizeof(char),
                                sequence_length, temp_mapping_output_file);
    return num_written_bytes;
  }

  size_t LoadFromFile(FILE *temp_mapping_output_file) {
    int num_read_bytes = 0;
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
    num_read_bytes +=
        fread(&num_dups_, sizeof(uint8_t), 1, temp_mapping_output_file);
    num_read_bytes +=
        fread(&pos_, sizeof(int64_t), 1, temp_mapping_output_file);
    num_read_bytes += fread(&rid_, sizeof(int), 1, temp_mapping_output_file);
    num_read_bytes +=
        fread(&mpos_, sizeof(int64_t), 1, temp_mapping_output_file);
    num_read_bytes += fread(&mrid_, sizeof(int), 1, temp_mapping_output_file);
    num_read_bytes += fread(&tlen_, sizeof(int), 1, temp_mapping_output_file);
    num_read_bytes += fread(&flag_, sizeof(int), 1, temp_mapping_output_file);
    uint32_t rev_alt_unique_mapq_NM = 0;
    num_read_bytes += fread(&rev_alt_unique_mapq_NM, sizeof(uint32_t), 1,
                            temp_mapping_output_file);
    is_rev_ = (rev_alt_unique_mapq_NM >> 31);
    is_alt_ = (rev_alt_unique_mapq_NM >> 30) & 1;
    is_unique_ = (rev_alt_unique_mapq_NM >> 29) & 1;
    mapq_ = ((rev_alt_unique_mapq_NM << 3) >> 25);
    NM_ = ((rev_alt_unique_mapq_NM << 10) >> 10);
    int previous_n_cigar_ = n_cigar_;
    num_read_bytes +=
        fread(&n_cigar_, sizeof(int), 1, temp_mapping_output_file);
    if (n_cigar_ > 0) {
      if (previous_n_cigar_ < n_cigar_) {
        cigar_.resize(n_cigar_);
      }
      num_read_bytes += fread(cigar_.data(), sizeof(uint32_t), n_cigar_,
                              temp_mapping_output_file);
    }
    uint16_t MD_length = 0;
    num_read_bytes +=
        fread(&MD_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    if (MD_length > 0) {
      MD_ = std::string(MD_length, '\0');
      num_read_bytes +=
          fread(&(MD_[0]), sizeof(char), MD_length, temp_mapping_output_file);
    }
    uint16_t sequence_length = 0;
    num_read_bytes +=
        fread(&sequence_length, sizeof(uint16_t), 1, temp_mapping_output_file);
    if (sequence_length > 0) {
      sequence_ = std::string(sequence_length, '\0');
      sequence_qual_ = std::string(sequence_length, '\0');
      num_read_bytes += fread(&(sequence_[0]), sizeof(char), sequence_length,
                              temp_mapping_output_file);
      num_read_bytes += fread(&(sequence_qual_[0]), sizeof(char),
                              sequence_length, temp_mapping_output_file);
    }
    return num_read_bytes;
  }
};

// TODO(Haowen) : Add PairedSAMMapping.
}  // namespace chromap

#endif  // SAMMAPPING_H_
