// Area C: the mapping record classes (bed_mapping.h, paf_mapping.h,
// sam_mapping.h) and the Mapping interface they implement.
//
// Every expected value is worked by hand in the comment above its check:
// sort keys and duplicate keys from the record definitions, Tn5 shift +4/-5,
// CIGAR strings and lengths from the SAM specification (M, I, S, =, X consume
// the query; M, D, N, =, X consume the reference), the BAM flag bits from the
// SAM specification, and round trips through the temp-file encoding.

#include "chromap_test.h"

#include <cstdio>

#include "bed_mapping.h"
#include "mapping.h"
#include "paf_mapping.h"
#include "sam_mapping.h"

using namespace chromap;

namespace {

// The SAMMapping constructor takes ownership of a malloc'ed CIGAR array and
// frees it, as chromap's ksw does.
uint32_t *MakeCigar(const std::string &ops, const std::vector<uint32_t> &lengths, int *n) {
  uint32_t *cigar = (uint32_t *)malloc(sizeof(uint32_t) * (lengths.size() + 1));
  for (size_t i = 0; i < lengths.size(); ++i) {
    cigar[i] = bam_cigar_gen(lengths[i], (uint32_t)bam_cigar_table[(unsigned char)ops[i]]);
  }
  *n = (int)lengths.size();
  return cigar;
}

SAMMapping MakeSam(uint32_t read_id, const std::string &name, uint64_t barcode, int64_t pos, int rid,
                   int64_t mpos, int mrid, int tlen, int flag, uint8_t positive, uint8_t mapq,
                   const std::string &ops, const std::vector<uint32_t> &lengths, const std::string &seq,
                   const std::string &qual, uint32_t nm = 0, const std::string &md = "50") {
  int n = 0;
  uint32_t *cigar = MakeCigar(ops, lengths, &n);
  return SAMMapping(read_id, name, barcode, /*num_dups=*/1, pos, rid, mpos, mrid, tlen, flag,
                    /*is_rev (1 = positive strand, as chromap passes it)=*/positive,
                    /*is_alt=*/0, /*is_unique=*/1, mapq, nm, n, cigar, md, seq, qual);
}

std::string Repeat(const std::string &unit, size_t n) {
  std::string s;
  while (s.size() < n) s += unit;
  return s.substr(0, n);
}

}  // namespace

// MappingWithoutBarcode sorts by start, then length, MAPQ, strand (0 before
// 1), uniqueness and read id; duplicates share a start whatever their strand
// or length.
CT_TEST(single_end_record_order_and_duplicates) {
  MappingWithoutBarcode a(/*read_id=*/9, /*start=*/100, /*length=*/50, /*mapq=*/60, /*direction=*/1, 1, 1);
  MappingWithoutBarcode b(/*read_id=*/1, /*start=*/101, 50, 60, 1, 1, 1);
  CT_CHECK(a < b);  // start dominates the read id
  CT_CHECK(!(b < a));
  MappingWithoutBarcode shorter(2, 100, 40, 60, 1, 1, 1);
  CT_CHECK(shorter < a);  // same start: shorter first
  MappingWithoutBarcode low_mapq(3, 100, 50, 10, 1, 1, 1);
  CT_CHECK(low_mapq < a);  // same interval: lower MAPQ first
  MappingWithoutBarcode minus(4, 100, 50, 60, 0, 1, 1);
  CT_CHECK(minus < a);  // same interval and MAPQ: strand 0 first
  MappingWithoutBarcode multi(5, 100, 50, 60, 1, 0, 1);
  CT_CHECK(multi < a);  // then non-unique first
  MappingWithoutBarcode later_id(10, 100, 50, 60, 1, 1, 1);
  CT_CHECK(a < later_id);  // then read id
  CT_CHECK(!(a < a));
  // Duplicates: same start, whatever the strand, length or MAPQ.
  CT_CHECK(a == minus);
  CT_CHECK(a == shorter);
  CT_CHECK(a == low_mapq);
  CT_CHECK(!(a == b));
  CT_CHECK(a.IsSamePosition(minus));
  CT_CHECK(!a.IsSamePosition(b));
  // Getters: [100, 150), strand from direction, no barcode.
  CT_CHECK_EQ(a.GetStartPosition(), 100u);
  CT_CHECK_EQ(a.GetEndPosition(), 150u);
  CT_CHECK(a.IsPositiveStrand());
  CT_CHECK(!minus.IsPositiveStrand());
  CT_CHECK_EQ(a.GetBarcode(), (uint64_t)0);
  CT_CHECK_EQ(a.num_dups_, (uint16_t)1);
  CT_CHECK_EQ((int)a.mapq_, 60);
}

// MappingWithBarcode sorts by start, length, then barcode before MAPQ;
// duplicates need the same barcode and start (cell level); IsSamePosition
// ignores the barcode (bulk level).
CT_TEST(single_end_barcoded_record) {
  MappingWithBarcode a(/*read_id=*/1, /*barcode=*/7, /*start=*/200, /*length=*/50, /*mapq=*/60,
                       /*direction=*/1, /*is_unique=*/1, /*num_dups=*/1);
  MappingWithBarcode other_bc(2, 3, 200, 50, 60, 1, 1, 1);
  MappingWithBarcode same_bc_minus(3, 7, 200, 60, 10, 0, 1, 1);
  MappingWithBarcode moved(4, 7, 201, 50, 60, 1, 1, 1);
  CT_CHECK(other_bc < a);           // barcode 3 before 7
  CT_CHECK(a < same_bc_minus);      // length 50 before 60
  CT_CHECK(a < moved);              // start 200 before 201
  MappingWithBarcode low_mapq(5, 7, 200, 50, 30, 1, 1, 1);
  CT_CHECK(low_mapq < a);           // same start, length, barcode: MAPQ
  CT_CHECK(a == same_bc_minus);     // same barcode and start
  CT_CHECK(!(a == other_bc));       // different barcode
  CT_CHECK(!(a == moved));
  CT_CHECK(a.IsSamePosition(other_bc));
  CT_CHECK(!a.IsSamePosition(moved));
  CT_CHECK_EQ(a.GetBarcode(), (uint64_t)7);
  CT_CHECK_EQ(a.GetEndPosition(), 250u);
  // Tn5: + start 200 -> 204 (end 254); - end 260 -> 255.
  a.Tn5Shift();
  CT_CHECK_EQ(a.GetStartPosition(), 204u);
  CT_CHECK_EQ(a.GetEndPosition(), 254u);
  same_bc_minus.Tn5Shift();
  CT_CHECK_EQ(same_bc_minus.GetStartPosition(), 200u);
  CT_CHECK_EQ(same_bc_minus.GetEndPosition(), 255u);
}

// Fragments: duplicates share start and length; strand and the per-read
// alignment lengths do not matter.
CT_TEST(paired_end_record_order_and_duplicates) {
  PairedEndMappingWithoutBarcode a(/*read_id=*/1, /*start=*/1000, /*length=*/300, /*mapq=*/60,
                                   /*direction=*/1, /*is_unique=*/1, /*num_dups=*/1, 50, 50);
  PairedEndMappingWithoutBarcode flipped(2, 1000, 300, 60, 0, 1, 1, 48, 50);
  PairedEndMappingWithoutBarcode longer(3, 1000, 301, 60, 1, 1, 1, 50, 50);
  PairedEndMappingWithoutBarcode later(4, 1001, 299, 60, 1, 1, 1, 50, 50);
  CT_CHECK(a == flipped);
  CT_CHECK(!(a == longer));
  CT_CHECK(!(a == later));
  CT_CHECK(a.IsSamePosition(flipped));
  CT_CHECK(!a.IsSamePosition(longer));
  CT_CHECK(flipped < a);  // same fragment and MAPQ: direction 0 first
  CT_CHECK(a < longer);
  CT_CHECK(longer < later);
  PairedEndMappingWithoutBarcode low_mapq(5, 1000, 300, 0, 1, 1, 1, 50, 50);
  CT_CHECK(low_mapq < flipped);
  CT_CHECK_EQ(a.GetEndPosition(), 1300u);
  CT_CHECK(!flipped.IsPositiveStrand());
  CT_CHECK_EQ(a.GetBarcode(), (uint64_t)0);
  // Tn5 on a fragment whatever read 1's strand: [1000, 1300) -> [1004,
  // 1295); the + read (50) loses 4 bases at its start, the - read (50) 5 at
  // its end.
  a.Tn5Shift();
  flipped.Tn5Shift();
  CT_CHECK_EQ(a.GetStartPosition(), 1004u);
  CT_CHECK_EQ(a.GetEndPosition(), 1295u);
  CT_CHECK_EQ(a.positive_alignment_length_, (uint16_t)46);
  CT_CHECK_EQ(a.negative_alignment_length_, (uint16_t)45);
  CT_CHECK_EQ(flipped.GetStartPosition(), 1004u);
  CT_CHECK_EQ(flipped.GetEndPosition(), 1295u);
  CT_CHECK_EQ(flipped.positive_alignment_length_, (uint16_t)44);
  CT_CHECK_EQ(flipped.negative_alignment_length_, (uint16_t)45);
}

// Barcoded fragments: cell-level duplicates need barcode, start and length;
// IsSamePosition (bulk level) only start and length.
CT_TEST(paired_end_barcoded_record) {
  PairedEndMappingWithBarcode a(/*read_id=*/1, /*barcode=*/11, /*start=*/500, /*length=*/200, /*mapq=*/60,
                                /*direction=*/0, /*is_unique=*/1, /*num_dups=*/1, 50, 50);
  PairedEndMappingWithBarcode other_bc(2, 12, 500, 200, 60, 0, 1, 1, 50, 50);
  PairedEndMappingWithBarcode other_len(3, 11, 500, 199, 60, 0, 1, 1, 50, 50);
  CT_CHECK(!(a == other_bc));
  CT_CHECK(a.IsSamePosition(other_bc));
  CT_CHECK(!(a == other_len));
  CT_CHECK(!a.IsSamePosition(other_len));
  CT_CHECK(a < other_bc);       // barcode 11 before 12
  CT_CHECK(other_len < a);      // length 199 before 200
  CT_CHECK_EQ(a.GetBarcode(), (uint64_t)11);
  // Tn5 on a fragment whatever read 1's strand: [500, 700) -> [504, 695);
  // the + read loses 4 bases at its start, the - read 5 at its end.
  a.Tn5Shift();
  CT_CHECK_EQ(a.GetStartPosition(), 504u);
  CT_CHECK_EQ(a.GetEndPosition(), 695u);
  CT_CHECK_EQ(a.positive_alignment_length_, (uint16_t)46);
  CT_CHECK_EQ(a.negative_alignment_length_, (uint16_t)45);
}

// Single-end PAF records sort and deduplicate like single-end BED records;
// paired PAF records like fragments.
CT_TEST(paf_records) {
  PAFMapping a(/*read_id=*/3, "q3", /*read_length=*/50, /*start=*/400, /*length=*/50, /*mapq=*/60,
               /*direction=*/0, /*is_unique=*/1, /*num_dups=*/1);
  PAFMapping b(4, "q4", 60, 400, 60, 60, 1, 1, 1);
  PAFMapping c(1, "q1", 50, 401, 50, 60, 1, 1, 1);
  CT_CHECK(a == b);
  CT_CHECK(!(a == c));
  CT_CHECK(a.IsSamePosition(b));
  CT_CHECK(a < b);
  CT_CHECK(b < c);
  CT_CHECK(!a.IsPositiveStrand());
  CT_CHECK_EQ(a.GetEndPosition(), 450u);
  CT_CHECK_EQ(a.GetBarcode(), (uint64_t)0);
  // - strand: end 450 -> 445; + strand: start 400 -> 404, end 460 -> 464.
  a.Tn5Shift();
  CT_CHECK_EQ(a.GetStartPosition(), 400u);
  CT_CHECK_EQ(a.GetEndPosition(), 445u);
  b.Tn5Shift();
  CT_CHECK_EQ(b.GetStartPosition(), 404u);
  CT_CHECK_EQ(b.GetEndPosition(), 464u);

  PairedPAFMapping p(/*read_id=*/1, "p1/1", "p1/2", 50, 50, /*start=*/1000, /*length=*/300, 50, 50,
                     /*mapq=*/60, /*mapq1=*/60, /*mapq2=*/60, /*direction=*/1, /*is_unique=*/1, 1);
  PairedPAFMapping q(2, "p2/1", "p2/2", 50, 50, 1000, 300, 50, 50, 60, 30, 60, 0, 1, 1);
  PairedPAFMapping r(3, "p3/1", "p3/2", 50, 50, 1000, 301, 50, 50, 60, 60, 60, 1, 1, 1);
  CT_CHECK(p == q);
  CT_CHECK(!(p == r));
  CT_CHECK(p.IsSamePosition(q));
  CT_CHECK(q < p);  // same fragment: mapq1 30 before 60
  CT_CHECK(p < r);
  CT_CHECK_EQ(p.GetEndPosition(), 1300u);
  p.Tn5Shift();
  CT_CHECK_EQ(p.GetStartPosition(), 1004u);
  CT_CHECK_EQ(p.GetEndPosition(), 1295u);
  CT_CHECK_EQ(p.positive_alignment_length_, (uint16_t)46);
  CT_CHECK_EQ(p.negative_alignment_length_, (uint16_t)45);
}

// Every record class answers the Mapping interface through a base reference.
CT_TEST(mapping_interface_dispatch) {
  MappingWithoutBarcode se(0, 10, 50, 60, 1, 1, 1);
  MappingWithBarcode seb(0, 5, 20, 50, 60, 0, 1, 1);
  PairedEndMappingWithoutBarcode pe(0, 30, 200, 60, 1, 1, 1, 50, 50);
  PairedEndMappingWithBarcode peb(0, 6, 40, 200, 60, 0, 1, 1, 50, 50);
  PAFMapping paf(0, "r", 50, 50, 50, 60, 1, 1, 1);
  PairedPAFMapping ppaf(0, "a", "b", 50, 50, 60, 200, 50, 50, 60, 60, 60, 0, 1, 1);
  SAMMapping sam = MakeSam(0, "s", 9, /*pos=*/70, 0, 0, -1, 0, 0, /*positive=*/1, 60, "M", {50},
                           Repeat("ACGT", 50), Repeat("I", 50));
  Mapping *records[] = {&se, &seb, &pe, &peb, &paf, &ppaf, &sam};
  const uint32_t starts[] = {10, 20, 30, 40, 50, 60, 71};
  const uint32_t ends[] = {60, 70, 230, 240, 100, 260, 120};
  const bool positive[] = {true, false, true, false, true, false, true};
  const uint64_t barcodes[] = {0, 5, 0, 6, 0, 0, 9};
  for (int i = 0; i < 7; ++i) {
    Mapping &m = *records[i];
    CT_CHECK_EQ(m.GetStartPosition(), starts[i]);
    CT_CHECK_EQ(m.GetEndPosition(), ends[i]);
    CT_CHECK_EQ(m.IsPositiveStrand(), positive[i]);
    CT_CHECK_EQ(m.GetBarcode(), barcodes[i]);
  }
  // Tn5 through the interface: se + start 10 -> 14; peb fragment [40, 240) ->
  // [44, 235); SAM records are not shifted (documented: --Tn5-shift does not
  // apply to --SAM).
  se.Tn5Shift();
  static_cast<Mapping &>(peb).Tn5Shift();
  static_cast<Mapping &>(sam).Tn5Shift();
  CT_CHECK_EQ(se.GetStartPosition(), 14u);
  CT_CHECK_EQ(peb.GetStartPosition(), 44u);
  CT_CHECK_EQ(peb.GetEndPosition(), 235u);
  CT_CHECK_EQ(sam.GetStartPosition(), 71u);
  CT_CHECK_EQ(sam.GetEndPosition(), 120u);
}

// CIGAR strings and lengths (SAM specification): M, I, S, =, X consume the
// query; M, D, N, =, X consume the reference.
CT_TEST(sam_cigar_strings_and_lengths) {
  const std::string seq50 = Repeat("ACGTTGCA", 50);
  const std::string qual50 = Repeat("ABCDEFGHIJ", 50);
  SAMMapping ins = MakeSam(0, "i", 0, 99, 0, 0, -1, 0, 0, 1, 60, "MIM", {20, 1, 29}, seq50, qual50);
  CT_CHECK_EQ(ins.GenerateCigarString(), std::string("20M1I29M"));
  CT_CHECK_EQ(ins.GetSequenceLength(), 50u);  // 20 + 1 + 29
  CT_CHECK_EQ(ins.GetAlignmentLength(), 49u);  // 20 + 29
  // POS is 1-based: pos_ 99 -> 100; end (1-based inclusive) 99 + 49 = 148.
  CT_CHECK_EQ(ins.GetStartPosition(), 100u);
  CT_CHECK_EQ(ins.GetEndPosition(), 148u);
  SAMMapping del = MakeSam(0, "d", 0, 0, 0, 0, -1, 0, 0, 1, 60, "MDM", {20, 1, 30}, seq50, qual50);
  CT_CHECK_EQ(del.GenerateCigarString(), std::string("20M1D30M"));
  CT_CHECK_EQ(del.GetSequenceLength(), 50u);
  CT_CHECK_EQ(del.GetAlignmentLength(), 51u);
  // Every operator once: query 1+2+5+8+9 = 25, reference 1+3+4+8+9 = 25.
  const std::string seq25 = Repeat("ACGT", 25);
  SAMMapping all = MakeSam(0, "a", 0, 0, 0, 0, -1, 0, 0, 1, 60, "MIDNSHP=X", {1, 2, 3, 4, 5, 6, 7, 8, 9},
                           seq25, Repeat("I", 25));
  CT_CHECK_EQ(all.GenerateCigarString(), std::string("1M2I3D4N5S6H7P8=9X"));
  CT_CHECK_EQ(all.GetSequenceLength(), 25u);
  CT_CHECK_EQ(all.GetAlignmentLength(), 25u);
  // Soft clip: 5S45M consumes 50 query and 45 reference bases.
  SAMMapping clip = MakeSam(0, "c", 0, 0, 0, 0, -1, 0, 0, 1, 60, "SM", {5, 45}, seq50, qual50);
  CT_CHECK_EQ(clip.GetSequenceLength(), 50u);
  CT_CHECK_EQ(clip.GetAlignmentLength(), 45u);
  // No CIGAR: "*".
  SAMMapping none = MakeSam(0, "n", 0, 0, 0, 0, -1, 0, 4, 1, 0, "", {}, "", "");
  CT_CHECK_EQ(none.GenerateCigarString(), std::string("*"));
  CT_CHECK_EQ(none.GenerateIntTagString("NM", 3), std::string("NM:i:3"));
  CT_CHECK_EQ(none.GenerateIntTagString("XS", -2), std::string("XS:i:-2"));
}

// The constructor reverses QUAL for a reverse-strand record (SEQ is passed
// already in reference orientation) and truncates SEQ and QUAL to the
// query length the CIGAR implies.
CT_TEST(sam_constructor_orients_and_truncates) {
  const std::string seq = "ACGTACGTAC";
  const std::string qual = "0123456789";
  SAMMapping fwd = MakeSam(0, "f", 0, 0, 0, 0, -1, 0, 0, /*positive=*/1, 60, "M", {10}, seq, qual);
  CT_CHECK_EQ(fwd.sequence_, seq);
  CT_CHECK_EQ(fwd.sequence_qual_, qual);
  SAMMapping rev = MakeSam(0, "r", 0, 0, 0, 0, -1, 0, 16, /*positive=*/0, 60, "M", {10}, seq, qual);
  CT_CHECK_EQ(rev.sequence_, seq);
  CT_CHECK_EQ(rev.sequence_qual_, std::string("9876543210"));
  CT_CHECK(!rev.IsPositiveStrand());
  // CIGAR 7M on a 10-base read keeps the first 7 bases and qualities (after
  // the reversal for a reverse-strand record).
  SAMMapping cut = MakeSam(0, "t", 0, 0, 0, 0, -1, 0, 0, 1, 60, "M", {7}, seq, qual);
  CT_CHECK_EQ(cut.sequence_, std::string("ACGTACG"));
  CT_CHECK_EQ(cut.sequence_qual_, std::string("0123456"));
  SAMMapping cut_rev = MakeSam(0, "u", 0, 0, 0, 0, -1, 0, 16, 0, 60, "M", {7}, seq, qual);
  CT_CHECK_EQ(cut_rev.sequence_, std::string("ACGTACG"));
  CT_CHECK_EQ(cut_rev.sequence_qual_, std::string("9876543"));
}

// SAM records sort by reference, position, barcode, mate reference, mate
// position, read-1 flag, MAPQ, read id; duplicates share position,
// reference, barcode, read-1 flag and the mate's reference and position.
CT_TEST(sam_record_order_and_duplicates) {
  const std::string s = Repeat("ACGT", 50), q = Repeat("I", 50);
  // A read-1 record of a pair on reference 1 at pos 100, mate at 300.
  SAMMapping a = MakeSam(1, "a", 0, 100, 1, 300, 1, 250, 99, 1, 60, "M", {50}, s, q);
  SAMMapping a_dup = MakeSam(2, "b", 0, 100, 1, 300, 1, 250, 99, 1, 10, "M", {50}, s, q);
  SAMMapping a_read2 = MakeSam(3, "c", 0, 100, 1, 300, 1, 250, 163, 1, 60, "M", {50}, s, q);
  SAMMapping other_mate = MakeSam(4, "d", 0, 100, 1, 301, 1, 251, 99, 1, 60, "M", {50}, s, q);
  SAMMapping other_bc = MakeSam(5, "e", 3, 100, 1, 300, 1, 250, 99, 1, 60, "M", {50}, s, q);
  SAMMapping on_ref0 = MakeSam(6, "f", 0, 900, 0, 950, 0, 100, 99, 1, 60, "M", {50}, s, q);
  CT_CHECK(a == a_dup);  // MAPQ and read id do not matter
  CT_CHECK(!(a == a_read2));
  CT_CHECK(!(a == other_mate));
  CT_CHECK(!(a == other_bc));
  CT_CHECK(on_ref0 < a);       // reference first, whatever the position
  CT_CHECK(a_dup < a);         // MAPQ 10 before 60
  CT_CHECK(a_read2 < a);       // read-2 record (flag bit 64 clear) first
  CT_CHECK(a < other_mate);    // mate position 300 before 301
  CT_CHECK(a < other_bc);      // barcode 0 before 3
  // Two paired-end records at the same place, mate on the same reference:
  // same position (bulk level), whatever the barcode.
  CT_CHECK(a.IsSamePosition(other_bc));
  CT_CHECK(!a.IsSamePosition(other_mate));
  CT_CHECK_EQ(other_bc.GetBarcode(), (uint64_t)3);
}

// The BAM flag bits and CIGAR codes of the SAM specification.
CT_TEST(sam_flag_bits_and_cigar_table) {
  CT_CHECK_EQ(BAM_FPAIRED, 1);
  CT_CHECK_EQ(BAM_FPROPER_PAIR, 2);
  CT_CHECK_EQ(BAM_FUNMAP, 4);
  CT_CHECK_EQ(BAM_FMUNMAP, 8);
  CT_CHECK_EQ(BAM_FREVERSE, 16);
  CT_CHECK_EQ(BAM_FMREVERSE, 32);
  CT_CHECK_EQ(BAM_FREAD1, 64);
  CT_CHECK_EQ(BAM_FREAD2, 128);
  CT_CHECK_EQ(BAM_FSECONDARY, 256);
  CT_CHECK_EQ(BAM_FQCFAIL, 512);
  CT_CHECK_EQ(BAM_FDUP, 1024);
  CT_CHECK_EQ(BAM_FSUPPLEMENTARY, 2048);
  const std::string ops = "MIDNSHP=XB";
  for (int code = 0; code < 10; ++code) {
    CT_CHECK_EQ((int)bam_cigar_table[(unsigned char)ops[code]], code);
    CT_CHECK_EQ(bam_cigar_opchr(code), ops[code]);
  }
  int unknown = 0;
  for (int c = 0; c < 256; ++c) {
    if (ops.find((char)c) == std::string::npos && c != 0 && bam_cigar_table[c] != -1) ++unknown;
  }
  CT_CHECK_EQ(unknown, 0);
  CT_CHECK_EQ((int)bam_cigar_table[(unsigned char)'m'], -1);
  // 123M packs as 123 << 4 | 0; 7D as 7 << 4 | 2.
  CT_CHECK_EQ((uint32_t)bam_cigar_gen(123u, BAM_CMATCH), (uint32_t)(123u << 4));
  CT_CHECK_EQ((uint32_t)bam_cigar_op(bam_cigar_gen(7u, BAM_CDEL)), (uint32_t)BAM_CDEL);
  CT_CHECK_EQ((uint32_t)bam_cigar_oplen(bam_cigar_gen(7u, BAM_CDEL)), 7u);
  // Query (bit 1) and reference (bit 2) consumption per operator.
  const int consumes[10] = {3, 1, 2, 2, 1, 0, 0, 3, 3, 0};
  for (int code = 0; code < 10; ++code) CT_CHECK_EQ((int)bam_cigar_type(code), consumes[code]);
}

// PAF records survive the temp-file encoding with every field.
CT_TEST(paf_write_load_round_trip) {
  FILE *f = tmpfile();
  CT_CHECK(f != NULL);
  if (f == NULL) return;
  PAFMapping a(/*read_id=*/123456, "read/1 with a long name", 151, 4000000000u, 151, 59, 0, 1, 7);
  PAFMapping b(7, "", 30, 0, 29, 0, 1, 0, 255);
  a.WriteToFile(f);
  b.WriteToFile(f);
  rewind(f);
  PAFMapping x, y;
  x.LoadFromFile(f);
  y.LoadFromFile(f);
  fclose(f);
  CT_CHECK_EQ(x.read_id_, 123456u);
  CT_CHECK_EQ(x.read_name_, std::string("read/1 with a long name"));
  CT_CHECK_EQ(x.read_length_, (uint16_t)151);
  CT_CHECK_EQ(x.fragment_start_position_, 4000000000u);
  CT_CHECK_EQ(x.fragment_length_, (uint16_t)151);
  CT_CHECK_EQ((int)x.mapq_, 59);
  CT_CHECK_EQ((int)x.direction_, 0);
  CT_CHECK_EQ((int)x.is_unique_, 1);
  CT_CHECK_EQ((int)x.num_dups_, 7);
  CT_CHECK_EQ(y.read_id_, 7u);
  CT_CHECK_EQ(y.read_name_, std::string(""));
  CT_CHECK_EQ(y.fragment_start_position_, 0u);
  CT_CHECK_EQ(y.fragment_length_, (uint16_t)29);
  CT_CHECK_EQ((int)y.mapq_, 0);
  CT_CHECK_EQ((int)y.direction_, 1);
  CT_CHECK_EQ((int)y.is_unique_, 0);
  CT_CHECK_EQ((int)y.num_dups_, 255);
}

CT_TEST(paired_paf_write_load_round_trip) {
  FILE *f = tmpfile();
  CT_CHECK(f != NULL);
  if (f == NULL) return;
  PairedPAFMapping a(/*read_id=*/99, "p/1", "p/2", 50, 76, /*start=*/123456, /*length=*/420, 50, 76,
                     /*mapq=*/42, /*mapq1=*/63, /*mapq2=*/17, /*direction=*/0, /*is_unique=*/1, 3);
  a.WriteToFile(f);
  rewind(f);
  PairedPAFMapping x;
  x.LoadFromFile(f);
  fclose(f);
  CT_CHECK_EQ(x.read_id_, 99u);
  CT_CHECK_EQ(x.read1_name_, std::string("p/1"));
  CT_CHECK_EQ(x.read2_name_, std::string("p/2"));
  CT_CHECK_EQ(x.read1_length_, (uint16_t)50);
  CT_CHECK_EQ(x.read2_length_, (uint16_t)76);
  CT_CHECK_EQ(x.fragment_start_position_, 123456u);
  CT_CHECK_EQ(x.fragment_length_, (uint16_t)420);
  CT_CHECK_EQ(x.positive_alignment_length_, (uint16_t)50);
  CT_CHECK_EQ(x.negative_alignment_length_, (uint16_t)76);
  CT_CHECK_EQ((int)x.mapq_, 42);
  CT_CHECK_EQ((int)x.mapq1_, 63);
  CT_CHECK_EQ((int)x.mapq2_, 17);
  CT_CHECK_EQ((int)x.direction_, 0);
  CT_CHECK_EQ((int)x.is_unique_, 1);
  CT_CHECK_EQ((int)x.num_dups_, 3);
}

CT_TEST(sam_write_load_round_trip) {
  FILE *f = tmpfile();
  CT_CHECK(f != NULL);
  if (f == NULL) return;
  const std::string seq = Repeat("ACGTTGCAAC", 51);
  const std::string qual = Repeat("!#%')+-/13579", 51);
  SAMMapping a = MakeSam(/*read_id=*/77, "pair17/2", /*barcode=*/0x123456789ull, /*pos=*/3999999, /*rid=*/3,
                         /*mpos=*/3999700, /*mrid=*/3, /*tlen=*/-350, /*flag=*/147, /*positive=*/0,
                         /*mapq=*/57, "MIM", {20, 1, 30}, seq, qual, /*nm=*/2, "20^A3C26");
  SAMMapping b = MakeSam(5, "x", 0, 0, 0, 0, -1, 0, 0, 1, 0, "", {}, "", "", 0, "");
  a.WriteToFile(f);
  b.WriteToFile(f);
  rewind(f);
  SAMMapping x, y;
  x.LoadFromFile(f);
  y.LoadFromFile(f);
  fclose(f);
  CT_CHECK_EQ(x.read_id_, 77u);
  CT_CHECK_EQ(x.read_name_, std::string("pair17/2"));
  CT_CHECK_EQ(x.cell_barcode_, (uint64_t)0x123456789ull);
  CT_CHECK_EQ((int)x.num_dups_, 1);
  CT_CHECK_EQ(x.pos_, (int64_t)3999999);
  CT_CHECK_EQ(x.rid_, 3);
  CT_CHECK_EQ(x.mpos_, (int64_t)3999700);
  CT_CHECK_EQ(x.mrid_, 3);
  CT_CHECK_EQ(x.tlen_, -350);
  CT_CHECK_EQ(x.flag_, 147);
  CT_CHECK(!x.IsPositiveStrand());
  CT_CHECK_EQ((int)x.is_unique_, 1);
  CT_CHECK_EQ((int)x.is_alt_, 0);
  CT_CHECK_EQ((int)x.mapq_, 57);
  CT_CHECK_EQ((int)x.NM_, 2);
  CT_CHECK_EQ(x.GenerateCigarString(), std::string("20M1I30M"));
  CT_CHECK_EQ(x.MD_, std::string("20^A3C26"));
  CT_CHECK_EQ(x.sequence_, seq);
  std::string reversed(qual.rbegin(), qual.rend());
  CT_CHECK_EQ(x.sequence_qual_, reversed);
  CT_CHECK_EQ(y.read_name_, std::string("x"));
  CT_CHECK_EQ(y.mrid_, -1);
  CT_CHECK_EQ(y.GenerateCigarString(), std::string("*"));
  CT_CHECK_EQ(y.sequence_, std::string(""));
  CT_CHECK(y.IsPositiveStrand());
}

CT_MAIN
