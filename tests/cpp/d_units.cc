// Area D C++ unit-test driver: pairs records, the barcode-correction
// candidate order, the read-format ranges, and the feature/barcode matrix
// classes (not reachable from the command line at 949043c: the peak options
// are never registered).
//
// CT_TEST cases check values worked by hand in their comments; CT_PROBE
// functions print what chromap computes so tests/test_d_unit.py can compare
// it with expectations computed in Python.

#include "chromap_test.h"

#include <cstdio>

#include "bed_mapping.h"
#include "feature_barcode_matrix.h"
#include "feature_barcode_matrix_writer.h"
#include "pairs_mapping.h"
#include "sequence_batch.h"
#include "sequence_effective_range.h"
#include "utils.h"

using namespace chromap;

// The pairs columns are 1-based: pos1_ = 999 (0-based) is written as 1000.
// strand 1 is '+', 0 is '-'. Start and end are pos1_ and pos2_ as stored.
CT_TEST(pairs_mapping_positions_and_strands) {
  PairsMapping m(/*read_id=*/7, "read7", /*cell_barcode=*/12345, /*rid1=*/1,
                 /*rid2=*/2, /*pos1=*/999, /*pos2=*/0, /*strand1=*/1,
                 /*strand2=*/0, /*mapq=*/60, /*is_unique=*/1, /*num_dups=*/1);
  CT_CHECK_EQ(m.GetPosition(1), 1000);
  CT_CHECK_EQ(m.GetPosition(2), 1);
  CT_CHECK_EQ(m.GetStrand(1), '+');
  CT_CHECK_EQ(m.GetStrand(2), '-');
  CT_CHECK(m.IsPositiveStrand());
  CT_CHECK_EQ(m.GetStartPosition(), 999u);
  CT_CHECK_EQ(m.GetEndPosition(), 0u);
  PairsMapping n(7, "read7", 0, 1, 2, 5, 6, /*strand1=*/0, /*strand2=*/1, 60, 1, 1);
  CT_CHECK(!n.IsPositiveStrand());
  CT_CHECK_EQ(n.GetStrand(1), '-');
  CT_CHECK_EQ(n.GetStrand(2), '+');
}

// Pins the current behaviour: a pairs record reports barcode key 0 whatever
// barcode it carries, and Tn5Shift leaves it unchanged (the class comment:
// no Tn5 shift for formats whose fields depend on the mapping position).
CT_TEST(pairs_mapping_barcode_and_tn5_shift) {
  PairsMapping m(3, "r", /*cell_barcode=*/987654321, 0, 1, 100, 200, 1, 0, 30, 1, 1);
  CT_CHECK_EQ(m.GetBarcode(), (uint64_t)0);
  m.Tn5Shift();
  CT_CHECK_EQ(m.pos1_, 100u);
  CT_CHECK_EQ(m.pos2_, 200u);
  CT_CHECK_EQ(m.GetPosition(1), 101);
  CT_CHECK_EQ(m.GetPosition(2), 201);
}

// operator< orders by rid1, rid2, pos1, pos2, mapq, then read id; strands
// and names take no part.
CT_TEST(pairs_mapping_order) {
  // rid1 decides before everything else.
  PairsMapping a(9, "a", 0, /*rid1=*/0, /*rid2=*/3, 500, 900, 1, 1, 60, 1, 1);
  PairsMapping b(1, "b", 0, /*rid1=*/1, /*rid2=*/0, 0, 0, 1, 1, 0, 1, 1);
  CT_CHECK(a < b);
  CT_CHECK(!(b < a));
  // Same rid1: rid2 before positions.
  PairsMapping c(1, "c", 0, 0, /*rid2=*/1, /*pos1=*/800, 900, 1, 1, 60, 1, 1);
  PairsMapping d(1, "d", 0, 0, /*rid2=*/2, /*pos1=*/10, 20, 1, 1, 60, 1, 1);
  CT_CHECK(c < d);
  // Same rids: pos1, then pos2.
  PairsMapping e(5, "e", 0, 0, 0, /*pos1=*/10, /*pos2=*/90, 1, 1, 60, 1, 1);
  PairsMapping f(4, "f", 0, 0, 0, /*pos1=*/11, /*pos2=*/20, 1, 1, 60, 1, 1);
  PairsMapping g(4, "g", 0, 0, 0, /*pos1=*/10, /*pos2=*/91, 1, 1, 60, 1, 1);
  CT_CHECK(e < f);
  CT_CHECK(e < g);
  // Same positions: lower MAPQ first, then lower read id.
  PairsMapping h(8, "h", 0, 0, 0, 10, 90, 1, 1, /*mapq=*/30, 1, 1);
  PairsMapping i(2, "i", 0, 0, 0, 10, 90, 1, 1, /*mapq=*/31, 1, 1);
  PairsMapping j(3, "j", 0, 0, 0, 10, 90, 1, 1, /*mapq=*/31, 1, 1);
  CT_CHECK(h < i);
  CT_CHECK(i < j);
  CT_CHECK(!(j < i));
  // Strands do not order records.
  PairsMapping k(3, "k", 0, 0, 0, 10, 90, /*strand1=*/0, /*strand2=*/0, 31, 1, 1);
  CT_CHECK(!(j < k));
  CT_CHECK(!(k < j));
}

// Duplicates are records with the same rid1, pos1, rid2 and pos2; strands,
// MAPQ, read id and name do not matter.
CT_TEST(pairs_mapping_duplicate_keys) {
  PairsMapping a(1, "a", 0, 0, 2, 10, 90, 1, 0, 60, 1, 1);
  PairsMapping b(2, "b", 77, 0, 2, 10, 90, 0, 1, 5, 0, 3);
  CT_CHECK(a == b);
  CT_CHECK(a.IsSamePosition(b));
  PairsMapping c(1, "a", 0, 0, 2, 10, /*pos2=*/91, 1, 0, 60, 1, 1);
  PairsMapping d(1, "a", 0, 0, /*rid2=*/1, 10, 90, 1, 0, 60, 1, 1);
  PairsMapping e(1, "a", 0, /*rid1=*/1, 2, 10, 90, 1, 0, 60, 1, 1);
  PairsMapping f(1, "a", 0, 0, 2, /*pos1=*/11, 90, 1, 0, 60, 1, 1);
  CT_CHECK(!(a == c));
  CT_CHECK(!(a == d));
  CT_CHECK(!(a == e));
  CT_CHECK(!(a == f));
  CT_CHECK(!a.IsSamePosition(c));
  CT_CHECK(!a.IsSamePosition(d));
  CT_CHECK(!a.IsSamePosition(e));
  CT_CHECK(!a.IsSamePosition(f));
}

static void CheckSamePairsRecord(const PairsMapping &x, const PairsMapping &y) {
  CT_CHECK_EQ(x.read_id_, y.read_id_);
  CT_CHECK_EQ(x.read_name_, y.read_name_);
  CT_CHECK_EQ(x.cell_barcode_, y.cell_barcode_);
  CT_CHECK_EQ(x.rid1_, y.rid1_);
  CT_CHECK_EQ(x.rid2_, y.rid2_);
  CT_CHECK_EQ(x.pos1_, y.pos1_);
  CT_CHECK_EQ(x.pos2_, y.pos2_);
  CT_CHECK_EQ(x.strand1_, y.strand1_);
  CT_CHECK_EQ(x.strand2_, y.strand2_);
  CT_CHECK_EQ((int)x.mapq_, (int)y.mapq_);
  CT_CHECK_EQ((int)x.is_unique_, (int)y.is_unique_);
  CT_CHECK_EQ((int)x.num_dups_, (int)y.num_dups_);
}

// WriteToFile then LoadFromFile restores every field, including the MAPQ (8
// bits), uniqueness (1 bit) and duplicate count (7 bits) packed together,
// at their extremes; records written back to back load back in order.
CT_TEST(pairs_mapping_write_load_round_trip) {
  std::vector<PairsMapping> records;
  records.push_back(PairsMapping(0, "p000001", 0, 0, 0, 10, 300, 1, 0, 60, 1, 1));
  records.push_back(PairsMapping(4294967295u, "x", 0xffffffffffffffffull, 3, 1,
                                 4000000000u, 7, 0, 1, 255, 0, 127));
  records.push_back(PairsMapping(42, "a_much_longer_read_name/1", 5, 2, 3, 0, 0, 0, 0,
                                 0, 1, 0));
  records.push_back(PairsMapping(43, "", 0, 1, 1, 123, 456, 1, 1, 1, 0, 64));
  FILE *f = tmpfile();
  CT_CHECK(f != NULL);
  if (f == NULL) return;
  for (size_t i = 0; i < records.size(); ++i) records[i].WriteToFile(f);
  rewind(f);
  for (size_t i = 0; i < records.size(); ++i) {
    PairsMapping loaded;
    loaded.LoadFromFile(f);
    CheckSamePairsRecord(loaded, records[i]);
  }
  // Nothing is left after the last record.
  CT_CHECK_EQ(fgetc(f), EOF);
  fclose(f);
}

// GetByteSize grows by one per read-name character.
CT_TEST(pairs_mapping_byte_size_counts_the_name) {
  PairsMapping a(1, "ab", 0, 0, 0, 1, 2, 1, 1, 60, 1, 1);
  PairsMapping b(1, "abcdefgh", 0, 0, 0, 1, 2, 1, 1, 60, 1, 1);
  CT_CHECK_EQ((int)b.GetByteSize() - (int)a.GetByteSize(), 6);
}

// BarcodeWithQual orders by score, then the first corrected position, its
// base, the second position and its base (std::tie).
CT_TEST(barcode_with_qual_order) {
  BarcodeWithQual a = {5, 'A', 0, 0, 0.5};
  BarcodeWithQual b = {9, 'T', 3, 'G', 0.25};
  BarcodeWithQual c = {7, 'A', 0, 0, 0.5};
  BarcodeWithQual d = {5, 'G', 0, 0, 0.5};
  BarcodeWithQual e = {5, 'A', 2, 'C', 0.5};
  BarcodeWithQual f = {5, 'A', 2, 'T', 0.5};
  // A higher score wins whatever the positions.
  CT_CHECK(a > b);
  CT_CHECK(!(b > a));
  // Equal scores: the larger first position, then the larger base.
  CT_CHECK(c > a);
  CT_CHECK(d > a);
  CT_CHECK(c > d);
  // Then the second position and its base.
  CT_CHECK(e > a);
  CT_CHECK(f > e);
  CT_CHECK(!(a > a));
  // Sorting with std::greater puts the best candidate first.
  std::vector<BarcodeWithQual> v;
  v.push_back(b);
  v.push_back(a);
  v.push_back(f);
  v.push_back(c);
  v.push_back(d);
  v.push_back(e);
  std::sort(v.begin(), v.end(), std::greater<BarcodeWithQual>());
  CT_CHECK_EQ(v[0].corrected_base_index1, 7u);
  CT_CHECK_EQ(v[1].correct_base1, 'G');
  CT_CHECK_EQ(v[2].correct_base2, 'T');
  CT_CHECK_EQ(v[3].correct_base2, 'C');
  CT_CHECK_EQ(v[4].corrected_base_index2, 0u);
  CT_CHECK_EQ(v[5].score, 0.25);
}

// effective_range FASTQ FIELD...: the fields ("r1:start:end[:strand]", one
// read's ranges) parsed into a SequenceEffectiveRange; prints "rejected" if
// one is refused, otherwise every record of FASTQ loaded through a
// SequenceBatch with that range as "SEQ<TAB>QUAL".
CT_PROBE(effective_range) {
  if (args.size() < 1) {
    std::cerr << "usage: effective_range FASTQ FIELD..." << std::endl;
    return 2;
  }
  SequenceEffectiveRange range;
  range.InitializeParsing();
  for (size_t i = 1; i < args.size(); ++i) {
    if (!range.ParseFormatStringAndAppendEffectiveRange(args[i].c_str(),
                                                        args[i].size())) {
      std::cout << "rejected\n";
      return 0;
    }
  }
  range.FinalizeParsing();
  SequenceBatch batch(1000, range);
  batch.InitializeLoading(args[0]);
  uint32_t n = batch.LoadBatch();
  for (uint32_t i = 0; i < n; ++i) {
    std::cout << std::string(batch.GetSequenceAt(i), batch.GetSequenceLengthAt(i))
              << "\t"
              << std::string(batch.GetSequenceQualAt(i), batch.GetSequenceLengthAt(i))
              << "\n";
  }
  batch.FinalizeLoading();
  return 0;
}

// feature_matrix FASTA PREFIX CELL_BY_BIN BIN_SIZE OVERLAP DEPTH: fragments
// "rid start length barcode" on stdin become PairedEndMappingWithBarcode
// records on their reference (in input order), and
// FeatureBarcodeMatrix::OutputFeatureMatrix writes PREFIX_peaks.bed,
// PREFIX_barcode.tsv and PREFIX_matrix.mtx.
CT_PROBE(feature_matrix) {
  if (args.size() != 6) {
    std::cerr << "usage: feature_matrix FASTA PREFIX CELL_BY_BIN BIN_SIZE OVERLAP DEPTH"
              << std::endl;
    return 2;
  }
  SequenceBatch reference;
  reference.InitializeLoading(args[0]);
  reference.LoadAllSequences();
  const uint32_t num_sequences = reference.GetNumSequences();
  std::vector<std::vector<PairedEndMappingWithBarcode> > mappings(num_sequences);
  uint32_t rid, start, length;
  std::string barcode;
  uint32_t read_id = 0;
  while (std::cin >> rid >> start >> length >> barcode) {
    uint64_t key = GenerateSeedFromSequence(barcode.c_str(), barcode.size(), 0,
                                            barcode.size());
    mappings[rid].push_back(PairedEndMappingWithBarcode(
        read_id++, key, start, (uint16_t)length, 60, 1, 1, 1, 50, 50));
  }
  FeatureBarcodeMatrix matrix(args[2] == "1", std::stoi(args[3]), std::stoi(args[4]),
                              (uint16_t)std::stoi(args[5]));
  matrix.OutputFeatureMatrix(num_sequences, reference, mappings, args[1]);
  reference.FinalizeLoading();
  return 0;
}

// barcode_output PREFIX LENGTH TABLE BARCODE...: a FeatureBarcodeMatrixWriter
// with barcode length LENGTH and, unless TABLE is "-", the translation table
// TABLE writes each barcode (as its 2-bit key) to PREFIX_barcode.tsv.
CT_PROBE(barcode_output) {
  if (args.size() < 3) {
    std::cerr << "usage: barcode_output PREFIX LENGTH TABLE BARCODE..." << std::endl;
    return 2;
  }
  FeatureBarcodeMatrixWriter writer;
  writer.InitializeMatrixOutput(args[0]);
  if (args[2] != "-") writer.SetBarcodeTranslateTable(args[2]);
  writer.SetBarcodeLength(std::stoi(args[1]));
  for (size_t i = 3; i < args.size(); ++i) {
    writer.AppendBarcodeOutput(GenerateSeedFromSequence(
        args[i].c_str(), args[i].size(), 0, args[i].size()));
  }
  writer.FinalizeMatrixOutput();
  return 0;
}

CT_MAIN
