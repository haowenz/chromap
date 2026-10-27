// Pilot C++ unit-test driver: the pattern for tests/cpp/<area>_*.cc.
//
// CT_TEST cases check values worked by hand in their comments. CT_PROBE
// functions expose a chromap function to pytest, which computes the expected
// values independently (tests/test_pilot_unit.py).

#include "chromap_test.h"

#include "bed_mapping.h"
#include "hit_utils.h"
#include "minimizer.h"
#include "minimizer_generator.h"
#include "sequence_batch.h"
#include "utils.h"

using namespace chromap;

// A=0, C=1, G=2, T=3, first base in the highest bits:
// ACGT = 00 01 10 11 = 27; TTTT = 255.
CT_TEST(seed_encoding_hand_worked) {
  CT_CHECK_EQ(GenerateSeedFromSequence("ACGT", 4, 0, 4), (uint64_t)27);
  CT_CHECK_EQ(GenerateSeedFromSequence("TTTT", 4, 0, 4), (uint64_t)255);
  // Lower case is the same base.
  CT_CHECK_EQ(GenerateSeedFromSequence("acgt", 4, 0, 4), (uint64_t)27);
  // N encodes as A: ANGT = 00 00 10 11 = 11.
  CT_CHECK_EQ(GenerateSeedFromSequence("ANGT", 4, 0, 4), (uint64_t)11);
  // Past the end pads with A: AC + AA = 00 01 00 00 = 16.
  CT_CHECK_EQ(GenerateSeedFromSequence("AC", 2, 0, 4), (uint64_t)16);
  // A start offset: GT of ACGT = 10 11 = 11.
  CT_CHECK_EQ(GenerateSeedFromSequence("ACGT", 4, 2, 2), (uint64_t)11);
}

// hit = ((sequence index << 32 | position) << 1) | strand (src/minimizer.h).
CT_TEST(hit_encoding_round_trip) {
  const uint64_t hit = GenerateMinimizer(3, 123456, kNegative);
  CT_CHECK_EQ(hit, ((((uint64_t)3) << 32 | 123456u) << 1) | 1u);
  CT_CHECK_EQ(HitToSequenceIndex(hit), 3u);
  CT_CHECK_EQ(HitToSequencePosition(hit), 123456u);
  CT_CHECK(HitToStrand(hit) == kNegative);
  const uint64_t hit2 = GenerateMinimizer(0, 7, kPositive);
  CT_CHECK_EQ(hit2, (uint64_t)14);
  CT_CHECK(HitToStrand(hit2) == kPositive);
  CT_CHECK(!AreTwoHitsOnTheSameStrand(hit, hit2));
  Minimizer m(42, hit);
  CT_CHECK_EQ(m.GetHash(), (uint64_t)42);
  CT_CHECK_EQ(m.GetSequenceIndex(), 3u);
  CT_CHECK_EQ(m.GetSequencePosition(), 123456u);
}

// Tn5 shift on a single-end record: + strand moves the start right by 4,
// - strand pulls the end left by 5 (BED end = start + length).
CT_TEST(tn5_shift_single_end_record) {
  MappingWithoutBarcode plus(/*read_id=*/0, /*start=*/100, /*length=*/50, /*mapq=*/60,
                             /*direction=*/1, /*is_unique=*/1, /*num_dups=*/1);
  plus.Tn5Shift();
  CT_CHECK_EQ(plus.GetStartPosition(), 104u);
  CT_CHECK_EQ(plus.GetEndPosition(), 154u);
  MappingWithoutBarcode minus(0, 100, 50, 60, /*direction=*/0, 1, 1);
  minus.Tn5Shift();
  CT_CHECK_EQ(minus.GetStartPosition(), 100u);
  CT_CHECK_EQ(minus.GetEndPosition(), 145u);
}

// Tn5 shift on a fragment: start + 4, end - 5, so the length drops by 9.
CT_TEST(tn5_shift_paired_end_record) {
  PairedEndMappingWithoutBarcode m(/*read_id=*/0, /*start=*/1000, /*length=*/300, /*mapq=*/60,
                                   /*direction=*/1, /*is_unique=*/1, /*num_dups=*/1,
                                   /*positive_alignment_length=*/50,
                                   /*negative_alignment_length=*/50);
  m.Tn5Shift();
  CT_CHECK_EQ(m.GetStartPosition(), 1004u);
  CT_CHECK_EQ(m.GetEndPosition(), 1295u);
  CT_CHECK_EQ(m.positive_alignment_length_, (uint16_t)46);
  CT_CHECK_EQ(m.negative_alignment_length_, (uint16_t)45);
}

// minimizers FASTA K W: every sequence of FASTA loaded with SequenceBatch,
// minimizers from MinimizerGenerator(K, W), one "index<TAB>hash<TAB>hit" line
// each, in the generator's order.
CT_PROBE(minimizers) {
  if (args.size() != 3) {
    std::cerr << "usage: minimizers FASTA K W" << std::endl;
    return 2;
  }
  SequenceBatch batch;
  batch.InitializeLoading(args[0]);
  batch.LoadAllSequences();
  MinimizerGenerator generator(std::stoi(args[1]), std::stoi(args[2]));
  for (uint32_t i = 0; i < batch.GetNumSequences(); ++i) {
    std::vector<Minimizer> minimizers;
    generator.GenerateMinimizers(batch, i, minimizers);
    for (size_t j = 0; j < minimizers.size(); ++j) {
      std::cout << i << "\t" << minimizers[j].GetHash() << "\t" << minimizers[j].GetHit() << "\n";
    }
  }
  batch.FinalizeLoading();
  return 0;
}

CT_MAIN
