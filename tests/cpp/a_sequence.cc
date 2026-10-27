// Area A C++ driver: base and seed encoding, Hash64, hit encoding,
// Minimizer, MinimizerGenerator and SequenceBatch (tests/test_a_sequence.py).
//
// CT_TEST cases check values worked by hand in their comments. CT_PROBE
// functions print what chromap returns so tests/test_a_sequence.py can compare
// it with values computed independently in Python.

#include "chromap_test.h"

#include "hit_utils.h"
#include "minimizer.h"
#include "minimizer_generator.h"
#include "sequence_batch.h"
#include "strand.h"
#include "utils.h"

using namespace chromap;

namespace {

// One record of a batch as a tab-separated line:
// rec slot id name name_length length sequence quality(or *).
void PrintRecord(const SequenceBatch &batch, uint32_t i) {
  const char *qual = batch.GetSequenceQualAt(i);
  std::cout << "rec\t" << i << "\t" << batch.GetSequenceIdAt(i) << "\t"
            << batch.GetSequenceNameAt(i) << "\t"
            << batch.GetSequenceNameLengthAt(i) << "\t"
            << batch.GetSequenceLengthAt(i) << "\t" << batch.GetSequenceAt(i)
            << "\t" << (qual == NULL ? "*" : qual) << "\n";
}

std::string JoinInts(const std::vector<int> &v) {
  std::ostringstream os;
  for (size_t i = 0; i < v.size(); ++i) os << (i ? "," : "") << v[i];
  return v.empty() ? "-" : os.str();
}

}  // namespace

// ACGT and acgt are 0-3, every other byte 4; 0-3 decode to ACGT, 4-7 to N.
// The tables and the functions agree byte for byte.
CT_TEST(base_tables_hand_worked) {
  const std::string upper = "ACGT", lower = "acgt";
  for (int b = 0; b < 4; ++b) {
    CT_CHECK_EQ((int)CharToUint8(upper[b]), b);
    CT_CHECK_EQ((int)CharToUint8(lower[b]), b);
    CT_CHECK_EQ(Uint8ToChar((uint8_t)b), upper[b]);
  }
  for (int i = 4; i < 8; ++i) CT_CHECK_EQ(Uint8ToChar((uint8_t)i), 'N');
  int num_fours = 0;
  for (int c = 0; c < 256; ++c) {
    CT_CHECK_EQ(CharToUint8((char)c), char_to_uint8_table_[c]);
    if (CharToUint8((char)c) == 4) ++num_fours;
  }
  CT_CHECK_EQ(num_fours, 248);
  for (int i = 0; i < 8; ++i) CT_CHECK_EQ(Uint8ToChar((uint8_t)i), uint8_to_char_table_[i]);
  // N, n, U, R, '-', '\0' and high bytes are ambiguous.
  CT_CHECK_EQ((int)CharToUint8('N'), 4);
  CT_CHECK_EQ((int)CharToUint8('n'), 4);
  CT_CHECK_EQ((int)CharToUint8('U'), 4);
  CT_CHECK_EQ((int)CharToUint8('R'), 4);
  CT_CHECK_EQ((int)CharToUint8('-'), 4);
  CT_CHECK_EQ((int)CharToUint8('\0'), 4);
  CT_CHECK_EQ((int)CharToUint8((char)0xC1), 4);
}

// Seeds of length 1 and 32: "T" = 3; 32 T's fill all 64 bits; "G" + 31 A's
// puts 10 in the top two bits (2 << 62); an empty range is 0.
CT_TEST(seed_lengths_one_and_32) {
  CT_CHECK_EQ(GenerateSeedFromSequence("T", 1, 0, 1), (uint64_t)3);
  CT_CHECK_EQ(GenerateSeedFromSequence("C", 1, 0, 1), (uint64_t)1);
  const std::string t32(32, 'T');
  CT_CHECK_EQ(GenerateSeedFromSequence(t32.c_str(), 32, 0, 32), UINT64_MAX);
  const std::string g31a = "G" + std::string(31, 'A');
  CT_CHECK_EQ(GenerateSeedFromSequence(g31a.c_str(), 32, 0, 32), ((uint64_t)2) << 62);
  // The last base sits in the lowest bits: 31 A's + "C" = 1.
  const std::string a31c = std::string(31, 'A') + "C";
  CT_CHECK_EQ(GenerateSeedFromSequence(a31c.c_str(), 32, 0, 32), (uint64_t)1);
  CT_CHECK_EQ(GenerateSeedFromSequence("ACGT", 4, 1, 0), (uint64_t)0);
  // A start past the end pads the whole seed with A: 0.
  CT_CHECK_EQ(GenerateSeedFromSequence("TTTT", 4, 4, 3), (uint64_t)0);
}

// hit = ((index << 32 | position) << 1) | strand at the extremes: index
// 2^31 - 1 fills bits 33-63, position 2^32 - 1 bits 1-32.
CT_TEST(hit_encoding_extremes) {
  const uint32_t max_index = 0x7FFFFFFFu, max_position = 0xFFFFFFFFu;
  const uint64_t hit = GenerateMinimizer(max_index, max_position, kNegative);
  CT_CHECK_EQ(hit, UINT64_MAX);
  CT_CHECK_EQ(HitToSequenceIndex(hit), max_index);
  CT_CHECK_EQ(HitToSequencePosition(hit), max_position);
  CT_CHECK(HitToStrand(hit) == kNegative);
  const uint64_t hit2 = GenerateMinimizer(max_index, 0, kPositive);
  CT_CHECK_EQ(hit2, ((uint64_t)max_index) << 33);
  CT_CHECK_EQ(HitToSequencePosition(hit2), 0u);
  CT_CHECK(HitToStrand(hit2) == kPositive);
  const uint64_t hit3 = GenerateMinimizer(0, max_position, kPositive);
  CT_CHECK_EQ(hit3, ((uint64_t)max_position) << 1);
  CT_CHECK_EQ(HitToSequenceIndex(hit3), 0u);
  CT_CHECK_EQ(GenerateMinimizer(0, 0, kNegative), (uint64_t)1);
  CT_CHECK_EQ(GenerateMinimizer(1, 0, kPositive), ((uint64_t)1) << 33);
  CT_CHECK_EQ(GenerateMinimizer(0, 1, kPositive), (uint64_t)2);
  CT_CHECK(AreTwoHitsOnTheSameStrand(hit, GenerateMinimizer(5, 9, kNegative)));
  CT_CHECK(AreTwoHitsOnTheSameStrand(hit2, hit3));
  CT_CHECK(!AreTwoHitsOnTheSameStrand(hit, hit2));
  CT_CHECK_EQ((int)kPositive, 0);
  CT_CHECK_EQ((int)kNegative, 1);
}

// Both constructors give the same getters; operator< orders by hash, then by
// hit, and is a strict weak order (irreflexive).
CT_TEST(minimizer_class) {
  const uint64_t hit = GenerateMinimizer(7, 1000, kNegative);
  Minimizer a(std::make_pair((uint64_t)99, hit));
  Minimizer b((uint64_t)99, hit);
  CT_CHECK_EQ(a.GetHash(), (uint64_t)99);
  CT_CHECK_EQ(a.GetHit(), hit);
  CT_CHECK_EQ(b.GetHit(), hit);
  CT_CHECK_EQ(a.GetSequenceIndex(), 7u);
  CT_CHECK_EQ(a.GetSequencePosition(), 1000u);
  CT_CHECK(a.GetSequenceStrand() == kNegative);
  CT_CHECK(!(a < b) && !(b < a));
  Minimizer smaller_hash(98, UINT64_MAX), larger_hit(99, hit + 1), smaller_hit(99, hit - 1);
  CT_CHECK(smaller_hash < a);
  CT_CHECK(!(a < smaller_hash));
  CT_CHECK(a < larger_hit);
  CT_CHECK(!(larger_hit < a));
  CT_CHECK(smaller_hit < a);
  Minimizer plus(5, GenerateMinimizer(0, 3, kPositive));
  CT_CHECK(plus.GetSequenceStrand() == kPositive);
  CT_CHECK_EQ(plus.GetSequencePosition(), 3u);
}

// seeds: stdin lines "SEQ START LENGTH" -> GenerateSeedFromSequence, one per
// line. A sequence of "." is the empty sequence.
CT_PROBE(seeds) {
  std::string seq;
  uint32_t start = 0, length = 0;
  while (std::cin >> seq >> start >> length) {
    if (seq == ".") seq.clear();
    std::cout << GenerateSeedFromSequence(seq.c_str(), seq.size(), start, length) << "\n";
  }
  return 0;
}

// hash64: stdin lines "KEY BITS" -> Hash64(KEY, 2^BITS - 1).
CT_PROBE(hash64) {
  uint64_t key = 0;
  int bits = 0;
  while (std::cin >> key >> bits) {
    const uint64_t mask = bits >= 64 ? UINT64_MAX : ((((uint64_t)1) << bits) - 1);
    std::cout << Hash64(key, mask) << "\n";
  }
  return 0;
}

// minimizers FASTA K W [append]: every sequence of FASTA loaded with
// LoadAllSequences, minimizers from MinimizerGenerator(K, W) as
// "index<TAB>hash<TAB>hit" lines in the generator's order. With "append" the
// output vector starts with one sentinel (hash 12345, hit 678), printed first.
CT_PROBE(minimizers) {
  if (args.size() < 3) {
    std::cerr << "usage: minimizers FASTA K W [append]" << std::endl;
    return 2;
  }
  SequenceBatch batch;
  batch.InitializeLoading(args[0]);
  batch.LoadAllSequences();
  MinimizerGenerator generator(std::stoi(args[1]), std::stoi(args[2]));
  const bool append = args.size() > 3 && args[3] == "append";
  for (uint32_t i = 0; i < batch.GetNumSequences(); ++i) {
    std::vector<Minimizer> minimizers;
    if (append) minimizers.emplace_back((uint64_t)12345, (uint64_t)678);
    generator.GenerateMinimizers(batch, i, minimizers);
    for (size_t j = 0; j < minimizers.size(); ++j) {
      std::cout << i << "\t" << minimizers[j].GetHash() << "\t" << minimizers[j].GetHit() << "\n";
    }
  }
  batch.FinalizeLoading();
  return 0;
}

// load all FILE: LoadAllSequences with the default constructor.
// load batch FILE N: SequenceBatch(N) and LoadBatch until it returns 0;
// "batch <count>" before each batch's records.
// load one FILE N: SequenceBatch(N), LoadOneSequenceAndSaveAt(0, 1, ...)
// until it returns true; "one <slot> <end>" for each call.
// Each record: rec slot id name name_length length sequence quality.
// Then num_sequences, num_bases and max_batch.
CT_PROBE(load) {
  if (args.size() < 2) {
    std::cerr << "usage: load all|batch|one FILE [N]" << std::endl;
    return 2;
  }
  const std::string mode = args[0];
  if (mode == "all") {
    SequenceBatch batch;
    batch.InitializeLoading(args[1]);
    batch.LoadAllSequences();
    for (uint32_t i = 0; i < batch.GetNumSequences(); ++i) PrintRecord(batch, i);
    std::cout << "num_sequences\t" << batch.GetNumSequences() << "\n";
    std::cout << "num_bases\t" << batch.GetNumBases() << "\n";
    std::cout << "max_batch\t" << batch.GetMaxBatchSize() << "\n";
    batch.FinalizeLoading();
    return 0;
  }
  const uint32_t n = std::stoul(args[2]);
  SequenceBatch batch(n, SequenceEffectiveRange());
  batch.InitializeLoading(args[1]);
  if (mode == "batch") {
    for (;;) {
      const uint32_t count = batch.LoadBatch();
      std::cout << "batch\t" << count << "\t" << batch.GetNumSequences() << "\n";
      if (count == 0) break;
      for (uint32_t i = 0; i < count; ++i) PrintRecord(batch, i);
    }
  } else if (mode == "one") {
    for (uint32_t slot = 0; slot < n; ++slot) {
      const bool end = batch.LoadOneSequenceAndSaveAt(slot);
      std::cout << "one\t" << slot << "\t" << end << "\t" << batch.GetNumSequences() << "\n";
      if (end) break;
      PrintRecord(batch, slot);
    }
  } else {
    std::cerr << "unknown mode " << mode << std::endl;
    return 2;
  }
  std::cout << "num_sequences\t" << batch.GetNumSequences() << "\n";
  std::cout << "num_bases\t" << batch.GetNumBases() << "\n";
  std::cout << "max_batch\t" << batch.GetMaxBatchSize() << "\n";
  batch.FinalizeLoading();
  return 0;
}

// slots FILE N S1 S2 ...: SequenceBatch(N), then LoadOneSequenceAndSaveAt(Si)
// for each slot in turn, printing "slot <Si> <end> <num_sequences>" and the
// record (when one was read). chromap exits when asked to overwrite a slot
// other than the last loaded one.
CT_PROBE(slots) {
  const uint32_t n = std::stoul(args[1]);
  SequenceBatch batch(n, SequenceEffectiveRange());
  batch.InitializeLoading(args[0]);
  for (size_t a = 2; a < args.size(); ++a) {
    const uint32_t slot = std::stoul(args[a]);
    const bool end = batch.LoadOneSequenceAndSaveAt(slot);
    std::cout << "slot\t" << slot << "\t" << end << "\t" << batch.GetNumSequences() << std::endl;
    if (!end) PrintRecord(batch, slot);
  }
  batch.FinalizeLoading();
  return 0;
}

// strands FILE N: SequenceBatch(N) and one LoadBatch; for each record the
// sequence, PrepareNegativeSequenceAt's reverse complement (also read back
// through GetNegativeSequenceBatch), IsNInSequenceAt and GetSequenceNsAt in
// sequence order and in seed (right-to-left) order.
CT_PROBE(strands) {
  const uint32_t n = std::stoul(args[1]);
  SequenceBatch batch(n, SequenceEffectiveRange());
  batch.InitializeLoading(args[0]);
  const uint32_t count = batch.LoadBatch();
  for (uint32_t i = 0; i < count; ++i) {
    batch.PrepareNegativeSequenceAt(i);
    std::vector<int> big, little;
    batch.GetSequenceNsAt(i, false, big);
    batch.GetSequenceNsAt(i, true, little);
    std::cout << batch.GetSequenceAt(i) << "\t" << batch.GetNegativeSequenceAt(i) << "\t"
              << batch.GetNegativeSequenceBatch()[i] << "\t" << batch.IsNInSequenceAt(i) << "\t"
              << JoinInts(big) << "\t" << JoinInts(little) << "\n";
  }
  batch.FinalizeLoading();
  return 0;
}

// trim FILE N LENGTH: SequenceBatch(N), one LoadBatch, PrepareNegativeSequenceAt
// then TrimSequenceAt(LENGTH) for every record; prints length, sequence,
// quality and reverse complement after trimming.
CT_PROBE(trim) {
  const uint32_t n = std::stoul(args[1]);
  const int length = std::stoi(args[2]);
  SequenceBatch batch(n, SequenceEffectiveRange());
  batch.InitializeLoading(args[0]);
  const uint32_t count = batch.LoadBatch();
  for (uint32_t i = 0; i < count; ++i) {
    batch.PrepareNegativeSequenceAt(i);
    batch.TrimSequenceAt(i, length);
    std::cout << batch.GetSequenceLengthAt(i) << "\t" << batch.GetSequenceAt(i) << "\t"
              << batch.GetSequenceQualAt(i) << "\t" << batch.GetNegativeSequenceAt(i) << "\n";
  }
  batch.FinalizeLoading();
  return 0;
}

// edit FILE N I POSITION BASE: SequenceBatch(N), one LoadBatch,
// CorrectBaseAt(I, POSITION, BASE); prints every sequence.
CT_PROBE(edit) {
  const uint32_t n = std::stoul(args[1]);
  SequenceBatch batch(n, SequenceEffectiveRange());
  batch.InitializeLoading(args[0]);
  const uint32_t count = batch.LoadBatch();
  batch.CorrectBaseAt(std::stoul(args[2]), std::stoul(args[3]), args[4][0]);
  for (uint32_t i = 0; i < count; ++i) std::cout << batch.GetSequenceAt(i) << "\n";
  batch.FinalizeLoading();
  return 0;
}

// reorder all FILE R0 R1 ...: LoadAllSequences, ReorderSequences(ranks).
// reorder batch FILE R0 R1 ...: SequenceBatch(number of ranks), one LoadBatch,
// PrepareNegativeSequenceAt for each, ReorderSequences(ranks).
// Prints "name<TAB>sequence<TAB>reverse complement or *" for slots 0, 1, ...
CT_PROBE(reorder) {
  std::vector<int> ranks;
  for (size_t a = 2; a < args.size(); ++a) ranks.push_back(std::stoi(args[a]));
  const bool all = args[0] == "all";
  SequenceBatch *batch = all ? new SequenceBatch()
                             : new SequenceBatch(ranks.size(), SequenceEffectiveRange());
  batch->InitializeLoading(args[1]);
  uint32_t count = 0;
  if (all) {
    batch->LoadAllSequences();
    count = batch->GetNumSequences();
  } else {
    count = batch->LoadBatch();
    for (uint32_t i = 0; i < count; ++i) batch->PrepareNegativeSequenceAt(i);
  }
  batch->ReorderSequences(ranks);
  for (uint32_t i = 0; i < count; ++i) {
    std::cout << batch->GetSequenceNameAt(i) << "\t" << batch->GetSequenceAt(i) << "\t"
              << (all ? std::string("*") : batch->GetNegativeSequenceAt(i)) << "\n";
  }
  batch->FinalizeLoading();
  delete batch;
  return 0;
}

// swap FILE1 FILE2 N: two SequenceBatch(N), one LoadBatch each and the reverse
// complements prepared; after first.SwapSequenceBatch(second), prints
// "first<TAB>sequence<TAB>reverse complement" for the first batch's slots and
// "second<TAB>..." for the second's (N slots each, empty slots skipped).
CT_PROBE(swap) {
  const uint32_t n = std::stoul(args[2]);
  SequenceBatch first(n, SequenceEffectiveRange()), second(n, SequenceEffectiveRange());
  first.InitializeLoading(args[0]);
  second.InitializeLoading(args[1]);
  const uint32_t c1 = first.LoadBatch(), c2 = second.LoadBatch();
  for (uint32_t i = 0; i < c1; ++i) first.PrepareNegativeSequenceAt(i);
  for (uint32_t i = 0; i < c2; ++i) second.PrepareNegativeSequenceAt(i);
  first.SwapSequenceBatch(second);
  for (uint32_t i = 0; i < n; ++i) {
    if (first.GetSequenceAt(i) != NULL)
      std::cout << "first\t" << first.GetSequenceAt(i) << "\t" << first.GetNegativeSequenceAt(i) << "\n";
  }
  for (uint32_t i = 0; i < n; ++i) {
    if (second.GetSequenceAt(i) != NULL)
      std::cout << "second\t" << second.GetSequenceAt(i) << "\t" << second.GetNegativeSequenceAt(i) << "\n";
  }
  first.FinalizeLoading();
  second.FinalizeLoading();
  return 0;
}

// batch_seeds FILE START LENGTH: LoadAllSequences, then
// GenerateSeedFromSequenceAt(i, START, LENGTH) for every sequence.
CT_PROBE(batch_seeds) {
  SequenceBatch batch;
  batch.InitializeLoading(args[0]);
  batch.LoadAllSequences();
  for (uint32_t i = 0; i < batch.GetNumSequences(); ++i) {
    std::cout << batch.GenerateSeedFromSequenceAt(i, std::stoul(args[1]), std::stoul(args[2])) << "\n";
  }
  batch.FinalizeLoading();
  return 0;
}

CT_MAIN
