// Area A C++ driver: index_utils.h bit layouts and heap merge,
// IndexParameters, and Index load, save, check and destroy
// (tests/test_a_index.py).
//
// CT_TEST cases check values worked by hand in their comments. CT_PROBE
// functions print what chromap returns so the Python tests can compare it
// with values computed independently.

#include "chromap_test.h"

#include "index.h"
#include "index_parameters.h"
#include "index_utils.h"
#include "sequence_batch.h"

using namespace chromap;

// The lookup key of a minimizer hash is the hash shifted left by one; bit 0
// marks singletons. Entry values pack offset << 32 | count.
CT_TEST(lookup_table_layout_hand_worked) {
  CT_CHECK_EQ(GenerateHashInLookupTable(0x5ULL), 0xAULL);
  // 56-bit hashes (k = 28) keep every bit after the shift.
  const uint64_t h56 = (((uint64_t)1) << 56) - 1;
  CT_CHECK_EQ(GenerateHashInLookupTable(h56) >> 1, h56);
  CT_CHECK(!IsSingletonLookupKey(GenerateHashInLookupTable(h56)));
  CT_CHECK(IsSingletonLookupKey(GenerateHashInLookupTable(h56) | 1));
  CT_CHECK(IsSingletonLookupKey(1));
  CT_CHECK(!IsSingletonLookupKey(2));
  // offset 3, count 7: 3 << 32 | 7.
  const uint64_t value = GenerateEntryValueInLookupTable(3, 7);
  CT_CHECK_EQ(value, (((uint64_t)3) << 32) | 7);
  CT_CHECK_EQ(GenerateOffsetInOccurrenceTable(value), 3u);
  CT_CHECK_EQ(GenerateNumOccurrenceInOccurrenceTable(value), 7u);
  // Extremes: offset and count of 2^32 - 1 each.
  const uint64_t full = GenerateEntryValueInLookupTable(0xFFFFFFFFu, 0xFFFFFFFFu);
  CT_CHECK_EQ(full, UINT64_MAX);
  CT_CHECK_EQ(GenerateOffsetInOccurrenceTable(full), 0xFFFFFFFFu);
  CT_CHECK_EQ(GenerateNumOccurrenceInOccurrenceTable(full), 0xFFFFFFFFu);
  CT_CHECK_EQ(GenerateNumOccurrenceInOccurrenceTable(GenerateEntryValueInLookupTable(9, 0)), 0u);
}

// Candidate positions: sequence id in the high 32 bits, position in the low
// 32; an occurrence-table entry (a hit) drops its strand bit: hit 2 * (5 << 32
// | 100) + 1 gives 5 << 32 | 100.
CT_TEST(candidate_position_layout_hand_worked) {
  CT_CHECK_EQ(SequenceIndexAndPositionToCandidatePosition(5, 100), (((uint64_t)5) << 32) | 100);
  CT_CHECK_EQ(SequenceIndexAndPositionToCandidatePosition(0, 0xFFFFFFFFu), (uint64_t)0xFFFFFFFFu);
  const uint64_t hit = ((((uint64_t)5) << 32 | 100) << 1) | 1;
  CT_CHECK_EQ(GenerateCandidatePositionFromOccurrenceTableEntry(hit), (((uint64_t)5) << 32) | 100);
  CT_CHECK_EQ(GenerateCandidatePositionFromOccurrenceTableEntry(hit - 1), (((uint64_t)5) << 32) | 100);
}

// The index's khash hashes and compares keys without their singleton bit.
CT_TEST(khash_functions_ignore_singleton_bit) {
  CT_CHECK_EQ(KHashFunctionForIndex((uint64_t)0x15), (uint64_t)0xA);
  CT_CHECK(KHashEqForIndex((uint64_t)0x14, (uint64_t)0x15));
  CT_CHECK(!KHashEqForIndex((uint64_t)0x14, (uint64_t)0x16));
}

// CandidatePositionWithListIndex compares in reverse so a priority_queue pops
// the smallest position; RepetitiveSeedStats starts empty.
CT_TEST(heap_entry_order_and_seed_stats) {
  CandidatePositionWithListIndex a(0, 10), b(1, 20);
  CT_CHECK(b < a);
  CT_CHECK(!(a < b));
  CT_CHECK_EQ(a.list_index, 0u);
  CT_CHECK_EQ(b.position, (uint64_t)20);
  std::priority_queue<CandidatePositionWithListIndex> heap;
  heap.push(b);
  heap.push(a);
  heap.emplace(2, 15);
  CT_CHECK_EQ(heap.top().position, (uint64_t)10);
  RepetitiveSeedStats stats;
  CT_CHECK_EQ(stats.repetitive_seed_length, 0u);
  CT_CHECK_EQ(stats.previous_repetitive_seed_position, std::numeric_limits<uint32_t>::max());
  CT_CHECK_EQ(stats.repetitive_seed_count, 0);
}

// chromap's index defaults: k 17, w 7, one thread, no paths.
CT_TEST(index_parameters_defaults) {
  IndexParameters p;
  CT_CHECK_EQ(p.kmer_size, 17);
  CT_CHECK_EQ(p.window_size, 7);
  CT_CHECK_EQ(p.num_threads, 1);
  CT_CHECK(p.reference_file_path.empty());
  CT_CHECK(p.index_output_file_path.empty());
}

// heap_merge: stdin lines of space-separated sorted positions (an empty line
// is an empty list; "-" alone too) -> HeapMergeCandidatePositionLists, the
// merged positions on one line. The output vector starts with a sentinel 7 when
// the first argument is "append".
CT_PROBE(heap_merge) {
  std::vector<std::vector<uint64_t> > lists;
  std::string line;
  while (std::getline(std::cin, line)) {
    std::istringstream is(line);
    std::vector<uint64_t> list;
    std::string token;
    while (is >> token) {
      if (token != "-") list.push_back(std::stoull(token));
    }
    lists.push_back(list);
  }
  std::vector<uint64_t> merged;
  if (!args.empty() && args[0] == "append") merged.push_back(7);
  HeapMergeCandidatePositionLists(lists, merged);
  for (size_t i = 0; i < merged.size(); ++i) std::cout << (i ? " " : "") << merged[i];
  std::cout << "\n";
  return 0;
}

// load_save INDEX: Index(INDEX).Load(), print "k w lookup_table_size", then
// Save() back to the same path, Destroy() twice (the destructor destroys a
// third time).
CT_PROBE(load_save) {
  Index index(args[0]);
  index.Load();
  std::cout << index.GetKmerSize() << " " << index.GetWindowSize() << " " << index.GetLookupTableSize() << "\n";
  index.Save();
  index.Destroy();
  index.Destroy();
  return 0;
}

// check INDEX FASTA: load the index and the reference and run
// Index::CheckIndex, which asserts that every minimizer of the reference is
// stored with its hit; prints "checked" when it returns.
CT_PROBE(check) {
  Index index(args[0]);
  index.Load();
  SequenceBatch reference;
  reference.InitializeLoading(args[1]);
  reference.LoadAllSequences();
  index.CheckIndex(reference.GetNumSequences(), reference);
  reference.FinalizeLoading();
  std::cout << "checked\n";
  return 0;
}

CT_MAIN
