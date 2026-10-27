// Area B unit tests of the candidate cache (src/mmcache.hpp). The command
// line reaches a filled cache only after a batch of 500,000 reads (one test
// in tests/test_b_cache.py does), so its rules are tested here directly.
// They are part of chromap's output: which reads hit the cache decides the
// summary's cachehit, fric, estfrip and numcacheslots columns (README), so a
// change to them changes the summary of any large input even where the
// mapping output stays the same:
//
// - Update(minimizers, candidates) counts updates per slot (slot = (first +
//   last minimizer hash) mod size) and per fingerprint ((first ^ last) mod
//   103). A slot stores an entry only once it has seen at least 10 updates
//   and the entry's fingerprint accounts for at least a fifth of them; a
//   stored entry is replaced when competing dominant updates drive its
//   weight below zero; after 100 updates the slot no longer changes.
// - Query returns the slot index for the same minimizers (hashes, strands and
//   spacing) or for the reverse-complemented read, with the stored
//   candidates moved to the query's coordinates, and -1 otherwise.
// Expected values are worked in the comments.

#include "chromap_test.h"

#include <omp.h>

#include "candidate.h"
#include "mapping_metadata.h"
#include "minimizer.h"
#include "mmcache.hpp"

using namespace chromap;

namespace {

class CacheMetadata : public MappingMetadata {
 public:
  using MappingMetadata::minimizers_;
  using MappingMetadata::negative_candidates_;
  using MappingMetadata::positive_candidates_;
  using MappingMetadata::repetitive_seed_length_;
};

const int kSlots = 1009;
const int kK = 17;
const uint32_t kReadLength = 50;

Candidate MakeCandidate(uint32_t rid, uint32_t position, uint8_t count) {
  Candidate c;
  c.position = (static_cast<uint64_t>(rid) << 32) | position;
  c.count = count;
  return c;
}

// hit = position << 1 | strand (read index 0).
std::vector<Minimizer> MakeMinimizers(const std::vector<uint64_t> &hashes, const std::vector<uint32_t> &positions,
                                      const std::vector<int> &strands) {
  std::vector<Minimizer> out;
  for (size_t i = 0; i < hashes.size(); ++i) {
    out.emplace_back(hashes[i], (static_cast<uint64_t>(positions[i]) << 1) | (uint64_t)strands[i]);
  }
  return out;
}

// Read A: four minimizers ending at read positions 20, 24, 31, 40. Its slot
// is (1001 + 4004) mod 1009 = 969.
std::vector<Minimizer> ReadA(uint32_t shift = 0) {
  return MakeMinimizers({1001, 2002, 3003, 4004}, {20 + shift, 24 + shift, 31 + shift, 40 + shift}, {0, 1, 0, 0});
}
const int kSlotA = 969;
std::vector<Candidate> PositiveA() { return {MakeCandidate(0, 5000, 4), MakeCandidate(2, 700, 4)}; }
std::vector<Candidate> NegativeA() { return {MakeCandidate(1, 9000, 3)}; }

// B: same slot as A ((1500 + 3505) mod 1009 = 969), another fingerprint.
std::vector<Minimizer> ReadB() { return MakeMinimizers({1500, 2000, 3505}, {18, 25, 33}, {1, 1, 0}); }
std::vector<Candidate> PositiveB() { return {MakeCandidate(3, 4000, 3)}; }

int Fingerprint(const std::vector<Minimizer> &m) {
  return (int)((m.front().GetHash() ^ m.back().GetHash()) % FINGER_PRINT_SIZE);
}

// `n` two-minimizer reads in A's slot with pairwise different fingerprints,
// none equal to A's or B's: updates that never dominate the slot.
std::vector<std::vector<Minimizer>> Noise(int n) {
  std::vector<std::vector<Minimizer>> out;
  std::vector<bool> used(FINGER_PRINT_SIZE, false);
  used[Fingerprint(ReadA())] = true;
  used[Fingerprint(ReadB())] = true;
  // (x + y) mod 1009 = 969 for y = 5005 + 1009 j - x.
  for (uint64_t j = 0; (int)out.size() < n && j < 64; ++j) {
    for (uint64_t x = 10; (int)out.size() < n && x < 5005; ++x) {
      std::vector<Minimizer> m = MakeMinimizers({x, 5005 + kSlots * j - x}, {20, 30}, {0, 0});
      int f = Fingerprint(m);
      if (used[f]) continue;
      used[f] = true;
      out.push_back(m);
    }
  }
  return out;
}

void UpdateTimes(mm_cache &cache, const std::vector<Minimizer> &m, std::vector<Candidate> positive,
                 std::vector<Candidate> negative, int times, uint32_t repetitive = 0) {
  for (int i = 0; i < times; ++i) cache.Update(m, positive, negative, repetitive);
}

void UpdateNoise(mm_cache &cache, int n) {
  std::vector<std::vector<Minimizer>> noise = Noise(n);
  CT_CHECK_EQ((int)noise.size(), n);
  std::vector<Candidate> positive = {MakeCandidate(0, 3000, 2)};
  std::vector<Candidate> negative;
  for (size_t i = 0; i < noise.size(); ++i) cache.Update(noise[i], positive, negative, 0);
}

int QueryWith(mm_cache &cache, const std::vector<Minimizer> &m, CacheMetadata &out) {
  out.PrepareForMappingNextRead(8);
  out.minimizers_ = m;
  return cache.Query(out, kReadLength);
}

bool SameCandidates(const std::vector<Candidate> &a, const std::vector<Candidate> &b) {
  if (a.size() != b.size()) return false;
  for (size_t i = 0; i < a.size(); ++i) {
    if (a[i].position != b[i].position || a[i].count != b[i].count) return false;
  }
  return true;
}

}  // namespace

// Nine updates leave the slot empty; the tenth stores A, and the same
// minimizers then return the slot with A's candidates and repetitive seed
// length unchanged.
CT_TEST(cache_stores_after_ten_updates) {
  mm_cache cache(kSlots);
  cache.SetKmerLength(kK);
  CacheMetadata q;
  CT_CHECK_EQ(QueryWith(cache, ReadA(), q), -1);
  UpdateTimes(cache, ReadA(), PositiveA(), NegativeA(), 9, 7);
  CT_CHECK_EQ(QueryWith(cache, ReadA(), q), -1);
  CT_CHECK(q.positive_candidates_.empty());
  UpdateTimes(cache, ReadA(), PositiveA(), NegativeA(), 1, 7);
  CT_CHECK_EQ(QueryWith(cache, ReadA(), q), kSlotA);
  CT_CHECK(SameCandidates(q.positive_candidates_, PositiveA()));
  CT_CHECK(SameCandidates(q.negative_candidates_, NegativeA()));
  CT_CHECK_EQ(q.repetitive_seed_length_, 7u);
}

// A read whose minimizers sit 3 positions later starts 3 bases earlier on
// the reference: forward candidates (alignment starts) move to 5000 - 3 and
// 700 - 3; the reverse candidate (alignment end, the read's 5' end) moves to
// 9000 + 3. A changed spacing, strand, hash or minimizer count misses.
CT_TEST(cache_query_moves_candidates_to_query_positions) {
  mm_cache cache(kSlots);
  cache.SetKmerLength(kK);
  UpdateTimes(cache, ReadA(), PositiveA(), NegativeA(), 10);
  CacheMetadata q;
  CT_CHECK_EQ(QueryWith(cache, ReadA(3), q), kSlotA);
  CT_CHECK_EQ(q.positive_candidates_.size(), (size_t)2);
  if (q.positive_candidates_.size() == 2) {
    CT_CHECK_EQ(q.positive_candidates_[0].position, (uint64_t)4997);
    CT_CHECK_EQ(q.positive_candidates_[1].position, (((uint64_t)2) << 32) | 697u);
    CT_CHECK_EQ((int)q.positive_candidates_[0].count, 4);
  }
  CT_CHECK_EQ(q.negative_candidates_.size(), (size_t)1);
  if (q.negative_candidates_.size() == 1) {
    CT_CHECK_EQ(q.negative_candidates_[0].position, (((uint64_t)1) << 32) | 9003u);
  }
  std::vector<Minimizer> spacing = MakeMinimizers({1001, 2002, 3003, 4004}, {20, 25, 31, 40}, {0, 1, 0, 0});
  CT_CHECK_EQ(QueryWith(cache, spacing, q), -1);
  std::vector<Minimizer> strand = MakeMinimizers({1001, 2002, 3003, 4004}, {20, 24, 31, 40}, {0, 0, 0, 0});
  CT_CHECK_EQ(QueryWith(cache, strand, q), -1);
  std::vector<Minimizer> hash = MakeMinimizers({1001, 2003, 3003, 4004}, {20, 24, 31, 40}, {0, 1, 0, 0});
  CT_CHECK_EQ(QueryWith(cache, hash, q), -1);
  std::vector<Minimizer> count = MakeMinimizers({1001, 2002, 4004}, {20, 24, 40}, {0, 1, 0});
  CT_CHECK_EQ(QueryWith(cache, count, q), -1);
  CT_CHECK_EQ(QueryWith(cache, ReadB(), q), -1);
}

// The reverse complement of a 50 bp read has the same minimizers in reverse
// order with flipped strands, a k-mer ending at p now ending at
// 50 - p + 17 - 2. It maps to the same locus on the other strand: its
// forward candidates are A's reverse candidates turned from end to start
// (9000 - 50 + 1 = 8951) and its reverse candidates A's forward candidates
// turned from start to end (5000 + 49, 700 + 49).
CT_TEST(cache_query_reverse_complement) {
  mm_cache cache(kSlots);
  cache.SetKmerLength(kK);
  UpdateTimes(cache, ReadA(), PositiveA(), NegativeA(), 10, 5);
  std::vector<Minimizer> a = ReadA();
  std::vector<uint64_t> hashes;
  std::vector<uint32_t> positions;
  std::vector<int> strands;
  for (int i = (int)a.size() - 1; i >= 0; --i) {
    hashes.push_back(a[i].GetHash());
    positions.push_back(kReadLength - a[i].GetSequencePosition() + kK - 2);
    strands.push_back(1 - (int)(a[i].GetHit() & 1));
  }
  CacheMetadata q;
  CT_CHECK_EQ(QueryWith(cache, MakeMinimizers(hashes, positions, strands), q), kSlotA);
  CT_CHECK_EQ(q.positive_candidates_.size(), (size_t)1);
  if (q.positive_candidates_.size() == 1) {
    CT_CHECK_EQ(q.positive_candidates_[0].position, (((uint64_t)1) << 32) | 8951u);
  }
  CT_CHECK_EQ(q.negative_candidates_.size(), (size_t)2);
  if (q.negative_candidates_.size() == 2) {
    CT_CHECK_EQ(q.negative_candidates_[0].position, (uint64_t)5049);
    CT_CHECK_EQ(q.negative_candidates_[1].position, (((uint64_t)2) << 32) | 749u);
  }
  CT_CHECK_EQ(q.repetitive_seed_length_, 5u);
}

// Updates of other reads count towards the 10: after 5 never-dominant
// updates, A is stored at its 5th update (sum 10, 5 * 5 >= 10). After 45 of
// them A needs 5 * n >= 45 + n, so n = 12: still absent after 11.
CT_TEST(cache_needs_a_dominant_fingerprint) {
  {
    mm_cache cache(kSlots);
    cache.SetKmerLength(kK);
    UpdateNoise(cache, 5);
    CacheMetadata q;
    UpdateTimes(cache, ReadA(), PositiveA(), NegativeA(), 4);
    CT_CHECK_EQ(QueryWith(cache, ReadA(), q), -1);
    UpdateTimes(cache, ReadA(), PositiveA(), NegativeA(), 1);
    CT_CHECK_EQ(QueryWith(cache, ReadA(), q), kSlotA);
  }
  {
    mm_cache cache(kSlots);
    cache.SetKmerLength(kK);
    UpdateNoise(cache, 45);
    CacheMetadata q;
    UpdateTimes(cache, ReadA(), PositiveA(), NegativeA(), 11);
    CT_CHECK_EQ(QueryWith(cache, ReadA(), q), -1);
    UpdateTimes(cache, ReadA(), PositiveA(), NegativeA(), 1);
    CT_CHECK_EQ(QueryWith(cache, ReadA(), q), kSlotA);
  }
}

// Replacement and saturation. A is stored after 10 updates (weight 1). With
// 60 more never-dominant updates (sum 70), B dominates from its 18th update
// (5 * 18 >= 88): weight 1 -> 0, and the 19th (sum 89) drops it below 0 and
// stores B. With 80 such updates (sum 90) B could only dominate from its
// 23rd update (sum 113), past the 100 after which a slot no longer changes,
// so A stays.
CT_TEST(cache_replacement_and_saturation) {
  {
    mm_cache cache(kSlots);
    cache.SetKmerLength(kK);
    UpdateTimes(cache, ReadA(), PositiveA(), NegativeA(), 10);
    UpdateNoise(cache, 60);
    CacheMetadata q;
    UpdateTimes(cache, ReadB(), PositiveB(), std::vector<Candidate>(), 18);
    CT_CHECK_EQ(QueryWith(cache, ReadB(), q), -1);
    CT_CHECK_EQ(QueryWith(cache, ReadA(), q), kSlotA);
    UpdateTimes(cache, ReadB(), PositiveB(), std::vector<Candidate>(), 1);
    CT_CHECK_EQ(QueryWith(cache, ReadB(), q), kSlotA);
    CT_CHECK(SameCandidates(q.positive_candidates_, PositiveB()));
    CT_CHECK_EQ(QueryWith(cache, ReadA(), q), -1);
  }
  {
    mm_cache cache(kSlots);
    cache.SetKmerLength(kK);
    UpdateTimes(cache, ReadA(), PositiveA(), NegativeA(), 10);
    UpdateNoise(cache, 80);
    CacheMetadata q;
    UpdateTimes(cache, ReadB(), PositiveB(), std::vector<Candidate>(), 40);
    CT_CHECK_EQ(QueryWith(cache, ReadB(), q), -1);
    CT_CHECK_EQ(QueryWith(cache, ReadA(), q), kSlotA);
  }
}

// Candidates that start within k + (first minimizer position) of the
// reference start (here 30 < 17 + 20) are not cached: the read keeps
// missing.
CT_TEST(cache_skips_candidates_near_reference_start) {
  mm_cache cache(kSlots);
  cache.SetKmerLength(kK);
  std::vector<Candidate> near = {MakeCandidate(0, 30, 4)};
  UpdateTimes(cache, ReadA(), near, std::vector<Candidate>(), 20);
  CacheMetadata q;
  CT_CHECK_EQ(QueryWith(cache, ReadA(), q), -1);
  // Reverse candidates: end - last minimizer position (40) must reach k + 20.
  std::vector<Candidate> near_end = {MakeCandidate(0, 76, 4)};
  mm_cache cache2(kSlots);
  cache2.SetKmerLength(kK);
  UpdateTimes(cache2, ReadA(), PositiveA(), near_end, 20);
  CT_CHECK_EQ(QueryWith(cache2, ReadA(), q), -1);
  std::vector<Candidate> far_end = {MakeCandidate(0, 77, 4)};
  mm_cache cache3(kSlots);
  cache3.SetKmerLength(kK);
  UpdateTimes(cache3, ReadA(), PositiveA(), far_end, 10);
  CT_CHECK_EQ(QueryWith(cache3, ReadA(), q), kSlotA);
}

// How many reads of a batch update the cache: all of them while the reads
// seen so far fit in one block (5,000,000 single-end reads, 2,500,000 for
// pairs), then n / (1 + param * (reads seen / block)), integer blocks.
CT_TEST(cache_update_threshold) {
  mm_cache cache(13);
  CT_CHECK_EQ(cache.GetUpdateThreshold(500000, 5000000, false, 0.01), 500000u);
  CT_CHECK_EQ(cache.GetUpdateThreshold(500000, 5000001, false, 0.01), 495049u);
  CT_CHECK_EQ(cache.GetUpdateThreshold(500000, 2500000, true, 0.01), 500000u);
  CT_CHECK_EQ(cache.GetUpdateThreshold(500000, 12500000, true, 0.01), 476190u);
  CT_CHECK_EQ(cache.GetUpdateThreshold(500000, 2500001, true, 1.0), 250000u);
  CT_CHECK_EQ(cache.GetUpdateThreshold(1000, 99999999, true, 0.0), 1000u);
}

// printstats SIZE: a cache of SIZE slots, read A updated 10 times with
// debug output, then PrintStats.
CT_PROBE(printstats) {
  mm_cache cache(std::stoi(args.at(0)));
  cache.SetKmerLength(kK);
  std::vector<Candidate> positive = PositiveA();
  std::vector<Candidate> negative = NegativeA();
  for (int i = 0; i < 10; ++i) cache.Update(ReadA(), positive, negative, 0, /*debug=*/true);
  std::cout.flush();
  cache.PrintStats();
  fflush(stdout);
  return 0;
}

CT_MAIN
