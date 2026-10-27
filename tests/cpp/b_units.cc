// Area B C++ unit tests worked by hand: the small value types of the
// mapping core (Candidate, DraftMapping, MappingInMemory, the generating
// config), MappingParameters' defaults, the per-read metadata containers and
// the ksw_semi_global wrappers. Run from tests/test_b_units.py.

#include "chromap_test.h"

#include <cstdio>

#include "candidate.h"
#include "candidate_position_generating_config.h"
#include "draft_mapping.h"
#include "ksw.h"
#include "mapping_in_memory.h"
#include "mapping_metadata.h"
#include "mapping_parameters.h"
#include "paired_end_mapping_metadata.h"
#include "strand.h"

using namespace chromap;

namespace {

// Exposes MappingMetadata's protected vectors to the tests (no access hacks:
// a subclass with using-declarations).
class TestMetadata : public MappingMetadata {
 public:
  using MappingMetadata::minimizers_;
  using MappingMetadata::negative_candidates_;
  using MappingMetadata::negative_candidates_buffer_;
  using MappingMetadata::negative_mappings_;
  using MappingMetadata::positive_candidates_;
  using MappingMetadata::positive_candidates_buffer_;
  using MappingMetadata::positive_mappings_;
  using MappingMetadata::repetitive_seed_length_;
};

class TestPairedMetadata : public PairedEndMappingMetadata {
 public:
  using PairedEndMappingMetadata::mapping_metadata1_;
  using PairedEndMappingMetadata::mapping_metadata2_;
};

Candidate MakeCandidate(uint32_t rid, uint32_t position, uint8_t count) {
  Candidate c;
  c.position = (static_cast<uint64_t>(rid) << 32) | position;
  c.count = count;
  return c;
}

void AddOneCandidate(std::vector<Candidate> &candidates) {
  candidates.push_back(MakeCandidate(1, 7, 3));
}

}  // namespace

// position = rid << 32 | position; operator< orders by count descending,
// then by position ascending.
CT_TEST(candidate_layout_and_order) {
  Candidate c = MakeCandidate(3, 123456, 9);
  CT_CHECK_EQ(c.position, (((uint64_t)3) << 32) + 123456u);
  CT_CHECK_EQ(c.GetReferenceSequenceIndex(), 3u);
  CT_CHECK_EQ(c.GetReferenceSequencePosition(), 123456u);
  CT_CHECK_EQ((int)c.GetCount(), 9);
  Candidate high = MakeCandidate(0, 900, 5);
  Candidate low = MakeCandidate(0, 100, 2);
  CT_CHECK(high < low);
  CT_CHECK(!(low < high));
  Candidate left = MakeCandidate(0, 100, 5);
  Candidate right = MakeCandidate(0, 101, 5);
  CT_CHECK(left < right);
  CT_CHECK(!(right < left));
  CT_CHECK(!(left < left));
  // A later reference sorts after an earlier one at equal counts.
  CT_CHECK(MakeCandidate(0, 4000000000u, 5) < MakeCandidate(1, 0, 5));
  std::vector<Candidate> v = {MakeCandidate(1, 5, 2), MakeCandidate(0, 9, 7), MakeCandidate(0, 3, 2),
                              MakeCandidate(0, 8, 7)};
  std::sort(v.begin(), v.end());
  CT_CHECK_EQ(v[0].position, (uint64_t)8);
  CT_CHECK_EQ(v[1].position, (uint64_t)9);
  CT_CHECK_EQ(v[2].position, (uint64_t)3);
  CT_CHECK_EQ(v[3].position, (((uint64_t)1) << 32) + 5u);
}

// Both thresholds are inclusive: a seed is frequent at frequency >= the max
// seed frequency and repetitive at >= the repetitive frequency.
CT_TEST(generating_config_thresholds) {
  CandidatePositionGeneratingConfig config(/*max_seed_frequency=*/500, /*repetitive_seed_frequency=*/300,
                                           /*use_heap_merge=*/true);
  CT_CHECK(!config.IsFrequentSeed(499));
  CT_CHECK(config.IsFrequentSeed(500));
  CT_CHECK(config.IsFrequentSeed(501));
  CT_CHECK(!config.IsRepetitiveSeed(299));
  CT_CHECK(config.IsRepetitiveSeed(300));
  CT_CHECK(config.UseHeapMerge());
  CT_CHECK_EQ(config.GetMaxSeedFrequency(), 500u);
  CandidatePositionGeneratingConfig other(2, 2, false);
  CT_CHECK(!other.UseHeapMerge());
  CT_CHECK(!other.IsFrequentSeed(1));
  CT_CHECK(other.IsFrequentSeed(2));
}

// DraftMapping keeps the error count and the end position rid << 32 | pos.
CT_TEST(draft_mapping_layout) {
  DraftMapping d(3, (((uint64_t)2) << 32) | 777u);
  CT_CHECK_EQ(d.GetNumErrors(), 3);
  CT_CHECK_EQ(d.GetReferenceSequenceIndex(), 2u);
  CT_CHECK_EQ(d.GetReferenceSequencePosition(), 777u);
}

// Closed intervals [start, end]: read 1 forward at [100, 139] (40 bp) and
// read 2 reverse at [260, 309] (50 bp) make the fragment [100, 309], length
// 210; with read 1 reverse at [260, 309] and read 2 forward at [100, 139]
// the fragment is the same and the strand flag 0.
CT_TEST(mapping_in_memory_fragments) {
  MappingInMemory single;
  single.ref_start_position = 100;
  single.ref_end_position = 149;
  single.strand = kPositive;
  CT_CHECK_EQ(single.GetFragmentStartPosition(), 100u);
  CT_CHECK_EQ((int)single.GetFragmentLength(), 50);
  CT_CHECK_EQ((int)single.GetAlignmentLength(), 50);
  CT_CHECK_EQ((int)single.GetStrand(), 1);
  single.strand = kNegative;
  CT_CHECK_EQ((int)single.GetStrand(), 0);

  PairedEndMappingInMemory pair;
  pair.mapping_in_memory1.read_id = 42;
  pair.mapping_in_memory1.barcode_key = 99;
  pair.mapping_in_memory1.strand = kPositive;
  pair.mapping_in_memory1.ref_start_position = 100;
  pair.mapping_in_memory1.ref_end_position = 139;
  pair.mapping_in_memory2.strand = kNegative;
  pair.mapping_in_memory2.ref_start_position = 260;
  pair.mapping_in_memory2.ref_end_position = 309;
  CT_CHECK_EQ(pair.GetReadId(), 42u);
  CT_CHECK_EQ(pair.GetBarcode(), (uint64_t)99);
  CT_CHECK_EQ((int)pair.GetStrand(), 1);
  CT_CHECK_EQ(pair.GetFragmentStartPosition(), 100u);
  CT_CHECK_EQ(pair.GetFragmentLength(), 210);
  CT_CHECK_EQ(pair.GetPositiveAlignmentLength(), 40u);
  CT_CHECK_EQ(pair.GetNegativeAlignmentLength(), 50u);

  PairedEndMappingInMemory flipped;
  flipped.mapping_in_memory1.strand = kNegative;
  flipped.mapping_in_memory1.ref_start_position = 260;
  flipped.mapping_in_memory1.ref_end_position = 309;
  flipped.mapping_in_memory2.strand = kPositive;
  flipped.mapping_in_memory2.ref_start_position = 100;
  flipped.mapping_in_memory2.ref_end_position = 139;
  CT_CHECK_EQ((int)flipped.GetStrand(), 0);
  CT_CHECK_EQ(flipped.GetFragmentStartPosition(), 100u);
  CT_CHECK_EQ(flipped.GetFragmentLength(), 210);
  CT_CHECK_EQ(flipped.GetPositiveAlignmentLength(), 40u);
  CT_CHECK_EQ(flipped.GetNegativeAlignmentLength(), 50u);
}

// The defaults chromap documents (chromap -h, README) and the SIMD lane
// rule: 8 lanes of 16 bits while the band 2e+1 fits in 16 bits (e < 8),
// 4 lanes of 32 bits for e < 16, none above.
CT_TEST(mapping_parameters_defaults) {
  MappingParameters p;
  CT_CHECK_EQ(p.error_threshold, 8);
  CT_CHECK_EQ(p.match_score, 1);
  CT_CHECK_EQ(p.mismatch_penalty, 4);
  CT_CHECK(p.gap_open_penalties == std::vector<int>({6, 6}));
  CT_CHECK(p.gap_extension_penalties == std::vector<int>({1, 1}));
  CT_CHECK_EQ(p.min_num_seeds_required_for_mapping, 2);
  CT_CHECK(p.max_seed_frequencies == std::vector<int>({500, 1000}));
  CT_CHECK_EQ(p.cache_update_param, 0.01);
  CT_CHECK_EQ(p.cache_size, 4000003);
  CT_CHECK(!p.debug_cache);
  CT_CHECK_EQ(p.frip_est_params, std::string("-1.0996;4.2391;3.0164e-05;-2.1087e-04;-5.5825e-05"));
  CT_CHECK(p.output_num_uniq_cache_slots);
  CT_CHECK_EQ(p.k_for_minhash, 250);
  CT_CHECK_EQ(p.max_num_best_mappings, 1);
  CT_CHECK_EQ(p.max_insert_size, 1000);
  CT_CHECK_EQ((int)p.mapq_threshold, 30);
  CT_CHECK_EQ(p.num_threads, 1);
  CT_CHECK_EQ(p.min_read_length, 30);
  CT_CHECK_EQ(p.barcode_correction_error_threshold, 1);
  CT_CHECK_EQ(p.barcode_correction_probability_threshold, 0.9);
  CT_CHECK_EQ(p.multi_mapping_allocation_distance, 0);
  CT_CHECK_EQ(p.multi_mapping_allocation_seed, 11);
  CT_CHECK_EQ(p.drop_repetitive_reads, 500000);
  CT_CHECK(!p.trim_adapters);
  CT_CHECK(!p.remove_pcr_duplicates);
  CT_CHECK(p.remove_pcr_duplicates_at_bulk_level);
  CT_CHECK(p.is_bulk_data);
  CT_CHECK(!p.allocate_multi_mappings);
  CT_CHECK(!p.Tn5_shift);
  CT_CHECK(!p.split_alignment);
  CT_CHECK(p.mapping_output_format == MAPPINGFORMAT_BED);
  CT_CHECK(!p.low_memory_mode);
  CT_CHECK(!p.skip_barcode_check);
  CT_CHECK(!p.output_mappings_not_in_whitelist);
  const int thresholds[] = {0, 1, 7, 8, 15, 16, 30};
  const int lanes[] = {8, 8, 8, 4, 4, 0, 0};
  for (int i = 0; i < 7; ++i) {
    p.error_threshold = thresholds[i];
    CT_CHECK_EQ(p.GetNumVPULanes(), lanes[i]);
  }
}

// Counters follow the vectors; PrepareForMappingNextRead leaves no
// minimizer, candidate or draft mapping and resets the repetitive seed
// length. (Which scratch vectors it clears is not asserted.)
CT_TEST(mapping_metadata_counters_and_reset) {
  TestMetadata m;
  m.PrepareForMappingNextRead(16);
  CT_CHECK_EQ(m.GetNumMinimizers(), (size_t)0);
  CT_CHECK_EQ(m.GetNumCandidates(), (size_t)0);
  m.minimizers_.emplace_back(1, 2);
  m.minimizers_.emplace_back(3, 4);
  m.minimizers_.emplace_back(5, 6);
  m.positive_candidates_.push_back(MakeCandidate(0, 10, 1));
  m.positive_candidates_.push_back(MakeCandidate(0, 20, 1));
  m.negative_candidates_.push_back(MakeCandidate(1, 30, 1));
  m.positive_mappings_.emplace_back(0, 5);
  m.negative_mappings_.emplace_back(1, 6);
  m.negative_mappings_.emplace_back(2, 7);
  m.repetitive_seed_length_ = 33;
  CT_CHECK_EQ(m.GetNumMinimizers(), (size_t)3);
  CT_CHECK_EQ(m.GetNumPositiveCandidates(), (size_t)2);
  CT_CHECK_EQ(m.GetNumNegativeCandidates(), (size_t)1);
  CT_CHECK_EQ(m.GetNumCandidates(), (size_t)3);
  CT_CHECK_EQ(m.GetNumDraftMappings(), (size_t)3);
  m.SetMinNumErrors(1);
  m.SetSecondMinNumErrors(4);
  m.SetNumBestMappings(2);
  m.SetNumSecondBestMappings(5);
  CT_CHECK_EQ(m.GetMinNumErrors(), 1);
  CT_CHECK_EQ(m.GetSecondMinNumErrors(), 4);
  CT_CHECK_EQ(m.GetNumBestMappings(), 2);
  CT_CHECK_EQ(m.GetNumSecondBestMappings(), 5);
  m.PrepareForMappingNextRead(16);
  CT_CHECK_EQ(m.GetNumMinimizers(), (size_t)0);
  CT_CHECK_EQ(m.GetNumCandidates(), (size_t)0);
  CT_CHECK_EQ(m.GetNumDraftMappings(), (size_t)0);
  CT_CHECK_EQ(m.repetitive_seed_length_, 0u);
}

// SortCandidates orders by count descending then position;
// SortMappingsByPositions by position; MoveCandidiatesToBuffer moves the
// candidates into the buffers and leaves the candidate vectors empty;
// UpdateCandidates runs the callback on both strands.
CT_TEST(mapping_metadata_sort_move_update) {
  TestMetadata m;
  m.PrepareForMappingNextRead(8);
  m.positive_candidates_ = {MakeCandidate(0, 50, 1), MakeCandidate(0, 10, 3), MakeCandidate(0, 30, 3)};
  m.negative_candidates_ = {MakeCandidate(2, 5, 2), MakeCandidate(1, 5, 2)};
  m.SortCandidates();
  CT_CHECK_EQ(m.positive_candidates_[0].position, (uint64_t)10);
  CT_CHECK_EQ(m.positive_candidates_[1].position, (uint64_t)30);
  CT_CHECK_EQ(m.positive_candidates_[2].position, (uint64_t)50);
  CT_CHECK_EQ(m.negative_candidates_[0].GetReferenceSequenceIndex(), 1u);
  m.positive_mappings_ = {DraftMapping(0, 90), DraftMapping(1, 20), DraftMapping(0, 55)};
  m.negative_mappings_ = {DraftMapping(3, (((uint64_t)1) << 32) | 1u), DraftMapping(2, 999)};
  m.SortMappingsByPositions();
  CT_CHECK_EQ(m.positive_mappings_[0].position, (uint64_t)20);
  CT_CHECK_EQ(m.positive_mappings_[1].position, (uint64_t)55);
  CT_CHECK_EQ(m.positive_mappings_[2].position, (uint64_t)90);
  CT_CHECK_EQ(m.negative_mappings_[0].position, (uint64_t)999);
  CT_CHECK_EQ(m.negative_mappings_[0].GetNumErrors(), 2);
  m.MoveCandidiatesToBuffer();
  CT_CHECK_EQ(m.GetNumCandidates(), (size_t)0);
  CT_CHECK_EQ(m.positive_candidates_buffer_.size(), (size_t)3);
  CT_CHECK_EQ(m.negative_candidates_buffer_.size(), (size_t)2);
  CT_CHECK_EQ(m.positive_candidates_buffer_[0].position, (uint64_t)10);
  m.UpdateCandidates(AddOneCandidate);
  CT_CHECK_EQ(m.GetNumPositiveCandidates(), (size_t)1);
  CT_CHECK_EQ(m.GetNumNegativeCandidates(), (size_t)1);
  CT_CHECK_EQ(m.negative_candidates_[0].position, (((uint64_t)1) << 32) | 7u);
}

// PrintCandidates writes "<strand> <index> <rid> <position> <count>" per
// candidate, positive strand first.
CT_TEST(mapping_metadata_print_candidates) {
  TestMetadata m;
  m.PrepareForMappingNextRead(4);
  m.positive_candidates_ = {MakeCandidate(0, 15, 4), MakeCandidate(2, 8, 1)};
  m.negative_candidates_ = {MakeCandidate(1, 300, 2)};
  FILE *fp = tmpfile();
  CT_CHECK(fp != NULL);
  if (fp == NULL) return;
  m.PrintCandidates(fp);
  rewind(fp);
  std::string text;
  char buffer[256];
  size_t n;
  while ((n = fread(buffer, 1, sizeof(buffer), fp)) > 0) text.append(buffer, n);
  fclose(fp);
  CT_CHECK_EQ(text, std::string("+ 0 0 15 4\n+ 1 2 8 1\n- 0 1 300 2\n"));
}

// GetBestMappings returns one vector per strand combination of (read 1,
// read 2): (+,-) F1R2, (-,+) F2R1, (+,+) F1F2, (-,-) R1R2, the same vector
// through the const and non-const overloads; PreparedForMappingNextReadPair
// empties all four and both mates' metadata.
CT_TEST(paired_metadata_best_mappings_by_strand) {
  TestPairedMetadata pe;
  pe.PreparedForMappingNextReadPair(8);
  const Strand strands[2] = {kPositive, kNegative};
  uint32_t tag = 1;
  for (int a = 0; a < 2; ++a) {
    for (int b = 0; b < 2; ++b) {
      pe.GetBestMappings(strands[a], strands[b]).emplace_back(tag, tag + 100);
      ++tag;
    }
  }
  const TestPairedMetadata &cpe = pe;
  tag = 1;
  for (int a = 0; a < 2; ++a) {
    for (int b = 0; b < 2; ++b) {
      const std::vector<std::pair<uint32_t, uint32_t>> &v = cpe.GetBestMappings(strands[a], strands[b]);
      CT_CHECK_EQ(v.size(), (size_t)1);
      if (!v.empty()) {
        CT_CHECK_EQ(v[0].first, tag);
        CT_CHECK_EQ(v[0].second, tag + 100);
      }
      CT_CHECK(&v == &pe.GetBestMappings(strands[a], strands[b]));
      ++tag;
    }
  }
  CT_CHECK(&cpe.GetBestMappings(kPositive, kNegative) != &cpe.GetBestMappings(kNegative, kPositive));
  CT_CHECK(&cpe.GetBestMappings(kPositive, kPositive) != &cpe.GetBestMappings(kNegative, kNegative));
  pe.SetMinSumErrors(3);
  pe.SetSecondMinSumErrors(5);
  pe.SetNumBestMappings(1);
  pe.SetNumSecondBestMappings(2);
  CT_CHECK_EQ(pe.GetMinSumErrors(), 3);
  CT_CHECK_EQ(pe.GetSecondMinSumErrors(), 5);
  CT_CHECK_EQ(pe.GetNumBestMappings(), 1);
  CT_CHECK_EQ(pe.GetNumSecondBestMappings(), 2);
  TestMetadata filled;
  filled.PrepareForMappingNextRead(4);
  filled.minimizers_.emplace_back(1, 2);
  filled.positive_candidates_.push_back(MakeCandidate(0, 1, 1));
  pe.mapping_metadata1_ = filled;
  pe.mapping_metadata2_ = filled;
  pe.PreparedForMappingNextReadPair(8);
  for (int a = 0; a < 2; ++a) {
    for (int b = 0; b < 2; ++b) CT_CHECK(cpe.GetBestMappings(strands[a], strands[b]).empty());
  }
  CT_CHECK_EQ(pe.mapping_metadata1_.GetNumMinimizers(), (size_t)0);
  CT_CHECK_EQ(pe.mapping_metadata2_.GetNumCandidates(), (size_t)0);
}

// The pair-level helpers act on both mates: BothEndsHaveMinimizers needs a
// minimizer in each; MoveCandidiatesToBuffer, UpdateCandidates and
// SortMappingsByPositions apply to mate 1 and mate 2.
CT_TEST(paired_metadata_acts_on_both_mates) {
  TestPairedMetadata pe;
  pe.PreparedForMappingNextReadPair(8);
  CT_CHECK(!pe.BothEndsHaveMinimizers());
  TestMetadata one;
  one.PrepareForMappingNextRead(4);
  one.minimizers_.emplace_back(1, 2);
  one.positive_candidates_.push_back(MakeCandidate(0, 11, 2));
  one.positive_mappings_ = {DraftMapping(0, 70), DraftMapping(0, 30)};
  pe.mapping_metadata1_ = one;
  CT_CHECK(!pe.BothEndsHaveMinimizers());
  pe.mapping_metadata2_ = one;
  CT_CHECK(pe.BothEndsHaveMinimizers());
  pe.UpdateCandidates(AddOneCandidate);
  CT_CHECK_EQ(pe.mapping_metadata1_.GetNumPositiveCandidates(), (size_t)2);
  CT_CHECK_EQ(pe.mapping_metadata2_.GetNumNegativeCandidates(), (size_t)1);
  pe.SortMappingsByPositions();
  TestMetadata out;
  static_cast<MappingMetadata &>(out) = pe.mapping_metadata2_;
  CT_CHECK_EQ(out.positive_mappings_[0].position, (uint64_t)30);
  CT_CHECK_EQ(out.positive_mappings_[1].position, (uint64_t)70);
  pe.MoveCandidiatesToBuffer();
  CT_CHECK_EQ(pe.mapping_metadata1_.GetNumCandidates(), (size_t)0);
  CT_CHECK_EQ(pe.mapping_metadata2_.GetNumCandidates(), (size_t)0);
  static_cast<MappingMetadata &>(out) = pe.mapping_metadata1_;
  CT_CHECK_EQ(out.positive_candidates_buffer_.size(), (size_t)2);
  CT_CHECK_EQ(out.negative_candidates_buffer_.size(), (size_t)1);
}

namespace {

int8_t g_mat[25];

const int8_t *Matrix(int match, int mismatch) {
  int k = 0;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) g_mat[k++] = i == j ? match : -mismatch;
    g_mat[k++] = 0;
  }
  for (int j = 0; j < 5; ++j) g_mat[k++] = 0;
  return g_mat;
}

std::string Cigar(int n_cigar, const uint32_t *cigar) {
  std::string s;
  for (int i = 0; i < n_cigar; ++i) {
    s += std::to_string(cigar[i] >> 4);
    s += "MIDNSHP=X"[cigar[i] & 0xf];
  }
  return s;
}

}  // namespace

// An exact 8 bp read inside a 12 bp window (2 bases either side, e = 2,
// band 5): score 8 * A = 8, CIGAR 8M at window [2, 10). The wrappers
// ksw_semi_global (one gap penalty pair for both gap types) and
// ksw_semi_global2 (separate pairs) give the same score and CIGAR as
// ksw_semi_global3, which also returns the start and end; without CIGAR
// output the score is the same.
CT_TEST(ksw_semi_global_wrappers_agree) {
  const char *window = "GGACGTTGCAGG";
  const char *read = "ACGTTGCA";
  const int8_t *mat = Matrix(1, 4);
  int n3 = 0, start = -1, end = -1;
  uint32_t *c3 = NULL;
  int s3 = ksw_semi_global3(12, window, 8, read, 5, mat, 6, 1, 6, 1, 5, &n3, &c3, &start, &end);
  CT_CHECK_EQ(s3, 8);
  CT_CHECK_EQ(Cigar(n3, c3), std::string("8M"));
  CT_CHECK_EQ(start, 2);
  CT_CHECK_EQ(end, 10);
  int n1 = 0, n2 = 0;
  uint32_t *c1 = NULL, *c2 = NULL;
  CT_CHECK_EQ(ksw_semi_global(12, window, 8, read, 5, mat, 6, 1, 5, &n1, &c1), 8);
  CT_CHECK_EQ(Cigar(n1, c1), std::string("8M"));
  CT_CHECK_EQ(ksw_semi_global2(12, window, 8, read, 5, mat, 6, 1, 6, 1, 5, &n2, &c2), 8);
  CT_CHECK_EQ(Cigar(n2, c2), std::string("8M"));
  int end_only = -1;
  CT_CHECK_EQ(ksw_semi_global3(12, window, 8, read, 5, mat, 6, 1, 6, 1, 5, NULL, NULL, NULL, &end_only), 8);
  CT_CHECK_EQ(end_only, 10);
  free(c1);
  free(c2);
  free(c3);
  // One mismatch (read base 3 T -> A): 7 * 1 - 4 = 3, still 8M at [2, 10).
  n3 = 0;
  c3 = NULL;
  CT_CHECK_EQ(ksw_semi_global3(12, window, 8, "ACGATGCA", 5, mat, 6, 1, 6, 1, 5, &n3, &c3, &start, &end), 3);
  CT_CHECK_EQ(Cigar(n3, c3), std::string("8M"));
  CT_CHECK_EQ(start, 2);
  free(c3);
}

CT_MAIN
