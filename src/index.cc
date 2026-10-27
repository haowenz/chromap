#include "index.h"

#include <assert.h>

#include <algorithm>
#include <atomic>
#include <cstring>
#include <fstream>
#include <iostream>
#include <thread>

#include "minimizer_generator.h"

namespace chromap {

namespace {

// A run of bytes of a file and the memory it is read into.
struct FileRange {
  uint64_t file_offset;
  char *destination;
  uint64_t num_bytes;
};

// Reads every range of the file at 'file_path' into its destination on up to
// 'num_threads' threads. Each thread opens the file and takes the next piece of
// at most 64 MiB, in file order, until none is left, so the first touch of the
// destination pages is spread over the threads too. Returns whether every byte
// was read.
bool ReadFileRanges(const std::string &file_path,
                    const std::vector<FileRange> &ranges, int num_threads) {
  const uint64_t kMaxPieceBytes = static_cast<uint64_t>(64) << 20;
  std::vector<FileRange> pieces;
  for (const FileRange &range : ranges) {
    for (uint64_t start = 0; start < range.num_bytes; start += kMaxPieceBytes) {
      pieces.push_back({range.file_offset + start, range.destination + start,
                        std::min(kMaxPieceBytes, range.num_bytes - start)});
    }
  }

  std::atomic<size_t> next_piece(0);
  std::atomic<bool> is_complete(true);
  auto read_pieces = [&file_path, &pieces, &next_piece, &is_complete]() {
    std::ifstream file(file_path, std::ios::in | std::ios::binary);
    for (size_t pi = next_piece++; pi < pieces.size(); pi = next_piece++) {
      const FileRange &piece = pieces[pi];
      file.seekg(static_cast<std::streamoff>(piece.file_offset));
      file.read(piece.destination,
                static_cast<std::streamsize>(piece.num_bytes));
      if (!file) {
        is_complete = false;
        return;
      }
    }
  };

  const size_t num_workers = std::max<size_t>(
      1, std::min<size_t>(std::max(num_threads, 1), pieces.size()));
  std::vector<std::thread> workers;
  for (size_t ti = 1; ti < num_workers; ++ti) {
    workers.emplace_back(read_pieces);
  }
  read_pieces();
  for (std::thread &worker : workers) {
    worker.join();
  }
  return is_complete;
}

// Merges the sorted runs [a, a_end) and [b, b_end) into 'out', from both ends
// at once: the front takes the smaller head (the head of a on a tie) and the
// back the larger tail (the tail of b on a tie), so the two halves are
// disjoint and each step's selection is a conditional move, not a branch.
inline void MergeTwoSortedRuns(const uint64_t *a, const uint64_t *a_end,
                               const uint64_t *b, const uint64_t *b_end,
                               uint64_t *out) {
  const size_t total = (a_end - a) + (b_end - b);
  // The back pointers point one past the next tail to take.
  const uint64_t *a_front = a, *b_front = b;
  const uint64_t *a_back = a_end, *b_back = b_end;
  uint64_t *out_front = out, *out_back = out + total;
  for (size_t step = total / 2; step > 0; --step) {
    if (a_front < a_end && b_front < b_end) {
      const uint64_t x = *a_front, y = *b_front;
      const bool take_b = y < x;
      *out_front++ = take_b ? y : x;
      a_front += !take_b;
      b_front += take_b;
    } else if (a_front < a_end) {
      *out_front++ = *a_front++;
    } else {
      *out_front++ = *b_front++;
    }
    if (a_back > a && b_back > b) {
      const uint64_t x = a_back[-1], y = b_back[-1];
      const bool take_a = y < x;
      *--out_back = take_a ? x : y;
      a_back -= take_a;
      b_back -= !take_a;
    } else if (a_back > a) {
      *--out_back = *--a_back;
    } else {
      *--out_back = *--b_back;
    }
  }
  if (total & 1) {
    // One element is left. A head that the back has already taken is not
    // smaller than it, so the smaller of the two heads is its value.
    if (a_front < a_end && b_front < b_end) {
      *out_front = *b_front < *a_front ? *b_front : *a_front;
    } else if (a_front < a_end) {
      *out_front = *a_front;
    } else {
      *out_front = *b_front;
    }
  }
}

// The most ascending runs SortHits merges; past this it calls std::sort.
const uint32_t kMaxHitRuns = 128;

// Merges the 'num_runs' (2 to kMaxHitRuns) non-descending runs of 'hits' that
// start at run_starts[0] = 0 < run_starts[1] < ... pass by pass through a
// per-thread scratch buffer, about log2(runs) passes over the hits. Overwrites
// 'run_starts', which has room for num_runs + 1 entries.
inline void MergeHitRuns(std::vector<uint64_t> &hits, uint32_t *run_starts,
                         uint32_t num_runs) {
  const uint32_t num_hits = hits.size();
  uint64_t *data = hits.data();
  run_starts[num_runs] = num_hits;

  static thread_local std::vector<uint64_t> scratch;
  if (scratch.size() < num_hits) {
    scratch.resize(num_hits);
  }
  uint64_t *source = data;
  uint64_t *destination = scratch.data();
  while (num_runs > 1) {
    uint32_t num_merged_runs = 0;
    uint32_t ri = 0;
    for (; ri + 1 < num_runs; ri += 2) {
      const uint32_t start = run_starts[ri];
      const uint32_t middle = run_starts[ri + 1];
      const uint32_t end = run_starts[ri + 2];
      MergeTwoSortedRuns(source + start, source + middle, source + middle,
                         source + end, destination + start);
      run_starts[num_merged_runs++] = start;
    }
    if (ri < num_runs) {
      const uint32_t start = run_starts[ri];
      const uint32_t end = run_starts[ri + 1];
      std::memcpy(destination + start, source + start,
                  (end - start) * sizeof(uint64_t));
      run_starts[num_merged_runs++] = start;
    }
    run_starts[num_merged_runs] = num_hits;
    num_runs = num_merged_runs;
    std::swap(source, destination);
  }
  if (source != data) {
    std::memcpy(data, source, num_hits * sizeof(uint64_t));
  }
}

// Sorts 'hits' in ascending order. The hits arrive as a few ascending runs,
// one per minimizer and strand (a minimizer's occurrences are sorted by
// reference position), so this finds the maximal non-descending runs and
// merges them with MergeHitRuns instead of std::sort's n log n compares. Plain
// uint64 keys, so the result equals std::sort's. Falls back to std::sort when
// there are more than kMaxHitRuns runs.
void SortHits(std::vector<uint64_t> &hits) {
  const uint32_t num_hits = hits.size();
  if (num_hits < 2) {
    return;
  }
  const uint64_t *data = hits.data();
  uint32_t run_starts[kMaxHitRuns + 1];
  uint32_t num_runs = 1;
  run_starts[0] = 0;
  for (uint32_t i = 1; i < num_hits; ++i) {
    if (data[i] < data[i - 1]) {
      if (num_runs == kMaxHitRuns) {
        std::sort(hits.begin(), hits.end());
        return;
      }
      run_starts[num_runs++] = i;
    }
  }
  if (num_runs == 1) {
    return;
  }
  MergeHitRuns(hits, run_starts, num_runs);
}

// Appends positions[0, num_positions) to 'hits' and sorts 'hits'. When
// 'are_runs_known' and 'hits' was empty, the positions are the
// 'num_runs' maximal non-descending runs starting at run_starts[0] = 0 <
// run_starts[1] < ... (run_starts has room for kMaxHitRuns + 1 entries and
// holds the first min(num_runs, kMaxHitRuns + 1) of them), so SortHits' pass
// that finds them is skipped; otherwise SortHits sorts the whole vector.
void AppendAndSortHits(const uint64_t *positions, uint32_t num_positions,
                       uint32_t *run_starts, uint32_t num_runs,
                       bool are_runs_known, std::vector<uint64_t> &hits) {
  const bool was_empty = hits.empty();
  hits.insert(hits.end(), positions, positions + num_positions);
  if (!was_empty || !are_runs_known) {
    SortHits(hits);
    return;
  }
  if (num_runs <= 1) {
    return;
  }
  if (num_runs > kMaxHitRuns) {
    std::sort(hits.begin(), hits.end());
    return;
  }
  MergeHitRuns(hits, run_starts, num_runs);
}

// The state of one repetitive minimizer's binary searches in
// GenerateCandidatePositionsFromRepetitiveReadWithMateInfoOnOneStrand: the
// search for boundary 'boundary_index', which starts at 'boundary', over the
// minimizer's occurrences [offset, offset + last].
struct OccurrenceSearch {
  uint64_t boundary;
  uint32_t offset;
  int32_t l;
  int32_t r;
  int32_t m;
  int32_t last;
  uint32_t boundary_index;
  // Where the final m of each boundary's search is stored.
  uint32_t result_start;
};

// A hint to load the cache line holding 'address'. It must be inlined: GCC 13
// at -O3 deduces a function whose only effect is __builtin_prefetch to be pure
// and deletes calls to it. Where the compiler has no prefetch builtin it does
// nothing.
#if defined(__GNUC__) || defined(__clang__)
__attribute__((always_inline)) inline void PrefetchOccurrence(
    const uint64_t *address) {
  __builtin_prefetch(address);
}
#else
inline void PrefetchOccurrence(const uint64_t * /*address*/) {}
#endif

// The middle of [l, r] that a binary search probes, (l + r) / 2, for
// 0 <= l <= r.
inline int32_t MiddleOfRange(int32_t l, int32_t r) {
  return static_cast<int32_t>(static_cast<uint32_t>(l + r) >> 1);
}

// Prefetches the two probes a binary search over [l, r] that probes m next
// may make after it: the middle of [l, m - 1] and of [m + 1, r], or m itself
// for a side that is empty. 0 <= l <= m <= r.
#if defined(__GNUC__) || defined(__clang__)
__attribute__((always_inline))
#endif
inline void PrefetchNextProbes(const uint64_t *occurrences, int32_t l,
                               int32_t m, int32_t r) {
  const int32_t left = MiddleOfRange(l, m - static_cast<int32_t>(l < m));
  const int32_t right = MiddleOfRange(m + static_cast<int32_t>(m < r), r);
  PrefetchOccurrence(&occurrences[left]);
  PrefetchOccurrence(&occurrences[right]);
}

// Prefetches the cache lines holding the 'num_occurrences' >= 1 occurrences
// that start at 'occurrences', up to eight lines (57 to 64 occurrences, by
// alignment). A hint only, assuming 64-byte lines. Eight lines mapped faster
// than two, four or 64 in a rotated sweep of the ATAC, ChIP and Hi-C cases.
#if defined(__GNUC__) || defined(__clang__)
__attribute__((always_inline))
#endif
inline void PrefetchOccurrenceListStart(const uint64_t *occurrences,
                                        uint32_t num_occurrences) {
  const uintptr_t kLineSize = 64;
  const uintptr_t kMaxLines = 8;
  const uintptr_t first_line =
      reinterpret_cast<uintptr_t>(occurrences) & ~(kLineSize - 1);
  const uintptr_t last_line =
      reinterpret_cast<uintptr_t>(occurrences + (num_occurrences - 1)) &
      ~(kLineSize - 1);
  PrefetchOccurrence(occurrences);
  for (uintptr_t line = first_line + kLineSize,
                 stop = first_line + kMaxLines * kLineSize;
       line <= last_line && line < stop; line += kLineSize) {
    PrefetchOccurrence(reinterpret_cast<const uint64_t *>(line));
  }
}

// The kind of lookup table entry a minimizer found.
enum MinimizerLookupKind : uint8_t {
  kNotInLookupTable = 0,
  kSingletonLookup = 1,
  kRepetitiveLookup = 2,
};

}  // namespace

void Index::Construct(uint32_t num_sequences, const SequenceBatch &reference) {
  const double real_start_time = GetRealTime();

  std::vector<Minimizer> minimizers;
  minimizers.reserve(reference.GetNumBases() / window_size_ * 2);
  std::cerr << "Collecting minimizers.\n";
  MinimizerGenerator minimizer_generator(kmer_size_, window_size_);
  for (uint32_t sequence_index = 0; sequence_index < num_sequences;
       ++sequence_index) {
    minimizer_generator.GenerateMinimizers(reference, sequence_index,
                                           minimizers);
  }
  std::cerr << "Collected " << minimizers.size() << " minimizers.\n";
  std::cerr << "Sorting minimizers.\n";
  std::stable_sort(minimizers.begin(), minimizers.end());
  std::cerr << "Sorted all minimizers.\n";
  const size_t num_minimizers = minimizers.size();
  assert(num_minimizers > 0);
  // TODO: check this assert!
  // Here I make sure the # minimizers is less than the limit of signed int32,
  // so that I can use int to store position later.
  assert(num_minimizers <= static_cast<size_t>(INT_MAX));

  occurrence_table_.reserve(num_minimizers);
  uint64_t previous_lookup_hash =
      GenerateHashInLookupTable(minimizers[0].GetHash());
  uint32_t num_previous_minimizer_occurrences = 0;
  uint64_t num_nonsingletons = 0;
  uint32_t num_singletons = 0;
  for (size_t mi = 0; mi <= num_minimizers; ++mi) {
    const bool is_last_iteration = mi == num_minimizers;
    const uint64_t current_lookup_hash =
        is_last_iteration ? previous_lookup_hash + 1
                          : GenerateHashInLookupTable(minimizers[mi].GetHash());

    if (current_lookup_hash != previous_lookup_hash) {
      int khash_return_code = 0;
      khiter_t khash_iterator =
          kh_put(k64, lookup_table_, previous_lookup_hash, &khash_return_code);
      assert(khash_return_code != -1 && khash_return_code != 0);

      if (num_previous_minimizer_occurrences == 1) {
        // We set the lowest bit of the key value to 1 if the minimizer only
        // occurs once. And the occurrence is directly saved in the lookup
        // table.
        kh_key(lookup_table_, khash_iterator) |= 1;
        kh_value(lookup_table_, khash_iterator) = occurrence_table_.back();
        occurrence_table_.pop_back();
        ++num_singletons;
      } else {
        kh_value(lookup_table_, khash_iterator) =
            GenerateEntryValueInLookupTable(num_nonsingletons,
                                            num_previous_minimizer_occurrences);
        num_nonsingletons += num_previous_minimizer_occurrences;
      }
      num_previous_minimizer_occurrences = 1;
    } else {
      num_previous_minimizer_occurrences++;
    }

    if (is_last_iteration) {
      break;
    }

    occurrence_table_.push_back(minimizers[mi].GetHit());
    previous_lookup_hash = current_lookup_hash;
  }
  assert(num_nonsingletons + num_singletons == num_minimizers);

  std::cerr << "Kmer size: " << kmer_size_ << ", window size: " << window_size_
            << ".\n";
  std::cerr << "Lookup table size: " << kh_size(lookup_table_)
            << ", # buckets: " << kh_n_buckets(lookup_table_)
            << ", occurrence table size: " << occurrence_table_.size()
            << ", # singletons: " << num_singletons << ".\n";
  std::cerr << "Built index successfully in " << GetRealTime() - real_start_time
            << "s.\n";
}

void Index::Save() const {
  const double real_start_time = GetRealTime();
  FILE *index_file = fopen(index_file_path_.c_str(), "wb");
  assert(index_file != nullptr);

  uint64_t num_bytes = 0;
  int err = 0;

  err = fwrite(&kmer_size_, sizeof(int), 1, index_file);
  num_bytes += sizeof(int);
  assert(err != 0);

  err = fwrite(&window_size_, sizeof(int), 1, index_file);
  num_bytes += sizeof(int);
  assert(err != 0);

  const uint32_t lookup_table_size = kh_size(lookup_table_);
  err = fwrite(&lookup_table_size, sizeof(uint32_t), 1, index_file);
  num_bytes += sizeof(uint32_t);
  assert(err != 0);

  kh_save(k64, lookup_table_, index_file);
  num_bytes += sizeof(uint64_t) * 2 * lookup_table_size;

  const uint32_t occurrence_table_size = occurrence_table_.size();
  err = fwrite(&occurrence_table_size, sizeof(uint32_t), 1, index_file);
  num_bytes += sizeof(uint32_t);
  assert(err != 0);

  if (occurrence_table_size > 0) {
    err = fwrite(occurrence_table_.data(), sizeof(uint64_t),
                 occurrence_table_size, index_file);
    num_bytes += sizeof(uint64_t) * occurrence_table_size;
    assert(err != 0);
  }

  fclose(index_file);
  // std::cerr << "Index size: " << num_bytes / (1024.0 * 1024 * 1024) << "GB,
  std::cerr << "Saved in " << GetRealTime() - real_start_time << "s.\n";
}

void Index::Load() { Load(std::cerr); }

void Index::Load(std::ostream &log, int num_threads) {
  const double real_start_time = GetRealTime();
  FILE *index_file = fopen(index_file_path_.c_str(), "rb");
  assert(index_file != nullptr);

  int err = 0;
  err = fread(&kmer_size_, sizeof(int), 1, index_file);
  assert(err != 0);

  err = fread(&window_size_, sizeof(int), 1, index_file);
  assert(err != 0);

  uint32_t lookup_table_size = 0;
  err = fread(&lookup_table_size, sizeof(uint32_t), 1, index_file);
  assert(err != 0);

  // The lookup table's counts, after which come its flags, keys and values,
  // then the occurrence count and the occurrences: every offset follows from
  // the counts, so the arrays are read by offset, on several threads.
  kh_load_header(k64, lookup_table_, index_file);
  uint64_t offset = static_cast<uint64_t>(ftell(index_file));
  fclose(index_file);

  std::vector<FileRange> ranges;
  if (lookup_table_->n_buckets > 0) {
    const uint64_t num_buckets = lookup_table_->n_buckets;
    const uint64_t num_flag_bytes =
        sizeof(khint32_t) * static_cast<uint64_t>(__ac_fsize(num_buckets));
    ranges.push_back({offset, reinterpret_cast<char *>(lookup_table_->flags),
                      num_flag_bytes});
    offset += num_flag_bytes;
    ranges.push_back({offset, reinterpret_cast<char *>(lookup_table_->keys),
                      sizeof(uint64_t) * num_buckets});
    offset += sizeof(uint64_t) * num_buckets;
    ranges.push_back({offset, reinterpret_cast<char *>(lookup_table_->vals),
                      sizeof(uint64_t) * num_buckets});
    offset += sizeof(uint64_t) * num_buckets;
  }

  uint32_t occurrence_table_size = 0;
  const bool has_occurrence_table_size = ReadFileRanges(
      index_file_path_,
      {{offset, reinterpret_cast<char *>(&occurrence_table_size),
        sizeof(uint32_t)}},
      /*num_threads=*/1);
  assert(has_occurrence_table_size);
  offset += sizeof(uint32_t);

  if (occurrence_table_size > 0) {
    occurrence_table_.resize(occurrence_table_size);
    // Huge pages, as for the lookup table's arrays: the occurrences are read
    // at random, and resize left the pages untouched.
    kmadvise_huge(occurrence_table_.data(),
                  sizeof(uint64_t) * occurrence_table_.size());
    ranges.push_back({offset,
                      reinterpret_cast<char *>(occurrence_table_.data()),
                      sizeof(uint64_t) * occurrence_table_size});
  }

  const bool is_complete =
      ReadFileRanges(index_file_path_, ranges, num_threads);
  assert(is_complete);

  log << "Kmer size: " << kmer_size_ << ", window size: " << window_size_
      << ".\n";
  log << "Lookup table size: " << kh_size(lookup_table_)
      << ", occurrence table size: " << occurrence_table_.size() << ".\n";
  log << "Loaded index successfully in " << GetRealTime() - real_start_time
      << "s.\n";
}

void Index::Statistics(uint32_t num_sequences,
                       const SequenceBatch &reference) const {
  double real_start_time = GetRealTime();
  int n = 0, n1 = 0;
  uint32_t i;
  uint64_t sum = 0, len = 0;
  fprintf(stderr, "[M::%s] kmer size: %d; skip: %d; #seq: %d\n", __func__,
          kmer_size_, window_size_, num_sequences);
  for (i = 0; i < num_sequences; ++i) {
    len += reference.GetSequenceLengthAt(i);
  }
  assert(len == reference.GetNumBases());
  if (lookup_table_) {
    n += kh_size(lookup_table_);
  }
  for (khint_t k = 0; k < kh_end(lookup_table_); ++k) {
    if (kh_exist(lookup_table_, k)) {
      sum +=
          kh_key(lookup_table_, k) & 1 ? 1 : (uint32_t)kh_val(lookup_table_, k);
      if (kh_key(lookup_table_, k) & 1) ++n1;
    }
  }
  fprintf(stderr,
          "[M::%s::%.3f] distinct minimizers: %d (%.2f%% are singletons); "
          "average occurrences: %.3lf; average spacing: %.3lf\n",
          __func__, GetRealTime() - real_start_time, n, 100.0 * n1 / n,
          (double)sum / n, (double)len / sum);
}

void Index::CheckIndex(uint32_t num_sequences,
                       const SequenceBatch &reference) const {
  std::vector<Minimizer> minimizers;
  minimizers.reserve(reference.GetNumBases() / window_size_ * 2);
  MinimizerGenerator minimizer_generator(kmer_size_, window_size_);
  for (uint32_t sequence_index = 0; sequence_index < num_sequences;
       ++sequence_index) {
    minimizer_generator.GenerateMinimizers(reference, sequence_index,
                                           minimizers);
  }
  std::cerr << "Collected " << minimizers.size() << " minimizers.\n";
  std::stable_sort(minimizers.begin(), minimizers.end());
  std::cerr << "Sorted minimizers.\n";

  uint32_t count = 0;
  for (uint32_t i = 0; i < minimizers.size(); ++i) {
    khiter_t khash_iterator = kh_get(
        k64, lookup_table_, GenerateHashInLookupTable(minimizers[i].GetHash()));
    assert(khash_iterator != kh_end(lookup_table_));
    uint64_t key = kh_key(lookup_table_, khash_iterator);
    uint64_t value = kh_value(lookup_table_, khash_iterator);
    if (key & 1) {  // singleton
      assert(minimizers[i].GetHit() == value);
      count = 0;
    } else {
      uint32_t offset = GenerateOffsetInOccurrenceTable(value);
      uint32_t num_occ = GenerateNumOccurrenceInOccurrenceTable(value);
      uint64_t value_in_index = occurrence_table_[offset + count];
      assert(value_in_index == minimizers[i].GetHit());
      ++count;
      if (count == num_occ) {
        count = 0;
      }
    }
  }
}

int Index::GenerateCandidatePositions(
    const CandidatePositionGeneratingConfig &generating_config,
    MappingMetadata &mapping_metadata) const {
  const uint32_t num_minimizers = mapping_metadata.GetNumMinimizers();
  const std::vector<Minimizer> &minimizers = mapping_metadata.minimizers_;

  std::vector<std::vector<uint64_t>> positive_candidate_position_lists;
  std::vector<std::vector<uint64_t>> negative_candidate_position_lists;
  if (generating_config.UseHeapMerge()) {
    for (uint32_t i = 0; i < num_minimizers; ++i) {
      positive_candidate_position_lists.emplace_back(std::vector<uint64_t>());
      negative_candidate_position_lists.emplace_back(std::vector<uint64_t>());
    }
  }
  bool is_candidate_position_list_sorted = true;

  mapping_metadata.positive_hits_.reserve(
      generating_config.GetMaxSeedFrequency() * 2);
  mapping_metadata.negative_hits_.reserve(
      generating_config.GetMaxSeedFrequency() * 2);

  for (uint32_t mi = 0; mi < num_minimizers; ++mi) {
    PrefetchLookupTableBucket(minimizers[mi].GetHash());
  }

  // Look up every minimizer first, in minimizer order, and prefetch the start
  // of each occurrence list copied below, so that the first misses on those
  // lists overlap instead of each stalling the copy loop in turn.
  static thread_local std::vector<uint64_t> lookup_values;
  static thread_local std::vector<uint8_t> lookup_kinds;
  lookup_values.resize(num_minimizers);
  lookup_kinds.resize(num_minimizers);
  const uint64_t *occurrences = occurrence_table_.data();
  // The number of candidate positions the copy below can generate.
  uint32_t max_num_candidate_positions = 0;
  for (uint32_t mi = 0; mi < num_minimizers; ++mi) {
    khiter_t khash_iterator =
        kh_get(k64, lookup_table_,
               GenerateHashInLookupTable(minimizers[mi].GetHash()));
    if (khash_iterator == kh_end(lookup_table_)) {
      // std::cerr << "The minimizer is not in reference!\n";
      lookup_kinds[mi] = kNotInLookupTable;
      continue;
    }

    const uint64_t lookup_key = kh_key(lookup_table_, khash_iterator);
    const uint64_t lookup_value = kh_value(lookup_table_, khash_iterator);
    lookup_values[mi] = lookup_value;
    if (IsSingletonLookupKey(lookup_key)) {
      lookup_kinds[mi] = kSingletonLookup;
      ++max_num_candidate_positions;
      continue;
    }

    lookup_kinds[mi] = kRepetitiveLookup;
    const uint32_t num_occurrences =
        GenerateNumOccurrenceInOccurrenceTable(lookup_value);
    if (num_occurrences > 0 &&
        !generating_config.IsFrequentSeed(num_occurrences)) {
      PrefetchOccurrenceListStart(
          occurrences + GenerateOffsetInOccurrenceTable(lookup_value),
          num_occurrences);
      max_num_candidate_positions += num_occurrences;
    }
  }

  // Every candidate position is written to the tails of both strands' buffers
  // and only the tail of its own strand advances, so the copy below has no
  // per-occurrence strand branch. Each buffer has room for every candidate
  // position of the call, and DefaultInitAllocator leaves that room
  // uninitialized. The positions are then appended to the hit vectors (or the
  // heap merge's per-minimizer lists) in the order the copy produced them.
  // Each buffer also has a slot before its first position, holding 0, and one
  // after its last, so the run bookkeeping below reads in bounds.
  static thread_local std::vector<uint64_t, DefaultInitAllocator<uint64_t>>
      positive_buffer;
  static thread_local std::vector<uint64_t, DefaultInitAllocator<uint64_t>>
      negative_buffer;
  if (positive_buffer.size() < max_num_candidate_positions + 2) {
    positive_buffer.resize(max_num_candidate_positions + 2);
    negative_buffer.resize(max_num_candidate_positions + 2);
  }
  positive_buffer[0] = 0;
  negative_buffer[0] = 0;
  // positive_before[i] is positive_out[i - 1], and the same for negative.
  const uint64_t *const positive_before = positive_buffer.data();
  const uint64_t *const negative_before = negative_buffer.data();
  uint64_t *const positive_out = positive_buffer.data() + 1;
  uint64_t *const negative_out = negative_buffer.data() + 1;
  uint32_t num_positive = 0;
  uint32_t num_negative = 0;
  // Nonzero above bit 31 when a candidate position of either strand, computed
  // in 64 bits, falls outside [0, 2^32), so that its 32-bit start wraps and
  // may break the ascending order of its minimizer's list.
  uint64_t has_wrapped_position = 0;
  // Each minimizer's positions on a strand form one non-descending run while
  // nothing wraps (its occurrences are sorted by reference id and position).
  // A run starts where a minimizer's first position is below the last
  // position before it, so these are the maximal runs SortHits would find.
  uint32_t positive_run_starts[kMaxHitRuns + 1];
  uint32_t negative_run_starts[kMaxHitRuns + 1];
  positive_run_starts[0] = 0;
  negative_run_starts[0] = 0;
  uint32_t num_positive_runs = 1;
  uint32_t num_negative_runs = 1;

  RepetitiveSeedStats repetitive_seed_stats;
  for (uint32_t mi = 0; mi < num_minimizers; ++mi) {
    const uint8_t lookup_kind = lookup_kinds[mi];
    if (lookup_kind == kNotInLookupTable) {
      continue;
    }

    const uint64_t lookup_value = lookup_values[mi];
    const uint64_t read_hit = minimizers[mi].GetHit();
    const uint32_t read_position = HitToSequencePosition(read_hit);
    // A candidate position is (reference id << 32) | start, where the start is
    // reference_position - read_position on the read's strand and
    // reference_position + read_position - kmer_size_ + 1 on the other, both
    // in 32-bit arithmetic (GenerateCandidatePositionFromHits).
    const uint32_t negative_start_offset = read_position - kmer_size_ + 1;
    // The same offsets in 64 bits, so a start outside [0, 2^32) shows above
    // bit 31; the low 32 bits are the 32-bit starts.
    const uint64_t positive_start_offset64 =
        static_cast<uint64_t>(-static_cast<int64_t>(read_position));
    const uint64_t negative_start_offset64 = static_cast<uint64_t>(
        static_cast<int64_t>(read_position) - kmer_size_ + 1);
    const uint32_t first_positive = num_positive;
    const uint32_t first_negative = num_negative;
    if (lookup_kind == kSingletonLookup) {
      const uint64_t reference_hit = lookup_value;
      const uint32_t same_strand =
          static_cast<uint32_t>(~(reference_hit ^ read_hit)) & 1;
      const uint32_t reference_position = HitToSequencePosition(reference_hit);
      const uint64_t reference_id_bits = (reference_hit >> 33) << 32;
      positive_out[num_positive] =
          reference_id_bits | (reference_position - read_position);
      negative_out[num_negative] =
          reference_id_bits | (reference_position + negative_start_offset);
      num_positive += same_strand;
      num_negative += same_strand ^ 1;
    } else {
      const uint32_t num_occurrences =
          GenerateNumOccurrenceInOccurrenceTable(lookup_value);
      if (!generating_config.IsFrequentSeed(num_occurrences)) {
        const uint64_t *const occurrence_list =
            occurrences + GenerateOffsetInOccurrenceTable(lookup_value);
        uint32_t positive_end = num_positive;
        uint32_t negative_end = num_negative;
        uint64_t wrapped = 0;
        for (uint32_t oi = 0; oi < num_occurrences; ++oi) {
          const uint64_t reference_hit = occurrence_list[oi];
          const uint32_t same_strand =
              static_cast<uint32_t>(~(reference_hit ^ read_hit)) & 1;
          const uint64_t reference_position =
              HitToSequencePosition(reference_hit);
          const uint64_t reference_id_bits = (reference_hit >> 33) << 32;
          // Both strands' starts are checked, whichever strand the occurrence
          // is on: a wrap the copy discards only costs SortHits' run search.
          const uint64_t positive_start =
              reference_position + positive_start_offset64;
          const uint64_t negative_start =
              reference_position + negative_start_offset64;
          positive_out[positive_end] =
              reference_id_bits | static_cast<uint32_t>(positive_start);
          negative_out[negative_end] =
              reference_id_bits | static_cast<uint32_t>(negative_start);
          wrapped |= positive_start | negative_start;
          positive_end += same_strand;
          negative_end += same_strand ^ 1;
        }
        num_positive = positive_end;
        num_negative = negative_end;
        has_wrapped_position |= wrapped;
      }

      if (generating_config.IsRepetitiveSeed(num_occurrences)) {
        UpdateRepetitiveSeedStats(read_position, repetitive_seed_stats);
      }
    }

    // Start a run at this minimizer's first position on a strand when it has
    // positions there and the first is below the one before it (the slot
    // before the buffer holds 0, below no position). The start is stored
    // either way and kept only by counting it; past kMaxHitRuns runs the
    // slot is reused, since SortHits then calls std::sort.
    const uint32_t positive_run_break =
        static_cast<uint32_t>(num_positive != first_positive) &
        static_cast<uint32_t>(positive_out[first_positive] <
                              positive_before[first_positive]);
    const uint32_t negative_run_break =
        static_cast<uint32_t>(num_negative != first_negative) &
        static_cast<uint32_t>(negative_out[first_negative] <
                              negative_before[first_negative]);
    positive_run_starts[std::min(num_positive_runs, kMaxHitRuns)] =
        first_positive;
    negative_run_starts[std::min(num_negative_runs, kMaxHitRuns)] =
        first_negative;
    num_positive_runs += positive_run_break;
    num_negative_runs += negative_run_break;

    if (generating_config.UseHeapMerge()) {
      positive_candidate_position_lists[mi].assign(
          positive_out + first_positive, positive_out + num_positive);
      negative_candidate_position_lists[mi].assign(
          negative_out + first_negative, negative_out + num_negative);
    }
  }

  if ((has_wrapped_position >> 32) != 0) {
    is_candidate_position_list_sorted = false;
  }

  if (generating_config.UseHeapMerge()) {
    // TODO: try to remove this sorting.
    if (!is_candidate_position_list_sorted) {
      for (uint32_t mi = 0; mi < num_minimizers; ++mi) {
        std::sort(positive_candidate_position_lists[mi].begin(),
                  positive_candidate_position_lists[mi].end());
      }
    }
    HeapMergeCandidatePositionLists(positive_candidate_position_lists,
                                    mapping_metadata.positive_hits_);
    HeapMergeCandidatePositionLists(negative_candidate_position_lists,
                                    mapping_metadata.negative_hits_);
  } else {
    // Where a start wrapped, a minimizer's positions may not ascend, so the
    // runs found above are not trusted.
    AppendAndSortHits(positive_out, num_positive, positive_run_starts,
                      num_positive_runs, is_candidate_position_list_sorted,
                      mapping_metadata.positive_hits_);
    AppendAndSortHits(negative_out, num_negative, negative_run_starts,
                      num_negative_runs, is_candidate_position_list_sorted,
                      mapping_metadata.negative_hits_);
  }

#ifdef LI_DEBUG
  for (uint32_t mi = 0; mi < mapping_metadata.positive_hits_.size(); ++mi)
    printf("+ %llu %d %d\n", mapping_metadata.positive_hits_[mi],
           (int)(mapping_metadata.positive_hits_[mi] >> 32), (int)(mapping_metadata.positive_hits_[mi]));

  for (uint32_t mi = 0; mi < mapping_metadata.negative_hits_.size(); ++mi)
    printf("- %llu %d %d\n", mapping_metadata.negative_hits_[mi],
           (int)(mapping_metadata.negative_hits_[mi] >> 32), (int)(mapping_metadata.negative_hits_[mi]));
#endif

  mapping_metadata.repetitive_seed_length_ =
      repetitive_seed_stats.repetitive_seed_length;
  return repetitive_seed_stats.repetitive_seed_count;
}

int Index::GenerateCandidatePositionsFromRepetitiveReadWithMateInfoOnOneStrand(
    const Strand strand, uint32_t search_range,
    int min_num_seeds_required_for_mapping, int max_seed_frequency0,
    int error_threshold, const std::vector<Minimizer> &minimizers,
    const std::vector<Candidate> &mate_candidates,
    uint32_t &repetitive_seed_length,
    std::vector<uint64_t> &candidate_positions) const {
  const uint32_t mate_candidates_size = mate_candidates.size();
  int max_minimizer_count = 0;
  int best_candidate_num = 0;
  for (uint32_t i = 0; i < mate_candidates_size; ++i) {
    int count = mate_candidates[i].count;
    if (count > max_minimizer_count) {
      max_minimizer_count = count;
      best_candidate_num = 1;
    } else if (count == max_minimizer_count) {
      ++best_candidate_num;
    }
  }

  const bool mate_has_too_many_candidates =
      best_candidate_num >= 300 ||
      mate_candidates_size > static_cast<uint32_t>(max_seed_frequency0);
  const bool mate_has_too_many_low_support_candidates =
      max_minimizer_count <= min_num_seeds_required_for_mapping &&
      best_candidate_num >= 200;
  if (mate_has_too_many_candidates ||
      mate_has_too_many_low_support_candidates) {
    return -max_minimizer_count;
  }

  // TODO: reduce the search range based on the strand.
  std::vector<std::pair<uint64_t, uint64_t>> boundaries;
  boundaries.reserve(best_candidate_num);
  for (uint32_t ci = 0; ci < mate_candidates_size; ++ci) {
    if (mate_candidates[ci].count == max_minimizer_count) {
      const uint64_t boundary_start =
          (mate_candidates[ci].position < search_range)
              ? 0
              : (mate_candidates[ci].position - search_range);
      const uint64_t boundary_end = mate_candidates[ci].position + search_range;
      boundaries.emplace_back(boundary_start, boundary_end);
    }
  }

  const uint32_t raw_boundary_size = boundaries.size();
  if (raw_boundary_size == 0) {
    return max_minimizer_count;
  }

  // Merge adjacent boundary point. Assume the candidates are sorted by
  // coordinate, and thus boundaries are also sorted.
  uint32_t boundary_size = 1;
  for (uint32_t bi = 1; bi < raw_boundary_size; ++bi) {
    if (boundaries[boundary_size - 1].second < boundaries[bi].first) {
      boundaries[boundary_size] = boundaries[bi];
      ++boundary_size;
    } else {
      boundaries[boundary_size - 1].second = boundaries[bi].second;
    }
  }
  boundaries.resize(boundary_size);

  for (uint32_t mi = 0; mi < minimizers.size(); ++mi) {
    PrefetchLookupTableBucket(minimizers[mi].GetHash());
  }

  // Look up every minimizer first, in minimizer order.
  const uint32_t num_minimizers = minimizers.size();
  static thread_local std::vector<uint64_t> lookup_values;
  static thread_local std::vector<uint8_t> lookup_kinds;
  static thread_local std::vector<OccurrenceSearch> searches;
  static thread_local std::vector<int32_t> search_results;
  lookup_values.resize(num_minimizers);
  lookup_kinds.resize(num_minimizers);
  searches.clear();
  uint32_t num_repetitive_minimizers = 0;
  for (uint32_t mi = 0; mi < num_minimizers; ++mi) {
    khiter_t khash_iterator =
        kh_get(k64, lookup_table_,
               GenerateHashInLookupTable(minimizers[mi].GetHash()));
    if (khash_iterator == kh_end(lookup_table_)) {
      // std::cerr << "The minimizer is not in reference!\n";
      lookup_kinds[mi] = kNotInLookupTable;
      continue;
    }

    const uint64_t lookup_key = kh_key(lookup_table_, khash_iterator);
    const uint64_t lookup_value = kh_value(lookup_table_, khash_iterator);
    lookup_values[mi] = lookup_value;
    if (IsSingletonLookupKey(lookup_key)) {
      lookup_kinds[mi] = kSingletonLookup;
      continue;
    }

    lookup_kinds[mi] = kRepetitiveLookup;
    const uint32_t offset = GenerateOffsetInOccurrenceTable(lookup_value);
    const uint32_t num_occurrences =
        GenerateNumOccurrenceInOccurrenceTable(lookup_value);
    const int32_t last = num_occurrences - 1;
    // With no occurrence to search, every boundary's m stays 0.
    if (last >= 0) {
      searches.push_back({boundaries[0].first, offset, /*l=*/0, /*r=*/last,
                          /*m=*/0, last, /*boundary_index=*/0,
                          num_repetitive_minimizers * boundary_size});
    }
    ++num_repetitive_minimizers;
  }

  // Use binary search to locate the coordinate near mate position, for each
  // boundary of each repetitive minimizer. A minimizer searches its boundaries
  // in order, and the search for the next boundary starts from the previous
  // one's m. The searches of different minimizers are independent, so they run
  // in lockstep: in each round every active search compares its probe and
  // prefetches both probes that may follow its next one, so that each probe is
  // prefetched two rounds before it is compared and the cache misses of the
  // searches overlap instead of running one after another. l and r are
  // updated with conditional moves, so only the end of a search is a branch.
  // Each search makes the same probes as when run alone, so its final m is
  // unchanged.
  search_results.resize(num_repetitive_minimizers * boundary_size);
  const uint64_t *occurrences = occurrence_table_.data();
  OccurrenceSearch *active_searches = searches.data();
  uint32_t num_active_searches = searches.size();
  for (uint32_t si = 0; si < num_active_searches; ++si) {
    OccurrenceSearch &search = active_searches[si];
    search.m = MiddleOfRange(search.l, search.r);
    PrefetchOccurrence(&occurrences[search.offset + search.m]);
    PrefetchNextProbes(occurrences + search.offset, search.l, search.m,
                       search.r);
  }
  while (num_active_searches > 0) {
    for (uint32_t si = 0; si < num_active_searches;) {
      OccurrenceSearch &search = active_searches[si];
      const int32_t m = search.m;
      const uint64_t candidate_position =
          GenerateCandidatePositionFromOccurrenceTableEntry(
              occurrences[search.offset + m]);
      const bool below = candidate_position < search.boundary;
      const bool above = candidate_position > search.boundary;
      // All ones where the condition holds, zero elsewhere.
      const int32_t below_mask = -static_cast<int32_t>(below);
      const int32_t above_mask = -static_cast<int32_t>(above);
      const int32_t l = ((m + 1) & below_mask) | (search.l & ~below_mask);
      const int32_t r = ((m - 1) & above_mask) | (search.r & ~above_mask);
      if ((below | above) & (l <= r)) {
        search.l = l;
        search.r = r;
        search.m = MiddleOfRange(l, r);
        PrefetchNextProbes(occurrences + search.offset, l, search.m, r);
        ++si;
        continue;
      }
      search_results[search.result_start + search.boundary_index] = m;
      ++search.boundary_index;
      if (search.boundary_index < boundary_size) {
        // For next boundary, we don't have to start from l=0.
        search.boundary = boundaries[search.boundary_index].first;
        search.l = m;
        search.r = search.last;
        search.m = MiddleOfRange(search.l, search.r);
        PrefetchOccurrence(&occurrences[search.offset + search.m]);
        PrefetchNextProbes(occurrences + search.offset, search.l, search.m,
                           search.r);
        ++si;
        continue;
      }
      active_searches[si] = active_searches[--num_active_searches];
    }
  }

  // Collect the candidate positions in minimizer order.
  RepetitiveSeedStats repetitive_seed_stats;
  uint32_t repetitive_minimizer_index = 0;
  for (uint32_t mi = 0; mi < num_minimizers; ++mi) {
    if (lookup_kinds[mi] == kNotInLookupTable) {
      continue;
    }

    const uint64_t lookup_value = lookup_values[mi];
    const uint64_t read_hit = minimizers[mi].GetHit();
    const uint32_t read_position = HitToSequencePosition(read_hit);
    if (lookup_kinds[mi] == kSingletonLookup) {
      const uint64_t candidate_position =
          GenerateCandidatePositionFromHits(lookup_value, read_hit);
      const bool on_same_strand =
          AreTwoHitsOnTheSameStrand(lookup_value, read_hit);
      if ((on_same_strand && strand == kPositive) ||
          (!on_same_strand && strand == kNegative)) {
        candidate_positions.push_back(candidate_position);
      }
      continue;
    }

    const uint32_t offset = GenerateOffsetInOccurrenceTable(lookup_value);
    const uint32_t num_occurrences =
        GenerateNumOccurrenceInOccurrenceTable(lookup_value);
    const int32_t *boundary_results =
        search_results.data() + repetitive_minimizer_index * boundary_size;
    const bool searched = static_cast<int32_t>(num_occurrences - 1) >= 0;
    ++repetitive_minimizer_index;
    for (uint32_t bi = 0; bi < boundary_size; ++bi) {
      const int32_t m = searched ? boundary_results[bi] : 0;
      for (uint32_t oi = m; oi < num_occurrences; ++oi) {
        const uint64_t reference_hit = occurrence_table_[offset + oi];
        if ((GenerateCandidatePositionFromOccurrenceTableEntry(reference_hit)) >
            boundaries[bi].second) {
          break;
        }
        const uint64_t candidate_position =
            GenerateCandidatePositionFromHits(reference_hit, read_hit);
        const bool on_same_strand =
            AreTwoHitsOnTheSameStrand(reference_hit, read_hit);
        if ((on_same_strand && strand == kPositive) ||
            (!on_same_strand && strand == kNegative)) {
          candidate_positions.push_back(candidate_position);
        }
      }
    }

    if (num_occurrences >= (uint32_t)max_seed_frequency0) {
      UpdateRepetitiveSeedStats(read_position, repetitive_seed_stats);
    }
  }

  SortHits(candidate_positions);
  repetitive_seed_length = repetitive_seed_stats.repetitive_seed_length;
  return max_minimizer_count;
}

uint64_t Index::GenerateCandidatePositionFromHits(uint64_t reference_hit,
                                                  uint64_t read_hit) const {
  const uint32_t reference_position = HitToSequencePosition(reference_hit);
  const uint32_t read_position = HitToSequencePosition(read_hit);
  // For now we can't see the reference here. So let us don't validate this
  // candidate position. Instead, we do it later some time when we check the
  // candidates.
  const uint32_t reference_start_position =
      AreTwoHitsOnTheSameStrand(reference_hit, read_hit)
          ? reference_position - read_position
          : reference_position + read_position - kmer_size_ + 1;
  const uint64_t reference_id = HitToSequenceIndex(reference_hit);
  return SequenceIndexAndPositionToCandidatePosition(reference_id,
                                                     reference_start_position);
}

void Index::UpdateRepetitiveSeedStats(uint32_t read_position,
                                      RepetitiveSeedStats &stats) const {
  if (stats.previous_repetitive_seed_position > read_position) {
    // First minimizer.
    stats.repetitive_seed_length += kmer_size_;
  } else {
    if (read_position < stats.previous_repetitive_seed_position + kmer_size_ +
                            window_size_ - 1) {
      stats.repetitive_seed_length +=
          read_position - stats.previous_repetitive_seed_position;
    } else {
      stats.repetitive_seed_length += kmer_size_;
    }
  }
  stats.previous_repetitive_seed_position = read_position;
  ++stats.repetitive_seed_count;
}

}  // namespace chromap
