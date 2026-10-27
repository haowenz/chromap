#ifndef INDEX_H_
#define INDEX_H_

#include <limits>
#include <memory>
#include <new>
#include <ostream>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#include "candidate_position_generating_config.h"
#include "index_parameters.h"
#include "index_utils.h"
#include "mapping_metadata.h"
#include "minimizer.h"
#include "sequence_batch.h"
#include "utils.h"

namespace chromap {

// An allocator whose value construction (as in std::vector::resize) leaves a
// trivial type uninitialized, so Index::Load can size the occurrence table
// without writing every page on one thread before the threads that read it
// from the file touch those pages.
template <typename T>
struct DefaultInitAllocator : std::allocator<T> {
  template <typename U>
  struct rebind {
    typedef DefaultInitAllocator<U> other;
  };

  DefaultInitAllocator() = default;

  template <typename U>
  DefaultInitAllocator(const DefaultInitAllocator<U> &) noexcept {}

  template <typename U>
  void construct(U *pointer) {
    ::new (static_cast<void *>(pointer)) U;
  }

  template <typename U, typename... Args>
  void construct(U *pointer, Args &&...args) {
    ::new (static_cast<void *>(pointer)) U(std::forward<Args>(args)...);
  }
};

class Index {
 public:
  Index() = delete;

  // For read mapping.
  Index(const std::string &index_file_path)
      : index_file_path_(index_file_path) {
    lookup_table_ = kh_init(k64);
  }

  // For index construction.
  Index(const IndexParameters &index_parameters)
      : kmer_size_(index_parameters.kmer_size),
        window_size_(index_parameters.window_size),
        num_threads_(index_parameters.num_threads),
        index_file_path_(index_parameters.index_output_file_path) {
    lookup_table_ = kh_init(k64);
  }

  ~Index() { Destroy(); }

  void Destroy() {
    if (lookup_table_ != nullptr) {
      kh_destroy(k64, lookup_table_);
      lookup_table_ = nullptr;
    }

    std::vector<uint64_t, DefaultInitAllocator<uint64_t>>().swap(
        occurrence_table_);
  }

  void Construct(uint32_t num_sequences, const SequenceBatch &reference);

  void Save() const;

  void Load();

  // As Load(), but writes its log lines to 'log' and reads the lookup table's
  // arrays and the occurrence table on up to 'num_threads' threads.
  void Load(std::ostream &log, int num_threads = 1);

  // Output index stats.
  void Statistics(uint32_t num_sequences, const SequenceBatch &reference) const;

  // Check the index for some reference genome. Only for debug.
  void CheckIndex(uint32_t num_sequences, const SequenceBatch &reference) const;

  // Return the number of repetitive seeds.
  int GenerateCandidatePositions(
      const CandidatePositionGeneratingConfig &generating_config,
      MappingMetadata &mapping_metadata) const;

  // Input a search range, for each best mate candidate, serach for candidate
  // positions for the read. Return the minimizer count of the best candidate if
  // it finishes normally. Or return a negative value if it stops early due to
  // too many candidates with low minimizer count.
  // 'strand' is the strand to generate (augment) candidates.
  int GenerateCandidatePositionsFromRepetitiveReadWithMateInfoOnOneStrand(
      const Strand strand, uint32_t search_range,
      int min_num_seeds_required_for_mapping, int max_seed_frequency0,
      int error_threshold, const std::vector<Minimizer> &minimizers,
      const std::vector<Candidate> &mate_candidates,
      uint32_t &repetitive_seed_length,
      std::vector<uint64_t> &candidate_positions) const;

  int GetKmerSize() const { return kmer_size_; }

  int GetWindowSize() const { return window_size_; }

  uint32_t GetLookupTableSize() const { return kh_size(lookup_table_); }

 private:
  // Prefetches the flags word, key and value of the lookup table bucket that
  // kh_get(k64, ...) probes first for this minimizer. Called for every
  // minimizer of a read before its lookups, so that their cache misses overlap
  // instead of running one after another. A hint only: the lookups and their
  // results are unchanged, and where the compiler has no prefetch builtin it
  // does nothing. It must be inlined: GCC 13 at -O3 deduces a function whose
  // only effect is __builtin_prefetch to be pure and deletes calls to it.
#if defined(__GNUC__) || defined(__clang__)
  __attribute__((always_inline)) void PrefetchLookupTableBucket(
      uint64_t minimizer_hash) const {
    if (lookup_table_->n_buckets == 0) {
      return;
    }
    const khint_t i =
        (khint_t)KHashFunctionForIndex(
            GenerateHashInLookupTable(minimizer_hash)) &
        (lookup_table_->n_buckets - 1);
    __builtin_prefetch(&lookup_table_->flags[i >> 4]);
    __builtin_prefetch(&lookup_table_->keys[i]);
    __builtin_prefetch(&lookup_table_->vals[i]);
  }
#else
  void PrefetchLookupTableBucket(uint64_t /*minimizer_hash*/) const {}
#endif

  uint64_t GenerateCandidatePositionFromHits(uint64_t reference_hit,
                                             uint64_t read_hit) const;

  void UpdateRepetitiveSeedStats(uint32_t read_position,
                                 RepetitiveSeedStats &stats) const;

  int kmer_size_ = 0;
  int window_size_ = 0;
  // Number of threads to build the index, which is not used right now.
  int num_threads_ = 1;
  const std::string index_file_path_;
  khash_t(k64) *lookup_table_ = nullptr;
  std::vector<uint64_t, DefaultInitAllocator<uint64_t>> occurrence_table_;
};

}  // namespace chromap

#endif  // INDEX_H_
