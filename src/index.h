#ifndef INDEX_H_
#define INDEX_H_

#include <limits>
#include <queue>
#include <string>
#include <vector>

#include "index_parameters.h"
#include "khash.h"
#include "mapping_metadata.h"
#include "sequence_batch.h"
#include "utils.h"

namespace chromap {

#define KHashFunctionForIndex(a) ((a) >> 1)
#define KHashEqForIndex(a, b) ((a) >> 1 == (b) >> 1)
KHASH_INIT(k64, uint64_t, uint64_t, 1, KHashFunctionForIndex, KHashEqForIndex);

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

    std::vector<uint64_t>().swap(occurrence_table_);
  }

  void Construct(uint32_t num_sequences, const SequenceBatch &reference);

  void Save() const;

  void Load();

  // Output index stats.
  void Statistics(uint32_t num_sequences, const SequenceBatch &reference) const;

  // Check the index for some reference genome. Only for debug.
  void CheckIndex(uint32_t num_sequences, const SequenceBatch &reference) const;

  // Return the number of repetitive seeds.
  int CollectSeedHits(
      int max_seed_frequency, int repetitive_seed_frequency,
      const std::vector<std::pair<uint64_t, uint64_t> > &minimizers,
      uint32_t &repetitive_seed_length, std::vector<uint64_t> &positive_hits,
      std::vector<uint64_t> &negative_hits, bool use_heap) const;

  // Input a search range, for each best mate candidate, serach for minimizer
  // hits. Return the minimizer count of the best candidate if it finishes
  // normally. Or return a negative value if it stops early due to too many
  // candidates with low minimizer count.
  int CollectSeedHitsFromRepetitiveReadWithMateInfo(
      int error_threshold,
      const std::vector<std::pair<uint64_t, uint64_t> > &minimizers,
      uint32_t &repetitive_seed_length, std::vector<uint64_t> &hits,
      const std::vector<Candidate> &mate_candidates, const Direction direction,
      uint32_t search_range, int min_num_seeds_required_for_mapping,
      int max_seed_frequency0) const;

  int GetKmerSize() const { return kmer_size_; }

  int GetWindowSize() const { return window_size_; }

  uint32_t GetLookupTableSize() const { return kh_size(lookup_table_); }

  // TODO(Haowen): move this out to form a minimizer class or struct.
  // One should always reserve space for minimizers in other functions.
  void GenerateMinimizerSketch(
      const SequenceBatch &sequence_batch, uint32_t sequence_index,
      std::vector<std::pair<uint64_t, uint64_t> > &minimizers) const;

 protected:
  inline static uint64_t Hash64(uint64_t key, const uint64_t mask) {
    key = (~key + (key << 21)) & mask;  // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask;  // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask;  // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
  }

  int kmer_size_ = 0;
  int window_size_ = 0;
  // Number of threads to build the index, which is not used right now.
  int num_threads_ = 1;
  const std::string index_file_path_;
  khash_t(k64) *lookup_table_ = nullptr;
  std::vector<uint64_t> occurrence_table_;
};

}  // namespace chromap

#endif  // INDEX_H_
