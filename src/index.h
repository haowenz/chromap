#ifndef INDEX_H_
#define INDEX_H_

#include <queue>
#include <string>
#include <vector>

#include "khash.h"
#include "mapping_metadata.h"
#include "sequence_batch.h"

namespace chromap {
#define KHashFunctionForIndex(a) ((a) >> 1)
#define KHashEqForIndex(a, b) ((a) >> 1 == (b) >> 1)
KHASH_INIT(k64, uint64_t, uint64_t, 1, KHashFunctionForIndex, KHashEqForIndex);

enum Direction {
  kPositive,
  kNegative,
};

struct mmHit {
  uint32_t mi;
  uint64_t position;

  bool operator<(const mmHit &h) const {
    // the inversed direction is to make a min-heap
    return position > h.position;
  }
};

class Index {
 public:
  Index(int min_num_seeds_required_for_mapping,
        const std::vector<int> &max_seed_frequencies,
        const std::string &index_file_path)
      : index_file_path_(index_file_path) {  // for read mapping
    lookup_table_ = kh_init(k64);
  }

  Index(int kmer_size, int window_size, int num_threads,
        const std::string &index_file_path)
      : kmer_size_(kmer_size),
        window_size_(window_size),
        num_threads_(num_threads),
        index_file_path_(index_file_path) {  // for index construction
    lookup_table_ = kh_init(k64);
  }

  ~Index() { Destroy(); }

  void Destroy() {
    if (lookup_table_ != NULL) {
      kh_destroy(k64, lookup_table_);
      lookup_table_ = NULL;
    }

    std::vector<uint64_t>().swap(occurrence_table_);
  }

  void Construct(uint32_t num_sequences, const SequenceBatch &reference);

  void Save() const;

  void Load();

  void Statistics(uint32_t num_sequences, const SequenceBatch &reference) const;

  void CheckIndex(uint32_t num_sequences, const SequenceBatch &reference) const;

  int CollectSeedHits(
      int max_seed_frequency, int repetitive_seed_frequency,
      const std::vector<std::pair<uint64_t, uint64_t> > &minimizers,
      uint32_t &repetitive_seed_length, std::vector<uint64_t> &positive_hits,
      std::vector<uint64_t> &negative_hits, bool use_heap) const;

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
  void GenerateMinimizerSketch(
      const SequenceBatch &sequence_batch, uint32_t sequence_index,
      std::vector<std::pair<uint64_t, uint64_t> > &minimizers) const;

  // void GenerateCandidatesOnOneDirection(
  //    int error_threshold, int num_seeds_required, uint32_t num_minimizers,
  //    std::vector<uint64_t> &hits, std::vector<Candidate> &candidates) const;

  // void GenerateCandidates(int error_threshold,
  //                        MappingMetadata &mapping_metadata) const;

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

 protected:
  int kmer_size_;
  int window_size_;
  // int min_num_seeds_required_for_mapping_;
  // Vector of size 2. The first element is the frequency threshold, and the
  // second element is the frequency threshold to run rescue. The second element
  // should always larger than the first one.
  // TODO(Haowen): add an error check.
  // std::vector<int> max_seed_frequencies_;
  // Number of threads to build the index, which is not used right now.
  int num_threads_;
  const std::string index_file_path_;
  khash_t(k64) *lookup_table_ = NULL;
  std::vector<uint64_t> occurrence_table_;
};

}  // namespace chromap

#endif  // INDEX_H_
