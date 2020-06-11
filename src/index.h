#ifndef INDEX_H_
#define INDEX_H_

#include <string>
#include <vector>

#include "khash.h"
#include "sequence_batch.h"

namespace chromap {
#define KHashFunctionForIndex(a) ((a)>>1)
#define KHashEqForIndex(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(k64, uint64_t, uint64_t, 1, KHashFunctionForIndex, KHashEqForIndex);

class Index {
 public:
  Index(int min_num_seeds_required_for_mapping, const std::vector<int> &max_seed_frequencies, const std::string &index_file_path) : min_num_seeds_required_for_mapping_(min_num_seeds_required_for_mapping), max_seed_frequencies_(max_seed_frequencies), index_file_path_(index_file_path) { // for read mapping
    lookup_table_ = kh_init(k64);
  }
  Index(int kmer_size, int window_size, int num_threads, const std::string &index_file_path) : kmer_size_(kmer_size), window_size_(window_size), num_threads_(num_threads), index_file_path_(index_file_path) { // for index construction
    lookup_table_ = kh_init(k64);
  }
  ~Index(){
    kh_destroy(k64, lookup_table_);
  }
  khash_t(k64) const * GetLookupTable() const {
    return lookup_table_;
  }
  int GetKmerSize() const {
    return kmer_size_;
  }
  int GetWindowSize() const {
    return window_size_;
  }
  uint32_t GetLookupTableSize() const {
    return kh_size(lookup_table_);
  }
  std::vector<uint64_t> const & GetOccurrenceTable() const {
    return occurrence_table_;
  }
  void Statistics(uint32_t num_sequences, const SequenceBatch &reference);
  void CheckIndex(uint32_t num_sequences, const SequenceBatch &reference);
  void GenerateMinimizerSketch(const SequenceBatch &sequence_batch, uint32_t sequence_index, std::vector<std::pair<uint64_t, uint64_t> > *minimizers);
  void Construct(uint32_t num_sequences, const SequenceBatch &reference);
  void Save();
  void Load();
  void GenerateCandidatesOnOneDirection(int error_threshold, std::vector<uint64_t> *hits, std::vector<uint64_t> *candidates);
  void GenerateCandidates(int error_threshold, const std::vector<std::pair<uint64_t, uint64_t> > &minimizers, std::vector<uint64_t> *positive_hits, std::vector<uint64_t> *negative_hits, std::vector<uint64_t> *positive_candidates, std::vector<uint64_t> *negative_candidates);
  void CollectCandidates(int max_seed_frequency, const std::vector<std::pair<uint64_t, uint64_t> > &minimizers, std::vector<uint64_t> *positive_hits, std::vector<uint64_t> *negative_hits);
  inline static uint64_t Hash64(uint64_t key, const uint64_t mask) {
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
  }

 protected:
  int kmer_size_;
  int window_size_;
  int min_num_seeds_required_for_mapping_;
  std::vector<int> max_seed_frequencies_;
  int num_threads_;
  std::string index_file_path_;
  khash_t(k64)* lookup_table_ = NULL;
  std::vector<uint64_t> occurrence_table_;
};
} // namespace chromap

#endif // INDEX_H_
