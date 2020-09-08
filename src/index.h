#ifndef INDEX_H_
#define INDEX_H_

#include <string>
#include <vector>
#include <queue>

#include "khash.h"
#include "sequence_batch.h"

//#define LI_DEBUG

namespace chromap {
#define KHashFunctionForIndex(a) ((a)>>1)
#define KHashEqForIndex(a, b) ((a)>>1 == (b)>>1)
KHASH_INIT(k64, uint64_t, uint64_t, 1, KHashFunctionForIndex, KHashEqForIndex);

struct Candidate {
  uint64_t position;
  uint8_t count;
  bool operator<(const Candidate &c) const {
    if (count != c.count)		
      return count > c.count;
    else
      return position < c.position;
  }
};

class Index {
 public:
  Index(int min_num_seeds_required_for_mapping, const std::vector<int> &max_seed_frequencies, const std::string &index_file_path) : min_num_seeds_required_for_mapping_(min_num_seeds_required_for_mapping), max_seed_frequencies_(max_seed_frequencies), index_file_path_(index_file_path) { // for read mapping
    lookup_table_ = kh_init(k64);
  }
  Index(int kmer_size, int window_size, int num_threads, const std::string &index_file_path) : kmer_size_(kmer_size), window_size_(window_size), num_threads_(num_threads), index_file_path_(index_file_path) { // for index construction
    lookup_table_ = kh_init(k64);
  }
  ~Index(){
    if (lookup_table_ != NULL) {
      kh_destroy(k64, lookup_table_);
    }
  }
  void Destroy() {
    kh_destroy(k64, lookup_table_);
    lookup_table_ = NULL;
    std::vector<uint64_t>().swap(occurrence_table_); 
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
  void GenerateCandidatesOnOneDirection(int error_threshold, int num_seeds_required, std::vector<uint64_t> *hits, std::vector<Candidate> *candidates);
  void GenerateCandidates(int error_threshold, const std::vector<std::pair<uint64_t, uint64_t> > &minimizers, uint32_t *repetitive_seed_length, std::vector<uint64_t> *positive_hits, std::vector<uint64_t> *negative_hits, std::vector<Candidate> *positive_candidates, std::vector<Candidate> *negative_candidates);
  void GenerateCandidatesFromRepetitiveReadWithMateInfo(int error_threshold, const std::vector<std::pair<uint64_t, uint64_t> > &minimizers, uint32_t *repetitive_seed_length, std::vector<uint64_t> *hits, std::vector<Candidate> *candidates, std::vector<Candidate> *mate_candidates, int direction, unsigned int range);
  int CollectCandidates(int max_seed_frequency, int repetitive_seed_frequency, const std::vector<std::pair<uint64_t, uint64_t> > &minimizers, uint32_t *repetitive_seed_length, std::vector<uint64_t> *positive_hits, std::vector<uint64_t> *negative_hits, bool use_heap);
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
