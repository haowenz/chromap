#ifndef INDEX_H_
#define INDEX_H_

#include <limits>
#include <queue>
#include <string>
#include <vector>

#include "candidate_position_generating_config.h"
#include "index_parameters.h"
#include "khash.h"
#include "mapping_metadata.h"
#include "minimizer.h"
#include "sequence_batch.h"
#include "utils.h"

namespace chromap {

// Note that the max kmer size is 28 and it is saved in a unsigned 64-bit
// integer. So it can right shift by one bit is lossless and safe. The minimizer
// only occurs once when the lowest bit of the key value is set to 1.
#define KHashFunctionForIndex(a) ((a) >> 1)
#define KHashEqForIndex(a, b) ((a) >> 1 == (b) >> 1)
KHASH_INIT(/*name=*/k64, /*khkey_t=*/uint64_t, /*khval_t=*/uint64_t,
           /*kh_is_map=*/1, /*__hash_func=*/KHashFunctionForIndex,
           /*__hash_equal=*/KHashEqForIndex);

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
  int GenerateCandidatePositions(
      const CandidatePositionGeneratingConfig &generating_config,
      MappingMetadata &mapping_metadata) const;

  // Input a search range, for each best mate candidate, serach for candidate
  // positions for the read. Return the minimizer count of the best candidate if
  // it finishes normally. Or return a negative value if it stops early due to
  // too many candidates with low minimizer count.
  int GenerateCandidatePositionsFromRepetitiveReadWithMateInfoOnOneStrand(
      const Strand strand, uint32_t search_range,
      int min_num_seeds_required_for_mapping, int max_seed_frequency0,
      int error_threshold, const std::vector<Minimizer> &minimizers,
      const std::vector<Candidate> &mate_candidates,
      uint32_t &repetitive_seed_length, std::vector<uint64_t> &hits) const;

  int GetKmerSize() const { return kmer_size_; }

  int GetWindowSize() const { return window_size_; }

  uint32_t GetLookupTableSize() const { return kh_size(lookup_table_); }

 private:
  void HeapMergeCandidatePositionLists(
      const std::vector<std::vector<uint64_t>> sorted_candidate_position_lists,
      std::vector<uint64_t> &candidate_positions) const;

  uint64_t GenerateCandidatePositionForSingleSeedHit(
      uint64_t reference_seed_hit, uint64_t read_seed_hit) const;

 protected:
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
