#ifndef SEQUENCEBATCH_H_
#define SEQUENCEBATCH_H_

#include <iostream>
#include <string>
#include <vector>
#include <unistd.h>
#include <zlib.h>

#include "kseq.h"

namespace chromap {
class SequenceBatch {
 public:
  KSEQ_INIT(gzFile, gzread);
  SequenceBatch(){}
  SequenceBatch(uint32_t max_num_sequences) : max_num_sequences_(max_num_sequences) {
    // Construct once and use update methods when loading each batch
    sequence_batch_.reserve(max_num_sequences_);
    for (uint32_t i = 0; i < max_num_sequences_; ++i) {
      sequence_batch_.emplace_back((kseq_t*)calloc(1, sizeof(kseq_t)));
    }
    negative_sequence_batch_.assign(max_num_sequences_, "");
  }
  ~SequenceBatch(){
    if (sequence_batch_.size() > 0) {
      for (uint32_t i = 0; i < sequence_batch_.size(); ++i) {
        free(sequence_batch_[i]);
      }
    }
  }
  inline uint32_t GetMaxBatchSize() const {
    return max_num_sequences_;
  }
  inline uint64_t GetNumBases() const {
    return num_bases_;
  }
  inline std::vector<kseq_t*> & GetSequenceBatch() {
    return sequence_batch_;
  }
  inline std::vector<std::string> & GetNegativeSequenceBatch() {
    return negative_sequence_batch_;
  }
  inline const char * GetSequenceAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->seq.s;
  }
  inline uint32_t GetSequenceLengthAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->seq.l;
  }
  inline const char * GetSequenceNameAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->name.s;
  }
  inline uint32_t GetSequenceNameLengthAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->name.l;
  }
  inline uint32_t GetSequenceIdAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->id;
  }
  inline const std::string & GetNegativeSequenceAt(uint32_t sequence_index) const {
    return negative_sequence_batch_[sequence_index];
  }
//  inline char GetReverseComplementBaseOfSequenceAt(uint32_t sequence_index, uint32_t position) {
//    kseq_t *sequence = sequence_batch_[sequence_index];
//    return Uint8ToChar(((uint8_t)3) ^ (CharToUint8((sequence->seq.s)[sequence->seq.l - position - 1])));
//  }
  inline void PrepareNegativeSequenceAt(uint32_t sequence_index) {
    kseq_t *sequence = sequence_batch_[sequence_index];
    uint32_t sequence_length = sequence->seq.l;
    std::string &negative_sequence = negative_sequence_batch_[sequence_index];
    negative_sequence.clear();
    negative_sequence.reserve(sequence_length);
    for (uint32_t i = 0; i < sequence_length; ++i) {
      negative_sequence.push_back(Uint8ToChar(((uint8_t)3) ^ (CharToUint8((sequence->seq.s)[sequence_length - i - 1]))));
    }
  }
  inline void TrimSequenceAt(uint32_t sequence_index, int length_after_trim) {
    kseq_t *sequence = sequence_batch_[sequence_index];
    negative_sequence_batch_[sequence_index].erase(negative_sequence_batch_[sequence_index].begin(), negative_sequence_batch_[sequence_index].begin() + sequence->seq.l - length_after_trim);
    sequence->seq.l = length_after_trim;
  }
  inline void SwapSequenceBatch(SequenceBatch &batch) {
    sequence_batch_.swap(batch.GetSequenceBatch());
    negative_sequence_batch_.swap(batch.GetNegativeSequenceBatch());
  }
  void InitializeLoading(const std::string &sequence_file_path);
  void FinalizeLoading();
  // Return the number of reads loaded into the batch
  // and return 0 if there is no more reads
  uint32_t LoadBatch();
  bool LoadOneSequenceAndSaveAt(uint32_t sequence_index);
  uint32_t LoadAllSequences();

  inline static uint8_t CharToUint8(const char c) {
    return char_to_uint8_table_[(uint8_t)c];
  }
  inline static char Uint8ToChar(const uint8_t i) {
    return uint8_to_char_table_[i];
  }

  inline uint64_t GenerateSeedFromSequenceAt(uint32_t sequence_index, uint32_t start_position, uint32_t seed_length) const {
    const char *sequence = GetSequenceAt(sequence_index);
    uint32_t sequence_length = GetSequenceLengthAt(sequence_index);
    uint64_t mask = (((uint64_t)1) << (2 * seed_length)) - 1;
    uint64_t seed = 0;
    for (uint32_t i = 0; i < seed_length; ++i) {
      if (start_position + i < sequence_length) {
        uint8_t current_base = SequenceBatch::CharToUint8(sequence[i + start_position]);
        if (current_base < 4) { // not an ambiguous base
          seed = ((seed << 2) | current_base) & mask; // forward k-mer
        } else {
          seed = (seed << 2) & mask; // N->A
        }
      } else {
        seed = (seed << 2) & mask; // Pad A
      }
    }
    return seed;
  }

 protected:
  uint32_t num_loaded_sequences_ = 0;
  uint32_t max_num_sequences_;
  uint64_t num_bases_;
  std::string sequence_file_path_;
  gzFile sequence_file_; 
  kseq_t *sequence_kseq_;
  std::vector<kseq_t*> sequence_batch_;
  std::vector<std::string> negative_sequence_batch_;
  static constexpr uint8_t char_to_uint8_table_[256] = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
  static constexpr char uint8_to_char_table_[8] = {'A', 'C', 'G', 'T', 'N', 'N', 'N', 'N'};
};
} // namespace chromap

#endif // SEQUENCEBATCH_H_
