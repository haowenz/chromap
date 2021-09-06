#ifndef SEQUENCEBATCH_H_
#define SEQUENCEBATCH_H_

#include <unistd.h>
#include <zlib.h>

#include <iostream>
#include <string>
#include <vector>

#include "kseq.h"

namespace chromap {
class SequenceBatch {
 public:
  KSEQ_INIT(gzFile, gzread);
  SequenceBatch() {
    effective_range_[0] = 0;
    effective_range_[1] = -1;
  }
  SequenceBatch(uint32_t max_num_sequences)
      : max_num_sequences_(max_num_sequences) {
    // Construct once and use update methods when loading each batch
    sequence_batch_.reserve(max_num_sequences_);
    for (uint32_t i = 0; i < max_num_sequences_; ++i) {
      sequence_batch_.emplace_back((kseq_t *)calloc(1, sizeof(kseq_t)));
      sequence_batch_.back()->f = NULL;
    }
    negative_sequence_batch_.assign(max_num_sequences_, "");
    effective_range_[0] = 0;
    effective_range_[1] = -1;
  }
  ~SequenceBatch() {
    if (sequence_batch_.size() > 0) {
      for (uint32_t i = 0; i < sequence_batch_.size(); ++i) {
        kseq_destroy(sequence_batch_[i]);
      }
    }
  }
  inline uint32_t GetMaxBatchSize() const { return max_num_sequences_; }
  inline uint64_t GetNumBases() const { return num_bases_; }
  inline std::vector<kseq_t *> &GetSequenceBatch() { return sequence_batch_; }
  inline std::vector<std::string> &GetNegativeSequenceBatch() {
    return negative_sequence_batch_;
  }
  inline const char *GetSequenceAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->seq.s;
  }
  inline uint32_t GetSequenceLengthAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->seq.l;
  }
  inline const char *GetSequenceNameAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->name.s;
  }
  inline uint32_t GetSequenceNameLengthAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->name.l;
  }
  inline const char *GetSequenceQualAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->qual.s;
  }
  inline uint32_t GetSequenceIdAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->id;
  }
  inline const std::string &GetNegativeSequenceAt(
      uint32_t sequence_index) const {
    return negative_sequence_batch_[sequence_index];
  }
  inline int GetSequenceBatchSize() const { return sequence_batch_.size(); }
  inline void SetSeqEffectiveRange(int start, int end) {
    effective_range_[0] = start;
    effective_range_[1] = end;
  }
  //  inline char GetReverseComplementBaseOfSequenceAt(uint32_t sequence_index,
  //  uint32_t position) {
  //    kseq_t *sequence = sequence_batch_[sequence_index];
  //    return Uint8ToChar(((uint8_t)3) ^
  //    (CharToUint8((sequence->seq.s)[sequence->seq.l - position - 1])));
  //  }
  inline void PrepareNegativeSequenceAt(uint32_t sequence_index) {
    kseq_t *sequence = sequence_batch_[sequence_index];
    uint32_t sequence_length = sequence->seq.l;
    std::string &negative_sequence = negative_sequence_batch_[sequence_index];
    negative_sequence.clear();
    negative_sequence.reserve(sequence_length);
    for (uint32_t i = 0; i < sequence_length; ++i) {
      negative_sequence.push_back(Uint8ToChar(
          ((uint8_t)3) ^
          (CharToUint8((sequence->seq.s)[sequence_length - i - 1]))));
    }
  }
  inline void TrimSequenceAt(uint32_t sequence_index, int length_after_trim) {
    kseq_t *sequence = sequence_batch_[sequence_index];
    negative_sequence_batch_[sequence_index].erase(
        negative_sequence_batch_[sequence_index].begin(),
        negative_sequence_batch_[sequence_index].begin() + sequence->seq.l -
            length_after_trim);
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
  inline void CorrectBaseAt(uint32_t sequence_index, uint32_t base_position,
                            char correct_base) {
    kseq_t *sequence = sequence_batch_[sequence_index];
    sequence->seq.s[base_position] = correct_base;
  }

  inline static uint8_t CharToUint8(const char c) {
    return char_to_uint8_table_[(uint8_t)c];
  }
  inline static char Uint8ToChar(const uint8_t i) {
    return uint8_to_char_table_[i];
  }

  inline uint64_t GenerateSeedFromSequenceAt(uint32_t sequence_index,
                                             uint32_t start_position,
                                             uint32_t seed_length) const {
    const char *sequence = GetSequenceAt(sequence_index);
    uint32_t sequence_length = GetSequenceLengthAt(sequence_index);
    uint64_t mask = (((uint64_t)1) << (2 * seed_length)) - 1;
    uint64_t seed = 0;
    for (uint32_t i = 0; i < seed_length; ++i) {
      if (start_position + i < sequence_length) {
        uint8_t current_base =
            SequenceBatch::CharToUint8(sequence[i + start_position]);
        if (current_base < 4) {                        // not an ambiguous base
          seed = ((seed << 2) | current_base) & mask;  // forward k-mer
        } else {
          seed = (seed << 2) & mask;  // N->A
        }
      } else {
        seed = (seed << 2) & mask;  // Pad A
      }
    }
    return seed;
  }

  inline static uint64_t GenerateSeedFromSequence(const char *sequence,
                                                  uint32_t sequence_length,
                                                  uint32_t start_position,
                                                  uint32_t seed_length) {
    uint64_t mask = (((uint64_t)1) << (2 * seed_length)) - 1;
    uint64_t seed = 0;
    for (uint32_t i = 0; i < seed_length; ++i) {
      if (start_position + i < sequence_length) {
        uint8_t current_base =
            SequenceBatch::CharToUint8(sequence[i + start_position]);
        if (current_base < 4) {                        // not an ambiguous base
          seed = ((seed << 2) | current_base) & mask;  // forward k-mer
        } else {
          seed = (seed << 2) & mask;  // N->A
        }
      } else {
        seed = (seed << 2) & mask;  // Pad A
      }
    }
    return seed;
  }

  inline void ReorderSequences(const std::vector<int> &rid_rank) {
    std::vector<kseq_t *> tmp_sequence_batch_ = sequence_batch_;
    std::vector<std::string> tmp_negative_sequence_batch_ =
        negative_sequence_batch_;
    int i;
    int sequence_size = sequence_batch_.size();
    for (i = 0; i < sequence_size; ++i) {
      sequence_batch_[rid_rank[i]] = tmp_sequence_batch_[i];
    }
    if (negative_sequence_batch_.size() > 0) {
      for (i = 0; i < sequence_size; ++i) {
        negative_sequence_batch_[rid_rank[i]] = tmp_negative_sequence_batch_[i];
      }
    }
  }

 protected:
  uint32_t num_loaded_sequences_ = 0;
  uint32_t max_num_sequences_;
  uint64_t num_bases_;
  std::string sequence_file_path_;
  gzFile sequence_file_;
  kseq_t *sequence_kseq_;
  std::vector<kseq_t *> sequence_batch_;
  std::vector<std::string> negative_sequence_batch_;
  int effective_range_[2];  // actual range within each sequence.
  static constexpr uint8_t char_to_uint8_table_[256] = {
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
      4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
  static constexpr char uint8_to_char_table_[8] = {'A', 'C', 'G', 'T',
                                                   'N', 'N', 'N', 'N'};
  void ReplaceByEffectiveRange(kstring_t &seq) {
    if (effective_range_[0] == 0 && effective_range_[1] == -1) return;
    int i;
    int start = effective_range_[0];
    int end = effective_range_[1];
    if (effective_range_[1] == -1) end = seq.l - 1;
    for (i = 0; i < end - start + 1; ++i) {
      seq.s[i] = seq.s[start + i];
    }
    seq.s[i] = '\0';
    seq.l = end - start + 1;
  }
};
}  // namespace chromap

#endif  // SEQUENCEBATCH_H_
