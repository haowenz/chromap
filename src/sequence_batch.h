#ifndef SEQUENCEBATCH_H_
#define SEQUENCEBATCH_H_

#include <unistd.h>
#include <zlib.h>

#include <iostream>
#include <string>
#include <vector>

#include "kseq.h"
#include "sequence_effective_range.h"
#include "utils.h"

namespace chromap {

class SequenceBatch {
 public:
  KSEQ_INIT(gzFile, gzread);
  SequenceBatch() : effective_range_(SequenceEffectiveRange()) {}

  // Construct once and use update sequences when loading each batch.
  SequenceBatch(uint32_t max_num_sequences,
                const SequenceEffectiveRange &effective_range)
      : max_num_sequences_(max_num_sequences),
        effective_range_(effective_range) {
    sequence_batch_.reserve(max_num_sequences_);
    for (uint32_t i = 0; i < max_num_sequences_; ++i) {
      sequence_batch_.emplace_back((kseq_t *)calloc(1, sizeof(kseq_t)));
      sequence_batch_.back()->f = NULL;
    }
    negative_sequence_batch_.assign(max_num_sequences_, "");
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

  // big_endian: N_pos is in the order of sequence
  // little_endian: N_pos is in the order from the sequence right side to left,
  //                this is the order of the GenerateSeed
  // e.g: If the sequence is "ACN", big endian returns N at 2,
  //      little endian returns N at 0.
  inline void GetSequenceNsAt(uint32_t sequence_index, bool little_endian,
                              std::vector<int> &N_pos) {
    const int l = sequence_batch_[sequence_index]->seq.l;
    const char *s = sequence_batch_[sequence_index]->seq.s;
    N_pos.clear();
    if (little_endian) {
      for (int i = l - 1; i >= 0; --i) {
        if (s[i] == 'N') N_pos.push_back(l - 1 - i);
      }
    } else {
      for (int i = 0; i < l; ++i) {
        if (s[i] == 'N') N_pos.push_back(i);
      }
    }
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

    if (length_after_trim >= (int)sequence->seq.l) {
      return;
    }

    negative_sequence_batch_[sequence_index].erase(
        negative_sequence_batch_[sequence_index].begin(),
        negative_sequence_batch_[sequence_index].begin() + sequence->seq.l -
            length_after_trim);

    sequence->seq.l = length_after_trim;
    sequence->seq.s[sequence->seq.l] = '\0';
    sequence->qual.l = length_after_trim;
    sequence->qual.s[sequence->qual.l] = '\0';
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

  inline uint64_t GenerateSeedFromSequenceAt(uint32_t sequence_index,
                                             uint32_t start_position,
                                             uint32_t seed_length) const {
    const char *sequence = GetSequenceAt(sequence_index);
    const uint32_t sequence_length = GetSequenceLengthAt(sequence_index);
    return GenerateSeedFromSequence(sequence, sequence_length, start_position,
                                    seed_length);
  }

  inline void ReorderSequences(const std::vector<int> &rid_rank) {
    std::vector<kseq_t *> tmp_sequence_batch_ = sequence_batch_;
    std::vector<std::string> tmp_negative_sequence_batch_ =
        negative_sequence_batch_;
    for (size_t i = 0; i < sequence_batch_.size(); ++i) {
      sequence_batch_[rid_rank[i]] = tmp_sequence_batch_[i];
    }

    if (negative_sequence_batch_.size() > 0) {
      for (size_t i = 0; i < sequence_batch_.size(); ++i) {
        negative_sequence_batch_[rid_rank[i]] = tmp_negative_sequence_batch_[i];
      }
    }
  }

 protected:
  uint32_t num_loaded_sequences_ = 0;
  uint32_t max_num_sequences_ = 0;
  uint64_t num_bases_ = 0;
  gzFile sequence_file_;
  kseq_t *sequence_kseq_ = nullptr;
  std::vector<kseq_t *> sequence_batch_;
  std::vector<std::string> negative_sequence_batch_;

  // Actual range within each sequence.
  const SequenceEffectiveRange effective_range_;
  void ReplaceByEffectiveRange(kstring_t &seq, bool is_seq);
};

}  // namespace chromap

#endif  // SEQUENCEBATCH_H_
