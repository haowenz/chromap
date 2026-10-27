#ifndef SEQUENCEBATCH_H_
#define SEQUENCEBATCH_H_

#include <unistd.h>
#include <zlib.h>

#include <atomic>
#include <iostream>
#include <string>
#include <vector>

#include "kseq.h"
#include "sequence_effective_range.h"
#include "utils.h"

namespace chromap {

// Sequence files are read through a zlib buffer of kGzBufferSize bytes (zlib's
// default is 8 KiB) into a kseq stream buffer of kKseqBufferSize bytes (kseq's
// default is 16 KiB): a plain file is read 512 KiB per read() rather than
// 16 KiB, and a gzip file is inflated from 128 KiB of input at a time rather
// than 8 KiB, into 512 KiB of output rather than 16 KiB. The stream buffer is
// four times zlib's: gzread copies a request smaller than twice zlib's buffer
// through zlib's own output buffer, and reads or inflates a larger one straight
// into the stream buffer.
constexpr unsigned kGzBufferSize = 128 * 1024;
constexpr int kKseqBufferSize = 4 * kGzBufferSize;

class SequenceBatch {
 public:
  KSEQ_INIT(gzFile, gzread);

  // When 'max_num_sequences' is not specified. This batch can be used to load
  // any number of sequences with a positive full effective range.
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
    ReleaseArena();
    if (sequence_batch_.size() > 0) {
      for (uint32_t i = 0; i < sequence_batch_.size(); ++i) {
        kseq_destroy(sequence_batch_[i]);
      }
    }
  }

  inline uint64_t GetNumSequences() const { return num_loaded_sequences_; }

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

  inline bool IsNInSequenceAt(uint32_t sequence_index) {
    const int l = sequence_batch_[sequence_index]->seq.l;
    const char *s = sequence_batch_[sequence_index]->seq.s;
    for (int i = 0 ; i < l ; ++i)
      if (s[i] == 'N')
        return true;
    return false;
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

  // Opens the file and reads it through the buffers above, or, with
  // default_buffers, through zlib's and kseq's default buffers. Index
  // construction uses the defaults: the index file holds the uninitialized
  // keys and values of empty hash buckets (kh_save), whose bytes depend on the
  // heap's layout, and so on the buffers allocated before the index is built.
  void InitializeLoading(const std::string &sequence_file_path,
                         bool default_buffers = false);

  void FinalizeLoading();

  // The func should never override other sequences rather than the last, which
  // means 'sequence_index' cannot be smaller than 'num_loaded_sequences_' - 1.
  // Return true when reaching the end of the file.
  bool LoadOneSequenceAndSaveAt(uint32_t sequence_index);

  // Return the number of sequences loaded into the batch and return 0 if there
  // is no more sequences. This func now is only used to load barcodes.
  uint32_t LoadBatch();

  // As LoadBatch(), but writes its log line to 'log'.
  uint32_t LoadBatch(std::ostream &log);

  // Load all sequences in a file. This function should only be used to load
  // reference. And once the reference is loaded, the batch should never be
  // updated. This func is slow when there are large number of sequences.
  // When the file was opened with the larger buffers (not default_buffers),
  // the sequences are read into an arena on huge pages (ReserveArena), which
  // the batch frees as a whole.
  void LoadAllSequences();

  // As LoadAllSequences(), but writes its log line to 'log' and, when the file
  // ends inside a record, returns false instead of exiting, so that it can run
  // on a thread other than the main one. The caller then exits with
  // IncompleteFileMessage(). Given num_helper_threads, a plain FASTA file read
  // into the arena is read on this thread and on as many more as
  // *num_helper_threads allows when it is checked (LoadFastaInParallel), so
  // another thread may raise it while the file loads.
  bool LoadAllSequences(std::ostream &log,
                        const std::atomic<int> *num_helper_threads = nullptr);

  static const char *IncompleteFileMessage() {
    return "Didn't reach the end of sequence file, which might be corrupted!";
  }

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
  // When 'is_seq' is set to true, this func will complement the base when
  // necessary. Otherwise, it will just reverse the sequence.
  void ReplaceByEffectiveRange(kstring_t &seq, bool is_seq);

  // On Linux, reserves the address space of the arena LoadAllSequences reads
  // the sequences into, its first half for their bases and its second for the
  // quality strings of a FASTQ file, and asks for huge pages over it
  // (MADV_HUGEPAGE): the reference's 3.2 GB (hg38) is then first touched in
  // 2 MiB pages rather than 4 KiB ones, read by the alignments with fewer TLB
  // misses and unmapped at once. Leaves arena_ null where that is not
  // available or the reservation fails: the sequences are then allocated one
  // by one, by kseq_read.
  void ReserveArena();

  // Points the loader's sequence and quality buffers at the free space of
  // the arena's two halves, or, when there is no arena or too little of it
  // left, at no buffer, which kseq_read then allocates.
  void PointLoaderAtArena();

  // Takes the loader's buffers that are in the arena out of it, leaving them
  // empty.
  void ReleaseLoaderFromArena();

  // Whether 's' lies in the arena.
  bool IsInArena(const char *s) const {
    return arena_ != nullptr && s >= arena_ && s < arena_ + arena_size_;
  }

  // Unmaps the arena, leaving no loaded sequence pointing into it.
  void ReleaseArena();

  // Loads every sequence of an uncompressed FASTA file into the arena as the
  // kseq_read loop of LoadAllSequences would, with the same records, names,
  // ids, sequence bytes and places in the arena, reading the file with
  // pread on this thread and on helper threads (see LoadAllSequences): a first
  // pass counts the '\n's of each piece of the file and finds the headers, a
  // second copies each piece's bases, but the '\n's, to their place. Returns
  // false, having loaded nothing, where the result could differ from
  // kseq_read's or the file cannot be read so: without num_helper_threads,
  // without the arena (on other systems than Linux), when the file is not a
  // regular file, does not start with '>' (gzip, FASTQ), holds a '\r', '\v'
  // or '\f', has a line starting with '+' or '@', or a sequence longer than
  // an int.
  bool LoadFastaInParallel(const std::atomic<int> *num_helper_threads);

  // This is the accumulated number of sequences that have ever been loaded into
  // the batch. It is useful for tracking read ids.
  uint32_t total_num_loaded_sequences_ = 0;

  // This is the number of sequences loaded into the current batch.
  uint32_t num_loaded_sequences_ = 0;

  // This is the number of bases loaded into the current batch. It is only
  // populated for the reference.
  uint64_t num_bases_ = 0;

  // This is the max number of sequences that can be loaded into the batch. It
  // is set to 0 when there is no such restriction.
  uint32_t max_num_sequences_ = 0;

  gzFile sequence_file_;
  kseq_t *sequence_kseq_ = nullptr;
  std::vector<kseq_t *> sequence_batch_;

  // Set by InitializeLoading when the file is read through the larger
  // buffers: LoadAllSequences then reads the sequences into the arena.
  bool use_arena_ = false;
  // The arena: arena_size_ bytes of reserved address space, or null. The
  // sequences in it are not freed one by one (kseq_destroy) but unmapped with
  // it, so a batch holding them is never swapped with another
  // (SwapSequenceBatch).
  char *arena_ = nullptr;
  size_t arena_size_ = 0;
  // The bytes the loaded bases take in the arena's first half and the loaded
  // quality strings in its second.
  size_t arena_bases_used_ = 0;
  size_t arena_quals_used_ = 0;

  // TODO: avoid constructing the negative sequence batch.
  std::vector<std::string> negative_sequence_batch_;

  // Actual range within each sequence.
  const SequenceEffectiveRange effective_range_;

  // The path InitializeLoading opened, unless with default_buffers.
  std::string sequence_file_path_;
};

}  // namespace chromap

#endif  // SEQUENCEBATCH_H_
