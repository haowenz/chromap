#include "sequence_batch.h"

#include <tuple>

#include "utils.h"

namespace chromap {

constexpr uint8_t SequenceBatch::char_to_uint8_table_[256];
constexpr char SequenceBatch::uint8_to_char_table_[8];

void SequenceBatch::InitializeLoading(const std::string &sequence_file_path) {
  sequence_file_ = gzopen(sequence_file_path.c_str(), "r");
  if (sequence_file_ == NULL) {
    ExitWithMessage("Cannot find sequence file" + sequence_file_path);
  }
  sequence_kseq_ = kseq_init(sequence_file_);
}

uint32_t SequenceBatch::LoadBatch() {
  double real_start_time = GetRealTime();
  uint32_t num_sequences = 0;
  num_bases_ = 0;
  for (uint32_t sequence_index = 0; sequence_index < max_num_sequences_;
       ++sequence_index) {
    int length = kseq_read(sequence_kseq_);
    while (length == 0) {  // Skip the sequences of length 0
      length = kseq_read(sequence_kseq_);
    }
    if (length > 0) {
      kseq_t *sequence = sequence_batch_[sequence_index];
      std::swap(sequence_kseq_->seq, sequence->seq);
      ReplaceByEffectiveRange(sequence->seq, /*is_seq=*/true);
      std::swap(sequence_kseq_->name, sequence->name);
      std::swap(sequence_kseq_->comment, sequence->comment);
      if (sequence_kseq_->qual.l != 0) {  // fastq file
        std::swap(sequence_kseq_->qual, sequence->qual);
        ReplaceByEffectiveRange(sequence->qual, /*is_seq=*/false);
      }
      sequence->id = num_loaded_sequences_;
      ++num_loaded_sequences_;
      ++num_sequences;
      num_bases_ += length;
    } else {
      if (length != -1) {
        ExitWithMessage(
            "Didn't reach the end of sequence file, which might be corrupted!");
      }
      // make sure to reach the end of file rather than meet an error
      break;
    }
  }
  if (num_sequences != 0) {
    std::cerr << "Loaded sequence batch successfully in "
              << GetRealTime() - real_start_time << "s, ";
    std::cerr << "number of sequences: " << num_sequences << ", ";
    std::cerr << "number of bases: " << num_bases_ << ".\n";
  } else {
    std::cerr << "No more sequences.\n";
  }
  return num_sequences;
}

bool SequenceBatch::LoadOneSequenceAndSaveAt(uint32_t sequence_index) {
  // double real_start_time = Chromap::GetRealTime();
  bool no_more_sequence = false;
  int length = kseq_read(sequence_kseq_);
  while (length == 0) {  // Skip the sequences of length 0
    length = kseq_read(sequence_kseq_);
  }
  if (length > 0) {
    kseq_t *sequence = sequence_batch_[sequence_index];
    std::swap(sequence_kseq_->seq, sequence->seq);
    ReplaceByEffectiveRange(sequence->seq, /*is_seq=*/true);
    std::swap(sequence_kseq_->name, sequence->name);
    std::swap(sequence_kseq_->comment, sequence->comment);
    sequence->id = num_loaded_sequences_;
    ++num_loaded_sequences_;
    if (sequence_kseq_->qual.l != 0) {  // fastq file
      std::swap(sequence_kseq_->qual, sequence->qual);
      ReplaceByEffectiveRange(sequence->qual, /*is_seq=*/false);
    }
  } else {
    if (length != -1) {
      ExitWithMessage(
          "Didn't reach the end of sequence file, which might be corrupted!");
    }
    // make sure to reach the end of file rather than meet an error
    no_more_sequence = true;
  }
  return no_more_sequence;
}

uint32_t SequenceBatch::LoadAllSequences() {
  double real_start_time = GetRealTime();
  sequence_batch_.reserve(200);
  uint32_t num_sequences = 0;
  num_bases_ = 0;
  int length = kseq_read(sequence_kseq_);
  while (length >= 0) {
    if (length == 0) {  // Skip the sequences of length 0
      continue;
    } else if (length > 0) {
      sequence_batch_.emplace_back((kseq_t *)calloc(1, sizeof(kseq_t)));
      kseq_t *sequence = sequence_batch_.back();
      std::swap(sequence_kseq_->seq, sequence->seq);
      ReplaceByEffectiveRange(sequence->seq, /*is_seq=*/true);
      std::swap(sequence_kseq_->name, sequence->name);
      std::swap(sequence_kseq_->comment, sequence->comment);
      if (sequence_kseq_->qual.l != 0) {  // fastq file
        std::swap(sequence_kseq_->qual, sequence->qual);
        ReplaceByEffectiveRange(sequence->qual, /*is_seq=*/false);
      }
      sequence->id = num_loaded_sequences_;
      ++num_loaded_sequences_;
      ++num_sequences;
      num_bases_ += length;
    } else {
      if (length != -1) {
        ExitWithMessage(
            "Didn't reach the end of sequence file, which might be corrupted!");
      }
      // make sure to reach the end of file rather than meet an error
      break;
    }
    length = kseq_read(sequence_kseq_);
  }
  std::cerr << "Loaded all sequences successfully in "
            << GetRealTime() - real_start_time << "s, ";
  std::cerr << "number of sequences: " << num_sequences << ", ";
  std::cerr << "number of bases: " << num_bases_ << ".\n";
  return num_sequences;
}

void SequenceBatch::FinalizeLoading() {
  kseq_destroy(sequence_kseq_);
  gzclose(sequence_file_);
}

void SequenceBatch::ReplaceByEffectiveRange(kstring_t &seq, bool is_seq) {
  if (effective_range_[0] == 0 && effective_range_[1] == -1 &&
      effective_range_[2] == 1) {
    return;
  }
  int i, j;
  int start = effective_range_[0];
  int end = effective_range_[1];
  if (effective_range_[1] == -1) end = seq.l - 1;
  for (i = 0; i < end - start + 1; ++i) {
    seq.s[i] = seq.s[start + i];
  }
  seq.s[i] = '\0';
  seq.l = end - start + 1;
  if (effective_range_[2] == -1) {
    if (is_seq) {
      for (i = 0; i < (int)seq.l; ++i) {
        seq.s[i] = Uint8ToChar(((uint8_t)3) ^ (CharToUint8(seq.s[i])));
      }
    }
    for (i = 0, j = seq.l - 1; i < j; ++i, --j) {
      char tmp = seq.s[i];
      seq.s[i] = seq.s[j];
      seq.s[j] = tmp;
    }
  }
}

}  // namespace chromap
