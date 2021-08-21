#include "sequence_batch.h"

#include <tuple>

#include "chromap.h"

namespace chromap {
constexpr uint8_t SequenceBatch::char_to_uint8_table_[256];
constexpr char SequenceBatch::uint8_to_char_table_[8];

void SequenceBatch::InitializeLoading(const std::string &sequence_file_path) {
  sequence_file_path_ = sequence_file_path;
  sequence_file_ = gzopen(sequence_file_path_.c_str(), "r");
  if (sequence_file_ == NULL) {
    Chromap<>::ExitWithMessage("Cannot find sequence file" + sequence_file_path);
  }
  sequence_kseq_ = kseq_init(sequence_file_);
}

uint32_t SequenceBatch::LoadBatch() {
  double real_start_time = Chromap<>::GetRealTime();
  uint32_t num_sequences = 0;
  num_bases_ = 0;
  for (uint32_t sequence_index = 0; sequence_index < max_num_sequences_; ++sequence_index) { 
    int length = kseq_read(sequence_kseq_);
    while (length == 0) { // Skip the sequences of length 0
      length = kseq_read(sequence_kseq_);
    }
    if (length > 0) {
      kseq_t *sequence = sequence_batch_[sequence_index];
      std::swap(sequence_kseq_->seq, sequence->seq);
      ReplaceByEffectiveRange(sequence->seq);
      std::swap(sequence_kseq_->name, sequence->name);
      std::swap(sequence_kseq_->comment, sequence->comment);
      if (sequence_kseq_->qual.l != 0) { // fastq file
        std::swap(sequence_kseq_->qual, sequence->qual);
        ReplaceByEffectiveRange(sequence->qual);
      }
      sequence->id = num_loaded_sequences_;
      ++num_loaded_sequences_;
      ++num_sequences;
      num_bases_ += length;
    } else {
      if (length != -1) {
        Chromap<>::ExitWithMessage("Didn't reach the end of sequence file, which might be corrupted!");
      }
      // make sure to reach the end of file rather than meet an error
      break;
    }
  }
  if (num_sequences != 0) {
    std::cerr << "Loaded sequence batch successfully in " << Chromap<>::GetRealTime() - real_start_time << "s, ";
    std::cerr << "number of sequences: " << num_sequences << ", ";
    std::cerr << "number of bases: " << num_bases_ << ".\n";
  } else {
    std::cerr << "No more sequences.\n";
  }
  return num_sequences;
}

bool SequenceBatch::LoadOneSequenceAndSaveAt(uint32_t sequence_index) {
  //double real_start_time = Chromap::GetRealTime();
  bool no_more_sequence = false;
  int length = kseq_read(sequence_kseq_);
  while (length == 0) { // Skip the sequences of length 0
    length = kseq_read(sequence_kseq_);
  }
  if (length > 0) {
    kseq_t *sequence = sequence_batch_[sequence_index];
    std::swap(sequence_kseq_->seq, sequence->seq);
    ReplaceByEffectiveRange(sequence->seq);
    std::swap(sequence_kseq_->name, sequence->name);
    std::swap(sequence_kseq_->comment, sequence->comment);
    sequence->id = num_loaded_sequences_;
    ++num_loaded_sequences_;
    if (sequence_kseq_->qual.l != 0) { // fastq file
      std::swap(sequence_kseq_->qual, sequence->qual);
        ReplaceByEffectiveRange(sequence->qual);
    } 
  } else {
    if (length != -1) {
      Chromap<>::ExitWithMessage("Didn't reach the end of sequence file, which might be corrupted!");
    }
    // make sure to reach the end of file rather than meet an error
    no_more_sequence = true;
  }
  return no_more_sequence;
}

uint32_t SequenceBatch::LoadAllSequences() {
  double real_start_time = Chromap<>::GetRealTime();
  sequence_batch_.reserve(200);
  uint32_t num_sequences = 0;
  num_bases_ = 0;
  int length = kseq_read(sequence_kseq_);
  while (length >= 0) { 
    if (length == 0) { // Skip the sequences of length 0
      continue;
    } else if (length > 0) {
      sequence_batch_.emplace_back((kseq_t*)calloc(1, sizeof(kseq_t)));
      kseq_t *sequence = sequence_batch_.back();
      std::swap(sequence_kseq_->seq, sequence->seq);
      ReplaceByEffectiveRange(sequence->seq);
      std::swap(sequence_kseq_->name, sequence->name);
      std::swap(sequence_kseq_->comment, sequence->comment);
      if (sequence_kseq_->qual.l != 0) { // fastq file
        std::swap(sequence_kseq_->qual, sequence->qual);
        ReplaceByEffectiveRange(sequence->qual);
      }
      sequence->id = num_loaded_sequences_;
      ++num_loaded_sequences_;
      ++num_sequences;
      num_bases_ += length;
    } else {
      if (length != -1) {
        Chromap<>::ExitWithMessage("Didn't reach the end of sequence file, which might be corrupted!");
      }
      // make sure to reach the end of file rather than meet an error
      break;
    }
    length = kseq_read(sequence_kseq_);
  }
  std::cerr << "Loaded all sequences successfully in " << Chromap<>::GetRealTime() - real_start_time << "s, ";
  std::cerr << "number of sequences: " << num_sequences << ", ";
  std::cerr << "number of bases: " << num_bases_ << ".\n";
  return num_sequences;
}

void SequenceBatch::FinalizeLoading() {
  kseq_destroy(sequence_kseq_);
  gzclose(sequence_file_);
}
} // namespace chromap
