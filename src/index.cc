#include "index.h"

#include <algorithm>
#include <assert.h>
#include <iostream>

#include "chromap.h"

namespace chromap {
void Index::Statistics(uint32_t num_sequences, const SequenceBatch &reference) {
  double real_start_time = Chromap<>::GetRealTime();
  int n = 0, n1 = 0;
  uint32_t i;
  uint64_t sum = 0, len = 0;
  fprintf(stderr, "[M::%s] kmer size: %d; skip: %d; #seq: %d\n", __func__, kmer_size_, window_size_, num_sequences);
  for (i = 0; i < num_sequences; ++i) {
    len += reference.GetSequenceLengthAt(i);
  }
  assert(len == reference.GetNumBases());
  if (lookup_table_) {
    n += kh_size(lookup_table_);
  }
  for (khint_t k = 0; k < kh_end(lookup_table_); ++k) {
    if (kh_exist(lookup_table_, k)) {
      sum += kh_key(lookup_table_, k) & 1 ? 1 : (uint32_t)kh_val(lookup_table_, k);
      if (kh_key(lookup_table_, k) & 1) 
        ++n1;
    }
  }
  fprintf(stderr, "[M::%s::%.3f] distinct minimizers: %d (%.2f%% are singletons); average occurrences: %.3lf; average spacing: %.3lf\n",
      __func__, Chromap<>::GetRealTime() - real_start_time, n, 100.0*n1/n, (double)sum / n, (double)len / sum);
}

// always reserve space for minimizers in other functions
void Index::GenerateMinimizerSketch(const SequenceBatch &sequence_batch, uint32_t sequence_index, std::vector<std::pair<uint64_t, uint64_t> > *minimizers) {
  uint64_t num_shifted_bits = 2 * (kmer_size_ - 1); 
  uint64_t mask = (((uint64_t)1) << (2 * kmer_size_)) - 1;
  uint64_t seeds_in_two_strands[2] = {0, 0};
  std::pair<uint64_t, uint64_t> buffer[256];
  std::pair<uint64_t, uint64_t> min_seed = {UINT64_MAX, UINT64_MAX};
  uint32_t sequence_length = sequence_batch.GetSequenceLengthAt(sequence_index);
  const char *sequence = sequence_batch.GetSequenceAt(sequence_index);
  assert(sequence_length > 0 && (window_size_ > 0 && window_size_ < 256) && (kmer_size_ > 0 && kmer_size_ <= 28)); // 56 bits for k-mer; could use long k-mers, but 28 enough in practice
  memset(buffer, 0xff, window_size_ * 16); // 2 uint64_t cost 16 bytes
  int unambiguous_length = 0;
  int position_in_buffer = 0;
  int min_position = 0;
  for (uint32_t position = 0; position < sequence_length; ++position) {
    uint8_t current_base = SequenceBatch::CharToUint8(sequence[position]);
    std::pair<uint64_t, uint64_t> current_seed = {UINT64_MAX, UINT64_MAX};
    if (current_base < 4) { // not an ambiguous base
      seeds_in_two_strands[0] = ((seeds_in_two_strands[0] << 2) | current_base) & mask; // forward k-mer
      seeds_in_two_strands[1] = (seeds_in_two_strands[1] >> 2) | (((uint64_t)(3 ^ current_base)) << num_shifted_bits); // reverse k-mer
      if (seeds_in_two_strands[0] == seeds_in_two_strands[1]) {
        continue; // skip "symmetric k-mers" as we don't know it strand
      }
      uint64_t hash_keys_for_two_seeds[2] = {Hash64(seeds_in_two_strands[0], mask), Hash64(seeds_in_two_strands[1], mask)};
      uint64_t strand = hash_keys_for_two_seeds[0] < hash_keys_for_two_seeds[1] ? 0 : 1; // strand
      //uint64_t strand = seeds_in_two_strands[0] < seeds_in_two_strands[1] ? 0 : 1; // strand
      ++unambiguous_length;
      if (unambiguous_length >= kmer_size_) {
        //current_seed.first = Hash64(seeds_in_two_strands[strand], mask);
        current_seed.first = Hash64(hash_keys_for_two_seeds[strand], mask);
        current_seed.second = ((((uint64_t)sequence_index) << 32 | (uint32_t)position) << 1) | strand;
      }
    } else {
      unambiguous_length = 0;
    }
    buffer[position_in_buffer] = current_seed; // need to do this here as appropriate position_in_buffer and buf[position_in_buffer] are needed below
    if (unambiguous_length == window_size_ + kmer_size_ - 1 && min_seed.first != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
      for (int j = position_in_buffer + 1; j < window_size_; ++j)
        if (min_seed.first == buffer[j].first && buffer[j].second != min_seed.second) 
          minimizers->push_back(buffer[j]);
      for (int j = 0; j < position_in_buffer; ++j)
        if (min_seed.first == buffer[j].first && buffer[j].second != min_seed.second) 
          minimizers->push_back(buffer[j]);
    }
    if (current_seed.first <= min_seed.first) { // a new minimum; then write the old min
      if (unambiguous_length >= window_size_ + kmer_size_ && min_seed.first != UINT64_MAX) {
        minimizers->push_back(min_seed);
      }
      min_seed = current_seed;
      min_position = position_in_buffer;
    } else if (position_in_buffer == min_position) { // old min has moved outside the window
      if (unambiguous_length >= window_size_ + kmer_size_ - 1 && min_seed.first != UINT64_MAX) {
        minimizers->push_back(min_seed);
      }
      min_seed.first = UINT64_MAX;
      for (int j = position_in_buffer + 1; j < window_size_; ++j) { // the two loops are necessary when there are identical k-mers
        if (min_seed.first >= buffer[j].first) {// >= is important s.t. min is always the closest k-mer
          min_seed = buffer[j];
          min_position = j; 
        }
      }
      for (int j = 0; j <= position_in_buffer; ++j) {
        if (min_seed.first >= buffer[j].first) {
          min_seed = buffer[j];
          min_position = j;
        }
      }
      if (unambiguous_length >= window_size_ + kmer_size_ - 1 && min_seed.first != UINT64_MAX) { // write identical k-mers
        for (int j = position_in_buffer + 1; j < window_size_; ++j) // these two loops make sure the output is sorted
          if (min_seed.first == buffer[j].first && min_seed.second != buffer[j].second) 
            minimizers->push_back(buffer[j]);
        for (int j = 0; j <= position_in_buffer; ++j)
          if (min_seed.first == buffer[j].first && min_seed.second != buffer[j].second) 
            minimizers->push_back(buffer[j]);
      }
    }
    ++position_in_buffer;
    if (position_in_buffer == window_size_) {
      position_in_buffer = 0;
    }
  }
  if (min_seed.first != UINT64_MAX) {
    minimizers->push_back(min_seed);
  }
}

void Index::Construct(uint32_t num_sequences, const SequenceBatch &reference) {
  double real_start_time = Chromap<>::GetRealTime();
  // tmp_table stores (minimizer, position)
  std::vector< std::pair<uint64_t, uint64_t> > tmp_table;
  tmp_table.reserve(reference.GetNumBases() / window_size_ * 2);
  for (uint32_t sequence_index = 0; sequence_index < num_sequences; ++sequence_index) {
    GenerateMinimizerSketch(reference, sequence_index, &tmp_table);
  }
  std::cerr << "Collected " << tmp_table.size() << " minimizers.\n";
  std::stable_sort(tmp_table.begin(), tmp_table.end());
  std::cerr << "Sorted minimizers.\n";
  uint32_t num_minimizers = tmp_table.size();
  assert(num_minimizers != 0 && num_minimizers <= INT_MAX); // Here I make sure the # minimizers is less than the limit of signed int32, so that I can use int to store position later.
  occurrence_table_.reserve(num_minimizers);
  uint64_t previous_key = tmp_table[0].first;
  uint32_t num_previous_minimizer_occurrences = 0;
  uint64_t num_nonsingletons = 0;
  uint32_t num_singletons = 0;
  for (uint32_t ti = 0; ti < num_minimizers; ++ti) {
    uint64_t current_key = tmp_table[ti].first;
    if (current_key != previous_key) {
      int khash_return_code;
      khiter_t khash_iterator = kh_put(k64, lookup_table_, previous_key << 1, &khash_return_code);
      assert(khash_return_code != -1 && khash_return_code != 0);
      if (num_previous_minimizer_occurrences == 1) { // singleton
        kh_key(lookup_table_, khash_iterator) |= 1;
        kh_value(lookup_table_, khash_iterator) = occurrence_table_.back();
        occurrence_table_.pop_back();
        ++num_singletons;
      } else {
        kh_value(lookup_table_, khash_iterator) = (num_nonsingletons << 32) | num_previous_minimizer_occurrences;
        num_nonsingletons += num_previous_minimizer_occurrences;
      }
      num_previous_minimizer_occurrences = 1;
    } else {
      num_previous_minimizer_occurrences++;
    }
    occurrence_table_.push_back(tmp_table[ti].second);
    previous_key = current_key;
  }
  int khash_return_code;
  khiter_t khash_iterator = kh_put(k64, lookup_table_, previous_key << 1, &khash_return_code);
  assert(khash_return_code != -1 && khash_return_code != 0);
  if (num_previous_minimizer_occurrences == 1) { // singleton
    kh_key(lookup_table_, khash_iterator) |= 1;
    kh_value(lookup_table_, khash_iterator) = occurrence_table_.back();
    occurrence_table_.pop_back();
    ++num_singletons;
  } else {
    kh_value(lookup_table_, khash_iterator) = (num_nonsingletons << 32) | num_previous_minimizer_occurrences;
    num_nonsingletons += num_previous_minimizer_occurrences;
  }
  assert(num_nonsingletons + num_singletons == num_minimizers);
  std::cerr << "Kmer size: " << kmer_size_ << ", window size: " << window_size_ << ".\n"; 
  std::cerr << "Lookup table size: " << kh_size(lookup_table_) << ", occurrence table size: " << occurrence_table_.size() << ", # singletons: " << num_singletons << ".\n";
  std::cerr << "Built index successfully in " << Chromap<>::GetRealTime() - real_start_time << "s.\n";
}

void Index::CheckIndex(uint32_t num_sequences, const SequenceBatch &reference) {
  std::vector< std::pair<uint64_t, uint64_t> > tmp_table;
  tmp_table.reserve(reference.GetNumBases() / window_size_ * 2);
  for (uint32_t sequence_index = 0; sequence_index < num_sequences; ++sequence_index) {
    GenerateMinimizerSketch(reference, sequence_index, &tmp_table);
  }
  std::cerr << "Collected " << tmp_table.size() << " minimizers.\n";
  std::stable_sort(tmp_table.begin(), tmp_table.end());
  std::cerr << "Sorted minimizers.\n";
  uint32_t count = 0;
  for (uint32_t i = 0; i < tmp_table.size(); ++i) {
    khiter_t khash_iterator = kh_get(k64, lookup_table_, tmp_table[i].first << 1);
    assert(khash_iterator != kh_end(lookup_table_));
    uint64_t key = kh_key(lookup_table_, khash_iterator);
    uint64_t value = kh_value(lookup_table_, khash_iterator);
    if (key & 1) { //singleton
      assert(tmp_table[i].second == value);
      count = 0;
    } else {
      uint32_t offset = value >> 32;
      uint32_t num_occ = value;
      uint64_t value_in_index = occurrence_table_[offset + count];
      assert(value_in_index == tmp_table[i].second);
      ++count;
      if (count == num_occ) {
        count = 0;
      }
    }
  }
}

void Index::Save() {
  double real_start_time = Chromap<>::GetRealTime();
  FILE *index_file = fopen(index_file_path_.c_str(), "wb");
  assert(index_file != NULL);
  uint64_t num_bytes = 0;
  int err = 0;
  err = fwrite(&kmer_size_, sizeof(int), 1, index_file);
  num_bytes += sizeof(int);
  assert(err != 0);
  err = fwrite(&window_size_, sizeof(int), 1, index_file);
  num_bytes += sizeof(int);
  assert(err != 0);
  uint32_t lookup_table_size = kh_size(lookup_table_);
  err = fwrite(&lookup_table_size, sizeof(uint32_t), 1, index_file);
  num_bytes += sizeof(uint32_t);
  assert(err != 0);
  kh_save(k64, lookup_table_, index_file);
  num_bytes += sizeof(uint64_t) * 2 * lookup_table_size;
  uint32_t occurrence_table_size = occurrence_table_.size();
  err = fwrite(&occurrence_table_size, sizeof(uint32_t), 1, index_file);
  num_bytes += sizeof(uint32_t);
  assert(err != 0);
  err = fwrite(occurrence_table_.data(), sizeof(uint64_t), occurrence_table_size, index_file);
  num_bytes += sizeof(uint64_t) * occurrence_table_size;
  assert(err != 0);
  fclose(index_file);
  //std::cerr << "Index size: " << num_bytes / (1024.0 * 1024 * 1024) << "GB, saved in " << Chromap<>::GetRealTime() - real_start_time << "s.\n";
  std::cerr << "Saved in " << Chromap<>::GetRealTime() - real_start_time << "s.\n";
}

void Index::Load() {
  double real_start_time = Chromap<>::GetRealTime();
  FILE *index_file = fopen(index_file_path_.c_str(), "rb");
  assert(index_file != NULL);
  int err = 0;
  err = fread(&kmer_size_, sizeof(int), 1, index_file);
  assert(err != 0);
  err = fread(&window_size_, sizeof(int), 1, index_file);
  assert(err != 0);
  uint32_t lookup_table_size = 0;
  err = fread(&lookup_table_size, sizeof(uint32_t), 1, index_file);
  assert(err != 0);
  kh_load(k64, lookup_table_, index_file);
  uint32_t occurrence_table_size = 0;
  err = fread(&occurrence_table_size, sizeof(uint32_t), 1, index_file);
  assert(err != 0);
  occurrence_table_.resize(occurrence_table_size); 
  err = fread(occurrence_table_.data(), sizeof(uint64_t), occurrence_table_size, index_file);
  assert(err != 0);
  fclose(index_file);
  std::cerr << "Kmer size: " << kmer_size_ << ", window size: " << window_size_ << ".\n";
  std::cerr << "Lookup table size: " << kh_size(lookup_table_) << ", occurrence table size: " << occurrence_table_.size() << ".\n";
  std::cerr << "Loaded index successfully in "<< Chromap<>::GetRealTime() - real_start_time << "s.\n";
}

void Index::GenerateCandidatesOnOneDirection(int error_threshold, std::vector<uint64_t> *hits, std::vector<uint64_t> *candidates) {
  hits->emplace_back(UINT64_MAX);
  if (hits->size() > 0) {
    std::sort(hits->begin(), hits->end());
    int count = 1;
    uint64_t previous_hit = (*hits)[0];
    uint32_t previous_reference_id = previous_hit >> 32;
    uint32_t previous_reference_position = previous_hit;
    for (uint32_t pi = 1; pi < hits->size(); ++pi) {
      uint32_t current_reference_id = (*hits)[pi] >> 32;
      uint32_t current_reference_position = (*hits)[pi];
      if (current_reference_id != previous_reference_id || current_reference_position > previous_reference_position + error_threshold) {
        if (count >= min_num_seeds_required_for_mapping_) {
          candidates->push_back(previous_hit);
        }
        count = 1;
      } else {
        ++count;
      }
      previous_hit = (*hits)[pi];
      previous_reference_id = current_reference_id;
      previous_reference_position = current_reference_position;
    }
  }
}

void Index::CollectCandiates(int max_seed_frequency, const std::vector<std::pair<uint64_t, uint64_t> > &minimizers, std::vector<uint64_t> *positive_hits, std::vector<uint64_t> *negative_hits) {
  uint32_t num_minimizers = minimizers.size();
  positive_hits->reserve(max_seed_frequencies_[0]);
  negative_hits->reserve(max_seed_frequencies_[0]);
  for (uint32_t mi = 0; mi < num_minimizers; ++mi) {
    khiter_t khash_iterator = kh_get(k64, lookup_table_, minimizers[mi].first << 1);
    if (khash_iterator == kh_end(lookup_table_)) {
      //std::cerr << "The minimizer is not in reference!\n";
      continue;
    }
    uint64_t value = kh_value(lookup_table_, khash_iterator);
    uint32_t read_position = minimizers[mi].second >> 1;
    if (kh_key(lookup_table_, khash_iterator) & 1) { // singleton
      uint64_t reference_id = value >> 33;
      uint32_t reference_position = value >> 1;
      // Check whether the strands of reference minimizer and read minimizer are the same
      // Later, we can play some tricks with 0,1 here to make it faster.
      if (((minimizers[mi].second & 1) ^ (value & 1)) == 0) { // same
        uint32_t candidate_position = reference_position - read_position;// > 0 ? reference_position - read_position : 0;
        // ok, for now we can't see the reference here. So let us don't do the check.
        // Instead, we do it later some time when we check the candidates.
        uint64_t candidate = (reference_id << 32) | candidate_position;
        positive_hits->push_back(candidate);
      } else {
        uint32_t candidate_position = reference_position + read_position - kmer_size_ + 1;// < reference_length ? reference_position - read_position : 0;
        uint64_t candidate = (reference_id << 32) | candidate_position;
        negative_hits->push_back(candidate);
      }
    } else {
      uint32_t offset = value >> 32;
      uint32_t num_occurrences = value;
      if (num_occurrences < (uint32_t)max_seed_frequency) {
        for (uint32_t oi = 0; oi < num_occurrences; ++oi) {
          uint64_t value = occurrence_table_[offset + oi];
          uint64_t reference_id = value >> 33;
          uint32_t reference_position = value >> 1;
          if (((minimizers[mi].second & 1) ^ (value & 1)) == 0) { // same
            uint32_t candidate_position = reference_position - read_position;
            uint64_t candidate = (reference_id << 32) | candidate_position;
            positive_hits->push_back(candidate);
          } else {
            uint32_t candidate_position = reference_position + read_position - kmer_size_ + 1;
            uint64_t candidate = (reference_id << 32) | candidate_position;
            negative_hits->push_back(candidate);
          }
        } 
      }
    }
  }
}

void Index::GenerateCandidates(int error_threshold, const std::vector<std::pair<uint64_t, uint64_t> > &minimizers, std::vector<uint64_t> *positive_hits, std::vector<uint64_t> *negative_hits, std::vector<uint64_t> *positive_candidates, std::vector<uint64_t> *negative_candidates) {
  CollectCandiates(max_seed_frequencies_[0], minimizers, positive_hits, negative_hits);
  // Now I can generate primer chain in candidates
  // Let me use sort for now, but I can use merge later.
  GenerateCandidatesOnOneDirection(error_threshold, positive_hits, positive_candidates);
  GenerateCandidatesOnOneDirection(error_threshold, negative_hits, negative_candidates);
  if (positive_candidates->size() + negative_candidates->size() == 0) {
    positive_hits->clear();
    negative_hits->clear();
    CollectCandiates(max_seed_frequencies_[1], minimizers, positive_hits, negative_hits);
    GenerateCandidatesOnOneDirection(error_threshold, positive_hits, positive_candidates);
    GenerateCandidatesOnOneDirection(error_threshold, negative_hits, negative_candidates);
    // TODO: if necessary, we can further improve the rescue. But the code below is not thread safe. We can think about this later
//    if (positive_candidates->size() + negative_candidates->size() == 0) {
//      --min_num_seeds_required_for_mapping_;
//      min_num_seeds_required_for_mapping_ = std::max(min_num_seeds_required_for_mapping_, 1);
//      GenerateCandidatesOnOneDirection(positive_hits, positive_candidates);
//      GenerateCandidatesOnOneDirection(negative_hits, negative_candidates);
//    }
  }
}
} // namespace chromap
