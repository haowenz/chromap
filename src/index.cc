#include "index.h"

#include <algorithm>
#include <assert.h>
#include <iostream>

#include "chromap.h"

namespace chromap {

struct mmHit {
	uint32_t mi ;
	uint64_t position ;
	bool operator<(const mmHit &h)  const {
		return position > h.position; // the inversed direction is to make a min-heap
	}
} ;

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
  //std::pair<uint64_t, uint64_t> pre_seed = {UINT64_MAX, UINT64_MAX};
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
    //if (unambiguous_length == window_size_ + kmer_size_ - 1 && min_seed.first != UINT64_MAX) { // special case for the first window - because identical k-mers are not stored yet
    //  //
    //  //bool found_first_min_in_first_loop = false;
    //  //bool found_first_min_in_second_loop = false; // For debug
    //  //int first_min_position_in_buffer = 0;
    //  //for (int j = position_in_buffer + 1; j < window_size_; ++j)
    //  //  if (min_seed.first == buffer[j].first) {
    //  //    found_first_min_in_first_loop = true;
    //  //    first_min_position_in_buffer = j;
    //  //    break;
    //  //  }
    //  //if (!found_first_min_in_first_loop) {
    //  //  for (int j = 0; j < position_in_buffer; ++j)
    //  //    if (min_seed.first == buffer[j].first) {
    //  //      found_first_min_in_second_loop = true;
    //  //      first_min_position_in_buffer = j;
    //  //      break;
    //  //    }
    //  //}
    //  //assert(found_first_min_in_first_loop || found_first_min_in_second_loop);
    //  //bool found_pre_seed = false;
    //  //if (found_first_min_in_first_loop) {
    //  //  for (int j = position_in_buffer + 1; j < first_min_position_in_buffer; ++j)
    //  //    if (buffer[j].first < pre_seed.first) {
    //  //      pre_seed = buffer[j];
    //  //      found_pre_seed = true;
    //  //    }
    //  //}
    //  //if (found_first_min_in_second_loop) {
    //  //  for (int j = 0; j < first_min_position_in_buffer; ++j)
    //  //    if (buffer[j].first < pre_seed.first) {
    //  //      pre_seed = buffer[j];
    //  //      found_pre_seed = true;
    //  //    }
    //  //}
    //  //if (found_pre_seed)
    //  //  minimizers->push_back(pre_seed);
    //  //
    //  for (int j = position_in_buffer + 1; j < window_size_; ++j)
    //    if (min_seed.first == buffer[j].first && buffer[j].second != min_seed.second) 
    //      minimizers->push_back(buffer[j]);
    //  for (int j = 0; j < position_in_buffer; ++j)
    //    if (min_seed.first == buffer[j].first && buffer[j].second != min_seed.second) 
    //      minimizers->push_back(buffer[j]);
    //}

    if (unambiguous_length == window_size_ + kmer_size_ - 1 && min_seed.first != UINT64_MAX && min_seed.first < current_seed.first) { // special case for the first window - because identical k-mers are not stored yet
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

void Index::GenerateCandidatesOnOneDirection(int error_threshold, int num_seeds_required, std::vector<uint64_t> *hits, std::vector<Candidate> *candidates) {
  //std::cerr << "Direction\n";
  hits->emplace_back(UINT64_MAX);
  if (hits->size() > 0) {
    //std::sort(hits->begin(), hits->end());
    int count = 1;
    int equal_count = 1 ; // the number of seeds with the exact same reference position.
    int best_equal_count = 1 ;

    uint64_t previous_hit = (*hits)[0];
    uint32_t previous_reference_id = previous_hit >> 32;
    uint32_t previous_reference_position = previous_hit;
    uint64_t best_local_hit = (*hits)[0];
    for (uint32_t pi = 1; pi < hits->size(); ++pi) {
      uint32_t current_reference_id = (*hits)[pi] >> 32;
      uint32_t current_reference_position = (*hits)[pi];
#ifdef LI_DEBUG
      printf("%s: %d %d\n", __func__, current_reference_id, current_reference_position);
#endif
      if (current_reference_id != previous_reference_id || current_reference_position > previous_reference_position + error_threshold) {
        if (count >= num_seeds_required) {
          Candidate nc;
          nc.position = best_local_hit;
          nc.count = count;
          candidates->push_back(nc);
          //std::cerr << count << " ";
        }
        count = 1;
	equal_count = 1;
	best_equal_count = 1;
	best_local_hit = (*hits)[pi];
      } else {
	//printf("%d %d %d: %d %d\n", (int)best_local_hit, (int)previous_hit, (int)(*hits)[pi], equal_count, best_equal_count);
        if ( (*hits)[pi] == best_local_hit ) { 
		++equal_count ;
		++best_equal_count ;
	} else if ( (*hits)[pi] == previous_hit ) {
		++equal_count ;
		if (equal_count > best_equal_count) {
			best_local_hit = previous_hit ;
			best_equal_count = equal_count ;
		}
	} else {
		equal_count = 1 ;
	}
        ++count;
      }
      
      previous_hit = (*hits)[pi];
      previous_reference_id = current_reference_id;
      previous_reference_position = current_reference_position;
    }
  }
}

// Return the number of repetitive seeds
int Index::CollectCandidates(int max_seed_frequency, int repetitive_seed_frequency, const std::vector<std::pair<uint64_t, uint64_t> > &minimizers, uint32_t *repetitive_seed_length, std::vector<uint64_t> *positive_hits, std::vector<uint64_t> *negative_hits, bool use_heap) {
  uint32_t num_minimizers = minimizers.size();
  int repetitive_seed_count = 0 ;
  std::vector<uint64_t> *mm_positive_hits = NULL, *mm_negative_hits = NULL;
  bool heap_resort = false ; // need to sort the elements of heap first
  if (use_heap) {
    mm_positive_hits = new std::vector<uint64_t>[num_minimizers] ;
    mm_negative_hits = new std::vector<uint64_t>[num_minimizers] ;
  }
  positive_hits->reserve(max_seed_frequencies_[0]);
  negative_hits->reserve(max_seed_frequencies_[0]);
  uint32_t previous_repetitive_seed_position = std::numeric_limits<uint32_t>::max();
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
        if (use_heap) {
	  mm_positive_hits[mi].push_back(candidate);
	}
	else
	  positive_hits->push_back(candidate);
      } else {
        uint32_t candidate_position = reference_position + read_position - kmer_size_ + 1;// < reference_length ? reference_position - read_position : 0;
        uint64_t candidate = (reference_id << 32) | candidate_position;
	if (use_heap)
          mm_negative_hits[mi].push_back(candidate);
	else
	  negative_hits->push_back(candidate);
      }
    } else {
      uint32_t offset = value >> 32;
      uint32_t num_occurrences = value;
      //printf("%s: %u %u\n", __func__, offset, num_occurrences) ;
      if (num_occurrences < (uint32_t)max_seed_frequency) {
        for (uint32_t oi = 0; oi < num_occurrences; ++oi) {
          uint64_t value = occurrence_table_[offset + oi];
          uint64_t reference_id = value >> 33;
          uint32_t reference_position = value >> 1;
          if (((minimizers[mi].second & 1) ^ (value & 1)) == 0) { // same
            uint32_t candidate_position = reference_position - read_position;
            uint64_t candidate = (reference_id << 32) | candidate_position;
	    if (use_heap) {
	      if (reference_position < read_position)
	        heap_resort = true ;
              mm_positive_hits[mi].push_back(candidate);
	    }
	    else
	      positive_hits->push_back(candidate);
          } else {
            uint32_t candidate_position = reference_position + read_position - kmer_size_ + 1;
            uint64_t candidate = (reference_id << 32) | candidate_position;
	    if (use_heap)
              mm_negative_hits[mi].push_back(candidate);
	    else
	      negative_hits->push_back(candidate);
          }
        } 
      } 
      
      if (num_occurrences >= (uint32_t)repetitive_seed_frequency){
        if (previous_repetitive_seed_position > read_position) { // first minimizer
          *repetitive_seed_length += kmer_size_;
        } else {
          if (read_position < previous_repetitive_seed_position + kmer_size_ + window_size_ - 1) {
            *repetitive_seed_length += read_position - previous_repetitive_seed_position;
          } else {
            *repetitive_seed_length += kmer_size_;
          }
        }
        previous_repetitive_seed_position = read_position;
        ++repetitive_seed_count ;
      }
    }
  }

  if (use_heap) {
	  std::priority_queue<struct mmHit> heap;
	  unsigned int *mm_pos = new unsigned int[num_minimizers];
	  positive_hits->clear();
	  for (uint32_t mi = 0; mi < num_minimizers; ++mi) {
		  if (mm_positive_hits[mi].size() == 0)
			  continue;
		  // only the positive part may have the underflow issue
		  if (heap_resort) 
		  	std::sort(mm_positive_hits[mi].begin(), mm_positive_hits[mi].end());
		  struct mmHit nh;
		  nh.mi = mi;
		  nh.position = mm_positive_hits[mi][0];
		  heap.push(nh);
		  mm_pos[mi] = 0;
	  }

	  while(!heap.empty()) {
		  struct mmHit top = heap.top();
		  heap.pop();
		  positive_hits->push_back(top.position) ;
		  ++mm_pos[top.mi];
		  if (mm_pos[top.mi] < mm_positive_hits[top.mi].size())
		  {
			  struct mmHit nh;
			  nh.mi = top.mi;
			  nh.position = mm_positive_hits[top.mi][mm_pos[top.mi]];
			  heap.push(nh);
		  }
	  }

	  negative_hits->clear();
	  for (uint32_t mi = 0; mi < num_minimizers; ++mi) {
		  if (mm_negative_hits[mi].size() == 0)
			  continue;
		  struct mmHit nh;
		  nh.mi = mi;
		  nh.position = mm_negative_hits[mi][0];
		  heap.push(nh);
		  mm_pos[mi] = 0;
	  }
	  while(!heap.empty()) {
		  struct mmHit top = heap.top();
		  heap.pop();
		  negative_hits->push_back(top.position) ;
		  ++mm_pos[top.mi];
		  if (mm_pos[top.mi] < mm_negative_hits[top.mi].size())
		  {
			  struct mmHit nh;
			  nh.mi = top.mi;
			  nh.position = mm_negative_hits[top.mi][mm_pos[top.mi]];
			  heap.push(nh);
		  }
	  }
	  delete[] mm_positive_hits;
	  delete[] mm_negative_hits;
	  delete[] mm_pos ;
  } else {
  	std::sort(positive_hits->begin(), positive_hits->end());
	std::sort(negative_hits->begin(), negative_hits->end());
  }
  /*for (uint32_t mi = 0 ; mi < positive_hits->size() ; ++mi)
  	printf("+ %llu %d %d\n", positive_hits->at(mi), (int)(positive_hits->at(mi)>>32), (int)(positive_hits->at(mi))) ;

  for (uint32_t mi = 0 ; mi < negative_hits->size() ; ++mi)
  	printf("- %llu %d %d\n", negative_hits->at(mi), (int)(negative_hits->at(mi)>>32), (int)(negative_hits->at(mi))) ;*/
  return repetitive_seed_count ;
}

void Index::GenerateCandidatesFromRepetitiveReadWithMateInfo(int error_threshold, const std::vector<std::pair<uint64_t, uint64_t> > &minimizers, uint32_t *repetitive_seed_length, std::vector<uint64_t> *hits, std::vector<Candidate> *candidates, std::vector<Candidate> *mate_candidates, int direction, unsigned int range) // directoin: +1: positive; -1: negative
{
  uint32_t num_minimizers = minimizers.size();
  hits->reserve(max_seed_frequencies_[0]);
  uint32_t previous_repetitive_seed_position = std::numeric_limits<uint32_t>::max();

  int best_candidate = -1;
  int max_count = 0;
  int best_candidate_num = 0;
  uint32_t mate_candidates_size = mate_candidates->size();
  for (uint32_t i = 0 ; i < mate_candidates_size; ++i) {
 	int count = mate_candidates->at(i).count;
  	if (count > max_count) {
		best_candidate = i;
		max_count = count;
		best_candidate_num = 1;
	}
	else if (count == max_count) {
		++best_candidate_num;
	}
  }
  if (best_candidate_num != 1 || max_count < min_num_seeds_required_for_mapping_) 
  	return;

  *repetitive_seed_length = 0 ;
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
        if (direction == 1) {
          uint32_t candidate_position = reference_position - read_position;// > 0 ? reference_position - read_position : 0;
          // ok, for now we can't see the reference here. So let us don't do the check.
          // Instead, we do it later some time when we check the candidates.
          uint64_t candidate = (reference_id << 32) | candidate_position;
          hits->push_back(candidate);
        }
      } else if (direction == -1) {
        uint32_t candidate_position = reference_position + read_position - kmer_size_ + 1;// < reference_length ? reference_position - read_position : 0;
        uint64_t candidate = (reference_id << 32) | candidate_position;
        hits->push_back(candidate);
      }
    } else {
      uint32_t offset = value >> 32;
      uint32_t num_occurrences = value;
      // use binary search to locate the coordinate near mate position
      int32_t l = 0, m = 0, r = num_occurrences - 1;
      uint64_t boundary = (mate_candidates->at(best_candidate).position < range) ? 0 : (mate_candidates->at(best_candidate).position - range) ;
      while (l <= r) {
        m = (l + r) / 2;
        uint64_t value = (occurrence_table_[offset + m])>>1;
        if (value <= boundary) {
          l = m + 1;
        } else if (value > boundary) {
          r = m - 1;
        } else {
          break ;
        }
      }
      //printf("%s: %d %d: %d %d\n", __func__, m, num_occurrences,
      //	(int)(boundary>>32), (int)boundary) ;
      for (uint32_t oi = m; oi < num_occurrences; ++oi) {
        uint64_t value = occurrence_table_[offset + oi];
        if ((value >> 1) > mate_candidates->at(best_candidate).position + range)
          break;
        uint64_t reference_id = value >> 33;
        uint32_t reference_position = value >> 1;
        if (((minimizers[mi].second & 1) ^ (value & 1)) == 0) { // same
          if (direction == 1) {
            uint32_t candidate_position = reference_position - read_position;
            uint64_t candidate = (reference_id << 32) | candidate_position;
            hits->push_back(candidate);
          }
        } else if (direction == -1){
          uint32_t candidate_position = reference_position + read_position - kmer_size_ + 1;
          uint64_t candidate = (reference_id << 32) | candidate_position;
          hits->push_back(candidate);
        }
      }  

      if (num_occurrences >= (uint32_t)max_seed_frequencies_[0]) {
        if (previous_repetitive_seed_position > read_position) { // first minimizer
          *repetitive_seed_length += kmer_size_;
        } else {
          if (read_position < previous_repetitive_seed_position + kmer_size_ + window_size_ - 1) {
            *repetitive_seed_length += read_position - previous_repetitive_seed_position;
          } else {
            *repetitive_seed_length += kmer_size_;
          }
        }
        previous_repetitive_seed_position = read_position;
      }
    }
  } // for mi

  std::sort(hits->begin(), hits->end());
  //for (uint32_t i = 0 ; i < hits->size(); ++i)
  //	  printf("%s: %d %d\n", __func__, (int)(hits->at(i)>>32),(int)hits->at(i));
  GenerateCandidatesOnOneDirection(error_threshold, 1, hits, candidates);
  //printf("%s: %d %d\n", __func__, hits->size(), candidates->size()) ;
}

void Index::GenerateCandidates(int error_threshold, const std::vector<std::pair<uint64_t, uint64_t> > &minimizers, uint32_t *repetitive_seed_length, std::vector<uint64_t> *positive_hits, std::vector<uint64_t> *negative_hits, std::vector<Candidate> *positive_candidates, std::vector<Candidate> *negative_candidates) {
  *repetitive_seed_length = 0;
  bool recollect = true;
  int repetitive_seed_count = CollectCandidates(max_seed_frequencies_[0], max_seed_frequencies_[0], minimizers, repetitive_seed_length, positive_hits, negative_hits, false);
  if (repetitive_seed_count > (int)minimizers.size() / 2 && minimizers.size() >= 10) {
    positive_hits->clear();
    negative_hits->clear();
    *repetitive_seed_length = 0;
    repetitive_seed_count = CollectCandidates(max_seed_frequencies_[1], max_seed_frequencies_[0], minimizers, repetitive_seed_length, positive_hits, negative_hits, true);
    recollect = false;
  }

  // Now I can generate primer chain in candidates
  // Let me use sort for now, but I can use merge later.
  //printf("p+n: %d. %d %d\n", positive_hits->size() + negative_hits->size(), repetitive_seed_count, minimizers.size()) ;
  GenerateCandidatesOnOneDirection(error_threshold, min_num_seeds_required_for_mapping_, positive_hits, positive_candidates);
  GenerateCandidatesOnOneDirection(error_threshold, min_num_seeds_required_for_mapping_, negative_hits, negative_candidates);
  //fprintf(stderr, "p+n: %d\n", positive_candidates->size() + negative_candidates->size()) ;
  if (positive_candidates->size() + negative_candidates->size() == 0 && recollect) {
    positive_hits->clear();
    negative_hits->clear();
    //printf("second round\n") ;
    *repetitive_seed_length = 0;
    CollectCandidates(max_seed_frequencies_[1], max_seed_frequencies_[0], minimizers, repetitive_seed_length, positive_hits, negative_hits, true);
    //printf("p+n2: %d\n", positive_hits->size() + negative_hits->size()) ;
    GenerateCandidatesOnOneDirection(error_threshold, min_num_seeds_required_for_mapping_, positive_hits, positive_candidates);
    GenerateCandidatesOnOneDirection(error_threshold, min_num_seeds_required_for_mapping_, negative_hits, negative_candidates);
    // TODO: if necessary, we can further improve the rescue. But the code below is not thread safe. We can think about this later
//    if (positive_candidates->size() + negative_candidates->size() == 0) {
//      --min_num_seeds_required_for_mapping_;
//      min_num_seeds_required_for_mapping_ = std::max(min_num_seeds_required_for_mapping_, 1);
//      GenerateCandidatesOnOneDirection(positive_hits, positive_candidates);
//      GenerateCandidatesOnOneDirection(negative_hits, negative_candidates);
//    }
  }
  //uint32_t i ;
  //uint32_t size = positive_candidates->size() ;
  //for (i = 0 ; i < size ; ++i)
  //	(*positive_candidates)[i].direction = kPositive ;
  //size = negative_candidates->size() ;
  //for (i = 0 ; i < size ; ++i)
  //	(*negative_candidates)[i].direction = kNegative ;
  //printf("p+n_candidates: %d\n", positive_candidates->size() + negative_candidates->size()) ;
}
} // namespace chromap
