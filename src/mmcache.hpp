#ifndef CHROMAP_CACHE_H_
#define CHROMAP_CACHE_H_

#include "index.h"

#define FINGER_PRINT_SIZE 103

namespace chromap {

class mm_cache_candidate_list {
private:
  uint32_t *positions;
  uint8_t *counts; // 0: the corresponding position is chr id. Otherwise, it is coordinate within the chr id
  uint32_t size;
  uint32_t actual_size;

  void Clear() {
    if (positions != NULL)
    {
      delete[] positions;
      positions = NULL;
    }
    if (counts != NULL)
    {
      delete[] counts;
      counts = NULL;
    }
    size = actual_size = 0;
  }
public:
  mm_cache_candidate_list(){
    positions = NULL;
    counts = NULL;
    size = actual_size = 0;
  }

  ~mm_cache_candidate_list() {
    Clear();
  }
  
  void Input(const std::vector<Candidate> &candidates) {
    Clear();
    int i, k;
    actual_size = candidates.size();
    if (actual_size == 0)
      return;
    size = actual_size + 1;
    
    // Collect the extra size introduced by the chr.
    for (i = 1; i < (int)actual_size; ++i) {
      if ((candidates[i].position >> 32) != (candidates[i - 1].position >> 32))
        ++size;
    }
    positions = new uint32_t[size];
    counts = new uint8_t[size];

    k = 0;
    for (i = 0; i < (int)actual_size; ++i) {
      if (i == 0 || (candidates[i].position >> 32) != (candidates[i - 1].position >> 32)) {
        counts[k] = 0;
	positions[k] = (uint32_t)(candidates[i].position >> 32);
	++k;
      }
      counts[k] = candidates[i].count;
      if (counts[k] == 0) // may happen on overflow for the tandem repeats read
        counts[k] = 1;
      positions[k] = (uint32_t)candidates[i].position;
      ++k;
    }
  }

  void Output(std::vector<Candidate> &candidates, int shift) {
    candidates.resize(actual_size);
    int i, k;
    k = 0;
    uint64_t rid = 0;
    for (i = 0; i < (int)size; ++i) {
     if (counts[i] == 0) {
       rid = (uint64_t)positions[i] << 32ull ;
     } else {
       candidates[k].position = rid | (uint32_t)(positions[i] + shift);
       candidates[k].count = counts[i];
       ++k;
     }
    }
  }

  void Shift(int offset) {
    int i;
    for (i = 0; i < (int)size; ++i) {
      if (counts[i] != 0)
        positions[i] += offset;
    }
  }

  uint32_t GetActualSize() {
    return actual_size;
  }
} ;

struct _mm_cache_entry {
  std::vector<uint64_t> minimizers;
  std::vector<int> offsets; // the distance to the next minimizer
  std::vector<uint8_t> strands;
  //std::vector<Candidate> positive_candidates;
  //std::vector<Candidate> negative_candidates;
 
  mm_cache_candidate_list positive_candidate_list ;
  mm_cache_candidate_list negative_candidate_list ;

  int weight;
  unsigned short finger_print_cnt[FINGER_PRINT_SIZE];
  int finger_print_cnt_sum;
  uint32_t repetitive_seed_length;
};

class mm_cache {
 private:
  int cache_size;
  struct _mm_cache_entry *cache;
  int kmer_length;
  int update_limit;

  // 0: not match. -1: opposite order. 1: same order
  int IsMinimizersMatchCache(const std::vector<std::pair<uint64_t, uint64_t> > &minimizers, const struct _mm_cache_entry &cache) {
    if (cache.minimizers.size() != minimizers.size())
      return 0;
    int size = minimizers.size();
    int i, j;
    int direction = 0;
    for (i = 0; i < size; ++i) {
      if (cache.minimizers[i] != minimizers[i].first || (minimizers[i].second & 1) != cache.strands[i])
        break;
    }
    if (i >= size) {
      for (i = 0; i < size - 1; ++i)	{
        if (cache.offsets[i] != ((int)minimizers[i + 1].second>>1) - ((int)minimizers[i].second>>1))
          break;
      }
      if (i >= size - 1)
        direction = 1;
    }

    if (direction == 1)
      return 1;

    for (i = 0, j = size - 1; i < size; ++i, --j) {
      if (cache.minimizers[i] != minimizers[j].first || (minimizers[j].second & 1) == cache.strands[i])
        break;
    }
    if (i >= size) {
      for (i = 0, j = size - 1; i < size - 1; ++i, --j) {
        if (cache.offsets[i] != ((int)minimizers[j].second>>1) - ((int)minimizers[j - 1].second>>1))
          break;
      }

      if (i >= size - 1) {
        direction = -1;
      }
    }
    return direction;
  }

 public:
  mm_cache(int size) {
    cache = new struct _mm_cache_entry[size];
    cache_size = size;
    //memset(cache, 0, sizeof(cache[0]) * size) ;
    for (int i = 0; i < size; ++i) {
      cache[i].weight = 0;
      memset(cache[i].finger_print_cnt, 0, sizeof(unsigned short) * FINGER_PRINT_SIZE);
      cache[i].finger_print_cnt_sum = 0;
    }
    update_limit = 10;
  }
  ~mm_cache() {
    delete[] cache;
  }

  void SetKmerLength(int kl) {
    kmer_length = kl;
  }

  // Return the hash entry index. -1 if failed.
  int Query(const std::vector<std::pair<uint64_t, uint64_t> > &minimizers, std::vector<Candidate> &pos_candidates, std::vector<Candidate> &neg_candidates, uint32_t &repetitive_seed_length, uint32_t read_len) {
    int i;
    int msize = minimizers.size();
    if (msize == 0)
      return -1;
    uint64_t h = 0;
    for (i = 0 ; i < msize; ++i)
      h += (minimizers[i].first);
    int hidx = h % cache_size;
    int direction = IsMinimizersMatchCache(minimizers, cache[hidx]);
    if (direction == 1) {
      int shift = (int)minimizers[0].second >> 1;
      cache[hidx].positive_candidate_list.Output(pos_candidates, -shift);
      cache[hidx].negative_candidate_list.Output(neg_candidates, shift);
      repetitive_seed_length = cache[hidx].repetitive_seed_length;
      return hidx;
    } else if (direction == -1) {// The "read" is on the other direction of the cached "read"
      // Start position of the last minimizer shoud equal the first minimizer's end position in rc "read".
      int shift = read_len - ((int)minimizers[msize - 1].second>>1) - 1 + kmer_length - 1; 
      
      cache[hidx].negative_candidate_list.Output(pos_candidates, shift - read_len + 1);
      cache[hidx].positive_candidate_list.Output(neg_candidates, -shift + read_len - 1);
      repetitive_seed_length = cache[hidx].repetitive_seed_length;
      return hidx;
    } else {
      return -1;
    }
  }

  void Update(const std::vector<std::pair<uint64_t, uint64_t> > &minimizers, std::vector<Candidate> &pos_candidates, std::vector<Candidate> &neg_candidates, uint32_t repetitive_seed_length) {
    int i;
    int msize = minimizers.size();
    if (msize == 0)
      return;
    uint64_t h = 0; // for hash
    uint64_t f = 0; // for finger printing
    for (i = 0; i < msize; ++i) {
      h += (minimizers[i].first);	
      f ^= (minimizers[i].first);
    }
    int hidx = h % cache_size;
    int finger_print = f % FINGER_PRINT_SIZE; 

    ++cache[hidx].finger_print_cnt[finger_print];
    ++cache[hidx].finger_print_cnt_sum;
    
    if (cache[hidx].finger_print_cnt_sum < 10 || (int)cache[hidx].finger_print_cnt[finger_print] * 5 < cache[hidx].finger_print_cnt_sum) {
      return;
    }

    int direction = IsMinimizersMatchCache(minimizers, cache[hidx]);
    if (direction != 0)
      ++cache[hidx].weight;
    else
      --cache[hidx].weight;

    // Renew the cache
    if (cache[hidx].weight < 0) {
      cache[hidx].weight = 1;
      cache[hidx].minimizers.resize(msize);
      if (msize == 0) {
        cache[hidx].offsets.resize(0);
	cache[hidx].strands.resize(0);
        return;
      }

      cache[hidx].offsets.resize(msize - 1);
      cache[hidx].strands.resize(msize);
      for (i = 0; i < msize; ++i)
      {
        cache[hidx].minimizers[i] = minimizers[i].first;
        cache[hidx].strands[i] = (minimizers[i].second & 1);
      }
      for (i = 0; i < msize - 1; ++i) {
        cache[hidx].offsets[i] = ((int)minimizers[i + 1].second>>1) - ((int)minimizers[i].second>>1);
      }
      /*std::vector<Candidate>().swap(cache[hidx].positive_candidates);
      std::vector<Candidate>().swap(cache[hidx].negative_candidates);
      cache[hidx].positive_candidates = pos_candidates;
      cache[hidx].negative_candidates = neg_candidates;*/
      cache[hidx].positive_candidate_list.Input(pos_candidates);
      cache[hidx].negative_candidate_list.Input(neg_candidates);
      cache[hidx].repetitive_seed_length = repetitive_seed_length;
      
      // adjust the candidate position.
      int shift = (int)minimizers[0].second>>1;
      cache[hidx].positive_candidate_list.Shift(shift);
      cache[hidx].negative_candidate_list.Shift(-shift);
    }
  }

  void DirectUpdateWeight(int idx, int weight) {
    cache[idx].weight += weight;
  }

  uint64_t GetMemoryBytes() {
    return 0;
    int i;
    uint64_t ret = 0;
    for (i = 0; i < cache_size; ++i) {
      ret += sizeof(cache[i]) + cache[i].minimizers.capacity() * sizeof(uint64_t) ;
        //+ cache[i].offsets.capacity() * sizeof(int)
        //+ cache[i].positive_candidates.capacity() * sizeof(Candidate) 
        //+ cache[i].negative_candidates.capacity() * sizeof(Candidate);
    }
    return ret;
  }

  void PrintStats() {
    for (int i = 0 ; i < cache_size ; ++i)
    {
      printf("%d %d %u ", cache[i].weight, cache[i].finger_print_cnt_sum, cache[i].positive_candidate_list.GetActualSize() + 
      			cache[i].negative_candidate_list.GetActualSize()) ;
      int tmp = 0;
      for (int j = 0 ; j < FINGER_PRINT_SIZE ; ++j)
        if (cache[i].finger_print_cnt[j] > tmp)
		tmp = cache[i].finger_print_cnt[j] ;
      printf("%d\n", tmp);
    }
  }
};
} // namespace chromap

#endif
