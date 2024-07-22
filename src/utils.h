#ifndef UTILS_H_
#define UTILS_H_

#include <sys/resource.h>
#include <sys/time.h>

#include <iostream>
#include <tuple>
#include <vector>

#include "candidate.h"
#include "khash.h"
#include "minimizer.h"
#include "strand.h"

namespace chromap {

struct uint128_t {
  uint64_t first;
  uint64_t second;
};

struct BarcodeWithQual {
  uint32_t corrected_base_index1;
  char correct_base1;
  uint32_t corrected_base_index2;
  char correct_base2;
  double score;
  bool operator>(const BarcodeWithQual &b) const {
    return std::tie(score, corrected_base_index1, correct_base1,
                    corrected_base_index2, correct_base2) >
           std::tie(b.score, b.corrected_base_index1, b.correct_base1,
                    b.corrected_base_index2, b.correct_base2);
  }
};

struct _mm_history {
  unsigned int timestamp = 0;
  std::vector<Minimizer> minimizers;
  std::vector<Candidate> positive_candidates;
  std::vector<Candidate> negative_candidates;
  uint32_t repetitive_seed_length;
};

KHASH_MAP_INIT_INT64(k128, uint128_t);
KHASH_MAP_INIT_INT64(k64_seq, uint64_t);
KHASH_SET_INIT_INT(k32_set);
KHASH_MAP_INIT_INT64(kmatrix, uint32_t);

struct StackCell {
  size_t x;  // node
  int k, w;  // k: level; w: 0 if left child hasn't been processed
  StackCell(){};
  StackCell(int k_, size_t x_, int w_) : x(x_), k(k_), w(w_){};
};

inline static double GetRealTime() {
  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp, &tzp);
  return tp.tv_sec + tp.tv_usec * 1e-6;
}

inline static double GetCPUTime() {
  struct rusage r;
  getrusage(RUSAGE_SELF, &r);
  return r.ru_utime.tv_sec + r.ru_stime.tv_sec +
         1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

inline static void ExitWithMessage(const std::string &message) {
  std::cerr << message << std::endl;
  exit(-1);
}

inline static uint64_t Hash64(uint64_t key, const uint64_t mask) {
  key = (~key + (key << 21)) & mask;  // key = (key << 21) - key - 1;
  key = key ^ key >> 24;
  key = ((key + (key << 3)) + (key << 8)) & mask;  // key * 265
  key = key ^ key >> 14;
  key = ((key + (key << 2)) + (key << 4)) & mask;  // key * 21
  key = key ^ key >> 28;
  key = (key + (key << 31)) & mask;
  return key;
}

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

inline static uint8_t CharToUint8(const char c) {
  return char_to_uint8_table_[(uint8_t)c];
}

inline static char Uint8ToChar(const uint8_t i) {
  return uint8_to_char_table_[i];
}

// Make sure the length is not greater than 32 before calling this function.
inline static uint64_t GenerateSeedFromSequence(const char *sequence,
                                                uint32_t sequence_length,
                                                uint32_t start_position,
                                                uint32_t seed_length) {
  uint64_t seed = 0;
  for (uint32_t i = 0; i < seed_length; ++i) {
    if (start_position + i < sequence_length) {
      uint8_t current_base = CharToUint8(sequence[i + start_position]);
      if (current_base < 4) {               // not an ambiguous base
        seed = (seed << 2) | current_base;  // forward k-mer
      } else {
        seed = seed << 2;  // N->A
      }
    } else {
      seed = seed << 2;  // Pad A
    }
  }
  return seed;
}

inline static uint64_t GenerateMinimizer(uint32_t sequence_index,
                                         uint32_t sequence_position,
                                         const Strand strand) {
  const uint64_t minimizer =
      (((uint64_t)sequence_index) << 32 | sequence_position) << 1;
  return minimizer | (strand == kPositive ? 0 : 1);
}

}  // namespace chromap

#endif  // UTILS_H_
