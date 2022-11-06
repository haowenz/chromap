#ifndef MINIMIZER_UTILS_H_
#define MINIMIZER_UTILS_H_

#include "strand.h"

namespace chromap {

inline static uint32_t GenerateSequenceIndex(uint64_t seed_hit) {
  return (seed_hit >> 33);
}

inline static uint32_t GenerateSequencePosition(uint64_t seed_hit) {
  return (seed_hit >> 1);
}

inline static Strand GenerateSequenceStrand(uint64_t seed_hit) {
  if ((seed_hit & 1) == 0) {
    return kPositive;
  }
  return kNegative;
}

inline static bool AreTwoHitsOnTheSameStrand(uint64_t seed_hit1,
                                             uint64_t seed_hit2) {
  return ((seed_hit1 & 1) == (seed_hit2 & 1));
}

}  // namespace chromap

#endif  // MINIMIZER_UTILS_H_
