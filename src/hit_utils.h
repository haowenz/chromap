#ifndef HIT_UTILS_H_
#define HIT_UTILS_H_

#include "strand.h"

namespace chromap {

inline static uint32_t HitToSequenceIndex(uint64_t hit) { return (hit >> 33); }

inline static uint32_t HitToSequencePosition(uint64_t hit) {
  return (hit >> 1);
}

inline static Strand HitToStrand(uint64_t hit) {
  if ((hit & 1) == 0) {
    return kPositive;
  }
  return kNegative;
}

inline static bool AreTwoHitsOnTheSameStrand(uint64_t hit1, uint64_t hit2) {
  return ((hit1 & 1) == (hit2 & 1));
}

}  // namespace chromap

#endif  // HIT_UTILS_H_
