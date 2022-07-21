#ifndef MINIMIZER_H_
#define MINIMIZER_H_

#include <utility>

#include "strand.h"

namespace chromap {

class Minimizer {
 public:
  Minimizer() = delete;

  Minimizer(std::pair<uint64_t, uint64_t> minimizer)
      : hash_key_(minimizer.first), minimizer_(minimizer.second) {}

  Minimizer(uint64_t hash_key, uint64_t minimizer)
      : hash_key_(hash_key), minimizer_(minimizer) {}

  ~Minimizer() = default;

  inline uint64_t GetHashKey() const { return hash_key_; }

  inline uint64_t GetMinimizer() const { return minimizer_; }

  inline uint32_t GetSequenceIndex() const { return (minimizer_ >> 33); }

  inline uint32_t GetSequencePosition() const { return (minimizer_ >> 1); }

  inline Strand GetSequenceStrand() const {
    if ((minimizer_ & 1) == 0) {
      return kPositive;
    }
    return kNegative;
  }

  inline bool operator<(const Minimizer &m) const {
    if (hash_key_ < m.hash_key_) {
      return true;
    }

    if (hash_key_ == m.hash_key_ && minimizer_ < m.minimizer_) {
      return true;
    }

    return false;
  }

 private:
  // The hash value of the kmer.
  uint64_t hash_key_ = 0;

  // The high 31 bits save the reference sequence index in the reference
  // sequence batch. The following 32 bits save the reference position on that
  // sequence. And the lowest bit encodes the strand (0 for positive).
  uint64_t minimizer_ = 0;
};

}  // namespace chromap

#endif  // MINIMIZER_H_
