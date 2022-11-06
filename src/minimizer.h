#ifndef MINIMIZER_H_
#define MINIMIZER_H_

#include <utility>

#include "hit.h"
#include "strand.h"

namespace chromap {

class Minimizer {
 public:
  Minimizer() = delete;

  Minimizer(std::pair<uint64_t, uint64_t> minimizer)
      : hash_key_(minimizer.first), hit_(minimizer.second) {}

  Minimizer(uint64_t hash_key, uint64_t hit) : hash_key_(hash_key), hit_(hit) {}

  ~Minimizer() = default;

  inline uint64_t GetHashKey() const { return hash_key_; }

  inline uint64_t GetHit() const { return hit_; }

  inline uint32_t GetSequenceIndex() const {
    return GenerateSequenceIndex(hit_);
  }

  inline uint32_t GetSequencePosition() const {
    return GenerateSequencePosition(hit_);
  }

  inline Strand GetSequenceStrand() const {
    return GenerateSequenceStrand(hit_);
  }

  inline bool operator<(const Minimizer &m) const {
    if (hash_key_ < m.hash_key_) {
      return true;
    }

    if (hash_key_ == m.hash_key_ && hit_ < m.hit_) {
      return true;
    }

    return false;
  }

 private:
  // The hash value of the kmer.
  uint64_t hash_key_ = 0;

  // The high 31 bits save the sequence index in the sequence batch. The
  // following 32 bits save the end position on that sequence. And the lowest
  // bit encodes the strand (0 for positive).
  uint64_t hit_ = 0;
};

}  // namespace chromap

#endif  // MINIMIZER_H_
