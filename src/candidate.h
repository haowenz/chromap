#ifndef CANDIDATE_H_
#define CANDIDATE_H_

namespace chromap {

struct Candidate {
  // The high 32 bits save the reference sequence index in the reference
  // sequence batch. The low 32 bits save the reference position on that
  // sequence.
  uint64_t position = 0;

  // The number of minimizers supports the position.
  uint8_t count = 0;

  inline uint32_t GetReferenceSequenceIndex() const { return (position >> 32); }

  inline uint32_t GetReferenceSequencePosition() const { return position; }

  inline bool operator<(const Candidate &c) const {
    if (count > c.count) {
      return true;
    }

    return position < c.position;
  }
};

}  // namespace chromap

#endif  // CANDIDATE_H_
