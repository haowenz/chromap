#ifndef CANDIDATE_H_
#define CANDIDATE_H_

namespace chromap {

struct Candidate {
  uint64_t position;
  uint8_t count;

  bool operator<(const Candidate &c) const {
    if (count != c.count) {
      return count > c.count;
    }

    return position < c.position;
  }
};

}  // namespace chromap

#endif  // CANDIDATE_H_
