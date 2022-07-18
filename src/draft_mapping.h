#ifndef DRAFT_MAPPING_H_
#define DRAFT_MAPPING_H_

#include <stdint.h>

namespace chromap {

struct DraftMapping {
  int num_errors = 0;

  // The high 32 bits save the reference sequence index in the reference
  // sequence batch. The low 32 bits save the mapping end position on the
  // reference sequence.
  uint64_t position = 0;

  DraftMapping(int num_errors, int position)
      : num_errors(num_errors), position(position) {}

  inline int GetNumErrors() const { return num_errors; }

  inline uint32_t GetReferenceSequenceIndex() const { return (position >> 32); }

  inline uint32_t GetReferenceSequencePosition() const { return position; }
};

}  // namespace chromap

#endif  // DRAFT_MAPPING_H_
