#ifndef INDEX_UTILS_H_
#define INDEX_UTILSH_

#include "khash.h"

// Note that the max kmer size is 28 and its hash value is always saved in the
// lowest 56 bits of an unsigned 64-bit integer. When an element is inserted
// into the hash table, its hash value is left shifted by 1 bit and the lowest
// bit of the key value is set to 1 when the minimizer only occurs once. So
// right shift by one bit is lossless and safe.
#define KHashFunctionForIndex(a) ((a) >> 1)
#define KHashEqForIndex(a, b) ((a) >> 1 == (b) >> 1)
KHASH_INIT(/*name=*/k64, /*khkey_t=*/uint64_t, /*khval_t=*/uint64_t,
           /*kh_is_map=*/1, /*__hash_func=*/KHashFunctionForIndex,
           /*__hash_equal=*/KHashEqForIndex);

namespace chromap {

inline static uint64_t GenerateHashKeyInLookupTable(
    uint64_t minimizer_hash_key) {
  return minimizer_hash_key << 1;
}

inline static uint64_t GenerateEntryValueInLookupTable(
    uint64_t occurrence_table_offset, uint32_t num_occurrences) {
  return (occurrence_table_offset << 32) | num_occurrences;
}

inline static uint32_t GenerateOffsetInOccurrenceTable(
    uint64_t lookup_table_entry_value) {
  return lookup_table_entry_value >> 32;
}

inline static uint32_t GenerateNumOccurrenceInOccurrenceTable(
    uint64_t lookup_table_entry_value) {
  return static_cast<uint32_t>(lookup_table_entry_value);
}

inline static uint64_t GenerateCandidatePosition(uint64_t sequence_id,
                                                 uint32_t sequence_position) {
  return (sequence_id << 32) | sequence_position;
}

inline static uint64_t GenerateCandidatePositionFromOccurrenceTableEntry(
    uint64_t entry) {
  return entry >> 1;
}

}  // namespace chromap

#endif  // INDEX_UTILS_H_
