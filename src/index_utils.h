#ifndef INDEX_UTILS_H_
#define INDEX_UTILS_H_

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

struct RepetitiveSeedStats {
  uint32_t repetitive_seed_length = 0;
  uint32_t previous_repetitive_seed_position =
      std::numeric_limits<uint32_t>::max();
  int repetitive_seed_count = 0;
};

inline static uint64_t GenerateHashInLookupTable(uint64_t minimizer_hash) {
  return minimizer_hash << 1;
}

inline static uint64_t GenerateEntryValueInLookupTable(
    uint64_t occurrence_table_offset, uint32_t num_occurrences) {
  return (occurrence_table_offset << 32) | num_occurrences;
}

inline static uint32_t GenerateOffsetInOccurrenceTable(uint64_t lookup_value) {
  return lookup_value >> 32;
}

inline static uint32_t GenerateNumOccurrenceInOccurrenceTable(
    uint64_t lookup_table_entry_value) {
  return static_cast<uint32_t>(lookup_table_entry_value);
}

inline static uint64_t SequenceIndexAndPositionToCandidatePosition(
    uint64_t sequence_id, uint32_t sequence_position) {
  return (sequence_id << 32) | sequence_position;
}

inline static uint64_t GenerateCandidatePositionFromOccurrenceTableEntry(
    uint64_t entry) {
  return entry >> 1;
}

inline static bool IsSingletonLookupKey(uint64_t lookup_key) {
  return (lookup_key & 1) > 0;
}

// Only used in Index to merge sorted candidate position lists using heap.
struct CandidatePositionWithListIndex {
  uint32_t list_index;
  uint64_t position;

  CandidatePositionWithListIndex(uint32_t list_index, uint64_t position)
      : list_index(list_index), position(position) {}

  bool operator<(const CandidatePositionWithListIndex &h) const {
    // The inversed direction is to make a min-heap.
    return position > h.position;
  }
};

inline static void HeapMergeCandidatePositionLists(
    const std::vector<std::vector<uint64_t>> sorted_candidate_position_lists,
    std::vector<uint64_t> &candidate_positions) {
  std::priority_queue<CandidatePositionWithListIndex> heap;
  std::vector<uint32_t> candidate_position_list_indices(
      sorted_candidate_position_lists.size(), 0);

  for (uint32_t li = 0; li < sorted_candidate_position_lists.size(); ++li) {
    if (sorted_candidate_position_lists[li].size() == 0) {
      continue;
    }
    heap.emplace(li, sorted_candidate_position_lists[li][0]);
  }

  while (!heap.empty()) {
    const CandidatePositionWithListIndex min_candidate_position = heap.top();
    heap.pop();
    candidate_positions.push_back(min_candidate_position.position);
    ++candidate_position_list_indices[min_candidate_position.list_index];

    const uint32_t min_candidate_position_list_index =
        candidate_position_list_indices[min_candidate_position.list_index];
    const std::vector<uint64_t> &min_sorted_candidate_position_list =
        sorted_candidate_position_lists[min_candidate_position.list_index];
    if (min_candidate_position_list_index <
        min_sorted_candidate_position_list.size()) {
      heap.emplace(min_candidate_position.list_index,
                   min_sorted_candidate_position_list
                       [min_candidate_position_list_index]);
    }
  }
}

}  // namespace chromap

#endif  // INDEX_UTILS_H_
