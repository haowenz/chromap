#ifndef MAPPING_PROCESSOR_H_
#define MAPPING_PROCESSOR_H_

#include <assert.h>

#include <cinttypes>
#include <cstring>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "bed_mapping.h"
#include "mapping.h"
#include "mapping_parameters.h"
#include "paf_mapping.h"
#include "pairs_mapping.h"
#include "sam_mapping.h"
#include "temp_mapping.h"
#include "utils.h"

namespace chromap {

template <typename MappingRecord>
bool ReadIdLess(const std::pair<uint32_t, MappingRecord> &a,
                const std::pair<uint32_t, MappingRecord> &b) {
  return a.second.read_id_ < b.second.read_id_;
}

// Class to process mappings. It supports multi-threadidng as only the
// parameters are owned by the class.
template <typename MappingRecord>
class MappingProcessor {
 public:
  MappingProcessor() = delete;
  MappingProcessor(const MappingParameters &mapping_parameters,
                   int min_unique_mapping_mapq)
      : min_unique_mapping_mapq_(min_unique_mapping_mapq),
        multi_mapping_allocation_seed_(
            mapping_parameters.multi_mapping_allocation_seed),
        multi_mapping_allocation_distance_(
            mapping_parameters.multi_mapping_allocation_distance),
        max_num_best_mappings_(mapping_parameters.max_num_best_mappings) {}

  ~MappingProcessor() = default;

  void SortOutputMappings(
      uint32_t num_reference_sequences,
      std::vector<std::vector<MappingRecord>> &mappings) const;
  
  void ParallelSortOutputMappings(
      uint32_t num_reference_sequences,
      std::vector<std::vector<MappingRecord>> &mappings,
      int num_threads) const;

  void RemovePCRDuplicate(
      uint32_t num_reference_sequences,
      std::vector<std::vector<MappingRecord>> &mappings,
      int num_threads) const;

  void AllocateMultiMappings(
      uint32_t num_reference_sequences, uint64_t num_multi_mappings,
      int multi_mapping_allocation_distance,
      std::vector<std::vector<MappingRecord>> &mappings) const;

  void ApplyTn5ShiftOnMappings(
      uint32_t num_reference_sequences,
      std::vector<std::vector<MappingRecord>> &mappings);

  uint32_t MoveMappingsInBuffersToMappingContainer(
      uint32_t num_reference_sequences,
      std::vector<std::vector<std::vector<MappingRecord>>>
          &mappings_on_diff_ref_seqs_for_diff_threads_for_saving,
      std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs);

  void OutputMappingStatistics(
      uint32_t num_reference_sequences,
      const std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs);

 private:
  void BuildAugmentedTree(
      uint32_t ref_id,
      std::vector<std::vector<MappingRecord>> &allocated_mappings,
      std::vector<std::pair<int, uint32_t>> &tree_info,
      std::vector<std::vector<uint32_t>> &tree_extras) const;

  uint32_t GetNumOverlappedMappings(
      uint32_t ref_id, int multi_mapping_allocation_distance,
      const MappingRecord &mapping,
      const std::vector<std::vector<MappingRecord>> &allocated_mappings,
      const std::vector<std::pair<int, uint32_t>> &tree_info,
      const std::vector<std::vector<uint32_t>> &tree_extras) const;

  const int min_unique_mapping_mapq_;
  const int multi_mapping_allocation_seed_;
  const int multi_mapping_allocation_distance_;
  const int max_num_best_mappings_;
};

template <typename MappingRecord>
void MappingProcessor<MappingRecord>::SortOutputMappings(
    uint32_t num_reference_sequences,
    std::vector<std::vector<MappingRecord>> &mappings) const {
  // double real_dedupe_start_time = Chromap<>::GetRealTime();
  uint32_t num_mappings = 0;
  for (uint32_t ri = 0; ri < num_reference_sequences; ++ri) {
    std::sort(mappings[ri].begin(), mappings[ri].end());
    num_mappings += mappings[ri].size();
  }
  // std::cerr << "Sorted " << num_mappings << " elements in " <<
  // Chromap<>::GetRealTime() - real_dedupe_start_time << "s.\n";
}

// If num_thread <= 0, then the number of thread is set externally
// Seems in this case omp task is much more efficient than omp parallel
template <typename MappingRecord>
void MappingProcessor<MappingRecord>::ParallelSortOutputMappings(
    uint32_t num_reference_sequences,
    std::vector<std::vector<MappingRecord>> &mappings,
    int num_threads) const {
  // double real_dedupe_start_time = Chromap<>::GetRealTime();
  if (num_threads <= 0)
  {
#pragma omp task shared(mappings)
    {
      for (uint32_t ri = 0; ri < num_reference_sequences; ri += 2) {
        std::sort(mappings[ri].begin(), mappings[ri].end());
      }
    }
    
#pragma omp task shared(mappings)
    {
      for (uint32_t ri = 1; ri < num_reference_sequences; ri += 2) {
        std::sort(mappings[ri].begin(), mappings[ri].end());
      }
    }
#pragma omp taskwait
  }
  else
  {
#pragma omp parallel shared(mappings) num_threads(num_threads)
    {
#pragma omp single
      {
#pragma omp taskloop
        for (uint32_t ri = 0; ri < num_reference_sequences; ++ri) {
          std::sort(mappings[ri].begin(), mappings[ri].end());
        }
      }
    }
  }
  
  //uint32_t num_mappings = 0;
  //for (uint32_t ri = 0; ri < num_reference_sequences; ++ri) {
  //  num_mappings += mappings[ri].size();
  //}
  // std::cerr << "Sorted " << num_mappings << " elements in " <<
  // Chromap<>::GetRealTime() - real_dedupe_start_time << "s.\n";
}

template <typename MappingRecord>
void MappingProcessor<MappingRecord>::RemovePCRDuplicate(
    uint32_t num_reference_sequences,
    std::vector<std::vector<MappingRecord>> &mappings,
    int num_threads) const {
  double real_dedupe_start_time = GetRealTime();
  ParallelSortOutputMappings(num_reference_sequences, mappings, num_threads);
  std::cerr << "Sorted in " << GetRealTime() - real_dedupe_start_time << "s.\n";

  std::vector<std::vector<MappingRecord>> deduped_mappings;
  uint32_t num_mappings = 0;
  for (uint32_t ri = 0; ri < num_reference_sequences; ++ri) {
    deduped_mappings.push_back(std::vector<MappingRecord>());
    if (mappings[ri].size() != 0) {
      // Haowen: Ideally I should output the last of the dups of first mappings.
      // Li: The mappings' mapq are sorted in increasing order, so we should put the last
      // map
      auto last_it = mappings[ri].begin();
      uint32_t last_dup_count = 1;

      for (auto it = ++(mappings[ri].begin()); it != mappings[ri].end(); ++it) {
        if (!((*it) == (*last_it))) {
          deduped_mappings[ri].emplace_back((*last_it));
          deduped_mappings[ri].back().num_dups_ = std::min(
              (uint32_t)std::numeric_limits<uint8_t>::max(), last_dup_count);
          last_dup_count = 1;
        } else {
          ++last_dup_count;
        }
        last_it = it;
      }

      deduped_mappings[ri].emplace_back((*last_it));
      deduped_mappings[ri].back().num_dups_ = std::min(
          (uint32_t)std::numeric_limits<uint8_t>::max(), last_dup_count);
      num_mappings += deduped_mappings[ri].size();
      deduped_mappings[ri].swap(mappings[ri]);
    }
  }
  std::cerr << num_mappings << " mappings left after deduplication in "
            << GetRealTime() - real_dedupe_start_time << "s.\n";
}

template <typename MappingRecord>
void MappingProcessor<MappingRecord>::BuildAugmentedTree(
    uint32_t ref_id,
    std::vector<std::vector<MappingRecord>> &allocated_mappings,
    std::vector<std::pair<int, uint32_t>> &tree_info,
    std::vector<std::vector<uint32_t>> &tree_extras) const {
  // std::sort(mappings.begin(), mappings.end(), IntervalLess());
  int max_level = 0;
  size_t i, last_i = 0;  // last_i points to the rightmost node in the tree
  uint32_t last = 0;     // last is the max value at node last_i
  int k;
  std::vector<MappingRecord> &mappings = allocated_mappings[ref_id];
  std::vector<uint32_t> &extras = tree_extras[ref_id];
  if (mappings.size() == 0) {
    max_level = -1;
  }
  for (i = 0; i < mappings.size(); i += 2) {
    last_i = i;
    // last = mappings[i].max = mappings[i].en; // leaves (i.e. at level 0)
    last = extras[i] =
        mappings[i].GetEndPosition();  // leaves (i.e. at level 0)
  }
  for (k = 1; 1LL << k <= (int64_t)mappings.size();
       ++k) {  // process internal nodes in the bottom-up order
    size_t x = 1LL << (k - 1);
    size_t i0 = (x << 1) - 1;
    size_t step = x << 2;  // i0 is the first node
    for (i = i0; i < mappings.size();
         i += step) {               // traverse all nodes at level k
      uint32_t el = extras[i - x];  // max value of the left child
      uint32_t er =
          i + x < mappings.size() ? extras[i + x] : last;  // of the right child
      uint32_t e = mappings[i].GetEndPosition();
      e = e > el ? e : el;
      e = e > er ? e : er;
      extras[i] = e;  // set the max value for node i
    }
    last_i =
        last_i >> k & 1
            ? last_i - x
            : last_i +
                  x;  // last_i now points to the parent of the original last_i
    if (last_i < mappings.size() &&
        extras[last_i] > last)  // update last accordingly
      last = extras[last_i];
  }
  max_level = k - 1;
  tree_info.emplace_back(max_level, mappings.size());
}

template <typename MappingRecord>
uint32_t MappingProcessor<MappingRecord>::GetNumOverlappedMappings(
    uint32_t ref_id, int multi_mapping_allocation_distance,
    const MappingRecord &mapping,
    const std::vector<std::vector<MappingRecord>> &allocated_mappings,
    const std::vector<std::pair<int, uint32_t>> &tree_info,
    const std::vector<std::vector<uint32_t>> &tree_extras) const {
  int t = 0;
  StackCell stack[64];
  // out.clear();
  int num_overlapped_mappings = 0;
  int max_level = tree_info[ref_id].first;
  uint32_t num_tree_nodes = tree_info[ref_id].second;
  const std::vector<MappingRecord> &mappings = allocated_mappings[ref_id];
  const std::vector<uint32_t> &extras = tree_extras[ref_id];
  // uint32_t interval_start = mapping.fragment_start_position;
  uint32_t interval_start =
      mapping.GetStartPosition() > (uint32_t)multi_mapping_allocation_distance
          ? mapping.GetStartPosition() - multi_mapping_allocation_distance
          : 0;
  uint32_t interval_end =
      mapping.GetEndPosition() + (uint32_t)multi_mapping_allocation_distance;
  // push the root; this is a top down traversal
  stack[t++] = StackCell(max_level, (1LL << max_level) - 1, 0);
  // the following guarantees that numbers in out[] are always sorted
  while (t) {
    StackCell z = stack[--t];
    // we are in a small subtree; traverse every node in this subtree
    if (z.k <= 3) {
      size_t i, i0 = z.x >> z.k << z.k, i1 = i0 + (1LL << (z.k + 1)) - 1;
      if (i1 >= num_tree_nodes) {
        i1 = num_tree_nodes;
      }
      for (i = i0; i < i1 && mappings[i].GetStartPosition() < interval_end;
           ++i) {
        if (interval_start <
            mappings[i].GetEndPosition()) {  // if overlap, append to out[]
          // out.push_back(i);
          ++num_overlapped_mappings;
        }
      }
    } else if (z.w == 0) {  // if left child not processed
      // the left child of z.x; NB: y may be out of range (i.e. y>=a.size())
      size_t y = z.x - (1LL << (z.k - 1));
      // re-add node z.x, but mark the left child having been processed
      stack[t++] = StackCell(z.k, z.x, 1);
      // push the left child if y is out of range or may overlap with the query
      if (y >= num_tree_nodes || extras[y] > interval_start)
        stack[t++] = StackCell(z.k - 1, y, 0);
    } else if (z.x < num_tree_nodes &&
               mappings[z.x].GetStartPosition() < interval_end) {
      // need to push the right child
      if (interval_start < mappings[z.x].GetEndPosition()) {
        // out.push_back(z.x);
        // test if z.x overlaps the query; if yes, append to out[]
        ++num_overlapped_mappings;
      }
      // push the right child
      stack[t++] = StackCell(z.k - 1, z.x + (1LL << (z.k - 1)), 0);
    }
  }
  return num_overlapped_mappings;
}

template <typename MappingRecord>
void MappingProcessor<MappingRecord>::AllocateMultiMappings(
    uint32_t num_reference_sequences, uint64_t num_multi_mappings,
    int multi_mapping_allocation_distance,
    std::vector<std::vector<MappingRecord>> &mappings) const {
  double real_start_time = GetRealTime();

  std::vector<std::pair<uint32_t, MappingRecord>> multi_mappings;
  multi_mappings.reserve(num_multi_mappings);

  std::vector<std::vector<MappingRecord>> allocated_mappings;
  allocated_mappings.reserve(num_reference_sequences);

  std::vector<std::pair<int, uint32_t>> tree_info;
  // max (max_level, # nodes)
  std::vector<std::vector<uint32_t>> tree_extras;
  tree_extras.reserve(num_reference_sequences);
  tree_info.reserve(num_reference_sequences);

  // two passes, one for memory pre-allocation, another to move the mappings.
  for (uint32_t ri = 0; ri < num_reference_sequences; ++ri) {
    allocated_mappings.emplace_back(std::vector<MappingRecord>());
    tree_extras.emplace_back(std::vector<uint32_t>());
    uint32_t num_uni_mappings = 0;
    uint32_t num_multi_mappings = 0;
    for (uint32_t mi = 0; mi < mappings[ri].size(); ++mi) {
      MappingRecord &mapping = mappings[ri][mi];
      if ((mapping.mapq_) <
          min_unique_mapping_mapq_) {  // we have to ensure that the mapq is
                                       // lower than this if and only if it is a
                                       // multi-read.
        ++num_multi_mappings;
      } else {
        ++num_uni_mappings;
      }
    }
    allocated_mappings[ri].reserve(num_uni_mappings);
    tree_extras[ri].reserve(num_uni_mappings);
    for (uint32_t mi = 0; mi < mappings[ri].size(); ++mi) {
      MappingRecord &mapping = mappings[ri][mi];
      if ((mapping.mapq_) < min_unique_mapping_mapq_) {
        multi_mappings.emplace_back(ri, mapping);
      } else {
        allocated_mappings[ri].emplace_back(mapping);
        tree_extras[ri].emplace_back(0);
      }
    }
    std::vector<MappingRecord>().swap(mappings[ri]);
    BuildAugmentedTree(ri, allocated_mappings, tree_info, tree_extras);
  }
  std::cerr << "Got all " << multi_mappings.size() << " multi-mappings!\n";

  std::stable_sort(multi_mappings.begin(), multi_mappings.end(),
                   ReadIdLess<MappingRecord>);
  std::vector<uint32_t> weights;
  weights.reserve(max_num_best_mappings_);
  uint32_t sum_weight = 0;
  assert(multi_mappings.size() > 0);
  uint32_t previous_read_id = multi_mappings[0].second.read_id_;
  uint32_t start_mapping_index = 0;
  // add a fake mapping at the end and make sure its id is different from the
  // last one
  assert(multi_mappings.size() != UINT32_MAX);
  std::pair<uint32_t, MappingRecord> foo_mapping = multi_mappings.back();
  foo_mapping.second.read_id_ = UINT32_MAX;
  multi_mappings.emplace_back(foo_mapping);
  std::mt19937 generator(multi_mapping_allocation_seed_);
  uint32_t current_read_id;  //, reference_id, mapping_index;
  // uint32_t allocated_read_id, allocated_reference_id,
  // allocated_mapping_index;
  uint32_t num_allocated_multi_mappings = 0;
  uint32_t num_multi_mappings_without_overlapping_unique_mappings = 0;
  for (uint32_t mi = 0; mi < multi_mappings.size(); ++mi) {
    std::pair<uint32_t, MappingRecord> &current_multi_mapping =
        multi_mappings[mi];  // mappings[reference_id][mapping_index];
    current_read_id = current_multi_mapping.second.read_id_;
    uint32_t num_overlaps = GetNumOverlappedMappings(
        current_multi_mapping.first, multi_mapping_allocation_distance,
        current_multi_mapping.second, allocated_mappings, tree_info,
        tree_extras);
    // std::cerr << mi << " " << current_read_id << " " << previous_read_id << "
    // " << reference_id << " " << mapping_index << " " << interval_start << " "
    // << num_overlaps << " " << sum_weight << "\n";
    if (current_read_id == previous_read_id) {
      weights.emplace_back(num_overlaps);
      sum_weight += num_overlaps;
    } else {
      // deal with the previous one.
      if (sum_weight == 0) {
        ++num_multi_mappings_without_overlapping_unique_mappings;
        // assert(weights.size() > 1); // After PCR dedupe, some multi-reads may
        // become uni-reads. For now, we just assign it to that unique mapping
        // positions. std::fill(weights.begin(), weights.end(), 1); // We drop
        // the multi-mappings that have no overlap with uni-mappings.
      } else {
        std::discrete_distribution<uint32_t> distribution(weights.begin(),
                                                          weights.end());
        uint32_t randomly_assigned_mapping_index = distribution(generator);
        allocated_mappings[multi_mappings[start_mapping_index +
                                          randomly_assigned_mapping_index]
                               .first]
            .emplace_back(multi_mappings[start_mapping_index +
                                         randomly_assigned_mapping_index]
                              .second);
        ++num_allocated_multi_mappings;
      }
      // update current
      weights.clear();
      weights.emplace_back(num_overlaps);
      sum_weight = num_overlaps;
      start_mapping_index = mi;
      previous_read_id = current_read_id;
    }
  }

  mappings.swap(allocated_mappings);

  std::cerr << "Allocated " << num_allocated_multi_mappings
            << " multi-mappings in " << GetRealTime() - real_start_time
            << "s.\n";
  std::cerr << "# multi-mappings that have no uni-mapping overlaps: "
            << num_multi_mappings_without_overlapping_unique_mappings << ".\n";
}

template <typename MappingRecord>
void MappingProcessor<MappingRecord>::ApplyTn5ShiftOnMappings(
    uint32_t num_reference_sequences,
    std::vector<std::vector<MappingRecord>> &mappings) {
  uint64_t num_shifted_mappings = 0;
  for (auto &mappings_on_one_ref_seq : mappings) {
    for (auto &mapping : mappings_on_one_ref_seq) {
      mapping.Tn5Shift();
      ++num_shifted_mappings;
    }
  }
  std::cerr << "# shifted mappings: " << num_shifted_mappings << ".\n";
}

template <typename MappingRecord>
uint32_t
MappingProcessor<MappingRecord>::MoveMappingsInBuffersToMappingContainer(
    uint32_t num_reference_sequences,
    std::vector<std::vector<std::vector<MappingRecord>>>
        &mappings_on_diff_ref_seqs_for_diff_threads_for_saving,
    std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs) {
  // double real_start_time = Chromap<>::GetRealTime();
  uint32_t num_moved_mappings = 0;
  for (size_t ti = 0;
       ti < mappings_on_diff_ref_seqs_for_diff_threads_for_saving.size();
       ++ti) {
    for (uint32_t i = 0; i < num_reference_sequences; ++i) {
      num_moved_mappings +=
          mappings_on_diff_ref_seqs_for_diff_threads_for_saving[ti][i].size();
      mappings_on_diff_ref_seqs[i].insert(
          mappings_on_diff_ref_seqs[i].end(),
          std::make_move_iterator(
              mappings_on_diff_ref_seqs_for_diff_threads_for_saving[ti][i]
                  .begin()),
          std::make_move_iterator(
              mappings_on_diff_ref_seqs_for_diff_threads_for_saving[ti][i]
                  .end()));
      mappings_on_diff_ref_seqs_for_diff_threads_for_saving[ti][i].clear();
    }
  }
  // std::cerr << "Moved mappings in " << Chromap<>::GetRealTime() -
  // real_start_time << "s.\n";
  return num_moved_mappings;
}

template <typename MappingRecord>
void MappingProcessor<MappingRecord>::OutputMappingStatistics(
    uint32_t num_reference_sequences,
    const std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs) {
  uint64_t num_uni_mappings = 0;
  uint64_t num_multi_mappings = 0;
  for (auto &mappings_on_one_ref_seq : mappings_on_diff_ref_seqs) {
    for (auto &mapping : mappings_on_one_ref_seq) {
      if ((mapping.is_unique_) == 1) {
        ++num_uni_mappings;
      } else {
        ++num_multi_mappings;
      }
    }
  }
  std::cerr << "# uni-mappings: " << num_uni_mappings
            << ", # multi-mappings: " << num_multi_mappings
            << ", total: " << num_uni_mappings + num_multi_mappings << ".\n";
}

}  // namespace chromap

#endif  // MAPPING_PROCESSOR_H_
