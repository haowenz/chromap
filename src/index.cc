#include "index.h"

#include <assert.h>

#include <algorithm>
#include <iostream>

#include "minimizer_generator.h"

namespace chromap {

void Index::Construct(uint32_t num_sequences, const SequenceBatch &reference) {
  const double real_start_time = GetRealTime();

  std::vector<Minimizer> minimizers;
  minimizers.reserve(reference.GetNumBases() / window_size_ * 2);

  std::cerr << "Collecting minimizers.\n";
  MinimizerGenerator minimizer_generator(kmer_size_, window_size_);
  for (uint32_t sequence_index = 0; sequence_index < num_sequences;
       ++sequence_index) {
    minimizer_generator.GenerateMinimizers(reference, sequence_index,
                                           minimizers);
  }
  std::cerr << "Collected " << minimizers.size() << " minimizers.\n";

  std::cerr << "Sorting minimizers.\n";
  std::stable_sort(minimizers.begin(), minimizers.end());
  std::cerr << "Sorted all minimizers.\n";

  const uint32_t num_minimizers = minimizers.size();

  assert(num_minimizers > 0);
  // TODO: check this assert!
  // Here I make sure the # minimizers is less than the limit of signed int32,
  // so that I can use int to store position later.
  assert(num_minimizers <= INT_MAX);

  occurrence_table_.reserve(num_minimizers);

  uint64_t previous_key = minimizers[0].GetHashKey();
  uint32_t num_previous_minimizer_occurrences = 0;
  uint64_t num_nonsingletons = 0;
  uint32_t num_singletons = 0;

  for (uint32_t ti = 0; ti <= num_minimizers; ++ti) {
    const bool is_last_iteration = ti == num_minimizers;
    const uint64_t current_key =
        is_last_iteration ? previous_key + 1 : minimizers[ti].GetHashKey();

    if (current_key != previous_key) {
      int khash_return_code = 0;
      khiter_t khash_iterator =
          kh_put(k64, lookup_table_, previous_key << 1, &khash_return_code);
      assert(khash_return_code != -1 && khash_return_code != 0);

      if (num_previous_minimizer_occurrences == 1) {
        // We set the lowest bit of the key value to 1 if the minimizer only
        // occurs once. And the occurrence is directly saved in the lookup
        // table.
        kh_key(lookup_table_, khash_iterator) |= 1;
        kh_value(lookup_table_, khash_iterator) = occurrence_table_.back();
        occurrence_table_.pop_back();
        ++num_singletons;
      } else {
        kh_value(lookup_table_, khash_iterator) =
            (num_nonsingletons << 32) | num_previous_minimizer_occurrences;
        num_nonsingletons += num_previous_minimizer_occurrences;
      }
      num_previous_minimizer_occurrences = 1;
    } else {
      num_previous_minimizer_occurrences++;
    }

    if (is_last_iteration) {
      break;
    }

    occurrence_table_.push_back(minimizers[ti].GetMinimizer());
    previous_key = current_key;
  }

  assert(num_nonsingletons + num_singletons == num_minimizers);

  std::cerr << "Kmer size: " << kmer_size_ << ", window size: " << window_size_
            << ".\n";
  std::cerr << "Lookup table size: " << kh_size(lookup_table_)
            << ", # buckets: " << kh_n_buckets(lookup_table_)
            << ", occurrence table size: " << occurrence_table_.size()
            << ", # singletons: " << num_singletons << ".\n";
  std::cerr << "Built index successfully in " << GetRealTime() - real_start_time
            << "s.\n";
}

void Index::Save() const {
  const double real_start_time = GetRealTime();
  FILE *index_file = fopen(index_file_path_.c_str(), "wb");
  assert(index_file != nullptr);

  uint64_t num_bytes = 0;
  int err = 0;

  err = fwrite(&kmer_size_, sizeof(int), 1, index_file);
  num_bytes += sizeof(int);
  assert(err != 0);

  err = fwrite(&window_size_, sizeof(int), 1, index_file);
  num_bytes += sizeof(int);
  assert(err != 0);

  const uint32_t lookup_table_size = kh_size(lookup_table_);
  err = fwrite(&lookup_table_size, sizeof(uint32_t), 1, index_file);
  num_bytes += sizeof(uint32_t);
  assert(err != 0);

  kh_save(k64, lookup_table_, index_file);
  num_bytes += sizeof(uint64_t) * 2 * lookup_table_size;

  const uint32_t occurrence_table_size = occurrence_table_.size();
  err = fwrite(&occurrence_table_size, sizeof(uint32_t), 1, index_file);
  num_bytes += sizeof(uint32_t);
  assert(err != 0);

  if (occurrence_table_size > 0) {
    err = fwrite(occurrence_table_.data(), sizeof(uint64_t),
                 occurrence_table_size, index_file);
    num_bytes += sizeof(uint64_t) * occurrence_table_size;
    assert(err != 0);
  }

  fclose(index_file);
  // std::cerr << "Index size: " << num_bytes / (1024.0 * 1024 * 1024) << "GB,
  std::cerr << "Saved in " << GetRealTime() - real_start_time << "s.\n";
}

void Index::Load() {
  const double real_start_time = GetRealTime();
  FILE *index_file = fopen(index_file_path_.c_str(), "rb");
  assert(index_file != nullptr);

  int err = 0;
  err = fread(&kmer_size_, sizeof(int), 1, index_file);
  assert(err != 0);

  err = fread(&window_size_, sizeof(int), 1, index_file);
  assert(err != 0);

  uint32_t lookup_table_size = 0;
  err = fread(&lookup_table_size, sizeof(uint32_t), 1, index_file);
  assert(err != 0);

  kh_load(k64, lookup_table_, index_file);

  uint32_t occurrence_table_size = 0;
  err = fread(&occurrence_table_size, sizeof(uint32_t), 1, index_file);
  assert(err != 0);

  if (occurrence_table_size > 0) {
    occurrence_table_.resize(occurrence_table_size);
    err = fread(occurrence_table_.data(), sizeof(uint64_t),
                occurrence_table_size, index_file);
    assert(err != 0);
  }

  fclose(index_file);

  std::cerr << "Kmer size: " << kmer_size_ << ", window size: " << window_size_
            << ".\n";
  std::cerr << "Lookup table size: " << kh_size(lookup_table_)
            << ", occurrence table size: " << occurrence_table_.size() << ".\n";
  std::cerr << "Loaded index successfully in "
            << GetRealTime() - real_start_time << "s.\n";
}

void Index::Statistics(uint32_t num_sequences,
                       const SequenceBatch &reference) const {
  double real_start_time = GetRealTime();
  int n = 0, n1 = 0;
  uint32_t i;
  uint64_t sum = 0, len = 0;
  fprintf(stderr, "[M::%s] kmer size: %d; skip: %d; #seq: %d\n", __func__,
          kmer_size_, window_size_, num_sequences);
  for (i = 0; i < num_sequences; ++i) {
    len += reference.GetSequenceLengthAt(i);
  }
  assert(len == reference.GetNumBases());
  if (lookup_table_) {
    n += kh_size(lookup_table_);
  }
  for (khint_t k = 0; k < kh_end(lookup_table_); ++k) {
    if (kh_exist(lookup_table_, k)) {
      sum +=
          kh_key(lookup_table_, k) & 1 ? 1 : (uint32_t)kh_val(lookup_table_, k);
      if (kh_key(lookup_table_, k) & 1) ++n1;
    }
  }
  fprintf(stderr,
          "[M::%s::%.3f] distinct minimizers: %d (%.2f%% are singletons); "
          "average occurrences: %.3lf; average spacing: %.3lf\n",
          __func__, GetRealTime() - real_start_time, n, 100.0 * n1 / n,
          (double)sum / n, (double)len / sum);
}

void Index::CheckIndex(uint32_t num_sequences,
                       const SequenceBatch &reference) const {
  std::vector<Minimizer> minimizers;
  minimizers.reserve(reference.GetNumBases() / window_size_ * 2);

  MinimizerGenerator minimizer_generator(kmer_size_, window_size_);
  for (uint32_t sequence_index = 0; sequence_index < num_sequences;
       ++sequence_index) {
    minimizer_generator.GenerateMinimizers(reference, sequence_index,
                                           minimizers);
  }
  std::cerr << "Collected " << minimizers.size() << " minimizers.\n";

  std::stable_sort(minimizers.begin(), minimizers.end());
  std::cerr << "Sorted minimizers.\n";

  uint32_t count = 0;
  for (uint32_t i = 0; i < minimizers.size(); ++i) {
    khiter_t khash_iterator =
        kh_get(k64, lookup_table_, minimizers[i].GetHashKey() << 1);
    assert(khash_iterator != kh_end(lookup_table_));
    uint64_t key = kh_key(lookup_table_, khash_iterator);
    uint64_t value = kh_value(lookup_table_, khash_iterator);
    if (key & 1) {  // singleton
      assert(minimizers[i].GetMinimizer() == value);
      count = 0;
    } else {
      uint32_t offset = value >> 32;
      uint32_t num_occ = value;
      uint64_t value_in_index = occurrence_table_[offset + count];
      assert(value_in_index == minimizers[i].GetMinimizer());
      ++count;
      if (count == num_occ) {
        count = 0;
      }
    }
  }
}

void Index::HeapMergeSeedHitLists(
    const std::vector<std::vector<uint64_t>> sorted_seed_hit_lists,
    std::vector<uint64_t> &seed_hits) const {
  std::priority_queue<SeedHitInList> heap;
  std::vector<uint32_t> seed_hit_list_indices(sorted_seed_hit_lists.size(), 0);

  for (uint32_t li = 0; li < sorted_seed_hit_lists.size(); ++li) {
    if (sorted_seed_hit_lists[li].size() == 0) {
      continue;
    }
    heap.emplace(li, sorted_seed_hit_lists[li][0]);
  }

  while (!heap.empty()) {
    const SeedHitInList min_seed_hit = heap.top();
    heap.pop();
    seed_hits.push_back(min_seed_hit.position);
    ++seed_hit_list_indices[min_seed_hit.list_index];

    const uint32_t min_seed_hit_list_index =
        seed_hit_list_indices[min_seed_hit.list_index];
    const std::vector<uint64_t> &min_sorted_seed_hit_list =
        sorted_seed_hit_lists[min_seed_hit.list_index];
    if (min_seed_hit_list_index < min_sorted_seed_hit_list.size()) {
      heap.emplace(min_seed_hit.list_index,
                   min_sorted_seed_hit_list[min_seed_hit_list_index]);
    }
  }
}

int Index::CollectSeedHits(int max_seed_frequency,
                           int repetitive_seed_frequency,
                           const std::vector<Minimizer> &minimizers,
                           uint32_t &repetitive_seed_length,
                           std::vector<uint64_t> &positive_hits,
                           std::vector<uint64_t> &negative_hits,
                           bool use_heap) const {
  const uint32_t num_minimizers = minimizers.size();

  std::vector<std::vector<uint64_t>> mm_positive_hits;
  std::vector<std::vector<uint64_t>> mm_negative_hits;

  if (use_heap) {
    for (uint32_t i = 0; i < num_minimizers; ++i) {
      mm_positive_hits.emplace_back(std::vector<uint64_t>());
      mm_negative_hits.emplace_back(std::vector<uint64_t>());
    }
  }

  bool heap_resort = false;  // need to sort the elements of heap first

  positive_hits.reserve(max_seed_frequency * 2);
  negative_hits.reserve(max_seed_frequency * 2);

  uint32_t previous_repetitive_seed_position =
      std::numeric_limits<uint32_t>::max();

  int repetitive_seed_count = 0;

  for (uint32_t mi = 0; mi < num_minimizers; ++mi) {
    khiter_t khash_iterator =
        kh_get(k64, lookup_table_, minimizers[mi].GetHashKey() << 1);
    if (khash_iterator == kh_end(lookup_table_)) {
      // std::cerr << "The minimizer is not in reference!\n";
      continue;
    }

    const uint64_t value = kh_value(lookup_table_, khash_iterator);

    const uint32_t read_position = minimizers[mi].GetSequencePosition();
    const Strand read_strand = minimizers[mi].GetSequenceStrand();

    const bool is_reference_minimizer_single =
        (kh_key(lookup_table_, khash_iterator) & 1) > 0;

    if (is_reference_minimizer_single) {
      const uint64_t reference_id = value >> 33;
      const uint32_t reference_position = value >> 1;
      const Strand reference_strand = (value & 1) == 0 ? kPositive : kNegative;

      // Check whether the strands of reference minimizer and read minimizer are
      // the same. Later, we can play some tricks with 0,1 here to make it
      // faster.
      if (read_strand == reference_strand) {
        const uint32_t candidate_position = reference_position - read_position;
        // Ok, for now we can't see the reference here. So let us don't validate
        // this candidate. Instead, we do it later some time when we check the
        // candidates.
        const uint64_t seed_hit = (reference_id << 32) | candidate_position;

        if (use_heap) {
          mm_positive_hits[mi].push_back(seed_hit);
        } else {
          positive_hits.push_back(seed_hit);
        }
      } else {
        const uint32_t candidate_position =
            reference_position + read_position - kmer_size_ + 1;
        const uint64_t seed_hit = (reference_id << 32) | candidate_position;

        if (use_heap) {
          mm_negative_hits[mi].push_back(seed_hit);
        } else {
          negative_hits.push_back(seed_hit);
        }
      }
      continue;
    }

    const uint32_t offset = value >> 32;
    const uint32_t num_occurrences = value;

    if (num_occurrences < (uint32_t)max_seed_frequency) {
      for (uint32_t oi = 0; oi < num_occurrences; ++oi) {
        const uint64_t value = occurrence_table_[offset + oi];
        const uint64_t reference_id = value >> 33;
        const uint32_t reference_position = value >> 1;
        const Strand reference_strand =
            (value & 1) == 0 ? kPositive : kNegative;

        if (read_strand == reference_strand) {
          const uint32_t candidate_position =
              reference_position - read_position;
          const uint64_t seed_hit = (reference_id << 32) | candidate_position;

          if (use_heap) {
            if (reference_position < read_position) {
              heap_resort = true;
            }
            mm_positive_hits[mi].push_back(seed_hit);
          } else {
            positive_hits.push_back(seed_hit);
          }
        } else {
          const uint32_t candidate_position =
              reference_position + read_position - kmer_size_ + 1;
          const uint64_t seed_hit = (reference_id << 32) | candidate_position;

          if (use_heap) {
            mm_negative_hits[mi].push_back(seed_hit);
          } else {
            negative_hits.push_back(seed_hit);
          }
        }
      }
    }

    if (num_occurrences >= (uint32_t)repetitive_seed_frequency) {
      if (previous_repetitive_seed_position > read_position) {
        // First minimizer.
        repetitive_seed_length += kmer_size_;
      } else {
        if (read_position <
            previous_repetitive_seed_position + kmer_size_ + window_size_ - 1) {
          repetitive_seed_length +=
              read_position - previous_repetitive_seed_position;
        } else {
          repetitive_seed_length += kmer_size_;
        }
      }
      previous_repetitive_seed_position = read_position;
      ++repetitive_seed_count;
    }
  }

  if (use_heap) {
    // TODO: try to remove this sorting.
    if (heap_resort) {
      for (uint32_t mi = 0; mi < num_minimizers; ++mi) {
        std::sort(mm_positive_hits[mi].begin(), mm_positive_hits[mi].end());
      }
    }
    HeapMergeSeedHitLists(mm_positive_hits, positive_hits);
    HeapMergeSeedHitLists(mm_negative_hits, negative_hits);
  } else {
    std::sort(positive_hits.begin(), positive_hits.end());
    std::sort(negative_hits.begin(), negative_hits.end());
  }

#ifdef LI_DEBUG
  for (uint32_t mi = 0; mi < positive_hits->size(); ++mi)
    printf("+ %llu %d %d\n", positive_hits->at(mi),
           (int)(positive_hits->at(mi) >> 32), (int)(positive_hits->at(mi)));

  for (uint32_t mi = 0; mi < negative_hits->size(); ++mi)
    printf("- %llu %d %d\n", negative_hits->at(mi),
           (int)(negative_hits->at(mi) >> 32), (int)(negative_hits->at(mi)));
#endif

  return repetitive_seed_count;
}

int Index::CollectSeedHitsFromRepetitiveReadWithMateInfo(
    int error_threshold, const std::vector<Minimizer> &minimizers,
    uint32_t &repetitive_seed_length, std::vector<uint64_t> &hits,
    const std::vector<Candidate> &mate_candidates, const Strand strand,
    uint32_t search_range, int min_num_seeds_required_for_mapping,
    int max_seed_frequency0) const {
  const uint32_t mate_candidates_size = mate_candidates.size();

  int max_minimizer_count = 0;
  int best_candidate_num = 0;

  for (uint32_t i = 0; i < mate_candidates_size; ++i) {
    int count = mate_candidates[i].count;
    if (count > max_minimizer_count) {
      max_minimizer_count = count;
      best_candidate_num = 1;
    } else if (count == max_minimizer_count) {
      ++best_candidate_num;
    }
  }

  const bool mate_has_too_many_candidates =
      best_candidate_num >= 300 ||
      mate_candidates_size > (uint32_t)max_seed_frequency0;

  const bool mate_has_too_many_low_support_candidates =
      max_minimizer_count <= min_num_seeds_required_for_mapping &&
      best_candidate_num >= 200;

  if (mate_has_too_many_candidates ||
      mate_has_too_many_low_support_candidates) {
    return -max_minimizer_count;
  }

  std::vector<std::pair<uint64_t, uint64_t>> boundaries;
  boundaries.reserve(300);

  for (uint32_t ci = 0; ci < mate_candidates_size; ++ci) {
    if (mate_candidates[ci].count == max_minimizer_count) {
      const uint64_t boundary_start =
          (mate_candidates[ci].position < search_range)
              ? 0
              : (mate_candidates[ci].position - search_range);
      const uint64_t boundary_end = mate_candidates[ci].position + search_range;
      boundaries.emplace_back(boundary_start, boundary_end);
    }
  }

  // Merge adjacent boundary point. Assume the candidates are sorted by
  // coordinate, and thus boundaries are also sorted.
  uint32_t raw_boundary_size = boundaries.size();
  if (raw_boundary_size == 0) {
    return max_minimizer_count;
  }

  uint32_t boundary_size = 1;
  for (uint32_t bi = 1; bi < raw_boundary_size; ++bi) {
    if (boundaries[boundary_size - 1].second < boundaries[bi].first) {
      boundaries[boundary_size] = boundaries[bi];
      ++boundary_size;
    } else {
      boundaries[boundary_size - 1].second = boundaries[bi].second;
    }
  }
  boundaries.resize(boundary_size);

  uint32_t previous_repetitive_seed_position =
      std::numeric_limits<uint32_t>::max();

  repetitive_seed_length = 0;

  for (uint32_t mi = 0; mi < minimizers.size(); ++mi) {
    khiter_t khash_iterator =
        kh_get(k64, lookup_table_, minimizers[mi].GetHashKey() << 1);
    if (khash_iterator == kh_end(lookup_table_)) {
      // std::cerr << "The minimizer is not in reference!\n";
      continue;
    }

    const uint64_t value = kh_value(lookup_table_, khash_iterator);
    const uint32_t read_position = minimizers[mi].GetSequencePosition();
    const Strand read_strand = minimizers[mi].GetSequenceStrand();

    const bool is_reference_minimizer_single =
        (kh_key(lookup_table_, khash_iterator) & 1) > 0;

    if (is_reference_minimizer_single) {
      const uint64_t reference_id = value >> 33;
      const uint32_t reference_position = value >> 1;
      const Strand reference_strand = (value & 1) == 0 ? kPositive : kNegative;

      if (read_strand == reference_strand) {
        if (strand == kPositive) {
          const uint32_t candidate_position =
              reference_position - read_position;
          const uint64_t seed_hit = (reference_id << 32) | candidate_position;
          hits.push_back(seed_hit);
        }
      } else if (strand == kNegative) {
        const uint32_t candidate_position =
            reference_position + read_position - kmer_size_ + 1;
        const uint64_t seed_hit = (reference_id << 32) | candidate_position;
        hits.push_back(seed_hit);
      }

      continue;
    }

    const uint32_t offset = value >> 32;
    const uint32_t num_occurrences = value;
    int32_t prev_l = 0;
    for (uint32_t bi = 0; bi < boundary_size; ++bi) {
      // use binary search to locate the coordinate near mate position
      int32_t l = prev_l, m = 0, r = num_occurrences - 1;
      uint64_t boundary = boundaries[bi].first;
      while (l <= r) {
        m = (l + r) / 2;

        uint64_t value = (occurrence_table_[offset + m]) >> 1;

        if (value < boundary) {
          l = m + 1;
        } else if (value > boundary) {
          r = m - 1;
        } else {
          break;
        }
      }

      prev_l = m;

      for (uint32_t oi = m; oi < num_occurrences; ++oi) {
        const uint64_t value = occurrence_table_[offset + oi];
        if ((value >> 1) > boundaries[bi].second) {
          break;
        }

        const uint64_t reference_id = value >> 33;
        const uint32_t reference_position = value >> 1;
        const Strand reference_strand =
            (value & 1) == 0 ? kPositive : kNegative;

        if (read_strand == reference_strand) {
          if (strand == kPositive) {
            const uint32_t candidate_position =
                reference_position - read_position;
            const uint64_t seed_hit = (reference_id << 32) | candidate_position;
            hits.push_back(seed_hit);
          }
        } else if (strand == kNegative) {
          const uint32_t candidate_position =
              reference_position + read_position - kmer_size_ + 1;
          const uint64_t seed_hit = (reference_id << 32) | candidate_position;
          hits.push_back(seed_hit);
        }
      }
    }  // for bi

    if (num_occurrences >= (uint32_t)max_seed_frequency0) {
      if (previous_repetitive_seed_position > read_position) {
        // First minimizer.
        repetitive_seed_length += kmer_size_;
      } else {
        if (read_position <
            previous_repetitive_seed_position + kmer_size_ + window_size_ - 1) {
          repetitive_seed_length +=
              read_position - previous_repetitive_seed_position;
        } else {
          repetitive_seed_length += kmer_size_;
        }
      }
      previous_repetitive_seed_position = read_position;
    }
  }  // for mi

  std::sort(hits.begin(), hits.end());

#ifdef LI_DEBUG
  for (uint32_t i = 0; i < hits->size(); ++i)
    printf("%s: %d %d\n", __func__, (int)(hits->at(i) >> 32), (int)hits->at(i));
  std::cerr << "Rescue gen on one dir\n ";
  printf("%s: %d\n", __func__, hits->size());
#endif

  return max_minimizer_count;
}

}  // namespace chromap
