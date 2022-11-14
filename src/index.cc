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
  const size_t num_minimizers = minimizers.size();
  assert(num_minimizers > 0);
  // TODO: check this assert!
  // Here I make sure the # minimizers is less than the limit of signed int32,
  // so that I can use int to store position later.
  assert(num_minimizers <= static_cast<size_t>(INT_MAX));

  occurrence_table_.reserve(num_minimizers);
  uint64_t previous_lookup_hash =
      GenerateHashInLookupTable(minimizers[0].GetHash());
  uint32_t num_previous_minimizer_occurrences = 0;
  uint64_t num_nonsingletons = 0;
  uint32_t num_singletons = 0;
  for (size_t mi = 0; mi <= num_minimizers; ++mi) {
    const bool is_last_iteration = mi == num_minimizers;
    const uint64_t current_lookup_hash =
        is_last_iteration ? previous_lookup_hash + 1
                          : GenerateHashInLookupTable(minimizers[mi].GetHash());

    if (current_lookup_hash != previous_lookup_hash) {
      int khash_return_code = 0;
      khiter_t khash_iterator =
          kh_put(k64, lookup_table_, previous_lookup_hash, &khash_return_code);
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
            GenerateEntryValueInLookupTable(num_nonsingletons,
                                            num_previous_minimizer_occurrences);
        num_nonsingletons += num_previous_minimizer_occurrences;
      }
      num_previous_minimizer_occurrences = 1;
    } else {
      num_previous_minimizer_occurrences++;
    }

    if (is_last_iteration) {
      break;
    }

    occurrence_table_.push_back(minimizers[mi].GetHit());
    previous_lookup_hash = current_lookup_hash;
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
    khiter_t khash_iterator = kh_get(
        k64, lookup_table_, GenerateHashInLookupTable(minimizers[i].GetHash()));
    assert(khash_iterator != kh_end(lookup_table_));
    uint64_t key = kh_key(lookup_table_, khash_iterator);
    uint64_t value = kh_value(lookup_table_, khash_iterator);
    if (key & 1) {  // singleton
      assert(minimizers[i].GetHit() == value);
      count = 0;
    } else {
      uint32_t offset = GenerateOffsetInOccurrenceTable(value);
      uint32_t num_occ = GenerateNumOccurrenceInOccurrenceTable(value);
      uint64_t value_in_index = occurrence_table_[offset + count];
      assert(value_in_index == minimizers[i].GetHit());
      ++count;
      if (count == num_occ) {
        count = 0;
      }
    }
  }
}

int Index::GenerateCandidatePositions(
    const CandidatePositionGeneratingConfig &generating_config,
    MappingMetadata &mapping_metadata) const {
  const uint32_t num_minimizers = mapping_metadata.GetNumMinimizers();
  const std::vector<Minimizer> &minimizers = mapping_metadata.minimizers_;

  std::vector<std::vector<uint64_t>> positive_candidate_position_lists;
  std::vector<std::vector<uint64_t>> negative_candidate_position_lists;
  if (generating_config.UseHeapMerge()) {
    for (uint32_t i = 0; i < num_minimizers; ++i) {
      positive_candidate_position_lists.emplace_back(std::vector<uint64_t>());
      negative_candidate_position_lists.emplace_back(std::vector<uint64_t>());
    }
  }
  bool is_candidate_position_list_sorted = true;

  mapping_metadata.positive_hits_.reserve(
      generating_config.GetMaxSeedFrequency() * 2);
  mapping_metadata.negative_hits_.reserve(
      generating_config.GetMaxSeedFrequency() * 2);

  RepetitiveSeedStats repetitive_seed_stats;
  for (uint32_t mi = 0; mi < num_minimizers; ++mi) {
    khiter_t khash_iterator =
        kh_get(k64, lookup_table_,
               GenerateHashInLookupTable(minimizers[mi].GetHash()));
    if (khash_iterator == kh_end(lookup_table_)) {
      // std::cerr << "The minimizer is not in reference!\n";
      continue;
    }

    std::vector<uint64_t> &positive_candidate_positions =
        generating_config.UseHeapMerge() ? positive_candidate_position_lists[mi]
                                         : mapping_metadata.positive_hits_;
    std::vector<uint64_t> &negative_candidate_positions =
        generating_config.UseHeapMerge() ? negative_candidate_position_lists[mi]
                                         : mapping_metadata.negative_hits_;

    const uint64_t lookup_key = kh_key(lookup_table_, khash_iterator);
    const uint64_t lookup_value = kh_value(lookup_table_, khash_iterator);
    const uint64_t read_hit = minimizers[mi].GetHit();
    if (IsSingletonLookupKey(lookup_key)) {
      const uint64_t candidate_position = GenerateCandidatePositionFromHits(
          /*reference_hit=*/lookup_value, read_hit);
      if (AreTwoHitsOnTheSameStrand(/*reference_hit=*/lookup_value, read_hit)) {
        positive_candidate_positions.push_back(candidate_position);
      } else {
        negative_candidate_positions.push_back(candidate_position);
      }
      continue;
    }

    const uint32_t num_occurrences =
        GenerateNumOccurrenceInOccurrenceTable(lookup_value);
    if (!generating_config.IsFrequentSeed(num_occurrences)) {
      const uint32_t read_position = HitToSequencePosition(read_hit);
      const uint32_t occ_offset = GenerateOffsetInOccurrenceTable(lookup_value);
      for (uint32_t oi = 0; oi < num_occurrences; ++oi) {
        const uint64_t reference_hit = occurrence_table_[occ_offset + oi];
        const uint64_t candidate_position =
            GenerateCandidatePositionFromHits(reference_hit, read_hit);
        if (AreTwoHitsOnTheSameStrand(reference_hit, read_hit)) {
          const uint32_t reference_position =
              HitToSequencePosition(reference_hit);
          if (reference_position < read_position) {
            is_candidate_position_list_sorted = false;
          }
          positive_candidate_positions.push_back(candidate_position);
        } else {
          negative_candidate_positions.push_back(candidate_position);
        }
      }
    }

    if (generating_config.IsRepetitiveSeed(num_occurrences)) {
      const uint32_t read_position = HitToSequencePosition(read_hit);
      UpdateRepetitiveSeedStats(read_position, repetitive_seed_stats);
    }
  }

  if (generating_config.UseHeapMerge()) {
    // TODO: try to remove this sorting.
    if (!is_candidate_position_list_sorted) {
      for (uint32_t mi = 0; mi < num_minimizers; ++mi) {
        std::sort(positive_candidate_position_lists[mi].begin(),
                  positive_candidate_position_lists[mi].end());
      }
    }
    HeapMergeCandidatePositionLists(positive_candidate_position_lists,
                                    mapping_metadata.positive_hits_);
    HeapMergeCandidatePositionLists(negative_candidate_position_lists,
                                    mapping_metadata.negative_hits_);
  } else {
    std::sort(mapping_metadata.positive_hits_.begin(),
              mapping_metadata.positive_hits_.end());
    std::sort(mapping_metadata.negative_hits_.begin(),
              mapping_metadata.negative_hits_.end());
  }

#ifdef LI_DEBUG
  for (uint32_t mi = 0; mi < positive_hits.size(); ++mi)
    printf("+ %llu %d %d\n", positive_hits[mi],
           (int)(positive_hits[mi] >> 32), (int)(positive_hits[mi]));

  for (uint32_t mi = 0; mi < negative_hits.size(); ++mi)
    printf("- %llu %d %d\n", negative_hits[mi],
           (int)(negative_hits[mi] >> 32), (int)(negative_hits[mi]));
#endif

  mapping_metadata.repetitive_seed_length_ =
      repetitive_seed_stats.repetitive_seed_length;
  return repetitive_seed_stats.repetitive_seed_count;
}

int Index::GenerateCandidatePositionsFromRepetitiveReadWithMateInfoOnOneStrand(
    const Strand strand, uint32_t search_range,
    int min_num_seeds_required_for_mapping, int max_seed_frequency0,
    int error_threshold, const std::vector<Minimizer> &minimizers,
    const std::vector<Candidate> &mate_candidates,
    uint32_t &repetitive_seed_length,
    std::vector<uint64_t> &candidate_positions) const {
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
      mate_candidates_size > static_cast<uint32_t>(max_seed_frequency0);
  const bool mate_has_too_many_low_support_candidates =
      max_minimizer_count <= min_num_seeds_required_for_mapping &&
      best_candidate_num >= 200;
  if (mate_has_too_many_candidates ||
      mate_has_too_many_low_support_candidates) {
    return -max_minimizer_count;
  }

  // TODO: reduce the search range based on the strand.
  std::vector<std::pair<uint64_t, uint64_t>> boundaries;
  boundaries.reserve(best_candidate_num);
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

  const uint32_t raw_boundary_size = boundaries.size();
  if (raw_boundary_size == 0) {
    return max_minimizer_count;
  }

  // Merge adjacent boundary point. Assume the candidates are sorted by
  // coordinate, and thus boundaries are also sorted.
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

  RepetitiveSeedStats repetitive_seed_stats;
  for (uint32_t mi = 0; mi < minimizers.size(); ++mi) {
    khiter_t khash_iterator =
        kh_get(k64, lookup_table_,
               GenerateHashInLookupTable(minimizers[mi].GetHash()));
    if (khash_iterator == kh_end(lookup_table_)) {
      // std::cerr << "The minimizer is not in reference!\n";
      continue;
    }

    const uint64_t lookup_key = kh_key(lookup_table_, khash_iterator);
    const uint64_t lookup_value = kh_value(lookup_table_, khash_iterator);
    const uint64_t read_hit = minimizers[mi].GetHit();
    const uint32_t read_position = HitToSequencePosition(read_hit);
    if (IsSingletonLookupKey(lookup_key)) {
      const uint64_t candidate_position =
          GenerateCandidatePositionFromHits(lookup_value, read_hit);
      const bool on_same_strand =
          AreTwoHitsOnTheSameStrand(lookup_value, read_hit);
      if ((on_same_strand && strand == kPositive) ||
          (!on_same_strand && strand == kNegative)) {
        candidate_positions.push_back(candidate_position);
      }
      continue;
    }

    const uint32_t offset = GenerateOffsetInOccurrenceTable(lookup_value);
    const uint32_t num_occurrences =
        GenerateNumOccurrenceInOccurrenceTable(lookup_value);
    int32_t prev_l = 0;
    for (uint32_t bi = 0; bi < boundary_size; ++bi) {
      // Use binary search to locate the coordinate near mate position.
      int32_t l = prev_l, m = 0, r = num_occurrences - 1;
      uint64_t boundary = boundaries[bi].first;
      while (l <= r) {
        m = (l + r) / 2;
        uint64_t candidate_position =
            GenerateCandidatePositionFromOccurrenceTableEntry(
                occurrence_table_[offset + m]);
        if (candidate_position < boundary) {
          l = m + 1;
        } else if (candidate_position > boundary) {
          r = m - 1;
        } else {
          break;
        }
      }
      // For next boundary, we don't have to start from l=0.
      prev_l = m;

      for (uint32_t oi = m; oi < num_occurrences; ++oi) {
        const uint64_t reference_hit = occurrence_table_[offset + oi];
        if ((GenerateCandidatePositionFromOccurrenceTableEntry(reference_hit)) >
            boundaries[bi].second) {
          break;
        }
        const uint64_t candidate_position =
            GenerateCandidatePositionFromHits(reference_hit, read_hit);
        const bool on_same_strand =
            AreTwoHitsOnTheSameStrand(reference_hit, read_hit);
        if ((on_same_strand && strand == kPositive) ||
            (!on_same_strand && strand == kNegative)) {
          candidate_positions.push_back(candidate_position);
        }
      }
    }

    if (num_occurrences >= (uint32_t)max_seed_frequency0) {
      UpdateRepetitiveSeedStats(read_position, repetitive_seed_stats);
    }
  }

  std::sort(candidate_positions.begin(), candidate_positions.end());
  repetitive_seed_length = repetitive_seed_stats.repetitive_seed_length;
  return max_minimizer_count;
}

uint64_t Index::GenerateCandidatePositionFromHits(uint64_t reference_hit,
                                                  uint64_t read_hit) const {
  const uint32_t reference_position = HitToSequencePosition(reference_hit);
  const uint32_t read_position = HitToSequencePosition(read_hit);
  // For now we can't see the reference here. So let us don't validate this
  // candidate position. Instead, we do it later some time when we check the
  // candidates.
  const uint32_t reference_start_position =
      AreTwoHitsOnTheSameStrand(reference_hit, read_hit)
          ? reference_position - read_position
          : reference_position + read_position - kmer_size_ + 1;
  const uint64_t reference_id = HitToSequenceIndex(reference_hit);
  return SequenceIndexAndPositionToCandidatePosition(reference_id,
                                                     reference_start_position);
}

void Index::UpdateRepetitiveSeedStats(uint32_t read_position,
                                      RepetitiveSeedStats &stats) const {
  if (stats.previous_repetitive_seed_position > read_position) {
    // First minimizer.
    stats.repetitive_seed_length += kmer_size_;
  } else {
    if (read_position < stats.previous_repetitive_seed_position + kmer_size_ +
                            window_size_ - 1) {
      stats.repetitive_seed_length +=
          read_position - stats.previous_repetitive_seed_position;
    } else {
      stats.repetitive_seed_length += kmer_size_;
    }
  }
  stats.previous_repetitive_seed_position = read_position;
  ++stats.repetitive_seed_count;
}

}  // namespace chromap
