#include "minimizer_generator.h"

#include "utils.h"

namespace chromap {

void MinimizerGenerator::GenerateMinimizers(
    const SequenceBatch &sequence_batch, uint32_t sequence_index,
    std::vector<Minimizer> &minimizers) const {
  const uint32_t sequence_length =
      sequence_batch.GetSequenceLengthAt(sequence_index);
  const char *sequence = sequence_batch.GetSequenceAt(sequence_index);

  const uint64_t num_shifted_bits = 2 * (kmer_size_ - 1);
  const uint64_t mask = (((uint64_t)1) << (2 * kmer_size_)) - 1;

  uint64_t seeds_in_two_strands[2] = {0, 0};
  std::pair<uint64_t, uint64_t> buffer[256];
  std::pair<uint64_t, uint64_t> min_seed = {UINT64_MAX, UINT64_MAX};

  // 2 uint64_t cost 16 bytes.
  memset(buffer, 0xff, window_size_ * 16);

  int unambiguous_length = 0;
  int position_in_buffer = 0;
  int min_position = 0;

  for (uint32_t position = 0; position < sequence_length; ++position) {
    const uint8_t current_base = CharToUint8(sequence[position]);
    std::pair<uint64_t, uint64_t> current_seed = {UINT64_MAX, UINT64_MAX};

    if (current_base < 4) {
      // Not an ambiguous base.
      // Forward k-mer.
      seeds_in_two_strands[0] =
          ((seeds_in_two_strands[0] << 2) | current_base) & mask;
      // Reverse k-mer.
      seeds_in_two_strands[1] =
          (seeds_in_two_strands[1] >> 2) |
          (((uint64_t)(3 ^ current_base)) << num_shifted_bits);

      if (seeds_in_two_strands[0] == seeds_in_two_strands[1]) {
        // Skip "symmetric k-mers" as we don't know it strand.
        continue;
      }

      uint64_t hash_keys_for_two_seeds[2] = {
          Hash64(seeds_in_two_strands[0], mask),
          Hash64(seeds_in_two_strands[1], mask)};

      uint64_t strand =
          hash_keys_for_two_seeds[0] < hash_keys_for_two_seeds[1] ? 0 : 1;

      ++unambiguous_length;

      if (unambiguous_length >= kmer_size_) {
        current_seed.first = Hash64(hash_keys_for_two_seeds[strand], mask);
        current_seed.second =
            ((((uint64_t)sequence_index) << 32 | (uint32_t)position) << 1) |
            strand;
      }
    } else {
      unambiguous_length = 0;
    }

    // Need to do this here as appropriate position_in_buffer and
    // buf[position_in_buffer] are needed below.
    buffer[position_in_buffer] = current_seed;
    if (unambiguous_length == window_size_ + kmer_size_ - 1 &&
        min_seed.first != UINT64_MAX && min_seed.first < current_seed.first) {
      // Special case for the first window - because identical k-mers are not
      // stored yet.
      for (int j = position_in_buffer + 1; j < window_size_; ++j)
        if (min_seed.first == buffer[j].first &&
            buffer[j].second != min_seed.second)
          minimizers.emplace_back(buffer[j]);
      for (int j = 0; j < position_in_buffer; ++j)
        if (min_seed.first == buffer[j].first &&
            buffer[j].second != min_seed.second)
          minimizers.emplace_back(buffer[j]);
    }

    if (current_seed.first <= min_seed.first) {
      // A new minimum; then write the old min.
      if (unambiguous_length >= window_size_ + kmer_size_ &&
          min_seed.first != UINT64_MAX) {
        minimizers.emplace_back(min_seed);
      }
      min_seed = current_seed;
      min_position = position_in_buffer;
    } else if (position_in_buffer == min_position) {
      // Old min has moved outside the window.
      if (unambiguous_length >= window_size_ + kmer_size_ - 1 &&
          min_seed.first != UINT64_MAX) {
        minimizers.emplace_back(min_seed);
      }

      min_seed.first = UINT64_MAX;
      for (int j = position_in_buffer + 1; j < window_size_; ++j) {
        // The two loops are necessary when there are identical k-mers.
        if (min_seed.first >= buffer[j].first) {
          // >= is important s.t. min is always the closest k-mer.
          min_seed = buffer[j];
          min_position = j;
        }
      }

      for (int j = 0; j <= position_in_buffer; ++j) {
        if (min_seed.first >= buffer[j].first) {
          min_seed = buffer[j];
          min_position = j;
        }
      }

      if (unambiguous_length >= window_size_ + kmer_size_ - 1 &&
          min_seed.first != UINT64_MAX) {
        // Write identical k-mers.
        // These two loops make sure the output is sorted.
        for (int j = position_in_buffer + 1; j < window_size_; ++j)
          if (min_seed.first == buffer[j].first &&
              min_seed.second != buffer[j].second)
            minimizers.emplace_back(buffer[j]);
        for (int j = 0; j <= position_in_buffer; ++j)
          if (min_seed.first == buffer[j].first &&
              min_seed.second != buffer[j].second)
            minimizers.emplace_back(buffer[j]);
      }
    }

    ++position_in_buffer;
    if (position_in_buffer == window_size_) {
      position_in_buffer = 0;
    }
  }

  if (min_seed.first != UINT64_MAX) {
    minimizers.emplace_back(min_seed);
  }
}

}  // namespace chromap
