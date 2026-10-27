#include "minimizer_generator.h"

#include <smmintrin.h>

#include <algorithm>

#include "utils.h"

namespace chromap {

namespace {

// GenerateMinimizers walks a sequence in blocks of this many positions: it
// builds the k-mers of the block, hashes them all in one loop and then runs
// the window logic over the stored hashes.
constexpr uint32_t kSeedBlockSize = 256;

// Hash64 on each 64-bit lane.
inline __m128i Hash64x2(__m128i key, const __m128i mask) {
  const __m128i ones = _mm_set1_epi64x(-1);
  key = _mm_and_si128(
      _mm_add_epi64(_mm_xor_si128(key, ones), _mm_slli_epi64(key, 21)), mask);
  key = _mm_xor_si128(key, _mm_srli_epi64(key, 24));
  key = _mm_and_si128(_mm_add_epi64(_mm_add_epi64(key, _mm_slli_epi64(key, 3)),
                                    _mm_slli_epi64(key, 8)),
                      mask);
  key = _mm_xor_si128(key, _mm_srli_epi64(key, 14));
  key = _mm_and_si128(_mm_add_epi64(_mm_add_epi64(key, _mm_slli_epi64(key, 2)),
                                    _mm_slli_epi64(key, 4)),
                      mask);
  key = _mm_xor_si128(key, _mm_srli_epi64(key, 28));
  key = _mm_and_si128(_mm_add_epi64(key, _mm_slli_epi64(key, 31)), mask);
  return key;
}

// Hash64 of one forward and one reverse k-mer, as in HashSeeds: returns the
// Hash64 of the smaller of the two hashes and sets *strand to 0 when the
// forward hash is the smaller one and to 1 otherwise.
inline uint64_t HashSeed(uint64_t forward_seed, uint64_t reverse_seed,
                         uint64_t mask, uint64_t *strand) {
  const uint64_t forward_hash = Hash64(forward_seed, mask);
  const uint64_t reverse_hash = Hash64(reverse_seed, mask);
  const uint64_t forward_is_smaller = (forward_hash - reverse_hash) >> 63;
  const uint64_t smaller_hash =
      reverse_hash ^ ((forward_hash ^ reverse_hash) & (0 - forward_is_smaller));
  *strand = forward_is_smaller ^ 1;
  return Hash64(smaller_hash, mask);
}

// HashSeed's choice on two lanes: smaller_hash takes the forward hash and
// strand is 0 in the lanes where the forward hash is the smaller; elsewhere
// they are the reverse hash and 1. blendv_pd selects by each lane's sign bit.
inline void ChooseStrands(__m128i forward_hash, __m128i reverse_hash,
                          __m128i *smaller_hash, __m128i *strand) {
  const __m128i difference = _mm_sub_epi64(forward_hash, reverse_hash);
  *smaller_hash = _mm_castpd_si128(
      _mm_blendv_pd(_mm_castsi128_pd(reverse_hash),
                    _mm_castsi128_pd(forward_hash),
                    _mm_castsi128_pd(difference)));
  *strand = _mm_srli_epi64(
      _mm_xor_si128(difference, _mm_set1_epi64x(-1)), 63);
}

// For every position i in [start, end), both multiples of 4: strands[i] is 0
// when Hash64(forward_seeds[i]) < Hash64(reverse_seeds[i]) and 1 otherwise,
// and hash_keys[i] is Hash64 of the hash of that strand. Hashes are at most
// mask < 2^56, so the sign bit of their difference is the comparison.
//
// Each step of eight positions hashes six in SSE registers, two lanes at a
// time, and two in general-purpose registers, so the vector and the integer
// units work at the same time; a last group of four is split two and two.
void HashSeeds(const uint64_t *forward_seeds, const uint64_t *reverse_seeds,
               uint32_t start, uint32_t end, uint64_t mask,
               uint64_t *hash_keys, uint64_t *strands) {
  const __m128i vmask = _mm_set1_epi64x((long long)mask);
  uint32_t i = start;
  for (; i + 8 <= end; i += 8) {
    __m128i forward_hashes[3];
    __m128i reverse_hashes[3];
    for (int j = 0; j < 3; ++j) {
      forward_hashes[j] = Hash64x2(
          _mm_load_si128((const __m128i *)(forward_seeds + i + 2 * j)), vmask);
      reverse_hashes[j] = Hash64x2(
          _mm_load_si128((const __m128i *)(reverse_seeds + i + 2 * j)), vmask);
    }
    uint64_t strand6, strand7;
    const uint64_t hash6 =
        HashSeed(forward_seeds[i + 6], reverse_seeds[i + 6], mask, &strand6);
    const uint64_t hash7 =
        HashSeed(forward_seeds[i + 7], reverse_seeds[i + 7], mask, &strand7);
    for (int j = 0; j < 3; ++j) {
      __m128i smaller_hash, strand;
      ChooseStrands(forward_hashes[j], reverse_hashes[j], &smaller_hash,
                    &strand);
      _mm_store_si128((__m128i *)(hash_keys + i + 2 * j),
                      Hash64x2(smaller_hash, vmask));
      _mm_store_si128((__m128i *)(strands + i + 2 * j), strand);
    }
    hash_keys[i + 6] = hash6;
    hash_keys[i + 7] = hash7;
    strands[i + 6] = strand6;
    strands[i + 7] = strand7;
  }
  if (i < end) {
    const __m128i forward_hash =
        Hash64x2(_mm_load_si128((const __m128i *)(forward_seeds + i)), vmask);
    const __m128i reverse_hash =
        Hash64x2(_mm_load_si128((const __m128i *)(reverse_seeds + i)), vmask);
    uint64_t strand2, strand3;
    const uint64_t hash2 =
        HashSeed(forward_seeds[i + 2], reverse_seeds[i + 2], mask, &strand2);
    const uint64_t hash3 =
        HashSeed(forward_seeds[i + 3], reverse_seeds[i + 3], mask, &strand3);
    __m128i smaller_hash, strand;
    ChooseStrands(forward_hash, reverse_hash, &smaller_hash, &strand);
    _mm_store_si128((__m128i *)(hash_keys + i), Hash64x2(smaller_hash, vmask));
    _mm_store_si128((__m128i *)(strands + i), strand);
    hash_keys[i + 2] = hash2;
    hash_keys[i + 3] = hash3;
    strands[i + 2] = strand2;
    strands[i + 3] = strand3;
  }
}

}  // namespace

void MinimizerGenerator::GenerateMinimizers(
    const SequenceBatch &sequence_batch, uint32_t sequence_index,
    std::vector<Minimizer> &minimizers) const {
  const uint32_t sequence_length =
      sequence_batch.GetSequenceLengthAt(sequence_index);
  const char *sequence = sequence_batch.GetSequenceAt(sequence_index);

  const uint64_t num_shifted_bits = 2 * (kmer_size_ - 1);
  const uint64_t mask = (((uint64_t)1) << (2 * kmer_size_)) - 1;

  uint64_t seeds_in_two_strands[2] = {0, 0};

  // The bases, k-mers, hashes and strands of the positions in one block.
  uint8_t bases[kSeedBlockSize];
  alignas(32) uint64_t forward_seeds[kSeedBlockSize];
  alignas(32) uint64_t reverse_seeds[kSeedBlockSize];
  alignas(32) uint64_t hash_keys[kSeedBlockSize];
  alignas(32) uint64_t strands[kSeedBlockSize];

  // The seeds in the window, UINT64_MAX when there is none.
  uint64_t buffer_hashes[256];
  uint64_t buffer_hits[256];
  std::fill(buffer_hashes, buffer_hashes + window_size_, UINT64_MAX);
  std::fill(buffer_hits, buffer_hits + window_size_, UINT64_MAX);
  uint64_t min_hash = UINT64_MAX;
  uint64_t min_hit = UINT64_MAX;

  int unambiguous_length = 0;
  int position_in_buffer = 0;
  int min_position = 0;

  for (uint32_t block_start = 0; block_start < sequence_length;
       block_start += kSeedBlockSize) {
    const uint32_t block_length =
        std::min(kSeedBlockSize, sequence_length - block_start);

    for (uint32_t i = 0; i < block_length; ++i) {
      const uint8_t current_base = CharToUint8(sequence[block_start + i]);
      bases[i] = current_base;
      if (current_base < 4) {
        // Not an ambiguous base.
        // Forward k-mer.
        seeds_in_two_strands[0] =
            ((seeds_in_two_strands[0] << 2) | current_base) & mask;
        // Reverse k-mer.
        seeds_in_two_strands[1] =
            (seeds_in_two_strands[1] >> 2) |
            (((uint64_t)(3 ^ current_base)) << num_shifted_bits);
      }
      forward_seeds[i] = seeds_in_two_strands[0];
      reverse_seeds[i] = seeds_in_two_strands[1];
    }

    // A position before k - 1 never ends k unambiguous bases, so its hash is
    // not used.
    const uint32_t hash_start =
        block_start == 0 ? ((uint32_t)(kmer_size_ - 1) & ~3u) : 0;
    const uint32_t hash_end = (block_length + 3) & ~3u;
    for (uint32_t i = block_length; i < hash_end; ++i) {
      forward_seeds[i] = 0;
      reverse_seeds[i] = 0;
    }
    if (hash_start < hash_end) {
      HashSeeds(forward_seeds, reverse_seeds, hash_start, hash_end, mask,
                hash_keys, strands);
    }

    for (uint32_t i = 0; i < block_length; ++i) {
      const uint32_t position = block_start + i;
      uint64_t current_hash = UINT64_MAX;
      uint64_t current_hit = UINT64_MAX;

      if (bases[i] < 4) {
        if (forward_seeds[i] == reverse_seeds[i]) {
          // Skip "symmetric k-mers" as we don't know it strand.
          continue;
        }

        ++unambiguous_length;

        if (unambiguous_length >= kmer_size_) {
          current_hash = hash_keys[i];
          current_hit =
              ((((uint64_t)sequence_index) << 32 | (uint32_t)position) << 1) |
              strands[i];
        }
      } else {
        unambiguous_length = 0;
      }

      // Need to do this here as appropriate position_in_buffer and
      // buf[position_in_buffer] are needed below.
      buffer_hashes[position_in_buffer] = current_hash;
      buffer_hits[position_in_buffer] = current_hit;
      if (unambiguous_length == window_size_ + kmer_size_ - 1 &&
          min_hash != UINT64_MAX && min_hash < current_hash) {
        // Special case for the first window - because identical k-mers are
        // not stored yet.
        for (int j = position_in_buffer + 1; j < window_size_; ++j)
          if (min_hash == buffer_hashes[j] && buffer_hits[j] != min_hit)
            minimizers.emplace_back(buffer_hashes[j], buffer_hits[j]);
        for (int j = 0; j < position_in_buffer; ++j)
          if (min_hash == buffer_hashes[j] && buffer_hits[j] != min_hit)
            minimizers.emplace_back(buffer_hashes[j], buffer_hits[j]);
      }

      if (current_hash <= min_hash) {
        // A new minimum; then write the old min.
        if (unambiguous_length >= window_size_ + kmer_size_ &&
            min_hash != UINT64_MAX) {
          minimizers.emplace_back(min_hash, min_hit);
        }
        min_hash = current_hash;
        min_hit = current_hit;
        min_position = position_in_buffer;
      } else if (position_in_buffer == min_position) {
        // Old min has moved outside the window.
        if (unambiguous_length >= window_size_ + kmer_size_ - 1 &&
            min_hash != UINT64_MAX) {
          minimizers.emplace_back(min_hash, min_hit);
        }

        min_hash = UINT64_MAX;
        for (int j = position_in_buffer + 1; j < window_size_; ++j) {
          // The two loops are necessary when there are identical k-mers.
          if (min_hash >= buffer_hashes[j]) {
            // >= is important s.t. min is always the closest k-mer.
            min_hash = buffer_hashes[j];
            min_hit = buffer_hits[j];
            min_position = j;
          }
        }

        for (int j = 0; j <= position_in_buffer; ++j) {
          if (min_hash >= buffer_hashes[j]) {
            min_hash = buffer_hashes[j];
            min_hit = buffer_hits[j];
            min_position = j;
          }
        }

        if (unambiguous_length >= window_size_ + kmer_size_ - 1 &&
            min_hash != UINT64_MAX) {
          // Write identical k-mers.
          // These two loops make sure the output is sorted.
          for (int j = position_in_buffer + 1; j < window_size_; ++j)
            if (min_hash == buffer_hashes[j] && min_hit != buffer_hits[j])
              minimizers.emplace_back(buffer_hashes[j], buffer_hits[j]);
          for (int j = 0; j <= position_in_buffer; ++j)
            if (min_hash == buffer_hashes[j] && min_hit != buffer_hits[j])
              minimizers.emplace_back(buffer_hashes[j], buffer_hits[j]);
        }
      }

      ++position_in_buffer;
      if (position_in_buffer == window_size_) {
        position_in_buffer = 0;
      }
    }
  }

  if (min_hash != UINT64_MAX) {
    minimizers.emplace_back(min_hash, min_hit);
  }
}

}  // namespace chromap
