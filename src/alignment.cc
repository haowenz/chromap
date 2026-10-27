#include "alignment.h"

#include <smmintrin.h>

namespace chromap {

int GetLongestMatchLength(const char *pattern, const char *text,
                          const int read_length) {
  int max_match = 0;
  int tmp = 0;
  for (int i = 0; i < read_length; ++i) {
    if (CharToUint8(pattern[i]) == CharToUint8(text[i])) {
      ++tmp;
    } else if (tmp > max_match) {
      max_match = tmp;
    }
  }
  if (tmp > max_match) {
    max_match = tmp;
  }
  return max_match;
}

int AdjustGapBeginning(const Strand mapping_strand, const char *ref,
                       const char *read, int *gap_beginning, int read_end,
                       int ref_start_position, int ref_end_position,
                       int *n_cigar, uint32_t **cigar) {
  int i, j;
  if (mapping_strand == kPositive) {
    if (*gap_beginning <= 0) {
      return ref_start_position;
    }

    // printf("%d\n", *gap_beginning);

    for (i = *gap_beginning - 1, j = ref_start_position - 1; i >= 0 && j >= 0;
         --i, --j) {
      // printf("%c %c\n", read[i], ref[j]);
      if (read[i] != ref[j] && read[i] != ref[j] - 'a' + 'A') {
        break;
      }
    }

    *gap_beginning = i + 1;
    // TODO: add soft clip in cigar
    if (n_cigar && *n_cigar > 0) {
      if (((*cigar)[0] & 0xf) == BAM_CMATCH) {
        (*cigar)[0] += (ref_start_position - 1 - j) << 4;
      }
    }

    return j + 1;
  }

  if (*gap_beginning <= 0) {
    return ref_end_position;
  }

  // printf("%d\n", *gap_beginning);
  /*char *tmp = new char[255] ;
  strncpy(tmp, ref + ref_start_position, ref_end_position - ref_start_position
  + 1 + 10) ; printf("%s %d. %d %d\n", tmp, strlen(tmp), ref_end_position -
  ref_start_position + 1 + 10, strlen(ref)) ; delete[] tmp;*/

  for (i = read_end + 1, j = ref_end_position + 1; read[i] && ref[j];
       ++i, ++j) {
    // printf("%c %c %c %c %c %c\n", read[i], ref[j - 1], ref[j], ref[j + 1],
    // ref[j + 2], ref[j + 3]);
    if (read[i] != ref[j] && read[i] != ref[j] - 'a' + 'A') {
      break;
    }
  }

  *gap_beginning = *gap_beginning + i - (read_end + 1);

  if (n_cigar && *n_cigar > 0) {
    if (((*cigar)[*n_cigar - 1] & 0xf) == BAM_CMATCH) {
      (*cigar)[*n_cigar - 1] += (j - (ref_end_position + 1)) << 4;
    }
  }

  return j - 1;
}

void GenerateNMAndMDTag(const char *pattern, const char *text,
                        int mapping_start_position,
                        MappingInMemory &mapping_in_memory) {
  const char *read = text;
  const char *reference = pattern + mapping_start_position;

  const uint32_t *cigar = mapping_in_memory.cigar;
  const int n_cigar = mapping_in_memory.n_cigar;
  mapping_in_memory.NM = 0;
  mapping_in_memory.MD_tag.clear();

  int num_matches = 0;
  int read_position = 0;
  int reference_position = 0;

  for (int ci = 0; ci < n_cigar; ++ci) {
    uint32_t current_cigar_uint = cigar[ci];
    uint8_t cigar_operation = bam_cigar_op(current_cigar_uint);
    int num_cigar_operations = bam_cigar_oplen(current_cigar_uint);
    if (cigar_operation == BAM_CMATCH) {
      for (int opi = 0; opi < num_cigar_operations; ++opi) {
        if (reference[reference_position] == read[read_position] ||
            reference[reference_position] - 'a' + 'A' == read[read_position]) {
          // a match
          ++num_matches;
        } else {
          // a mismatch
          ++mapping_in_memory.NM;
          
          mapping_in_memory.MD_tag.append(std::to_string(num_matches));
          num_matches = 0;
          mapping_in_memory.MD_tag.push_back(reference[reference_position]);
        }
        ++reference_position;
        ++read_position;
      }
    } else if (cigar_operation == BAM_CINS) {
      mapping_in_memory.NM += num_cigar_operations;
      read_position += num_cigar_operations;
    } else if (cigar_operation == BAM_CDEL) {
      mapping_in_memory.NM += num_cigar_operations;
      
      mapping_in_memory.MD_tag.append(std::to_string(num_matches));
      num_matches = 0;
      mapping_in_memory.MD_tag.push_back('^');
      for (int opi = 0; opi < num_cigar_operations; ++opi) {
        mapping_in_memory.MD_tag.push_back(reference[reference_position]);
        ++reference_position;
      }
    } else {
      std::cerr << "Unexpected cigar op: " << (int)cigar_operation << "\n";
    }
  }
  mapping_in_memory.MD_tag.append(std::to_string(num_matches));
}

namespace {

// The banded kernels below read Peq[c], the mask whose bit b is set when
// base c is at pattern position i + b, b in [0, 2 * error_threshold], at
// iteration i (mirrored for the 3' end). Two sources give it: the original
// five words updated every iteration, and bitmasks of the whole pattern window
// built once, used when the band fits in 32 bits and the window in the masks.
constexpr int kPatternMaskWords = 8;  // windows of up to 512 bases
constexpr int kPatternMaskStride = kPatternMaskWords + 1;

inline bool UsePatternMasks(int error_threshold, int read_length) {
  if (error_threshold < 0 || error_threshold > 15 || read_length <= 0) {
    return false;
  }
  const int window_length = read_length + 2 * error_threshold;
  return window_length >= 16 && window_length <= 64 * kPatternMaskWords;
}

// Peq in five words, shifted every iteration. The window is
// pattern[0, read_length + 2 * error_threshold), read from the 3' end when
// kFrom3End is true.
template <bool kFrom3End>
class PeqArray {
 public:
  PeqArray(int error_threshold, const char *pattern, int read_length)
      : pattern_(pattern),
        read_length_(read_length),
        error_threshold_(error_threshold),
        highest_bit_in_band_mask_(1 << (2 * error_threshold)) {
    for (int i = 0; i < 2 * error_threshold; i++) {
      uint8_t base =
          kFrom3End
              ? CharToUint8(pattern[read_length + 2 * error_threshold - 1 - i])
              : CharToUint8(pattern[i]);
      peq_[base] = peq_[base] | (1 << i);
    }
  }

  uint32_t Get(int i, uint8_t text_base) {
    uint8_t pattern_base =
        kFrom3End ? CharToUint8(pattern_[read_length_ - 1 - i])
                  : CharToUint8(pattern_[i + 2 * error_threshold_]);
    peq_[pattern_base] = peq_[pattern_base] | highest_bit_in_band_mask_;
    const uint32_t peq = peq_[text_base];
    for (int ai = 0; ai < 5; ai++) {
      peq_[ai] >>= 1;
    }
    return peq;
  }

 private:
  const char *pattern_;
  const int read_length_;
  const int error_threshold_;
  const uint32_t highest_bit_in_band_mask_;
  uint32_t peq_[5] = {0, 0, 0, 0, 0};
};

// Peq as bits [i, i + 2 * error_threshold] of one bitmask per base over the
// window: bit j of base c's mask is set when CharToUint8(window[j]) == c (or
// window[window_length - 1 - j] when kFrom3End is true). Built with one SSE
// compare per base per 16 bytes; every load lies inside the window. Needs
// UsePatternMasks.
template <bool kFrom3End>
class PeqFromMasks {
 public:
  PeqFromMasks(int error_threshold, const char *pattern, int read_length)
      : band_mask_(
            (uint32_t)(((uint64_t)2 << (2 * error_threshold)) - 1)) {
    const int window_length = read_length + 2 * error_threshold;
    const __m128i lower_case_bit = _mm_set1_epi8(0x20);
    const __m128i a = _mm_set1_epi8('a');
    const __m128i c = _mm_set1_epi8('c');
    const __m128i g = _mm_set1_epi8('g');
    const __m128i t = _mm_set1_epi8('t');
    const __m128i reversed_lanes =
        _mm_setr_epi8(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    const int num_words = (window_length + 63) >> 6;
    for (int w = 0; w < num_words; ++w) {
      uint64_t mask_a = 0, mask_c = 0, mask_g = 0, mask_t = 0, mask_n = 0;
      for (int k = 0; k < 4; ++k) {
        const int start = (w << 6) + (k << 4);
        const int remaining = window_length - start;
        if (remaining <= 0) {
          break;
        }
        // A last chunk of fewer than 16 bases is read as the window's
        // outermost 16 bytes, and its first (16 - remaining) lanes dropped.
        const char *chunk;
        int dropped_lanes = 0;
        if (remaining >= 16) {
          chunk = kFrom3End ? pattern + window_length - start - 16
                            : pattern + start;
        } else {
          chunk = kFrom3End ? pattern : pattern + window_length - 16;
          dropped_lanes = 16 - remaining;
        }
        __m128i bytes = _mm_loadu_si128((const __m128i *)chunk);
        if (kFrom3End) {
          bytes = _mm_shuffle_epi8(bytes, reversed_lanes);
        }
        // Only 'A' and 'a' become 'a' with bit 5 set, and so on.
        bytes = _mm_or_si128(bytes, lower_case_bit);
        const uint32_t lanes_a =
            (uint32_t)_mm_movemask_epi8(_mm_cmpeq_epi8(bytes, a));
        const uint32_t lanes_c =
            (uint32_t)_mm_movemask_epi8(_mm_cmpeq_epi8(bytes, c));
        const uint32_t lanes_g =
            (uint32_t)_mm_movemask_epi8(_mm_cmpeq_epi8(bytes, g));
        const uint32_t lanes_t =
            (uint32_t)_mm_movemask_epi8(_mm_cmpeq_epi8(bytes, t));
        const uint32_t lanes_n =
            ~(lanes_a | lanes_c | lanes_g | lanes_t) & 0xffff;
        const int shift = k << 4;
        mask_a |= (uint64_t)(lanes_a >> dropped_lanes) << shift;
        mask_c |= (uint64_t)(lanes_c >> dropped_lanes) << shift;
        mask_g |= (uint64_t)(lanes_g >> dropped_lanes) << shift;
        mask_t |= (uint64_t)(lanes_t >> dropped_lanes) << shift;
        mask_n |= (uint64_t)(lanes_n >> dropped_lanes) << shift;
      }
      masks_[0 * kPatternMaskStride + w] = mask_a;
      masks_[1 * kPatternMaskStride + w] = mask_c;
      masks_[2 * kPatternMaskStride + w] = mask_g;
      masks_[3 * kPatternMaskStride + w] = mask_t;
      masks_[4 * kPatternMaskStride + w] = mask_n;
    }
    // Get reads the word after the one holding bit i.
    for (int base = 0; base < 5; ++base) {
      masks_[base * kPatternMaskStride + num_words] = 0;
    }
  }

  uint32_t Get(int i, uint8_t text_base) const {
    const uint64_t *word =
        masks_ + text_base * kPatternMaskStride + (i >> 6);
    const int shift = i & 63;
    const uint64_t bits = (word[0] >> shift) | ((word[1] << 1) << (63 - shift));
    return (uint32_t)bits & band_mask_;
  }

 private:
  const uint32_t band_mask_;
  uint64_t masks_[5 * kPatternMaskStride];
};

template <typename PeqSource>
int AlignPatternToText(int error_threshold, const char *text,
                       const int read_length, int *mapping_end_position,
                       PeqSource &peq) {
  uint32_t lowest_bit_in_band_mask = 1;
  uint32_t VP = 0;
  uint32_t VN = 0;
  uint32_t X = 0;
  uint32_t D0 = 0;
  uint32_t HN = 0;
  int num_errors_at_band_start_position = 0;
  for (int i = 0; i < read_length; i++) {
    // The recurrence D0 = ((VP + (X & VP)) ^ VP) | X, HP = VN | ~(VP | D0),
    // VN' = X' & HP and VP' = HN | ~(X' | HP), with X = Peq | VN and
    // X' = D0 >> 1, in two dependent operations fewer per base: VP and VN
    // never share a bit, so X & VP = Peq & VP, and not_HP = ~HP =
    // (VP | D0) & ~VN.
    const uint32_t peq_bits = peq.Get(i, CharToUint8(text[i]));
    const uint32_t not_VN = ~VN;
    X = peq_bits | VN;
    D0 = ((VP + (peq_bits & VP)) ^ VP) | X;
    HN = VP & D0;
    const uint32_t not_HP = (VP | D0) & not_VN;
    X = D0 >> 1;
    VN = X & ~not_HP;
    VP = HN | (~X & not_HP);
    num_errors_at_band_start_position += 1 - (D0 & lowest_bit_in_band_mask);
    if (num_errors_at_band_start_position > 3 * error_threshold) {
      return error_threshold + 1;
    }
  }
  int band_start_position = read_length - 1;
  int min_num_errors = num_errors_at_band_start_position;
  *mapping_end_position = band_start_position;
  for (int i = 0; i < 2 * error_threshold; i++) {
    num_errors_at_band_start_position =
        num_errors_at_band_start_position + ((VP >> i) & (uint32_t)1);
    num_errors_at_band_start_position =
        num_errors_at_band_start_position - ((VN >> i) & (uint32_t)1);
    if (num_errors_at_band_start_position < min_num_errors ||
        (num_errors_at_band_start_position == min_num_errors &&
         i + 1 == error_threshold)) {
      min_num_errors = num_errors_at_band_start_position;
      *mapping_end_position = band_start_position + 1 + i;
    }
  }
  return min_num_errors;
}

// The drop-off kernels; kFrom3End aligns from the 3' end of the text.
template <bool kFrom3End, typename PeqSource>
int AlignPatternToTextWithDropOff(int error_threshold, const char *text,
                                  const int read_length,
                                  int *mapping_end_position,
                                  int *read_mapping_length, PeqSource &peq) {
  uint32_t lowest_bit_in_band_mask = 1;
  uint32_t VP = 0;
  uint32_t VN = 0;
  uint32_t X = 0;
  uint32_t D0 = 0;
  uint32_t HN = 0;
  uint32_t prev_VP = 0;
  uint32_t prev_VN = 0;
  int num_errors_at_band_start_position = 0;
  int i = 0;
  int fail_beginning = 0;  // the alignment failed at the beginning part
  int prev_num_errors_at_band_start_position = 0;
  for (; i < read_length; i++) {
    const char text_char = kFrom3End ? text[read_length - 1 - i] : text[i];
    // The recurrence of AlignPatternToText, written as there.
    const uint32_t peq_bits = peq.Get(i, CharToUint8(text_char));
    const uint32_t not_VN = ~VN;
    X = peq_bits | VN;
    D0 = ((VP + (peq_bits & VP)) ^ VP) | X;
    HN = VP & D0;
    const uint32_t not_HP = (VP | D0) & not_VN;
    X = D0 >> 1;
    prev_VN = VN;
    prev_VP = VP;
    VN = X & ~not_HP;
    VP = HN | (~X & not_HP);
    prev_num_errors_at_band_start_position = num_errors_at_band_start_position;
    num_errors_at_band_start_position += 1 - (D0 & lowest_bit_in_band_mask);
    if (num_errors_at_band_start_position > 2 * error_threshold) {
      // return error_threshold + 1;
      // the min error in this band could be still less than the
      // error_threshold, and could but this should be fine since it does not
      // affect the 5' end of the read.
      if (i < 4 * error_threshold && i < read_length / 2) {
        fail_beginning = 1;
      }
      break;
    }
  }

  if (i < read_length) {
    num_errors_at_band_start_position = prev_num_errors_at_band_start_position;
    VN = prev_VN;
    VP = prev_VP;
  }
  int band_start_position = i - 1;
  int min_num_errors = num_errors_at_band_start_position;
  *read_mapping_length = i;
  *mapping_end_position = band_start_position;

  for (i = 0; i < 2 * error_threshold; i++) {
    num_errors_at_band_start_position =
        num_errors_at_band_start_position + ((VP >> i) & (uint32_t)1);
    num_errors_at_band_start_position =
        num_errors_at_band_start_position - ((VN >> i) & (uint32_t)1);
    if (num_errors_at_band_start_position < min_num_errors ||
        (num_errors_at_band_start_position == min_num_errors &&
         i + 1 == error_threshold)) {
      min_num_errors = num_errors_at_band_start_position;
      *mapping_end_position = band_start_position + 1 + i;
    }
  }
  if (fail_beginning ||
      (read_length > 60 &&
       *mapping_end_position + 1 - error_threshold - min_num_errors < 30)) {
    *mapping_end_position = -*mapping_end_position;
  }
  return min_num_errors;
}

}  // namespace

int BandedAlignPatternToText(int error_threshold, const char *pattern,
                             const char *text, const int read_length,
                             int *mapping_end_position) {
  if (UsePatternMasks(error_threshold, read_length)) {
    PeqFromMasks<false> peq(error_threshold, pattern, read_length);
    return AlignPatternToText(error_threshold, text, read_length,
                              mapping_end_position, peq);
  }
  PeqArray<false> peq(error_threshold, pattern, read_length);
  return AlignPatternToText(error_threshold, text, read_length,
                            mapping_end_position, peq);
}

// Return negative number if the termination are deemed at the beginning of the
// read mappping_end_position is relative to pattern (reference)
// read_mapping_length is for text (read)
int BandedAlignPatternToTextWithDropOff(int error_threshold,
                                        const char *pattern, const char *text,
                                        const int read_length,
                                        int *mapping_end_position,
                                        int *read_mapping_length) {
  if (UsePatternMasks(error_threshold, read_length)) {
    PeqFromMasks<false> peq(error_threshold, pattern, read_length);
    return AlignPatternToTextWithDropOff<false>(
        error_threshold, text, read_length, mapping_end_position,
        read_mapping_length, peq);
  }
  PeqArray<false> peq(error_threshold, pattern, read_length);
  return AlignPatternToTextWithDropOff<false>(error_threshold, text,
                                              read_length, mapping_end_position,
                                              read_mapping_length, peq);
}

int BandedAlignPatternToTextWithDropOffFrom3End(int error_threshold,
                                                const char *pattern,
                                                const char *text,
                                                const int read_length,
                                                int *mapping_end_position,
                                                int *read_mapping_length) {
  if (UsePatternMasks(error_threshold, read_length)) {
    PeqFromMasks<true> peq(error_threshold, pattern, read_length);
    return AlignPatternToTextWithDropOff<true>(
        error_threshold, text, read_length, mapping_end_position,
        read_mapping_length, peq);
  }
  PeqArray<true> peq(error_threshold, pattern, read_length);
  return AlignPatternToTextWithDropOff<true>(error_threshold, text,
                                             read_length, mapping_end_position,
                                             read_mapping_length, peq);
}

namespace {

// BandedAlign4PatternsToText's Peq, one 32-bit lane per pattern: bit b of lane
// l of Get(i, c) is set when base c is at patterns[l][i + b], b in
// [0, 2 * error_threshold].

// The five vectors updated every iteration, indexed by base.
class LanePeqArray {
 public:
  LanePeqArray(int error_threshold, const char **patterns)
      : error_threshold_(error_threshold) {
    const uint32_t highest_bit_in_band_mask = 1 << (2 * error_threshold);
    for (int l = 0; l < 4; ++l) {
      patterns_[l] = patterns[l];
    }
    highest_bit_in_band_mask_vpu_[0] =
        _mm_set_epi32(0, 0, 0, highest_bit_in_band_mask);
    highest_bit_in_band_mask_vpu_[1] =
        _mm_set_epi32(0, 0, highest_bit_in_band_mask, 0);
    highest_bit_in_band_mask_vpu_[2] =
        _mm_set_epi32(0, highest_bit_in_band_mask, 0, 0);
    highest_bit_in_band_mask_vpu_[3] =
        _mm_set_epi32(highest_bit_in_band_mask, 0, 0, 0);
    for (int ai = 0; ai < 5; ai++) {
      peq_[ai] = _mm_setzero_si128();
    }
    for (int i = 0; i < 2 * error_threshold; i++) {
      Add(i);
      Shift();
    }
  }

  __m128i Get(int i, uint8_t text_base) {
    Add(i + 2 * error_threshold_);
    const __m128i peq = peq_[text_base];
    Shift();
    return peq;
  }

 private:
  void Add(int pattern_position) {
    for (int l = 0; l < 4; ++l) {
      const uint8_t base = CharToUint8(patterns_[l][pattern_position]);
      peq_[base] = _mm_or_si128(highest_bit_in_band_mask_vpu_[l], peq_[base]);
    }
  }

  void Shift() {
    for (int ai = 0; ai < 5; ai++) {
      peq_[ai] = _mm_srli_epi32(peq_[ai], 1);
    }
  }

  const char *patterns_[4];
  const int error_threshold_;
  __m128i highest_bit_in_band_mask_vpu_[4];
  __m128i peq_[5];
};

constexpr int kLaneMaskBlocks = 2 * kPatternMaskWords;  // 32-bit blocks

// Peq as bits [i, i + 2 * error_threshold] of one bitmask per base over each
// window, built once as in PeqFromMasks<false>: bit j of lane l of block k of
// base c is set when CharToUint8(patterns[l][32 * k + j]) == c. Needs
// UsePatternMasks.
class LanePeqFromMasks {
 public:
  LanePeqFromMasks(int error_threshold, const char **patterns,
                   int read_length)
      : band_mask_(_mm_set1_epi32((int)(uint32_t)(
            ((uint64_t)2 << (2 * error_threshold)) - 1))) {
    const int window_length = read_length + 2 * error_threshold;
    const int num_chunks = (window_length + 15) >> 4;
    const int num_blocks = (num_chunks + 1) >> 1;
    // Get reads the block after the one holding bit i, and an odd number of
    // chunks leaves the upper half of the last block unwritten.
    for (int base = 0; base < 5; ++base) {
      _mm_store_si128((__m128i *)masks_[base][num_blocks - 1],
                      _mm_setzero_si128());
      _mm_store_si128((__m128i *)masks_[base][num_blocks],
                      _mm_setzero_si128());
    }
    // Chunk by chunk, so the first blocks are complete first.
    for (int q = 0; q < num_chunks; ++q) {
      const int start = q << 4;
      const int remaining = window_length - start;
      if (remaining >= 16) {
        for (int l = 0; l < 4; ++l) {
          AddChunk(patterns[l] + start, 0, q, l);
        }
      } else {
        // A last chunk of fewer than 16 bases is read as the window's last
        // 16 bytes, and its first (16 - remaining) lanes dropped.
        for (int l = 0; l < 4; ++l) {
          AddChunk(patterns[l] + window_length - 16, 16 - remaining, q, l);
        }
      }
    }
  }

  __m128i Get(int i, uint8_t text_base) const {
    const __m128i *block = (const __m128i *)masks_[text_base][i >> 5];
    const int shift = i & 31;
    // A shift by 32 gives zero.
    const __m128i bits = _mm_or_si128(
        _mm_srl_epi32(_mm_load_si128(block), _mm_cvtsi32_si128(shift)),
        _mm_sll_epi32(_mm_load_si128(block + 1),
                      _mm_cvtsi32_si128(32 - shift)));
    return _mm_and_si128(bits, band_mask_);
  }

 private:
  // Bits [16 * q, 16 * q + 16) of lane l's masks from the 16 bytes at chunk,
  // shifted right by dropped_lanes.
  inline void AddChunk(const char *chunk, int dropped_lanes, int q, int l) {
    // Only 'A' and 'a' become 'a' with bit 5 set, and so on.
    const __m128i bytes = _mm_or_si128(
        _mm_loadu_si128((const __m128i *)chunk), _mm_set1_epi8(0x20));
    const uint32_t lanes_a = (uint32_t)_mm_movemask_epi8(
        _mm_cmpeq_epi8(bytes, _mm_set1_epi8('a')));
    const uint32_t lanes_c = (uint32_t)_mm_movemask_epi8(
        _mm_cmpeq_epi8(bytes, _mm_set1_epi8('c')));
    const uint32_t lanes_g = (uint32_t)_mm_movemask_epi8(
        _mm_cmpeq_epi8(bytes, _mm_set1_epi8('g')));
    const uint32_t lanes_t = (uint32_t)_mm_movemask_epi8(
        _mm_cmpeq_epi8(bytes, _mm_set1_epi8('t')));
    const uint32_t lanes_n = ~(lanes_a | lanes_c | lanes_g | lanes_t) & 0xffff;
    const int k = q >> 1;
    const int h = q & 1;
    masks_[0][k][l][h] = (uint16_t)(lanes_a >> dropped_lanes);
    masks_[1][k][l][h] = (uint16_t)(lanes_c >> dropped_lanes);
    masks_[2][k][l][h] = (uint16_t)(lanes_g >> dropped_lanes);
    masks_[3][k][l][h] = (uint16_t)(lanes_t >> dropped_lanes);
    masks_[4][k][l][h] = (uint16_t)(lanes_n >> dropped_lanes);
  }

  const __m128i band_mask_;
  // masks_[c][k][l][h]: bits [32 * k + 16 * h, 32 * k + 16 * h + 16) of lane
  // l's mask for base c, so masks_[c][k] is block k of the four lanes.
  alignas(16) uint16_t masks_[5][kLaneMaskBlocks + 1][4][2];
};

template <typename LanePeqSource>
void Align4PatternsToText(int error_threshold, const char *text,
                          int read_length, int32_t *mapping_edit_distances,
                          int32_t *mapping_end_positions,
                          LanePeqSource &peq) {
  uint32_t lowest_bit_in_band_mask = 1;
  __m128i lowest_bit_in_band_mask_vpu = _mm_set1_epi32(lowest_bit_in_band_mask);
  __m128i VP = _mm_setzero_si128();
  __m128i VN = _mm_setzero_si128();
  __m128i X = _mm_setzero_si128();
  __m128i D0 = _mm_setzero_si128();
  __m128i HN = _mm_setzero_si128();
  __m128i num_errors_at_band_start_position_vpu = _mm_setzero_si128();
  __m128i early_stop_threshold_vpu = _mm_set1_epi32(error_threshold * 3);
  for (int i = 0; i < read_length; i++) {
    X = _mm_or_si128(peq.Get(i, CharToUint8(text[i])), VN);
    D0 = _mm_and_si128(X, VP);
    D0 = _mm_add_epi32(D0, VP);
    D0 = _mm_xor_si128(D0, VP);
    D0 = _mm_or_si128(D0, X);
    HN = _mm_and_si128(VP, D0);
    // With not_HP = ~HP = (VP | D0) & ~VN: VN' = X' & HP and
    // VP' = HN | ~(X' | HP), two dependent operations fewer than computing HP.
    const __m128i not_HP = _mm_andnot_si128(VN, _mm_or_si128(VP, D0));
    X = _mm_srli_epi32(D0, 1);
    VN = _mm_andnot_si128(not_HP, X);
    VP = _mm_or_si128(HN, _mm_andnot_si128(X, not_HP));
    __m128i E = _mm_and_si128(D0, lowest_bit_in_band_mask_vpu);
    E = _mm_xor_si128(E, lowest_bit_in_band_mask_vpu);
    num_errors_at_band_start_position_vpu =
        _mm_add_epi32(num_errors_at_band_start_position_vpu, E);
    __m128i early_stop = _mm_cmpgt_epi32(num_errors_at_band_start_position_vpu,
                                         early_stop_threshold_vpu);
    int tmp = _mm_movemask_epi8(early_stop);
    if (tmp == 0xffff) {
      _mm_store_si128((__m128i *)mapping_edit_distances,
                      num_errors_at_band_start_position_vpu);
      return;
    }
  }
  int band_start_position = read_length - 1;
  __m128i min_num_errors_vpu = num_errors_at_band_start_position_vpu;
  for (int i = 0; i < 2 * error_threshold; i++) {
    __m128i lowest_bit_in_VP_vpu =
        _mm_and_si128(VP, lowest_bit_in_band_mask_vpu);
    __m128i lowest_bit_in_VN_vpu =
        _mm_and_si128(VN, lowest_bit_in_band_mask_vpu);
    num_errors_at_band_start_position_vpu = _mm_add_epi32(
        num_errors_at_band_start_position_vpu, lowest_bit_in_VP_vpu);
    num_errors_at_band_start_position_vpu = _mm_sub_epi32(
        num_errors_at_band_start_position_vpu, lowest_bit_in_VN_vpu);
    __m128i mapping_end_positions_update_mask_vpu = _mm_cmplt_epi32(
        num_errors_at_band_start_position_vpu, min_num_errors_vpu);
    __m128i mapping_end_positions_update_mask_vpu1 = _mm_cmpeq_epi32(
        num_errors_at_band_start_position_vpu, min_num_errors_vpu);
    int mapping_end_positions_update_mask =
        _mm_movemask_epi8(mapping_end_positions_update_mask_vpu);
    int mapping_end_positions_update_mask1 =
        _mm_movemask_epi8(mapping_end_positions_update_mask_vpu1);
    for (int li = 0; li < 4; ++li) {
      if ((mapping_end_positions_update_mask & 1) == 1 ||
          ((mapping_end_positions_update_mask1 & 1) == 1 &&
           i + 1 == error_threshold)) {
        mapping_end_positions[li] = band_start_position + 1 + i;
      }
      mapping_end_positions_update_mask =
          mapping_end_positions_update_mask >> 4;
      mapping_end_positions_update_mask1 =
          mapping_end_positions_update_mask1 >> 4;
    }
    min_num_errors_vpu = _mm_min_epi32(min_num_errors_vpu,
                                       num_errors_at_band_start_position_vpu);
    VP = _mm_srli_epi32(VP, 1);
    VN = _mm_srli_epi32(VN, 1);
  }
  _mm_store_si128((__m128i *)mapping_edit_distances, min_num_errors_vpu);
}

}  // namespace

void BandedAlign4PatternsToText(int error_threshold, const char **patterns,
                                const char *text, int read_length,
                                int32_t *mapping_edit_distances,
                                int32_t *mapping_end_positions) {
  if (UsePatternMasks(error_threshold, read_length)) {
    LanePeqFromMasks peq(error_threshold, patterns, read_length);
    Align4PatternsToText(error_threshold, text, read_length,
                         mapping_edit_distances, mapping_end_positions, peq);
    return;
  }
  LanePeqArray peq(error_threshold, patterns);
  Align4PatternsToText(error_threshold, text, read_length,
                       mapping_edit_distances, mapping_end_positions, peq);
}

void BandedAlign8PatternsToText(int error_threshold, const char **patterns,
                                const char *text, int read_length,
                                int16_t *mapping_edit_distances,
                                int16_t *mapping_end_positions) {
  int ALPHABET_SIZE = 5;
  const char *reference_sequence0 = patterns[0];
  const char *reference_sequence1 = patterns[1];
  const char *reference_sequence2 = patterns[2];
  const char *reference_sequence3 = patterns[3];
  const char *reference_sequence4 = patterns[4];
  const char *reference_sequence5 = patterns[5];
  const char *reference_sequence6 = patterns[6];
  const char *reference_sequence7 = patterns[7];
  uint16_t highest_bit_in_band_mask = 1 << (2 * error_threshold);
  __m128i highest_bit_in_band_mask_vpu0 =
      _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, highest_bit_in_band_mask);
  __m128i highest_bit_in_band_mask_vpu1 =
      _mm_set_epi16(0, 0, 0, 0, 0, 0, highest_bit_in_band_mask, 0);
  __m128i highest_bit_in_band_mask_vpu2 =
      _mm_set_epi16(0, 0, 0, 0, 0, highest_bit_in_band_mask, 0, 0);
  __m128i highest_bit_in_band_mask_vpu3 =
      _mm_set_epi16(0, 0, 0, 0, highest_bit_in_band_mask, 0, 0, 0);
  __m128i highest_bit_in_band_mask_vpu4 =
      _mm_set_epi16(0, 0, 0, highest_bit_in_band_mask, 0, 0, 0, 0);
  __m128i highest_bit_in_band_mask_vpu5 =
      _mm_set_epi16(0, 0, highest_bit_in_band_mask, 0, 0, 0, 0, 0);
  __m128i highest_bit_in_band_mask_vpu6 =
      _mm_set_epi16(0, highest_bit_in_band_mask, 0, 0, 0, 0, 0, 0);
  __m128i highest_bit_in_band_mask_vpu7 =
      _mm_set_epi16(highest_bit_in_band_mask, 0, 0, 0, 0, 0, 0, 0);
  // Init Peq
  __m128i Peq[ALPHABET_SIZE];
  for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
    Peq[ai] = _mm_setzero_si128();
  }
  for (int i = 0; i < 2 * error_threshold; i++) {
    uint8_t base0 = CharToUint8(reference_sequence0[i]);
    uint8_t base1 = CharToUint8(reference_sequence1[i]);
    uint8_t base2 = CharToUint8(reference_sequence2[i]);
    uint8_t base3 = CharToUint8(reference_sequence3[i]);
    uint8_t base4 = CharToUint8(reference_sequence4[i]);
    uint8_t base5 = CharToUint8(reference_sequence5[i]);
    uint8_t base6 = CharToUint8(reference_sequence6[i]);
    uint8_t base7 = CharToUint8(reference_sequence7[i]);
    Peq[base0] = _mm_or_si128(highest_bit_in_band_mask_vpu0, Peq[base0]);
    Peq[base1] = _mm_or_si128(highest_bit_in_band_mask_vpu1, Peq[base1]);
    Peq[base2] = _mm_or_si128(highest_bit_in_band_mask_vpu2, Peq[base2]);
    Peq[base3] = _mm_or_si128(highest_bit_in_band_mask_vpu3, Peq[base3]);
    Peq[base4] = _mm_or_si128(highest_bit_in_band_mask_vpu4, Peq[base4]);
    Peq[base5] = _mm_or_si128(highest_bit_in_band_mask_vpu5, Peq[base5]);
    Peq[base6] = _mm_or_si128(highest_bit_in_band_mask_vpu6, Peq[base6]);
    Peq[base7] = _mm_or_si128(highest_bit_in_band_mask_vpu7, Peq[base7]);
    for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
      Peq[ai] = _mm_srli_epi16(Peq[ai], 1);
    }
  }

  uint16_t lowest_bit_in_band_mask = 1;
  __m128i lowest_bit_in_band_mask_vpu = _mm_set1_epi16(lowest_bit_in_band_mask);
  __m128i VP = _mm_setzero_si128();
  __m128i VN = _mm_setzero_si128();
  __m128i X = _mm_setzero_si128();
  __m128i D0 = _mm_setzero_si128();
  __m128i HN = _mm_setzero_si128();
  __m128i HP = _mm_setzero_si128();
  __m128i max_mask_vpu = _mm_set1_epi16(0xffff);
  __m128i num_errors_at_band_start_position_vpu = _mm_setzero_si128();
  __m128i early_stop_threshold_vpu = _mm_set1_epi16(error_threshold * 3);
  for (int i = 0; i < read_length; i++) {
    uint8_t base0 = CharToUint8(reference_sequence0[i + 2 * error_threshold]);
    uint8_t base1 = CharToUint8(reference_sequence1[i + 2 * error_threshold]);
    uint8_t base2 = CharToUint8(reference_sequence2[i + 2 * error_threshold]);
    uint8_t base3 = CharToUint8(reference_sequence3[i + 2 * error_threshold]);
    uint8_t base4 = CharToUint8(reference_sequence4[i + 2 * error_threshold]);
    uint8_t base5 = CharToUint8(reference_sequence5[i + 2 * error_threshold]);
    uint8_t base6 = CharToUint8(reference_sequence6[i + 2 * error_threshold]);
    uint8_t base7 = CharToUint8(reference_sequence7[i + 2 * error_threshold]);
    Peq[base0] = _mm_or_si128(highest_bit_in_band_mask_vpu0, Peq[base0]);
    Peq[base1] = _mm_or_si128(highest_bit_in_band_mask_vpu1, Peq[base1]);
    Peq[base2] = _mm_or_si128(highest_bit_in_band_mask_vpu2, Peq[base2]);
    Peq[base3] = _mm_or_si128(highest_bit_in_band_mask_vpu3, Peq[base3]);
    Peq[base4] = _mm_or_si128(highest_bit_in_band_mask_vpu4, Peq[base4]);
    Peq[base5] = _mm_or_si128(highest_bit_in_band_mask_vpu5, Peq[base5]);
    Peq[base6] = _mm_or_si128(highest_bit_in_band_mask_vpu6, Peq[base6]);
    Peq[base7] = _mm_or_si128(highest_bit_in_band_mask_vpu7, Peq[base7]);
    X = _mm_or_si128(Peq[CharToUint8(text[i])], VN);
    D0 = _mm_and_si128(X, VP);
    D0 = _mm_add_epi16(D0, VP);
    D0 = _mm_xor_si128(D0, VP);
    D0 = _mm_or_si128(D0, X);
    HN = _mm_and_si128(VP, D0);
    HP = _mm_or_si128(VP, D0);
    HP = _mm_xor_si128(HP, max_mask_vpu);
    HP = _mm_or_si128(HP, VN);
    X = _mm_srli_epi16(D0, 1);
    VN = _mm_and_si128(X, HP);
    VP = _mm_or_si128(X, HP);
    VP = _mm_xor_si128(VP, max_mask_vpu);
    VP = _mm_or_si128(VP, HN);
    __m128i E = _mm_and_si128(D0, lowest_bit_in_band_mask_vpu);
    E = _mm_xor_si128(E, lowest_bit_in_band_mask_vpu);
    num_errors_at_band_start_position_vpu =
        _mm_add_epi16(num_errors_at_band_start_position_vpu, E);
    __m128i early_stop = _mm_cmpgt_epi16(num_errors_at_band_start_position_vpu,
                                         early_stop_threshold_vpu);
    int tmp = _mm_movemask_epi8(early_stop);
    if (tmp == 0xffff) {
      _mm_store_si128((__m128i *)mapping_edit_distances,
                      num_errors_at_band_start_position_vpu);
      return;
    }
    for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
      Peq[ai] = _mm_srli_epi16(Peq[ai], 1);
    }
  }
  int band_start_position = read_length - 1;
  __m128i min_num_errors_vpu = num_errors_at_band_start_position_vpu;
  for (int i = 0; i < 2 * error_threshold; i++) {
    __m128i lowest_bit_in_VP_vpu =
        _mm_and_si128(VP, lowest_bit_in_band_mask_vpu);
    __m128i lowest_bit_in_VN_vpu =
        _mm_and_si128(VN, lowest_bit_in_band_mask_vpu);
    num_errors_at_band_start_position_vpu = _mm_add_epi16(
        num_errors_at_band_start_position_vpu, lowest_bit_in_VP_vpu);
    num_errors_at_band_start_position_vpu = _mm_sub_epi16(
        num_errors_at_band_start_position_vpu, lowest_bit_in_VN_vpu);
    __m128i mapping_end_positions_update_mask_vpu = _mm_cmplt_epi16(
        num_errors_at_band_start_position_vpu, min_num_errors_vpu);
    __m128i mapping_end_positions_update_mask_vpu1 = _mm_cmpeq_epi16(
        num_errors_at_band_start_position_vpu, min_num_errors_vpu);
    int mapping_end_positions_update_mask =
        _mm_movemask_epi8(mapping_end_positions_update_mask_vpu);
    int mapping_end_positions_update_mask1 =
        _mm_movemask_epi8(mapping_end_positions_update_mask_vpu1);
    for (int li = 0; li < 8; ++li) {
      if ((mapping_end_positions_update_mask & 1) == 1 ||
          ((mapping_end_positions_update_mask1 & 1) == 1 &&
           i + 1 == error_threshold)) {
        mapping_end_positions[li] = band_start_position + 1 + i;
      }
      mapping_end_positions_update_mask =
          mapping_end_positions_update_mask >> 2;
      mapping_end_positions_update_mask1 =
          mapping_end_positions_update_mask1 >> 2;
    }
    min_num_errors_vpu = _mm_min_epi16(min_num_errors_vpu,
                                       num_errors_at_band_start_position_vpu);
    VP = _mm_srli_epi16(VP, 1);
    VN = _mm_srli_epi16(VN, 1);
  }
  _mm_store_si128((__m128i *)mapping_edit_distances, min_num_errors_vpu);
}

void BandedTraceback(int error_threshold, int min_num_errors,
                     const char *pattern, const char *text,
                     const int read_length, int *mapping_start_position) {
  // fisrt calculate the hamming distance and see whether it's equal to # errors
  if (min_num_errors == 0) {
    *mapping_start_position = error_threshold;
    return;
  }
  int error_count = 0;
  for (int i = 0; i < read_length; ++i) {
    if (pattern[i + error_threshold] != text[i]) {
      ++error_count;
    }
  }
  if (error_count == min_num_errors) {
    *mapping_start_position = error_threshold;
    return;
  }
  // if not then there are gaps so that we have to traceback with edit distance.
  uint32_t Peq[5] = {0, 0, 0, 0, 0};
  for (int i = 0; i < 2 * error_threshold; i++) {
    uint8_t base =
        CharToUint8(pattern[read_length - 1 + 2 * error_threshold - i]);
    Peq[base] = Peq[base] | (1 << i);
  }
  uint32_t highest_bit_in_band_mask = 1 << (2 * error_threshold);
  uint32_t lowest_bit_in_band_mask = 1;
  uint32_t VP = 0;
  uint32_t VN = 0;
  uint32_t X = 0;
  uint32_t D0 = 0;
  uint32_t HN = 0;
  uint32_t HP = 0;
  int num_errors_at_band_start_position = 0;
  for (int i = 0; i < read_length; i++) {
    uint8_t pattern_base = CharToUint8(pattern[read_length - 1 - i]);
    Peq[pattern_base] = Peq[pattern_base] | highest_bit_in_band_mask;
    X = Peq[CharToUint8(text[read_length - 1 - i])] | VN;
    D0 = ((VP + (X & VP)) ^ VP) | X;
    HN = VP & D0;
    HP = VN | ~(VP | D0);
    X = D0 >> 1;
    VN = X & HP;
    VP = HN | ~(X | HP);
    num_errors_at_band_start_position += 1 - (D0 & lowest_bit_in_band_mask);
    for (int ai = 0; ai < 5; ai++) {
      Peq[ai] >>= 1;
    }
  }
  *mapping_start_position = 2 * error_threshold;
  for (int i = 0; i < 2 * error_threshold; i++) {
    num_errors_at_band_start_position =
        num_errors_at_band_start_position + ((VP >> i) & (uint32_t)1);
    num_errors_at_band_start_position =
        num_errors_at_band_start_position - ((VN >> i) & (uint32_t)1);
    if (num_errors_at_band_start_position == min_num_errors) {
      *mapping_start_position = 2 * error_threshold - (1 + i);
      if (i + 1 == error_threshold) {
        return;
      }
    }
  }
}

void BandedTracebackToEnd(int error_threshold, int min_num_errors,
                          const char *pattern, const char *text,
                          const int read_length, int *mapping_end_position) {
  // fisrt calculate the hamming distance and see whether it's equal to # errors
  if (min_num_errors == 0) {
    *mapping_end_position = read_length + error_threshold;
    return;
  }
  int error_count = 0;
  for (int i = 0; i < read_length; ++i) {
    if (pattern[i + error_threshold] != text[i]) {
      ++error_count;
    }
  }
  if (error_count == min_num_errors) {
    *mapping_end_position = read_length + error_threshold;
    return;
  }
  // if not then there are gaps so that we have to traceback with edit distance.
  uint32_t Peq[5] = {0, 0, 0, 0, 0};
  for (int i = 0; i < 2 * error_threshold; i++) {
    uint8_t base = CharToUint8(pattern[i]);
    Peq[base] = Peq[base] | (1 << i);
  }
  uint32_t highest_bit_in_band_mask = 1 << (2 * error_threshold);
  uint32_t lowest_bit_in_band_mask = 1;
  uint32_t VP = 0;
  uint32_t VN = 0;
  uint32_t X = 0;
  uint32_t D0 = 0;
  uint32_t HN = 0;
  uint32_t HP = 0;
  int num_errors_at_band_start_position = 0;
  for (int i = 0; i < read_length; i++) {
    // printf("=>%d %d %c %c\n", i, num_errors_at_band_start_position, pattern[i
    // + 2 * error_threshold], text[i]) ;
    uint8_t pattern_base = CharToUint8(pattern[i + 2 * error_threshold]);
    Peq[pattern_base] = Peq[pattern_base] | highest_bit_in_band_mask;
    X = Peq[CharToUint8(text[i])] | VN;
    D0 = ((VP + (X & VP)) ^ VP) | X;
    HN = VP & D0;
    HP = VN | ~(VP | D0);
    X = D0 >> 1;
    VN = X & HP;
    VP = HN | ~(X | HP);
    num_errors_at_band_start_position += 1 - (D0 & lowest_bit_in_band_mask);
    for (int ai = 0; ai < 5; ai++) {
      Peq[ai] >>= 1;
    }
  }
  int band_start_position = read_length;
  *mapping_end_position = band_start_position + 1;
  for (int i = 0; i < 2 * error_threshold; i++) {
    num_errors_at_band_start_position =
        num_errors_at_band_start_position + ((VP >> i) & (uint32_t)1);
    num_errors_at_band_start_position =
        num_errors_at_band_start_position - ((VN >> i) & (uint32_t)1);
    if (num_errors_at_band_start_position == min_num_errors) {
      *mapping_end_position = band_start_position + (i + 1);
      if (i + 1 == error_threshold) {
        return;
      }
    }
  }
}

}  // namespace chromap
