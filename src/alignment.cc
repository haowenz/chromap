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

int BandedAlignPatternToText(int error_threshold, const char *pattern,
                             const char *text, const int read_length,
                             int *mapping_end_position) {
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
    if (num_errors_at_band_start_position > 3 * error_threshold) {
      return error_threshold + 1;
    }
    for (int ai = 0; ai < 5; ai++) {
      Peq[ai] >>= 1;
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

// Return negative number if the termination are deemed at the beginning of the
// read mappping_end_position is relative to pattern (reference)
// read_mapping_length is for text (read)
int BandedAlignPatternToTextWithDropOff(int error_threshold,
                                        const char *pattern, const char *text,
                                        const int read_length,
                                        int *mapping_end_position,
                                        int *read_mapping_length) {
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
  uint32_t prev_VP = 0;
  uint32_t prev_VN = 0;
  int num_errors_at_band_start_position = 0;
  int i = 0;
  int fail_beginning = 0;  // the alignment failed at the beginning part
  int prev_num_errors_at_band_start_position = 0;
  for (; i < read_length; i++) {
    uint8_t pattern_base = CharToUint8(pattern[i + 2 * error_threshold]);
    Peq[pattern_base] = Peq[pattern_base] | highest_bit_in_band_mask;
    X = Peq[CharToUint8(text[i])] | VN;
    D0 = ((VP + (X & VP)) ^ VP) | X;
    HN = VP & D0;
    HP = VN | ~(VP | D0);
    X = D0 >> 1;
    prev_VN = VN;
    prev_VP = VP;
    VN = X & HP;
    VP = HN | ~(X | HP);
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
    for (int ai = 0; ai < 5; ai++) {
      Peq[ai] >>= 1;
    }
  }

  /*char tmp[255] ;
  strncpy(tmp, pattern, read_length + 2 * error_threshold);
  printf("%s\n%s\n", tmp, text);
  printf("%d\n", i) ;
  fflush(stdout);*/
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

int BandedAlignPatternToTextWithDropOffFrom3End(int error_threshold,
                                                const char *pattern,
                                                const char *text,
                                                const int read_length,
                                                int *mapping_end_position,
                                                int *read_mapping_length) {
  uint32_t Peq[5] = {0, 0, 0, 0, 0};
  for (int i = 0; i < 2 * error_threshold; i++) {
    uint8_t base =
        CharToUint8(pattern[read_length + 2 * error_threshold - 1 - i]);
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
  uint32_t prev_VP = 0;
  uint32_t prev_VN = 0;
  int num_errors_at_band_start_position = 0;
  int i = 0;
  int fail_beginning = 0;  // the alignment failed at the beginning part
  int prev_num_errors_at_band_start_position = 0;
  for (; i < read_length; i++) {
    // printf("%c %c %d\n", pattern[read_length - 1 - i], pattern[read_length -
    // 1 - i + error_threshold], text[read_length - 1 - i]);
    uint8_t pattern_base = CharToUint8(pattern[read_length - 1 - i]);
    Peq[pattern_base] = Peq[pattern_base] | highest_bit_in_band_mask;
    X = Peq[CharToUint8(text[read_length - 1 - i])] | VN;
    D0 = ((VP + (X & VP)) ^ VP) | X;
    HN = VP & D0;
    HP = VN | ~(VP | D0);
    X = D0 >> 1;
    prev_VN = VN;
    prev_VP = VP;
    VN = X & HP;
    VP = HN | ~(X | HP);
    prev_num_errors_at_band_start_position = num_errors_at_band_start_position;
    num_errors_at_band_start_position += 1 - (D0 & lowest_bit_in_band_mask);
    /*printf("->%d %d %c %c", i, num_errors_at_band_start_position,
    pattern[read_length - 1 - i], text[read_length - 1 - i]) ; int tmp =
    num_errors_at_band_start_position; for (int j = 0; j < 2 * error_threshold;
    j++) { tmp = tmp + ((VP >> j) & (uint32_t) 1); tmp = tmp - ((VN >> j) &
    (uint32_t) 1); printf(" %d", tmp);
    }
    printf("\n");*/
    if (num_errors_at_band_start_position > 2 * error_threshold) {
      // return error_threshold + 1;
      if (i < 4 * error_threshold && i < read_length / 2) {
        fail_beginning = 1;
      }
      break;
    }
    for (int ai = 0; ai < 5; ai++) {
      Peq[ai] >>= 1;
    }
  }
  // printf("li %d: %d %d %d\n", fail_beginning, i, error_threshold,
  // read_length);
  if (i < read_length) {
    num_errors_at_band_start_position = prev_num_errors_at_band_start_position;
    VN = prev_VN;
    VP = prev_VP;
  }
  int band_start_position = i - 1;
  int min_num_errors = num_errors_at_band_start_position;
  *read_mapping_length = i;
  *mapping_end_position = band_start_position;
  // printf("-1: %d\n", num_errors_at_band_start_position);
  for (i = 0; i < 2 * error_threshold; i++) {
    num_errors_at_band_start_position =
        num_errors_at_band_start_position + ((VP >> i) & (uint32_t)1);
    num_errors_at_band_start_position =
        num_errors_at_band_start_position - ((VN >> i) & (uint32_t)1);
    // printf("%d: %d\n", i, num_errors_at_band_start_position);
    if (num_errors_at_band_start_position < min_num_errors ||
        (num_errors_at_band_start_position == min_num_errors &&
         i + 1 == error_threshold)) {
      min_num_errors = num_errors_at_band_start_position;
      *mapping_end_position = band_start_position + (1 + i);
    }
  }
  if (fail_beginning ||
      (read_length > 60 &&
       *mapping_end_position + 1 - error_threshold - min_num_errors < 30)) {
    *mapping_end_position = -*mapping_end_position;
  }
  return min_num_errors;
}

void BandedAlign4PatternsToText(int error_threshold, const char **patterns,
                                const char *text, int read_length,
                                int32_t *mapping_edit_distances,
                                int32_t *mapping_end_positions) {
  int ALPHABET_SIZE = 5;
  const char *reference_sequence0 = patterns[0];
  const char *reference_sequence1 = patterns[1];
  const char *reference_sequence2 = patterns[2];
  const char *reference_sequence3 = patterns[3];
  uint32_t highest_bit_in_band_mask = 1 << (2 * error_threshold);
  __m128i highest_bit_in_band_mask_vpu0 =
      _mm_set_epi32(0, 0, 0, highest_bit_in_band_mask);
  __m128i highest_bit_in_band_mask_vpu1 =
      _mm_set_epi32(0, 0, highest_bit_in_band_mask, 0);
  __m128i highest_bit_in_band_mask_vpu2 =
      _mm_set_epi32(0, highest_bit_in_band_mask, 0, 0);
  __m128i highest_bit_in_band_mask_vpu3 =
      _mm_set_epi32(highest_bit_in_band_mask, 0, 0, 0);
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
    Peq[base0] = _mm_or_si128(highest_bit_in_band_mask_vpu0, Peq[base0]);
    Peq[base1] = _mm_or_si128(highest_bit_in_band_mask_vpu1, Peq[base1]);
    Peq[base2] = _mm_or_si128(highest_bit_in_band_mask_vpu2, Peq[base2]);
    Peq[base3] = _mm_or_si128(highest_bit_in_band_mask_vpu3, Peq[base3]);
    for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
      Peq[ai] = _mm_srli_epi32(Peq[ai], 1);
    }
  }

  uint32_t lowest_bit_in_band_mask = 1;
  __m128i lowest_bit_in_band_mask_vpu = _mm_set1_epi32(lowest_bit_in_band_mask);
  __m128i VP = _mm_setzero_si128();
  __m128i VN = _mm_setzero_si128();
  __m128i X = _mm_setzero_si128();
  __m128i D0 = _mm_setzero_si128();
  __m128i HN = _mm_setzero_si128();
  __m128i HP = _mm_setzero_si128();
  __m128i max_mask_vpu = _mm_set1_epi32(0xffffffff);
  __m128i num_errors_at_band_start_position_vpu = _mm_setzero_si128();
  __m128i early_stop_threshold_vpu = _mm_set1_epi32(error_threshold * 3);
  for (int i = 0; i < read_length; i++) {
    uint8_t base0 = CharToUint8(reference_sequence0[i + 2 * error_threshold]);
    uint8_t base1 = CharToUint8(reference_sequence1[i + 2 * error_threshold]);
    uint8_t base2 = CharToUint8(reference_sequence2[i + 2 * error_threshold]);
    uint8_t base3 = CharToUint8(reference_sequence3[i + 2 * error_threshold]);
    Peq[base0] = _mm_or_si128(highest_bit_in_band_mask_vpu0, Peq[base0]);
    Peq[base1] = _mm_or_si128(highest_bit_in_band_mask_vpu1, Peq[base1]);
    Peq[base2] = _mm_or_si128(highest_bit_in_band_mask_vpu2, Peq[base2]);
    Peq[base3] = _mm_or_si128(highest_bit_in_band_mask_vpu3, Peq[base3]);
    X = _mm_or_si128(Peq[CharToUint8(text[i])], VN);
    D0 = _mm_and_si128(X, VP);
    D0 = _mm_add_epi32(D0, VP);
    D0 = _mm_xor_si128(D0, VP);
    D0 = _mm_or_si128(D0, X);
    HN = _mm_and_si128(VP, D0);
    HP = _mm_or_si128(VP, D0);
    HP = _mm_xor_si128(HP, max_mask_vpu);
    HP = _mm_or_si128(HP, VN);
    X = _mm_srli_epi32(D0, 1);
    VN = _mm_and_si128(X, HP);
    VP = _mm_or_si128(X, HP);
    VP = _mm_xor_si128(VP, max_mask_vpu);
    VP = _mm_or_si128(VP, HN);
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
    for (int ai = 0; ai < ALPHABET_SIZE; ai++) {
      Peq[ai] = _mm_srli_epi32(Peq[ai], 1);
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
