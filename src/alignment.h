#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include "sam_mapping.h"
#include "sequence_batch.h"
#include "utils.h"

namespace chromap {

int GetLongestMatchLength(const char *pattern, const char *text,
                          const int read_length);

// return: newly adjust reference start/end position (kPositive for start,
// kNegative for end)
int AdjustGapBeginning(Direction mapping_direction, const char *ref,
                       const char *read, int *gap_beginning, int read_end,
                       int ref_start_position, int ref_end_position,
                       int *n_cigar, uint32_t **cigar);

void GenerateMDTag(const char *pattern, const char *text,
                   int mapping_start_position, int n_cigar,
                   const uint32_t *cigar, int &NM, std::string &MD_tag);

int BandedAlignPatternToText(int error_threshold, const char *pattern,
                             const char *text, const int read_length,
                             int *mapping_end_position);

// Return negative number if the termination are deemed at the beginning of the
// read mappping_end_position is relative to pattern (reference)
// read_mapping_length is for text (read)
int BandedAlignPatternToTextWithDropOff(int error_threshold,
                                        const char *pattern, const char *text,
                                        const int read_length,
                                        int *mapping_end_position,
                                        int *read_mapping_length);

int BandedAlignPatternToTextWithDropOffFrom3End(
    int error_threshold, const char *pattern, const char *text,
    const int read_length, int *mapping_end_position, int *read_mapping_length);

void BandedAlign4PatternsToText(int error_threshold, const char **patterns,
                                const char *text, int read_length,
                                int32_t *mapping_edit_distances,
                                int32_t *mapping_end_positions);

void BandedAlign8PatternsToText(int error_threshold, const char **patterns,
                                const char *text, int read_length,
                                int16_t *mapping_edit_distances,
                                int16_t *mapping_end_positions);

void BandedTraceback(int error_threshold, int min_num_errors,
                     const char *pattern, const char *text,
                     const int read_length, int *mapping_start_position);

void BandedTracebackToEnd(int error_threshold, int min_num_errors,
                          const char *pattern, const char *text,
                          const int read_length, int *mapping_end_position);

}  // namespace chromap

#endif  // ALIGNMENT_H_
