#ifndef SEQUENCE_EFFECTIVE_RANGE_H_
#define SEQUENCE_EFFECTIVE_RANGE_H_

#include <stdlib.h>

#include <algorithm>
#include <vector>

#include "utils.h"

namespace chromap {

// The class handles the custom read format indicating the effective range on a
// sequence.
class SequenceEffectiveRange {
 public:
  SequenceEffectiveRange() {}
  ~SequenceEffectiveRange() {}

  void Init() {
    starts.push_back(0);
    ends.push_back(-1);
    strand = '+';
    default_range = true;
  }

  // Return false if it fails to parse the format string.
  bool ParseEffectiveRange(const char *s, int len) {
    int i;
    int j = 0;  // start, end, strand section
    char buffer[20];
    int blen = 0;

    if (default_range) {
      starts.clear();
      ends.clear();
      strand = '+';
      default_range = false;
    }

    for (i = 3; i <= len; ++i) {
      if (i == len || s[i] == ':') {
        buffer[blen] = '\0';
        if (j == 0) {
          starts.push_back(atoi(buffer));
        } else if (j == 1) {
          ends.push_back(atoi(buffer));
        } else {
          strand = buffer[0];
        }

        blen = 0;
        if (i < len && s[i] == ':') ++j;
      } else {
        buffer[blen] = s[i];
        ++blen;
      }
    }

    if (j >= 3 || starts.size() != ends.size()) {
      return false;
    }

    std::sort(starts.begin(), starts.end());
    std::sort(ends.begin(), ends.end());

    const int num_ranges = starts.size();
    if (ends[0] == -1) {
      for (i = 0; i < num_ranges - 1; ++i) {
        ends[i] = ends[i + 1];
      }
      ends[i] = -1;
    }

    return true;
  }

  // Replace by the range specified in the starts, ends section, but does not
  // apply the strand operation. Return new length.
  int Replace(char *s, int len, bool need_complement) {
    if (IsFullRangeAndPositiveStrand()) {
      return len;
    }

    int i, j, k;
    i = 0;
    const int num_ranges = starts.size();
    for (k = 0; k < num_ranges; ++k) {
      int start = starts[k];
      int end = ends[k];

      if (end == -1) {
        end = len - 1;
      }

      for (j = start; j <= end; ++i, ++j) {
        s[i] = s[j];
      }
    }

    s[i] = '\0';
    len = i;

    if (strand == '-') {
      if (need_complement) {
        for (i = 0; i < len; ++i) {
          s[i] = Uint8ToChar(((uint8_t)3) ^ (CharToUint8(s[i])));
        }
      }

      for (i = 0, j = len - 1; i < j; ++i, --j) {
        char tmp = s[i];
        s[i] = s[j];
        s[j] = tmp;
      }
    }
    return len;
  }

 private:
  bool IsFullRangeAndPositiveStrand() {
    if (strand == '+' && starts[0] == 0 && ends[0] == -1) {
      return true;
    }

    return false;
  }

  std::vector<int> starts;
  std::vector<int> ends;
  // Strand is either '+' or '-'. The barcode will be reverse-complemented after
  // extraction if strand is '-'.
  char strand;
  // Whether the range has been modified by new input.
  bool default_range;
};

}  // namespace chromap

#endif
