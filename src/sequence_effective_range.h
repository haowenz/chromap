#ifndef SEQUENCE_EFFECTIVE_RANGE_H_
#define SEQUENCE_EFFECTIVE_RANGE_H_

#include <stdlib.h>

#include <algorithm>
#include <vector>

#include "utils.h"

namespace chromap {

// The class handles the custom read format indicating
// the effective range on a sequence.
class SequenceEffectiveRange {
 public:
  SequenceEffectiveRange() {}
  ~SequenceEffectiveRange() {}

  void Init() {
    starts.push_back(0);
    ends.push_back(-1);
    strand = 1;
    range_num = 1;
  }

  // Return: false - fail to parse
  bool ParseEffectiveRange(const char *s, int len) {
    int i;
    int j = 0;  // start, end, strand section
    char buffer[20];
    int blen = 0;
    starts.clear();
    ends.clear();
    strand = 1;
    for (i = 3; i <= len; ++i) {
      if (i == len || s[i] == '/' || s[i] == ':') {
        buffer[blen] = '\0';
        if (j == 0)
          starts.push_back(atoi(buffer));
        else if (j == 1)
          ends.push_back(atoi(buffer));
        else
          strand = buffer[0] == '+' ? 1 : -1;
        blen = 0;
        if (i < len && s[i] == ':') ++j;
      } else {
        buffer[blen] = s[i];
        ++blen;
      }
    }
    if (j >= 3 || starts.size() != ends.size()) return false;
    std::sort(starts.begin(), starts.end());
    std::sort(ends.begin(), ends.end());
    range_num = starts.size();
    if (ends[0] == -1) {
      for (i = 0; i < range_num - 1; ++i) ends[i] = ends[i + 1];
      ends[i] = -1;
    }
    return true;
  }

  // Check whether the effective range is the full range
  bool IsFullRange() {
    if (strand == 1 && starts[0] == 0 && ends[0] == -1) return true;
    return false;
  }

  // Replace by the range specified in the starts,ends section, but does not
  // apply the strand operation. return new length
  int Replace(char *s, int len, bool need_complement) {
    if (IsFullRange()) return len;
    int i, j, k;
    i = 0;
    for (k = 0; k < range_num; ++k) {
      int start = starts[k];
      int end = ends[k];
      if (end == -1) end = len - 1;
      for (j = start; j <= end; ++i, ++j) s[i] = s[j];
    }
    s[i] = '\0';
    len = i;
    if (strand == -1) {
      if (need_complement) {
        for (i = 0; i < len; ++i)
          s[i] = Uint8ToChar(((uint8_t)3) ^ (CharToUint8(s[i])));
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
  std::vector<int> starts;
  std::vector<int> ends;
  int range_num;
  int strand;
};

}  // namespace chromap

#endif
