#ifndef BARCODETRANSLATOR_H_
#define BARCODETRANSLATOR_H_

#include <cinttypes>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "khash.h"
#include "sequence_batch.h"

namespace chromap {

KHASH_INIT(k64_str, uint64_t, char *, 1, kh_int64_hash_func,
           kh_int64_hash_equal);

// The class for handling barcode convertion.
class BarcodeTranslator {
 public:
  BarcodeTranslator() {
    barcode_translate_table_ = NULL;
    from_bc_length_ = -1;
  }

  ~BarcodeTranslator() {
    if (barcode_translate_table_ != NULL) {
      khiter_t k;
      for (k = kh_begin(barcode_translate_table_);
           k != kh_end(barcode_translate_table_); ++k) {
        if (kh_exist(barcode_translate_table_, k))
          free(kh_value(barcode_translate_table_, k));
      }
      kh_destroy(k64_str, barcode_translate_table_);
    }
  }

  void SetTranslateTable(const std::string &file) {
    barcode_translate_table_ = kh_init(k64_str);
    std::ifstream file_stream(file);
    std::string file_line;
    while (getline(file_stream, file_line)) {
      ProcessTranslateFileLine(file_line);
    }

    mask_ = (1ull << (2 * from_bc_length_)) - 1;
    /*for (int i = 0; i < from_bc_length_; ++i)
    {
      mask_ |= (3ull << (2*i));
    }*/
  }

  std::string Translate(uint64_t bc, uint32_t bc_length) {
    if (barcode_translate_table_ == NULL) {
      return Seed2Sequence(bc, bc_length);
    }

    std::string ret;
    uint64_t i;
    for (i = 0; i < bc_length / from_bc_length_; ++i) {
      uint64_t seed = (bc << (2 * i * from_bc_length_)) >>
                      (2 * (bc_length / from_bc_length_ - 1) * from_bc_length_);
      seed &= mask_;
      khiter_t barcode_translate_table_iter =
          kh_get(k64_str, barcode_translate_table_, seed);
      if (barcode_translate_table_iter == kh_end(barcode_translate_table_)) {
        std::cerr << "Barcode does not exist in the translation table."
                  << std::endl;
        exit(-1);
      }
      std::string bc_to(
          kh_value(barcode_translate_table_, barcode_translate_table_iter));
      if (i == 0) {
        ret = bc_to;
      } else {
        ret += "-" + bc_to;
      }
    }
    return ret;
  }

 private:
  khash_t(k64_str) * barcode_translate_table_;
  int from_bc_length_;
  uint64_t mask_;

  std::string Seed2Sequence(uint64_t seed, uint32_t seed_length) const {
    std::string sequence;
    sequence.reserve(seed_length);
    uint64_t mask_ = 3;
    for (uint32_t i = 0; i < seed_length; ++i) {
      sequence.push_back(SequenceBatch::Uint8ToChar(
          (seed >> ((seed_length - 1 - i) * 2)) & mask_));
    }
    return sequence;
  }

  void ProcessTranslateFileLine(std::string &line) {
    int i;
    int len = line.length();
    std::string to;
    for (i = 0; i < len; ++i) {
      if (line[i] == ',' || line[i] == '\t') break;
    }

    to = line.substr(0, i);
    // from = line.substr(i + 1, len - i - 1);
    from_bc_length_ = len - i - 1;
    uint64_t from_seed = SequenceBatch::GenerateSeedFromSequence(
        line.c_str(), len, i + 1, from_bc_length_);

    int khash_return_code;
    khiter_t barcode_translate_table_iter = kh_put(
        k64_str, barcode_translate_table_, from_seed, &khash_return_code);
    kh_value(barcode_translate_table_, barcode_translate_table_iter) =
        strdup(to.c_str());
  }
};

}  // namespace chromap
#endif
