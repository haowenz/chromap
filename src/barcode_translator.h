#include <cinttypes>
#include <cstring>
#include <functional>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "khash.h"

namespace chromap {

KHASH_INIT(k64_str, uint64_t, char *, 1, kh_int64_hash_func, kh_int64_hash_equal);

// The class for handling barcode convertion.
class BarcodeTranslator
{
private:
  khash_t(k64_str) *barcode_translate_table;
  int from_bc_length;
  uint64_t mask;

  std::string Seed2Sequence(uint64_t seed, uint32_t seed_length) const {
    std::string sequence;
    sequence.reserve(seed_length);
    uint64_t mask = 3;
    for (uint32_t i = 0; i < seed_length; ++i) {
      sequence.push_back(SequenceBatch::Uint8ToChar(
          (seed >> ((seed_length - 1 - i) * 2)) & mask));
    }
    return sequence;
  }
  
  void ProcessTranslateFileLine(std::string &line) {
    int i;
    int len = line.length();
    std::string to;
    for (i = 0; i < len; ++i) {
      if (line[i] == ',' || line[i] == '\t')
        break;
    }

    to = line.substr(0, i);
    //from = line.substr(i + 1, len - i - 1);
    from_bc_length = len - i - 1;
    uint64_t from_seed = SequenceBatch::GenerateSeedFromSequence(line.c_str(), len, i + 1, from_bc_length);
    
    int khash_return_code;
    khiter_t barcode_translate_table_iter = kh_put(k64_str, barcode_translate_table, from_seed, &khash_return_code);
    kh_value(barcode_translate_table, barcode_translate_table_iter) = strdup(to.c_str());
  }

public:
  BarcodeTranslator() {
    barcode_translate_table = NULL;
    from_bc_length = -1;
  }
  
  ~BarcodeTranslator() {
    if (barcode_translate_table != NULL) {
      khiter_t k ;
      for (k = kh_begin(barcode_translate_table) ; k != kh_end(barcode_translate_table); ++k)
      {
        if (kh_exist(barcode_translate_table, k))
          free(kh_value(barcode_translate_table, k)) ;
      }
      kh_destroy(k64_str, barcode_translate_table);
    }
  }

  void SetTranslateTable(const std::string &file) {
    barcode_translate_table = kh_init(k64_str);
    std::ifstream file_stream(file);
    std::string file_line;
    while (getline(file_stream, file_line)) {
      ProcessTranslateFileLine(file_line); 
    }
    
    mask = (1ull<<(2*from_bc_length)) - 1;
    /*for (int i = 0; i < from_bc_length; ++i)
    {
      mask |= (3ull << (2*i));
    }*/
  }
  
  std::string Translate(uint64_t bc, uint32_t bc_length) {
    if (barcode_translate_table == NULL) {
      return Seed2Sequence(bc, bc_length);
    } 

    std::string ret;  
    uint64_t i;
    for (i = 0; i < bc_length / from_bc_length; ++i) {
      uint64_t seed = (bc << (2 * i * from_bc_length)) >> (2 * (bc_length / from_bc_length - 1) * from_bc_length);
      seed &= mask;
      khiter_t barcode_translate_table_iter = kh_get(k64_str, barcode_translate_table, seed);
      if (barcode_translate_table_iter == kh_end(barcode_translate_table)) {
        std::cerr << "Barcode does not exist in the translation table.\n" << std::endl;
        exit(-1);
      }
      std::string bc_to(kh_value(barcode_translate_table, barcode_translate_table_iter));
      if (i == 0) {
        ret = bc_to;
      } else {
        ret += "-" + bc_to;
      }
    }
    return ret;
  }
};
}
