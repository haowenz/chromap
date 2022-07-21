#ifndef MINIMIZER_GENERATOR_H_
#define MINIMIZER_GENERATOR_H_

#include <cassert>
#include <cstdint>
#include <vector>

#include "minimizer.h"
#include "sequence_batch.h"

namespace chromap {

class MinimizerGenerator {
 public:
  MinimizerGenerator() = delete;

  MinimizerGenerator(int kmer_size, int window_size)
      : kmer_size_(kmer_size), window_size_(window_size) {
    // 56 bits for a k-mer. So the max kmer size is 28.
    assert(kmer_size_ > 0 && kmer_size_ <= 28);
    assert(window_size_ > 0 && window_size_ < 256);
  }

  ~MinimizerGenerator() = default;

  void GenerateMinimizers(const SequenceBatch &sequence_batch,
                          uint32_t sequence_index,
                          std::vector<Minimizer> &minimizers) const;

 private:
  const int kmer_size_;
  const int window_size_;
};

}  // namespace chromap

#endif  // MINIMIZER_GENERATOR_H_
