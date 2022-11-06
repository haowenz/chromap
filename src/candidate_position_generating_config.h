#ifndef CANDIDATE_POSITION_GENERATING_CONFIG_H_
#define CANDIDATE_POSITION_GENERATING_CONFIG_H_

namespace chromap {

// This class holds the parameters to generate candidate position. Using the
// parameters, it can check whether a seed is frequent or repetitive. Frequent
// seeds are usually skipped when generating candidate positions. Repetitive
// seeds are the seeds considered when computing repetitive seed length.
class CandidatePositionGeneratingConfig {
 public:
  CandidatePositionGeneratingConfig() = delete;

  CandidatePositionGeneratingConfig(uint32_t max_seed_frequency,
                                    uint32_t repetitive_seed_frequency,
                                    bool use_heap_merge)
      : max_seed_frequency_(max_seed_frequency),
        repetitive_seed_frequency_(repetitive_seed_frequency),
        use_heap_merge_(use_heap_merge) {}

  ~CandidatePositionGeneratingConfig() = default;

  inline bool IsFrequentSeed(uint32_t seed_frequency) const {
    return seed_frequency >= max_seed_frequency_;
  }

  inline bool IsRepetitiveSeed(uint32_t seed_frequency) const {
    return seed_frequency >= repetitive_seed_frequency_;
  }

  inline bool UseHeapMerge() const { return use_heap_merge_; }

  inline uint32_t GetMaxSeedFrequency() const { return max_seed_frequency_; }

 private:
  // Only seeds with frequency less than this threshold will be used.
  const uint32_t max_seed_frequency_;

  // Seeds with frequency greater than or equal to this threshold will be
  // considered as repetitive seeds.
  const uint32_t repetitive_seed_frequency_;

  // When the number of candidate positions is really large, use heap merge to
  // merge sorted candidate lists.
  const bool use_heap_merge_;
};

}  // namespace chromap

#endif  // CANDIDATE_POSITION_GENERATING_CONFIG_H_
