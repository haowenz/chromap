#ifndef DRAFT_MAPPING_GENERATOR_H_
#define DRAFT_MAPPING_GENERATOR_H_

#include <cstdint>
#include <vector>

#include "draft_mapping.h"
#include "mapping_metadata.h"
#include "mapping_parameters.h"
#include "sequence_batch.h"
#include "utils.h"

namespace chromap {

class DraftMappingGenerator {
 public:
  DraftMappingGenerator() = delete;

  DraftMappingGenerator(const MappingParameters &mapping_parameters)
      : error_threshold_(mapping_parameters.error_threshold),
        split_alignment_(mapping_parameters.split_alignment),
        num_vpu_lanes_(mapping_parameters.GetNumVPULanes()),
        mapping_output_format_(mapping_parameters.mapping_output_format) {}

  ~DraftMappingGenerator() = default;

  void GenerateDraftMappings(const SequenceBatch &read_batch,
                             uint32_t read_index,
                             const SequenceBatch &reference,
                             MappingMetadata &mapping_metadata);

 private:
  // Return true if the candidate position is valid on the reference with rid.
  // Note only the position is checked and the input rid is not checked in this
  // function. So the input rid must be valid.
  bool IsValidCandidate(uint32_t rid, uint32_t position, uint32_t read_length,
                        const SequenceBatch &reference);

  // How many candidates ahead of the one being verified the reference window
  // is prefetched.
  static constexpr uint32_t kReferenceWindowPrefetchDistance = 8;

  // Prefetches, one cache line at a time, the reference window that the
  // banded alignment of a valid candidate reads, [position - e, position +
  // read_length + e), so that the DRAM and TLB misses of a read's candidates
  // overlap with each other and with the alignments before them. Does nothing
  // for an invalid candidate. A hint only: no result depends on it, and where
  // the compiler has no prefetch builtin it does nothing. It must be inlined:
  // GCC 13 at -O3 deduces a function whose only effect is __builtin_prefetch
  // to be pure and deletes calls to it.
#if defined(__GNUC__) || defined(__clang__)
  __attribute__((always_inline)) void PrefetchCandidateReferenceWindow(
      const Candidate &candidate, const Strand candidate_strand,
      uint32_t read_length, const SequenceBatch &reference) {
    const uint32_t rid = candidate.GetReferenceSequenceIndex();
    uint32_t position = candidate.GetReferenceSequencePosition();
    if (candidate_strand == kNegative) {
      position = position - read_length + 1;
    }
    if (!IsValidCandidate(rid, position, read_length, reference)) {
      return;
    }
    const uintptr_t window_start = reinterpret_cast<uintptr_t>(
        reference.GetSequenceAt(rid) + position - error_threshold_);
    const uintptr_t window_end =
        window_start + read_length + 2 * error_threshold_;
    for (uintptr_t line = window_start & ~(uintptr_t)63; line < window_end;
         line += 64) {
      __builtin_prefetch(reinterpret_cast<const void *>(line));
    }
  }
#else
  void PrefetchCandidateReferenceWindow(const Candidate & /*candidate*/,
                                        const Strand /*candidate_strand*/,
                                        uint32_t /*read_length*/,
                                        const SequenceBatch & /*reference*/) {}
#endif

  // Prefetches the reference windows of the first
  // kReferenceWindowPrefetchDistance candidates on the strand. Inlined for the
  // same reason as PrefetchCandidateReferenceWindow.
#if defined(__GNUC__) || defined(__clang__)
  __attribute__((always_inline))
#endif
  void PrefetchFirstCandidateReferenceWindows(
      const std::vector<Candidate> &candidates, const Strand candidate_strand,
      uint32_t read_length, const SequenceBatch &reference) {
    const size_t num_candidates =
        candidates.size() < kReferenceWindowPrefetchDistance
            ? candidates.size()
            : kReferenceWindowPrefetchDistance;
    for (size_t ci = 0; ci < num_candidates; ++ci) {
      PrefetchCandidateReferenceWindow(candidates[ci], candidate_strand,
                                       read_length, reference);
    }
  }

  // Return true when there is one non-split mapping generated and the mapping
  // is supported by all the minimizers.
  bool GenerateNonSplitDraftMappingSupportedByAllMinimizers(
      const SequenceBatch &read_batch, uint32_t read_index,
      const SequenceBatch &reference, MappingMetadata &mapping_metadata);

  void GenerateDraftMappingsOnOneStrandUsingSIMD(
      const Strand candidate_strand, uint32_t read_index,
      const SequenceBatch &read_batch, const SequenceBatch &reference,
      MappingMetadata &mapping_metadata);

  void GenerateDraftMappingsOnOneStrand(const Strand candidate_strand,
                                        uint32_t read_index,
                                        const SequenceBatch &read_batch,
                                        const SequenceBatch &reference,
                                        MappingMetadata &mapping_metadata);

  const int error_threshold_;
  const bool split_alignment_;
  const int num_vpu_lanes_;
  const MappingOutputFormat mapping_output_format_;
};

}  // namespace chromap

#endif  // DRAFT_MAPPING_GENERATOR_H_
