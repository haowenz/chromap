#ifndef DRAFT_MAPPING_GENERATOR_H_
#define DRAFT_MAPPING_GENERATOR_H_

#include <cstdint>

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
