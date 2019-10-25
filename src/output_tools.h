#ifndef OUTPUTTOOLS_H_
#define OUTPUTTOOLS_H_

#include <assert.h>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

namespace chromap {
class OutputTools {
 public:
  OutputTools(const std::string &mapping_output_file_path) : mapping_output_file_path_(mapping_output_file_path) {}
  ~OutputTools() {}
  inline void InitializeMappingOutput() {
    mapping_output_file_ = fopen(mapping_output_file_path_.c_str(), "w");
    assert(mapping_output_file_ != NULL);
  }
  inline void AppendMappingOutput(const std::string &line) {
    fprintf(mapping_output_file_, line.data());
  }
  inline void FinalizeMappingOutput() {
    fclose(mapping_output_file_);
  }
  inline std::string GeneratePAFLine(const SequenceBatch &query_batch, uint32_t query_index, const int query_start, const int query_end, const char relative_strand, const SequenceBatch &target_batch, uint32_t target_index, const int target_start, const int target_end, const int num_matches, const int alignment_length, const int mapping_quality) {
    return std::string(query_batch.GetSequenceNameAt(query_index)) + "\t" + std::to_string(query_batch.GetSequenceLengthAt(query_index)) + "\t" + std::to_string(query_start) + "\t" + std::to_string(query_end) + "\t" + relative_strand + "\t" + std::string(target_batch.GetSequenceNameAt(target_index)) + "\t" + std::to_string(target_batch.GetSequenceLengthAt(target_index)) + "\t" + std::to_string(target_start) + "\t" + std::to_string(target_end) + "\t" + std::to_string(num_matches) + "\t" + std::to_string(alignment_length) + "\t" + std::to_string(mapping_quality) + "\n";
  }
  inline std::string GenerateBEDPELine(const char *reference_sequence_name, uint32_t fragment_start_position, uint16_t fragment_length){
    return std::string(reference_sequence_name) + "\t" + std::to_string(fragment_start_position) + "\t" + std::to_string(fragment_start_position + fragment_length) + "\n";
  }
  inline std::string GeneratePairedEndTagAlignLine(uint8_t strand, const char *reference_sequence_name1, uint32_t positive_read_start, uint32_t positive_read_end, const char *reference_sequence_name2, uint32_t negative_read_start, uint32_t negative_read_end) {
    if (strand == 1) {
      return std::string(reference_sequence_name1) + "\t" + std::to_string(positive_read_start) + "\t" + std::to_string(positive_read_end) + "\tN\t1000\t+\n" + std::string(reference_sequence_name2) + "\t" + std::to_string(negative_read_start) + "\t" + std::to_string(negative_read_end) + "\tN\t1000\t-\n";
    } else {
      return std::string(reference_sequence_name1) + "\t" + std::to_string(negative_read_start) + "\t" + std::to_string(negative_read_end) + "\tN\t1000\t-\n" + std::string(reference_sequence_name2) + "\t" + std::to_string(positive_read_start) + "\t" + std::to_string(positive_read_end) + "\tN\t1000\t+\n";
    }
  }
  inline uint32_t GetNumMappings() const {
    return num_mappings_;
  }

 protected:
  std::string mapping_output_file_path_; 
  FILE *mapping_output_file_;
  uint32_t num_mappings_;
};
} // namespace chromap

#endif // OUTPUTTOOLS_H_
