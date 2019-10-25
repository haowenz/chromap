#include "output_tools.h"

#include <unistd.h>

namespace fem {
void OutputTools::InitializeMappingOutput() {
  mapping_output_file_ = fopen(mapping_output_file_path_.c_str(), "w");
  assert(mapping_output_file_ != NULL);
}

void OutputTools::AppendMappingOutput(const std::string &line) {
  fprintf(mapping_output_file_, line.data());
}

void OutputTools::FinalizeMappingOutput() {
  fclose(mapping_output_file_);
}

//void OutputTools::GeneratePAFLine(const std::string &query_name, const int query_length, const int query_start, const int query_end, const char relative_strand, const std::string &target_name, const int target_length, const int target_start, const int target_end, const int num_matches, const int alignment_length, const int mapping_quality) {
//}

std::string OutputTools::GeneratePAFLine(const SequenceBatch &query_batch, uint32_t query_index, const int query_start, const int query_end, const char relative_strand, const SequenceBatch &target_batch, uint32_t target_index, const int target_start, const int target_end, const int num_matches, const int alignment_length, const int mapping_quality) {
  std::string paf_string = std::string(query_batch.GetSequenceNameAt(query_index)) + "\t" + std::to_string(query_batch.GetSequenceLengthAt(query_index)) + "\t" + std::to_string(query_start) + "\t" + std::to_string(query_end) + "\t" + relative_strand + "\t" + std::string(target_batch.GetSequenceNameAt(target_index)) + "\t" + std::to_string(target_batch.GetSequenceLengthAt(target_index)) + "\t" + std::to_string(target_start) + "\t" + std::to_string(target_end) + "\t" + std::to_string(num_matches) + "\t" + std::to_string(alignment_length) + "\t" + std::to_string(mapping_quality) + "\n";
    return paf_string;
}
} // namespace fem
