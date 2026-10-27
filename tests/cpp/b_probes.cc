// Area B probes: each prints what one public entry point of the mapping
// core returns, so tests/test_b_units.py can compare it with values computed
// in Python from definitions (helpers_b.py).
//
//   ksw WINDOW READ A B O_DEL E_DEL O_INS E_INS W [nocigar]
//       ksw_semi_global3 -> "score start end cigar" (ksw's own naming: the
//       o_del/e_del pair prices gaps in READ, o_ins/e_ins gaps in WINDOW).
//   draft REF_FASTA READS_FASTA E SPLIT   (stdin: read and candidate lines)
//       DraftMappingGenerator::GenerateDraftMappings.
//   candidates INDEX READS_FASTA E S F0 F1
//       CandidateProcessor::GenerateCandidates.
//   positions INDEX READS_FASTA MAX_FREQ REP_FREQ HEAP
//       Index::GenerateCandidatePositions.
//   matewindow INDEX READS_FASTA MIN_SEEDS F0 E  (stdin: one query per line)
//       Index::GenerateCandidatePositionsFromRepetitiveReadWithMateInfoOnOneStrand.
//   reduce DISTANCE   (stdin: four candidate lists)
//       CandidateProcessor::ReduceCandidatesForPairedEndRead.

#include "chromap_test.h"

#include "candidate.h"
#include "candidate_position_generating_config.h"
#include "candidate_processor.h"
#include "draft_mapping.h"
#include "draft_mapping_generator.h"
#include "index.h"
#include "ksw.h"
#include "mapping_metadata.h"
#include "mapping_parameters.h"
#include "minimizer_generator.h"
#include "paired_end_mapping_metadata.h"
#include "sequence_batch.h"
#include "sequence_effective_range.h"

using namespace chromap;

namespace {

class ProbeMetadata : public MappingMetadata {
 public:
  using MappingMetadata::minimizers_;
  using MappingMetadata::negative_candidates_;
  using MappingMetadata::negative_candidates_buffer_;
  using MappingMetadata::negative_hits_;
  using MappingMetadata::negative_mappings_;
  using MappingMetadata::negative_split_sites_;
  using MappingMetadata::positive_candidates_;
  using MappingMetadata::positive_candidates_buffer_;
  using MappingMetadata::positive_hits_;
  using MappingMetadata::positive_mappings_;
  using MappingMetadata::positive_split_sites_;
  using MappingMetadata::repetitive_seed_length_;
};

class ProbePairedMetadata : public PairedEndMappingMetadata {
 public:
  using PairedEndMappingMetadata::mapping_metadata1_;
  using PairedEndMappingMetadata::mapping_metadata2_;
};

// Reads (FASTA/FASTQ) in a batch that also holds their reverse complements.
uint32_t LoadReads(const std::string &path, SequenceBatch &batch) {
  batch.InitializeLoading(path);
  uint32_t n = batch.LoadBatch();
  batch.FinalizeLoading();
  for (uint32_t i = 0; i < n; ++i) batch.PrepareNegativeSequenceAt(i);
  return n;
}

void PrintCandidates(const char *tag, const std::vector<Candidate> &candidates) {
  for (size_t i = 0; i < candidates.size(); ++i) {
    std::cout << tag << " " << candidates[i].position << " " << (int)candidates[i].count << "\n";
  }
}

void PrintHits(const char *tag, const std::vector<uint64_t> &hits) {
  for (size_t i = 0; i < hits.size(); ++i) {
    if (hits[i] != UINT64_MAX) std::cout << tag << " " << hits[i] << "\n";
  }
}

std::vector<Candidate> ParseCandidates(std::istringstream &in) {
  std::vector<Candidate> out;
  std::string item;
  while (in >> item) {
    size_t colon = item.find(':');
    Candidate c;
    c.position = std::stoull(item.substr(0, colon));
    c.count = (uint8_t)std::stoi(item.substr(colon + 1));
    out.push_back(c);
  }
  return out;
}

}  // namespace

CT_PROBE(ksw) {
  if (args.size() < 9) {
    std::cerr << "usage: ksw WINDOW READ A B O_DEL E_DEL O_INS E_INS W [nocigar]" << std::endl;
    return 2;
  }
  const std::string &window = args[0];
  const std::string &read = args[1];
  int match = std::stoi(args[2]), mismatch = std::stoi(args[3]);
  int8_t mat[25];
  int k = 0;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) mat[k++] = i == j ? match : -mismatch;
    mat[k++] = 0;
  }
  for (int j = 0; j < 5; ++j) mat[k++] = 0;
  const bool want_cigar = args.size() < 10;
  int n_cigar = 0, start = -1, end = -1;
  uint32_t *cigar = NULL;
  int score = ksw_semi_global3((int)window.size(), window.c_str(), (int)read.size(), read.c_str(), 5, mat,
                               std::stoi(args[4]), std::stoi(args[5]), std::stoi(args[6]), std::stoi(args[7]),
                               std::stoi(args[8]), want_cigar ? &n_cigar : NULL, want_cigar ? &cigar : NULL,
                               &start, &end);
  std::cout << score << " " << start << " " << end << " ";
  for (int i = 0; i < n_cigar; ++i) std::cout << (cigar[i] >> 4) << "MIDNSHP=X"[cigar[i] & 0xf];
  std::cout << "\n";
  free(cigar);
  return 0;
}

// stdin: "R <read index> <number of minimizers>" starts a read, "C <+|-> <rid>
// <position> <count>" adds a candidate, "E" maps it. Output per read: "S
// <min errors> <second min errors> <#best> <#second best>", then "M <+|->
// <errors> <rid> <end position>" per draft mapping (+ first, each in the
// generator's order) and, with SPLIT, "X <+|-> <split site>".
CT_PROBE(draft) {
  if (args.size() != 4) {
    std::cerr << "usage: draft REF_FASTA READS_FASTA E SPLIT" << std::endl;
    return 2;
  }
  SequenceBatch reference;
  reference.InitializeLoading(args[0]);
  reference.LoadAllSequences();
  SequenceBatch reads(100000, SequenceEffectiveRange());
  LoadReads(args[1], reads);
  MappingParameters parameters;
  parameters.error_threshold = std::stoi(args[2]);
  parameters.split_alignment = std::stoi(args[3]) != 0;
  DraftMappingGenerator generator(parameters);
  ProbeMetadata metadata;
  uint32_t read_index = 0;
  std::string line;
  while (std::getline(std::cin, line)) {
    std::istringstream in(line);
    std::string tag;
    in >> tag;
    if (tag == "R") {
      int num_minimizers = 0;
      in >> read_index >> num_minimizers;
      metadata.PrepareForMappingNextRead(64);
      for (int i = 0; i < num_minimizers; ++i) metadata.minimizers_.emplace_back(0, 0);
    } else if (tag == "C") {
      std::string strand;
      uint64_t rid = 0, position = 0;
      int count = 0;
      in >> strand >> rid >> position >> count;
      Candidate c;
      c.position = (rid << 32) | position;
      c.count = (uint8_t)count;
      (strand == "+" ? metadata.positive_candidates_ : metadata.negative_candidates_).push_back(c);
    } else if (tag == "E") {
      generator.GenerateDraftMappings(reads, read_index, reference, metadata);
      std::cout << "S " << metadata.GetMinNumErrors() << " " << metadata.GetSecondMinNumErrors() << " "
                << metadata.GetNumBestMappings() << " " << metadata.GetNumSecondBestMappings() << "\n";
      for (int s = 0; s < 2; ++s) {
        const std::vector<DraftMapping> &mappings = s == 0 ? metadata.positive_mappings_ : metadata.negative_mappings_;
        for (size_t i = 0; i < mappings.size(); ++i) {
          std::cout << "M " << (s == 0 ? "+" : "-") << " " << mappings[i].GetNumErrors() << " "
                    << mappings[i].GetReferenceSequenceIndex() << " " << mappings[i].GetReferenceSequencePosition()
                    << "\n";
        }
        const std::vector<int> &sites = s == 0 ? metadata.positive_split_sites_ : metadata.negative_split_sites_;
        for (size_t i = 0; i < sites.size(); ++i) std::cout << "X " << (s == 0 ? "+" : "-") << " " << sites[i] << "\n";
      }
      std::cout << "END\n";
    }
  }
  return 0;
}

// Per read: "R <index> <#minimizers> <repetitive seed length>", then "P"/"N"
// <position> <count> per positive/negative candidate and "HP"/"HN" <hit> for
// the sorted candidate positions (hits) they were built from.
CT_PROBE(candidates) {
  if (args.size() != 6) {
    std::cerr << "usage: candidates INDEX READS_FASTA E S F0 F1" << std::endl;
    return 2;
  }
  Index index(args[0]);
  index.Load();
  SequenceBatch reads(100000, SequenceEffectiveRange());
  uint32_t n = LoadReads(args[1], reads);
  MinimizerGenerator minimizer_generator(index.GetKmerSize(), index.GetWindowSize());
  std::vector<int> frequencies = {std::stoi(args[4]), std::stoi(args[5])};
  CandidateProcessor processor(std::stoi(args[3]), frequencies);
  ProbeMetadata metadata;
  for (uint32_t i = 0; i < n; ++i) {
    metadata.PrepareForMappingNextRead(frequencies[0]);
    minimizer_generator.GenerateMinimizers(reads, i, metadata.minimizers_);
    processor.GenerateCandidates(std::stoi(args[2]), index, metadata);
    std::cout << "R " << i << " " << metadata.GetNumMinimizers() << " " << metadata.repetitive_seed_length_ << "\n";
    PrintCandidates("P", metadata.positive_candidates_);
    PrintCandidates("N", metadata.negative_candidates_);
    PrintHits("HP", metadata.positive_hits_);
    PrintHits("HN", metadata.negative_hits_);
  }
  return 0;
}

// Per read: "R <index> <#repetitive seeds> <repetitive seed length>" then
// "HP"/"HN" <hit> in the order returned.
CT_PROBE(positions) {
  if (args.size() != 5) {
    std::cerr << "usage: positions INDEX READS_FASTA MAX_FREQ REP_FREQ HEAP" << std::endl;
    return 2;
  }
  Index index(args[0]);
  index.Load();
  SequenceBatch reads(100000, SequenceEffectiveRange());
  uint32_t n = LoadReads(args[1], reads);
  MinimizerGenerator minimizer_generator(index.GetKmerSize(), index.GetWindowSize());
  CandidatePositionGeneratingConfig config(std::stoi(args[2]), std::stoi(args[3]), std::stoi(args[4]) != 0);
  ProbeMetadata metadata;
  for (uint32_t i = 0; i < n; ++i) {
    metadata.PrepareForMappingNextRead(64);
    minimizer_generator.GenerateMinimizers(reads, i, metadata.minimizers_);
    int repetitive = index.GenerateCandidatePositions(config, metadata);
    std::cout << "R " << i << " " << repetitive << " " << metadata.repetitive_seed_length_ << "\n";
    PrintHits("HP", metadata.positive_hits_);
    PrintHits("HN", metadata.negative_hits_);
  }
  return 0;
}

// stdin: one query per line, "<read index> <+|-> <search range>
// <position>:<count> ..." (the mate's candidates). Output per query: "RET
// <return value> <repetitive seed length>", then "H <position>" per
// candidate position found, then "END".
CT_PROBE(matewindow) {
  if (args.size() != 5) {
    std::cerr << "usage: matewindow INDEX READS_FASTA MIN_SEEDS F0 E" << std::endl;
    return 2;
  }
  Index index(args[0]);
  index.Load();
  SequenceBatch reads(100000, SequenceEffectiveRange());
  LoadReads(args[1], reads);
  MinimizerGenerator minimizer_generator(index.GetKmerSize(), index.GetWindowSize());
  std::string line;
  while (std::getline(std::cin, line)) {
    std::istringstream in(line);
    uint32_t read_index = 0, search_range = 0;
    std::string strand;
    if (!(in >> read_index >> strand >> search_range)) continue;
    std::vector<Candidate> mate_candidates = ParseCandidates(in);
    std::vector<Minimizer> minimizers;
    minimizer_generator.GenerateMinimizers(reads, read_index, minimizers);
    uint32_t repetitive_seed_length = 12345;
    std::vector<uint64_t> positions;
    int ret = index.GenerateCandidatePositionsFromRepetitiveReadWithMateInfoOnOneStrand(
        strand == "+" ? kPositive : kNegative, search_range, std::stoi(args[2]), std::stoi(args[3]),
        std::stoi(args[4]), minimizers, mate_candidates, repetitive_seed_length, positions);
    std::cout << "RET " << ret << " " << repetitive_seed_length << "\n";
    for (size_t i = 0; i < positions.size(); ++i) std::cout << "H " << positions[i] << "\n";
    std::cout << "END\n";
  }
  return 0;
}

// stdin: four lines "p1|n1|p2|n2 <position>:<count> ...", the candidate
// buffers of read 1 and read 2 on each strand. Output: the filtered lists in
// the same format.
CT_PROBE(reduce) {
  if (args.size() != 1) {
    std::cerr << "usage: reduce DISTANCE" << std::endl;
    return 2;
  }
  ProbeMetadata m1, m2;
  m1.PrepareForMappingNextRead(16);
  m2.PrepareForMappingNextRead(16);
  std::string line;
  while (std::getline(std::cin, line)) {
    std::istringstream in(line);
    std::string tag;
    if (!(in >> tag)) continue;
    std::vector<Candidate> list = ParseCandidates(in);
    if (tag == "p1") m1.positive_candidates_buffer_ = list;
    if (tag == "n1") m1.negative_candidates_buffer_ = list;
    if (tag == "p2") m2.positive_candidates_buffer_ = list;
    if (tag == "n2") m2.negative_candidates_buffer_ = list;
  }
  ProbePairedMetadata pe;
  pe.PreparedForMappingNextReadPair(16);
  pe.mapping_metadata1_ = m1;
  pe.mapping_metadata2_ = m2;
  CandidateProcessor processor(2, std::vector<int>({500, 1000}));
  processor.ReduceCandidatesForPairedEndRead((uint32_t)std::stoul(args[0]), pe);
  ProbeMetadata out1, out2;
  static_cast<MappingMetadata &>(out1) = pe.mapping_metadata1_;
  static_cast<MappingMetadata &>(out2) = pe.mapping_metadata2_;
  const char *tags[4] = {"p1", "n1", "p2", "n2"};
  const std::vector<Candidate> *lists[4] = {&out1.positive_candidates_, &out1.negative_candidates_,
                                            &out2.positive_candidates_, &out2.negative_candidates_};
  for (int i = 0; i < 4; ++i) {
    std::cout << tags[i];
    for (size_t j = 0; j < lists[i]->size(); ++j) {
      std::cout << " " << (*lists[i])[j].position << ":" << (int)(*lists[i])[j].count;
    }
    std::cout << "\n";
  }
  return 0;
}

CT_MAIN
