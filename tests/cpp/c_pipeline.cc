// Area C: the post-processing pipeline driven through the classes chromap
// itself calls (MappingProcessor, MappingWriter, TempMappingFileHandle).
//
// Probes read mapping records from stdin, one per line, and print or write
// what chromap produces; tests/test_c_pipeline.py computes the expected
// values with the reference functions in tests/helpers_c.py.
//
// Record lines (whitespace separated, the last field is the duplicate count):
//   se   rid read_id start length mapq direction unique dups
//   seb  rid read_id barcode start length mapq direction unique dups
//   pe   rid read_id start length mapq direction unique plen nlen dups
//   peb  rid read_id barcode start length mapq direction unique plen nlen dups
//   paf  rid read_id name read_length start length mapq direction unique dups
//   ppaf rid read_id name1 name2 rlen1 rlen2 start length plen nlen mapq mapq1
//        mapq2 direction unique dups
//   sam  rid read_id name barcode pos mrid mpos tlen flag positive mapq nm
//        cigar md seq qual dups      ("*" for no CIGAR, "-" for an empty MD,
//                                     SEQ or QUAL)

#include "chromap_test.h"

#include <fstream>

#include "bed_mapping.h"
#include "khash.h"
#include "mapping_parameters.h"
#include "mapping_processor.h"
#include "mapping_writer.h"
#include "paf_mapping.h"
#include "sam_mapping.h"
#include "sequence_batch.h"
#include "summary_metadata.h"
#include "temp_mapping.h"
#include "utils.h"

using namespace chromap;

namespace {

std::string Dash(const std::string &s) { return s == "-" ? std::string() : s; }
std::string Undash(const std::string &s) { return s.empty() ? std::string("-") : s; }

bool Read(std::istream &in, MappingWithoutBarcode &m, uint32_t &rid) {
  uint32_t id, start, len, mapq, dir, uniq, dups;
  if (!(in >> rid >> id >> start >> len >> mapq >> dir >> uniq >> dups)) return false;
  m = MappingWithoutBarcode(id, start, len, mapq, dir, uniq, 1);
  m.num_dups_ = dups;
  return true;
}
void Write(std::ostream &out, const MappingWithoutBarcode &m) {
  out << m.read_id_ << " " << m.fragment_start_position_ << " " << m.fragment_length_ << " " << (int)m.mapq_
      << " " << (int)m.direction_ << " " << (int)m.is_unique_ << " " << (int)m.num_dups_;
}

bool Read(std::istream &in, MappingWithBarcode &m, uint32_t &rid) {
  uint32_t id, start, len, mapq, dir, uniq, dups;
  uint64_t bc;
  if (!(in >> rid >> id >> bc >> start >> len >> mapq >> dir >> uniq >> dups)) return false;
  m = MappingWithBarcode(id, bc, start, len, mapq, dir, uniq, dups);
  return true;
}
void Write(std::ostream &out, const MappingWithBarcode &m) {
  out << m.read_id_ << " " << m.cell_barcode_ << " " << m.fragment_start_position_ << " " << m.fragment_length_
      << " " << (int)m.mapq_ << " " << (int)m.direction_ << " " << (int)m.is_unique_ << " " << (int)m.num_dups_;
}

bool Read(std::istream &in, PairedEndMappingWithoutBarcode &m, uint32_t &rid) {
  uint32_t id, start, len, mapq, dir, uniq, plen, nlen, dups;
  if (!(in >> rid >> id >> start >> len >> mapq >> dir >> uniq >> plen >> nlen >> dups)) return false;
  m = PairedEndMappingWithoutBarcode(id, start, len, mapq, dir, uniq, dups, plen, nlen);
  return true;
}
void Write(std::ostream &out, const PairedEndMappingWithoutBarcode &m) {
  out << m.read_id_ << " " << m.fragment_start_position_ << " " << m.fragment_length_ << " " << (int)m.mapq_
      << " " << (int)m.direction_ << " " << (int)m.is_unique_ << " " << m.positive_alignment_length_ << " "
      << m.negative_alignment_length_ << " " << (int)m.num_dups_;
}

bool Read(std::istream &in, PairedEndMappingWithBarcode &m, uint32_t &rid) {
  uint32_t id, start, len, mapq, dir, uniq, plen, nlen, dups;
  uint64_t bc;
  if (!(in >> rid >> id >> bc >> start >> len >> mapq >> dir >> uniq >> plen >> nlen >> dups)) return false;
  m = PairedEndMappingWithBarcode(id, bc, start, len, mapq, dir, uniq, dups, plen, nlen);
  return true;
}
void Write(std::ostream &out, const PairedEndMappingWithBarcode &m) {
  out << m.read_id_ << " " << m.cell_barcode_ << " " << m.fragment_start_position_ << " " << m.fragment_length_
      << " " << (int)m.mapq_ << " " << (int)m.direction_ << " " << (int)m.is_unique_ << " "
      << m.positive_alignment_length_ << " " << m.negative_alignment_length_ << " " << (int)m.num_dups_;
}

bool Read(std::istream &in, PAFMapping &m, uint32_t &rid) {
  uint32_t id, rlen, start, len, mapq, dir, uniq, dups;
  std::string name;
  if (!(in >> rid >> id >> name >> rlen >> start >> len >> mapq >> dir >> uniq >> dups)) return false;
  m = PAFMapping(id, name, rlen, start, len, mapq, dir, uniq, dups);
  return true;
}
void Write(std::ostream &out, const PAFMapping &m) {
  out << m.read_id_ << " " << m.read_name_ << " " << m.read_length_ << " " << m.fragment_start_position_ << " "
      << m.fragment_length_ << " " << (int)m.mapq_ << " " << (int)m.direction_ << " " << (int)m.is_unique_ << " "
      << (int)m.num_dups_;
}

bool Read(std::istream &in, PairedPAFMapping &m, uint32_t &rid) {
  uint32_t id, rlen1, rlen2, start, len, plen, nlen, mapq, mapq1, mapq2, dir, uniq, dups;
  std::string name1, name2;
  if (!(in >> rid >> id >> name1 >> name2 >> rlen1 >> rlen2 >> start >> len >> plen >> nlen >> mapq >> mapq1 >>
        mapq2 >> dir >> uniq >> dups))
    return false;
  m = PairedPAFMapping(id, name1, name2, rlen1, rlen2, start, len, plen, nlen, mapq, mapq1, mapq2, dir, uniq, dups);
  return true;
}
void Write(std::ostream &out, const PairedPAFMapping &m) {
  out << m.read_id_ << " " << m.read1_name_ << " " << m.read2_name_ << " " << m.read1_length_ << " "
      << m.read2_length_ << " " << m.fragment_start_position_ << " " << m.fragment_length_ << " "
      << m.positive_alignment_length_ << " " << m.negative_alignment_length_ << " " << (int)m.mapq_ << " "
      << (int)m.mapq1_ << " " << (int)m.mapq2_ << " " << (int)m.direction_ << " " << (int)m.is_unique_ << " "
      << (int)m.num_dups_;
}

// "20M1I29M" -> a malloc'ed BAM CIGAR array (the SAMMapping constructor frees
// it).
uint32_t *ParseCigar(const std::string &text, int *n) {
  std::vector<uint32_t> ops;
  uint32_t len = 0;
  if (text != "*") {
    for (size_t i = 0; i < text.size(); ++i) {
      if (text[i] >= '0' && text[i] <= '9') {
        len = len * 10 + (text[i] - '0');
      } else {
        ops.push_back(bam_cigar_gen(len, (uint32_t)bam_cigar_table[(unsigned char)text[i]]));
        len = 0;
      }
    }
  }
  uint32_t *cigar = (uint32_t *)malloc(sizeof(uint32_t) * (ops.size() + 1));
  for (size_t i = 0; i < ops.size(); ++i) cigar[i] = ops[i];
  *n = (int)ops.size();
  return cigar;
}

bool Read(std::istream &in, SAMMapping &m, uint32_t &rid) {
  uint32_t id, flag, positive, mapq, nm, dups;
  int64_t pos, mpos;
  int mrid, tlen;
  uint64_t bc;
  std::string name, cigar, md, seq, qual;
  if (!(in >> rid >> id >> name >> bc >> pos >> mrid >> mpos >> tlen >> flag >> positive >> mapq >> nm >> cigar >>
        md >> seq >> qual >> dups))
    return false;
  int n = 0;
  uint32_t *ops = ParseCigar(cigar, &n);
  m = SAMMapping(id, name, bc, dups, pos, rid, mpos, mrid, tlen, flag, positive, 0, 1, mapq, nm, n, ops, Dash(md),
                 Dash(seq), Dash(qual));
  return true;
}
void Write(std::ostream &out, const SAMMapping &m) {
  out << m.read_id_ << " " << m.read_name_ << " " << m.cell_barcode_ << " " << m.pos_ << " " << m.mrid_ << " "
      << m.mpos_ << " " << m.tlen_ << " " << m.flag_ << " " << (m.IsPositiveStrand() ? 1 : 0) << " " << (int)m.mapq_
      << " " << (int)m.NM_ << " " << m.GenerateCigarString() << " " << Undash(m.MD_) << " " << Undash(m.sequence_)
      << " " << Undash(m.sequence_qual_) << " " << (int)m.num_dups_;
}

template <typename T>
uint32_t ReadRecords(std::istream &in, uint32_t num_refs, std::vector<std::vector<T> > &mappings) {
  mappings.assign(num_refs, std::vector<T>());
  uint32_t n = 0;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::istringstream fields(line);
    T m;
    uint32_t rid = 0;
    if (!Read(fields, m, rid) || rid >= num_refs) {
      std::cerr << "bad record line: " << line << std::endl;
      exit(2);
    }
    mappings[rid].push_back(m);
    ++n;
  }
  return n;
}

MappingOutputFormat ParseFormat(const std::string &name) {
  if (name == "BED") return MAPPINGFORMAT_BED;
  if (name == "TagAlign") return MAPPINGFORMAT_TAGALIGN;
  if (name == "PAF") return MAPPINGFORMAT_PAF;
  if (name == "SAM") return MAPPINGFORMAT_SAM;
  std::cerr << "unknown format " << name << std::endl;
  exit(2);
}

// inmemory TYPE FORMAT FASTA OUT [OPS...]: records from stdin into
// per-reference vectors, the OPS applied in order through MappingProcessor
// (tn5, sort, psort:THREADS, dedup:THREADS), then MappingWriter writes the
// header and the mappings to OUT. Options: mapq:Q (MAPQ threshold),
// bc:LENGTH (barcode length given to the writer).
template <typename T>
int RunInMemory(const std::vector<std::string> &args) {
  MappingParameters params;
  params.mapping_output_format = ParseFormat(args[1]);
  params.mapping_output_file_path = args[3];
  uint32_t barcode_length = 0;
  std::vector<std::string> ops;
  for (size_t i = 4; i < args.size(); ++i) {
    if (args[i].compare(0, 5, "mapq:") == 0) {
      params.mapq_threshold = (uint8_t)std::stoi(args[i].substr(5));
    } else if (args[i].compare(0, 3, "bc:") == 0) {
      barcode_length = std::stoi(args[i].substr(3));
    } else {
      ops.push_back(args[i]);
    }
  }
  SequenceBatch reference;
  reference.InitializeLoading(args[2]);
  reference.LoadAllSequences();
  const uint32_t num_refs = reference.GetNumSequences();
  std::vector<std::vector<T> > mappings;
  ReadRecords(std::cin, num_refs, mappings);
  MappingProcessor<T> processor(params, /*min_unique_mapping_mapq=*/4);
  for (size_t i = 0; i < ops.size(); ++i) {
    const std::string &op = ops[i];
    if (op == "tn5") {
      processor.ApplyTn5ShiftOnMappings(num_refs, mappings);
    } else if (op == "sort") {
      processor.SortOutputMappings(num_refs, mappings);
    } else if (op.compare(0, 6, "psort:") == 0) {
      processor.ParallelSortOutputMappings(num_refs, mappings, std::stoi(op.substr(6)));
    } else if (op.compare(0, 6, "dedup:") == 0) {
      processor.RemovePCRDuplicate(num_refs, mappings, std::stoi(op.substr(6)));
    } else {
      std::cerr << "unknown op " << op << std::endl;
      return 2;
    }
  }
  {
    MappingWriter<T> writer(params, barcode_length, std::vector<int>());
    writer.OutputHeader(num_refs, reference);
    writer.OutputMappings(num_refs, reference, mappings);
  }
  reference.FinalizeLoading();
  return 0;
}

CT_PROBE(inmemory) {
  if (args.size() < 4) {
    std::cerr << "usage: inmemory TYPE FORMAT FASTA OUT [OPS...]" << std::endl;
    return 2;
  }
  const std::string &type = args[0];
  if (type == "se") return RunInMemory<MappingWithoutBarcode>(args);
  if (type == "seb") return RunInMemory<MappingWithBarcode>(args);
  if (type == "pe") return RunInMemory<PairedEndMappingWithoutBarcode>(args);
  if (type == "peb") return RunInMemory<PairedEndMappingWithBarcode>(args);
  if (type == "paf") return RunInMemory<PAFMapping>(args);
  if (type == "ppaf") return RunInMemory<PairedPAFMapping>(args);
  if (type == "sam") return RunInMemory<SAMMapping>(args);
  std::cerr << "unknown type " << type << std::endl;
  return 2;
}

// lowmem TYPE FORMAT FASTA OUT MAPQ TN5 DEDUP NUM_FILES: stdin lines are
// "FILE_INDEX <record>". Each file's records are sorted
// (MappingProcessor::SortOutputMappings) and written as one temp file
// (MappingWriter::OutputTempMappings), then
// ProcessAndOutputMappingsInLowMemory merges, deduplicates and writes OUT.
template <typename T>
int RunLowMem(const std::vector<std::string> &args) {
  MappingParameters params;
  params.mapping_output_format = ParseFormat(args[1]);
  params.mapping_output_file_path = args[3];
  params.mapq_threshold = (uint8_t)std::stoi(args[4]);
  params.Tn5_shift = std::stoi(args[5]) != 0;
  params.remove_pcr_duplicates = std::stoi(args[6]) != 0;
  const int num_files = std::stoi(args[7]);
  SequenceBatch reference;
  reference.InitializeLoading(args[2]);
  reference.LoadAllSequences();
  const uint32_t num_refs = reference.GetNumSequences();
  std::vector<std::string> per_file(num_files);
  std::string line;
  while (std::getline(std::cin, line)) {
    if (line.empty()) continue;
    std::istringstream fields(line);
    int file_index = -1;
    fields >> file_index;
    if (file_index < 0 || file_index >= num_files) {
      std::cerr << "bad file index: " << line << std::endl;
      return 2;
    }
    std::string rest;
    std::getline(fields, rest);
    per_file[file_index] += rest + "\n";
  }
  khash_t(k64_seq) *whitelist = kh_init(k64_seq);
  {
    MappingProcessor<T> processor(params, /*min_unique_mapping_mapq=*/4);
    MappingWriter<T> writer(params, 0, std::vector<int>());
    writer.OutputHeader(num_refs, reference);
    std::vector<TempMappingFileHandle<T> > handles;
    for (int f = 0; f < num_files; ++f) {
      std::istringstream in(per_file[f]);
      std::vector<std::vector<T> > mappings;
      ReadRecords(in, num_refs, mappings);
      processor.SortOutputMappings(num_refs, mappings);
      writer.OutputTempMappings(num_refs, mappings, handles);
      for (uint32_t r = 0; r < num_refs; ++r) {
        if (!mappings[r].empty()) {
          std::cerr << "OutputTempMappings left mappings behind" << std::endl;
          return 3;
        }
      }
    }
    writer.ProcessAndOutputMappingsInLowMemory(0, num_refs, reference, whitelist, handles);
  }
  kh_destroy(k64_seq, whitelist);
  reference.FinalizeLoading();
  return 0;
}

CT_PROBE(lowmem) {
  if (args.size() != 8) {
    std::cerr << "usage: lowmem TYPE FORMAT FASTA OUT MAPQ TN5 DEDUP NUM_FILES" << std::endl;
    return 2;
  }
  // Only record types whose low-memory merge buffers total 1 GiB (PAF, SAM);
  // BED records get 10 GiB of buffers (see tests/README.md, Runtime).
  if (args[0] == "paf") return RunLowMem<PAFMapping>(args);
  if (args[0] == "sam") return RunLowMem<SAMMapping>(args);
  std::cerr << "unsupported type " << args[0] << std::endl;
  return 2;
}

// temp TYPE BLOCK_SIZE NUM_REFS OUT: records from stdin written as one temp
// file with MappingWriter::OutputTempMappings, read back with
// TempMappingFileHandle in blocks of BLOCK_SIZE; prints "path <temp file>",
// then "rid <record>" per mapping in the order the handle yields them, then
// "end", and removes the file.
template <typename T>
int RunTemp(const std::vector<std::string> &args) {
  const uint32_t block = std::stoi(args[1]);
  const uint32_t num_refs = std::stoi(args[2]);
  MappingParameters params;
  params.mapping_output_format = MAPPINGFORMAT_PAF;
  params.mapping_output_file_path = args[3];
  std::vector<std::vector<T> > mappings;
  ReadRecords(std::cin, num_refs, mappings);
  std::vector<TempMappingFileHandle<T> > handles;
  {
    MappingWriter<T> writer(params, 0, std::vector<int>());
    writer.OutputTempMappings(num_refs, mappings, handles);
  }
  if (handles.size() != 1) {
    std::cerr << "expected one temp file handle" << std::endl;
    return 3;
  }
  TempMappingFileHandle<T> &h = handles[0];
  std::cout << "path " << h.file_path << "\n";
  h.InitializeTempMappingLoading(block);
  h.LoadTempMappingBlock(num_refs);
  uint64_t steps = 0;
  while (h.HasMappings()) {
    std::cout << h.current_rid << " ";
    Write(std::cout, h.GetCurrentMapping());
    std::cout << "\n";
    h.Next(num_refs);
    if (++steps > 10000000) return 4;
  }
  h.FinalizeTempMappingLoading();
  remove(h.file_path.c_str());
  std::cout << "end\n";
  return 0;
}

CT_PROBE(temp) {
  if (args.size() != 4) {
    std::cerr << "usage: temp TYPE BLOCK_SIZE NUM_REFS OUT" << std::endl;
    return 2;
  }
  const std::string &type = args[0];
  if (type == "se") return RunTemp<MappingWithoutBarcode>(args);
  if (type == "seb") return RunTemp<MappingWithBarcode>(args);
  if (type == "pe") return RunTemp<PairedEndMappingWithoutBarcode>(args);
  if (type == "peb") return RunTemp<PairedEndMappingWithBarcode>(args);
  if (type == "paf") return RunTemp<PAFMapping>(args);
  if (type == "ppaf") return RunTemp<PairedPAFMapping>(args);
  if (type == "sam") return RunTemp<SAMMapping>(args);
  std::cerr << "unknown type " << type << std::endl;
  return 2;
}

// allocate NUM_REFS SEED DISTANCE: single-end records ("se" lines) from
// stdin, sorted (SortOutputMappings) and passed to AllocateMultiMappings
// with chromap's minimum unique MAPQ of 4; prints "rid <record>" for every
// mapping left, per reference in order.
CT_PROBE(allocate) {
  if (args.size() != 3) {
    std::cerr << "usage: allocate NUM_REFS SEED DISTANCE" << std::endl;
    return 2;
  }
  const uint32_t num_refs = std::stoi(args[0]);
  MappingParameters params;
  params.multi_mapping_allocation_seed = std::stoi(args[1]);
  params.multi_mapping_allocation_distance = std::stoi(args[2]);
  params.max_num_best_mappings = 8;
  std::vector<std::vector<MappingWithoutBarcode> > mappings;
  ReadRecords(std::cin, num_refs, mappings);
  uint64_t num_multi = 0;
  for (uint32_t r = 0; r < num_refs; ++r)
    for (size_t i = 0; i < mappings[r].size(); ++i) num_multi += mappings[r][i].mapq_ < 4;
  MappingProcessor<MappingWithoutBarcode> processor(params, /*min_unique_mapping_mapq=*/4);
  processor.SortOutputMappings(num_refs, mappings);
  processor.AllocateMultiMappings(num_refs, num_multi, params.multi_mapping_allocation_distance, mappings);
  for (uint32_t r = 0; r < mappings.size(); ++r) {
    for (size_t i = 0; i < mappings[r].size(); ++i) {
      std::cout << r << " ";
      Write(std::cout, mappings[r][i]);
      std::cout << "\n";
    }
  }
  return 0;
}

int ParseField(const std::string &name) {
  if (name == "total") return SUMMARY_METADATA_TOTAL;
  if (name == "dup") return SUMMARY_METADATA_DUP;
  if (name == "mapped") return SUMMARY_METADATA_MAPPED;
  if (name == "lowmapq") return SUMMARY_METADATA_LOWMAPQ;
  if (name == "cachehit") return SUMMARY_METADATA_CACHEHIT;
  if (name == "cardinality") return SUMMARY_METADATA_CARDINALITY;
  std::cerr << "unknown summary field " << name << std::endl;
  exit(2);
}

// summary FORMAT BARCODE_LENGTH WHITELIST(0|1) NOT_IN_WHITELIST(0|1) OUT CSV:
// a MappingWriter whose summary goes to CSV, driven by stdin commands
//   count BARCODE_KEY FIELD CHANGE   UpdateSummaryMetadata
//   nonwl FIELD CHANGE               UpdateSpeicalCategorySummaryMetadata(0, ...)
//   adjust                           AdjustSummaryPairedEndOverCount
//   output C0 C1 C2 C3 C4 SLOTS      OutputSummaryMetadata({C0..C4}, SLOTS)
CT_PROBE(summary) {
  if (args.size() != 6) {
    std::cerr << "usage: summary FORMAT BARCODE_LENGTH WHITELIST NOT_IN_WHITELIST OUT CSV" << std::endl;
    return 2;
  }
  MappingParameters params;
  params.mapping_output_format = ParseFormat(args[0]);
  if (std::stoi(args[2])) params.barcode_whitelist_file_path = "whitelist.txt";
  params.output_mappings_not_in_whitelist = std::stoi(args[3]) != 0;
  params.mapping_output_file_path = args[4];
  params.summary_metadata_file_path = args[5];
  MappingWriter<MappingWithBarcode> writer(params, std::stoi(args[1]), std::vector<int>());
  std::string cmd;
  while (std::cin >> cmd) {
    if (cmd == "count") {
      uint64_t key;
      std::string field;
      int change;
      std::cin >> key >> field >> change;
      writer.UpdateSummaryMetadata(key, ParseField(field), change);
    } else if (cmd == "nonwl") {
      std::string field;
      int change;
      std::cin >> field >> change;
      writer.UpdateSpeicalCategorySummaryMetadata(0, ParseField(field), change);
    } else if (cmd == "adjust") {
      writer.AdjustSummaryPairedEndOverCount();
    } else if (cmd == "output") {
      std::vector<double> coeffs(5);
      int slots;
      for (int i = 0; i < 5; ++i) std::cin >> coeffs[i];
      std::cin >> slots;
      writer.OutputSummaryMetadata(coeffs, slots != 0);
    } else {
      std::cerr << "unknown command " << cmd << std::endl;
      return 2;
    }
  }
  return 0;
}

}  // namespace

// MoveMappingsInBuffersToMappingContainer adds every thread's buffer to the
// container per reference, keeping what the container held, empties the
// buffers and returns how many mappings moved (2 threads x 2 references, 4
// moved). The order within a reference is not checked: chromap sorts every
// reference before it writes or deduplicates.
CT_TEST(move_buffers_to_container) {
  MappingParameters params;
  MappingProcessor<MappingWithoutBarcode> processor(params, 4);
  std::vector<std::vector<MappingWithoutBarcode> > container(2);
  container[0].push_back(MappingWithoutBarcode(0, 5, 50, 60, 1, 1, 1));
  std::vector<std::vector<std::vector<MappingWithoutBarcode> > > buffers(
      2, std::vector<std::vector<MappingWithoutBarcode> >(2));
  buffers[0][0].push_back(MappingWithoutBarcode(1, 10, 50, 60, 1, 1, 1));
  buffers[0][1].push_back(MappingWithoutBarcode(2, 20, 50, 60, 1, 1, 1));
  buffers[1][0].push_back(MappingWithoutBarcode(3, 1, 50, 60, 1, 1, 1));
  buffers[1][1].push_back(MappingWithoutBarcode(4, 30, 50, 60, 0, 1, 1));
  const uint32_t moved = processor.MoveMappingsInBuffersToMappingContainer(2, buffers, container);
  CT_CHECK_EQ(moved, 4u);
  CT_CHECK_EQ(container[0].size(), (size_t)3);
  CT_CHECK_EQ(container[1].size(), (size_t)2);
  std::vector<uint32_t> ids0, ids1;
  for (size_t i = 0; i < container[0].size(); ++i) ids0.push_back(container[0][i].read_id_);
  for (size_t i = 0; i < container[1].size(); ++i) ids1.push_back(container[1][i].read_id_);
  std::sort(ids0.begin(), ids0.end());
  std::sort(ids1.begin(), ids1.end());
  CT_CHECK(ids0 == std::vector<uint32_t>({0, 1, 3}));
  CT_CHECK(ids1 == std::vector<uint32_t>({2, 4}));
  for (int t = 0; t < 2; ++t)
    for (int r = 0; r < 2; ++r) CT_CHECK(buffers[t][r].empty());
}

// OutputMappingStatistics reports the unique (is_unique 1) and multi
// mappings over every reference and their total, in that order: here 2, 1
// and 3. Only the numbers are checked, not the wording around them.
CT_TEST(mapping_statistics_counts) {
  MappingParameters params;
  MappingProcessor<MappingWithoutBarcode> processor(params, 4);
  std::vector<std::vector<MappingWithoutBarcode> > mappings(3);
  mappings[0].push_back(MappingWithoutBarcode(0, 5, 50, 60, 1, 1, 1));
  mappings[2].push_back(MappingWithoutBarcode(1, 5, 50, 0, 1, 0, 1));
  mappings[2].push_back(MappingWithoutBarcode(2, 9, 50, 60, 0, 1, 1));
  std::stringstream captured;
  std::streambuf *old = std::cerr.rdbuf(captured.rdbuf());
  processor.OutputMappingStatistics(3, mappings);
  std::cerr.rdbuf(old);
  std::vector<long> numbers;
  const std::string text = captured.str();
  for (size_t i = 0; i < text.size();) {
    if (text[i] >= '0' && text[i] <= '9') {
      size_t j = i;
      while (j < text.size() && text[j] >= '0' && text[j] <= '9') ++j;
      numbers.push_back(std::stol(text.substr(i, j - i)));
      i = j;
    } else {
      ++i;
    }
  }
  CT_CHECK(numbers == std::vector<long>({2, 1, 3}));
}

CT_MAIN
