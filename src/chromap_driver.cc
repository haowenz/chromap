#include "chromap_driver.h"

#include <glob.h>

#include <cassert>
#include <iomanip>
#include <string>
#include <vector>

#include "chromap.h"
#include "cxxopts.hpp"

namespace chromap {
namespace {

void AddIndexingOptions(cxxopts::Options &options) {
  options.add_options("Indexing")("i,build-index", "Build index")(
      "min-frag-length",
      "Min fragment length for choosing k and w automatically [30]",
      cxxopts::value<int>(),
      "INT")("k,kmer", "Kmer length [17]", cxxopts::value<int>(), "INT")(
      "w,window", "Window size [7]", cxxopts::value<int>(), "INT");
}

void AddMappingOptions(cxxopts::Options &options) {
  options.set_width(120).add_options("Mapping")(
      "preset",
      "Preset parameters for mapping reads (always applied before other "
      "options) []\natac: mapping ATAC-seq/scATAC-seq reads\nchip: mapping "
      "ChIP-seq reads\nhic: mapping Hi-C reads",
      cxxopts::value<std::string>(),
      "STR")("split-alignment", "Allow split alignments")(
      "e,error-threshold", "Max # errors allowed to map a read [8]",
      cxxopts::value<int>(), "INT")
      //("A,match-score", "Match score [1]", cxxopts::value<int>(), "INT")
      //("B,mismatch-penalty", "Mismatch penalty [4]", cxxopts::value<int>(),
      //"INT")
      //("O,gap-open-penalties", "Gap open penalty [6,6]",
      // cxxopts::value<std::vector<int>>(), "INT[,INT]")
      //("E,gap-extension-penalties", "Gap extension penalty [1,1]",
      // cxxopts::value<std::vector<int>>(), "INT[,INT]")
      ("s,min-num-seeds", "Min # seeds to try to map a read [2]",
       cxxopts::value<int>(),
       "INT")("f,max-seed-frequencies",
              "Max seed frequencies for a seed to be selected [500,1000]",
              cxxopts::value<std::vector<int>>(), "INT[,INT]")
      //("n,max-num-best-mappings", "Only report n best mappings [1]",
      // cxxopts::value<int>(), "INT")
      ("l,max-insert-size",
       "Max insert size, only for paired-end read mapping [1000]",
       cxxopts::value<int>(),
       "INT")("q,MAPQ-threshold",
              "Min MAPQ in range [0, 60] for mappings to be output [30]",
              cxxopts::value<uint8_t>(),
              "INT")("min-read-length", "Min read length [30]",
                     cxxopts::value<int>(), "INT")
      //("multi-mapping-allocation-distance", "Uni-mappings within this distance
      // from any end of multi-mappings are used for allocation [0]",
      // cxxopts::value<int>(), "INT")
      //("multi-mapping-allocation-seed", "Seed for random number generator in
      // multi-mapping allocation [11]", cxxopts::value<int>(), "INT")
      //("drop-repetitive-reads", "Drop reads with too many best mappings
      //[500000]", cxxopts::value<int>(), "INT")
      ("trim-adapters", "Try to trim adapters on 3'")("remove-pcr-duplicates",
                                                      "Remove PCR duplicates")(
          "remove-pcr-duplicates-at-bulk-level",
          "Remove PCR duplicates at bulk level for single cell data")(
          "remove-pcr-duplicates-at-cell-level",
          "Remove PCR duplicates at cell level for single cell data")
      //("allocate-multi-mappings", "Allocate multi-mappings")
      ("Tn5-shift", "Perform Tn5 shift")("low-mem", "Use low memory mode")(
          "bc-error-threshold",
          "Max Hamming distance allowed to correct a barcode [1]",
          cxxopts::value<int>(),
          "INT")("bc-probability-threshold",
                 "Min probability to correct a barcode [0.9]",
                 cxxopts::value<double>(),
                 "FLT")("t,num-threads", "# threads for mapping [1]",
                        cxxopts::value<int>(), "INT")
      ("cache-size", "number of cache entries [4000003]", cxxopts::value<int>(), "INT")
      ("cache-update-param", "value used to control number of reads sampled [0.01]", cxxopts::value<double>(), "FLT")
      ("use-all-reads", "use all reads for cache")
      ("debug-cache", "verbose output for debugging cache used in chromap");
}

void AddInputOptions(cxxopts::Options &options) {
  options.add_options("Input")("r,ref", "Reference file",
                               cxxopts::value<std::string>(), "FILE")(
      "x,index", "Index file", cxxopts::value<std::string>(), "FILE")(
      "1,read1", "Single-end read files or paired-end read files 1",
      cxxopts::value<std::vector<std::string>>(),
      "FILE")("2,read2", "Paired-end read files 2",
              cxxopts::value<std::vector<std::string>>(),
              "FILE")("b,barcode", "Cell barcode files",
                      cxxopts::value<std::vector<std::string>>(), "FILE")(
      "barcode-whitelist", "Cell barcode whitelist file",
      cxxopts::value<std::string>(),
      "FILE")("read-format",
              "Format for read files and barcode files  [\"r1:0:-1,bc:0:-1\" "
              "as 10x Genomics single-end format]",
              cxxopts::value<std::string>(), "STR");
}

void AddOutputOptions(cxxopts::Options &options) {
  options.add_options("Output")("o,output", "Output file",
                                cxxopts::value<std::string>(), "FILE")
      //("p,matrix-output-prefix", "Prefix of matrix output files",
      // cxxopts::value<std::string>(), "FILE")
      ("output-mappings-not-in-whitelist",
       "Output mappings with barcode not in the whitelist")(
          "chr-order",
          "Custom chromosome order file. If not specified, the order of "
          "reference sequences will be used",
          cxxopts::value<std::string>(),
          "FILE")("BED", "Output mappings in BED/BEDPE format")(
          "TagAlign", "Output mappings in TagAlign/PairedTagAlign format")(
          "SAM", "Output mappings in SAM format")(
          "pairs",
          "Output mappings in pairs format (defined by 4DN for HiC data)")(
          "pairs-natural-chr-order",
          "Custom chromosome order file for pairs flipping. If not specified, "
          "the custom chromosome order will be used",
          cxxopts::value<std::string>(),
          "FILE")("barcode-translate",
                  "Convert barcode to the specified sequences during output",
                  cxxopts::value<std::string>(), "FILE")(
          "summary",
          "Summarize the mapping statistics at bulk or barcode level",
          cxxopts::value<std::string>(), "FILE");
  //("PAF", "Output mappings in PAF format (only for test)");
}

void AddDevelopmentOptions(cxxopts::Options &options) {
  options.add_options("Development options")("A,match-score", "Match score [1]",
                                             cxxopts::value<int>(), "INT")(
      "B,mismatch-penalty", "Mismatch penalty [4]", cxxopts::value<int>(),
      "INT")("O,gap-open-penalties", "Gap open penalty [6,6]",
             cxxopts::value<std::vector<int>>(), "INT[,INT]")(
      "E,gap-extension-penalties", "Gap extension penalty [1,1]",
      cxxopts::value<std::vector<int>>(),
      "INT[,INT]")("n,max-num-best-mappings", "Only report n best mappings [1]",
                   cxxopts::value<int>(),
                   "INT")("multi-mapping-allocation-distance",
                          "Uni-mappings within this distance from any end of "
                          "multi-mappings are used for allocation [0]",
                          cxxopts::value<int>(), "INT")(
      "multi-mapping-allocation-seed",
      "Seed for random number generator in multi-mapping allocation [11]",
      cxxopts::value<int>(), "INT")(
      "drop-repetitive-reads",
      "Drop reads with too many best mappings [500000]", cxxopts::value<int>(),
      "INT")("allocate-multi-mappings", "Allocate multi-mappings")(
      "PAF", "Output mappings in PAF format (only for test)")(
      "skip-barcode-check",
      "Do not check whether too few barcodes are in the whitelist");
}

void AddPeakOptions(cxxopts::Options &options) {
  options.add_options("Peak")("cell-by-bin", "Generate cell-by-bin matrix")(
      "bin-size", "Bin size to generate cell-by-bin matrix [5000]",
      cxxopts::value<int>(),
      "INT")("depth-cutoff", "Depth cutoff for peak calling [3]",
             cxxopts::value<int>(),
             "INT")("peak-min-length", "Min length of peaks to report [30]",
                    cxxopts::value<int>(), "INT")(
      "peak-merge-max-length", "Peaks within this length will be merged [30]",
      cxxopts::value<int>(), "INT");
}

// Return all file paths that match the input pattern.
std::vector<std::string> GetMatchedFilePaths(const std::string &pattern) {
  glob_t glob_result;
  memset(&glob_result, 0, sizeof(glob_result));

  const int return_value =
      glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);

  if (return_value != 0) {
    globfree(&glob_result);
    chromap::ExitWithMessage("glob() failed with return value " +
                             std::to_string(return_value) + "\n");
  }

  std::vector<std::string> matched_file_paths;
  matched_file_paths.reserve(glob_result.gl_pathc);
  for (size_t i = 0; i < glob_result.gl_pathc; ++i) {
    matched_file_paths.push_back(std::string(glob_result.gl_pathv[i]));
    std::cerr << matched_file_paths.back() << "\n";
  }
  globfree(&glob_result);

  return matched_file_paths;
}

// Return all file paths that match the input patterns.
std::vector<std::string> GetMatchedFilePaths(
    const std::vector<std::string> &patterns) {
  std::vector<std::string> all_matched_file_paths;
  for (const auto &pattern : patterns) {
    std::vector<std::string> matched_file_paths = GetMatchedFilePaths(pattern);
    all_matched_file_paths.reserve(all_matched_file_paths.size() +
                                   matched_file_paths.size());
    all_matched_file_paths.insert(
        std::end(all_matched_file_paths),
        std::make_move_iterator(std::begin(matched_file_paths)),
        std::make_move_iterator(std::end(matched_file_paths)));
  }
  return all_matched_file_paths;
}

}  // namespace

void ChromapDriver::ParseArgsAndRun(int argc, char *argv[]) {
  cxxopts::Options options(
      "chromap", "Fast alignment and preprocessing of chromatin profiles");

  options.add_options()("v,version", "Print version")("h,help", "Print help");

  AddIndexingOptions(options);
  AddMappingOptions(options);

  // We don't support peak options for now.
  // AddPeakOptions(options);

  AddInputOptions(options);
  AddOutputOptions(options);

  AddDevelopmentOptions(options);

  auto result = options.parse(argc, argv);
  if (result.count("h")) {
    std::cerr << options.help(
        {"", "Indexing", "Mapping", "Peak", "Input", "Output"});
    return;
  }
  if (result.count("v")) {
    std::cerr << CHROMAP_VERSION << "\n";
    return;
  }
  // Parameters and their default
  IndexParameters index_parameters;
  MappingParameters mapping_parameters;

  if (result.count("preset")) {
    std::string read_type = result["preset"].as<std::string>();
    if (read_type == "atac") {
      std::cerr << "Preset parameters for ATAC-seq/scATAC-seq are used.\n";
      mapping_parameters.max_insert_size = 2000;
      mapping_parameters.trim_adapters = true;
      mapping_parameters.remove_pcr_duplicates = true;
      mapping_parameters.remove_pcr_duplicates_at_bulk_level = false;
      mapping_parameters.Tn5_shift = true;
      mapping_parameters.mapping_output_format = MAPPINGFORMAT_BED;
      mapping_parameters.low_memory_mode = true;
    } else if (read_type == "chip") {
      std::cerr << "Preset parameters for ChIP-seq are used.\n";
      mapping_parameters.max_insert_size = 2000;
      mapping_parameters.remove_pcr_duplicates = true;
      mapping_parameters.low_memory_mode = true;
      mapping_parameters.mapping_output_format = MAPPINGFORMAT_BED;
    } else if (read_type == "hic") {
      std::cerr << "Preset parameters for Hi-C are used.\n";
      mapping_parameters.error_threshold = 4;
      mapping_parameters.mapq_threshold = 1;
      mapping_parameters.split_alignment = true;
      mapping_parameters.low_memory_mode = true;
      mapping_parameters.mapping_output_format = MAPPINGFORMAT_PAIRS;
    } else {
      chromap::ExitWithMessage("Unrecognized preset parameters " + read_type +
                               "\n");
    }
  }
  // Optional parameters
  if (result.count("min-frag-length")) {
    int min_fragment_length = result["min-frag-length"].as<int>();
    if (min_fragment_length <= 60) {
      index_parameters.kmer_size = 17;
      index_parameters.window_size = 7;
    } else if (min_fragment_length <= 80) {
      index_parameters.kmer_size = 19;
      index_parameters.window_size = 10;
    } else {
      index_parameters.kmer_size = 23;
      index_parameters.window_size = 11;
    }
  }
  if (result.count("k")) {
    index_parameters.kmer_size = result["kmer"].as<int>();
  }
  if (result.count("w")) {
    index_parameters.window_size = result["window"].as<int>();
  }
  if (result.count("e")) {
    mapping_parameters.error_threshold = result["error-threshold"].as<int>();
  }
  if (result.count("A")) {
    mapping_parameters.match_score = result["match-score"].as<int>();
  }
  if (result.count("B")) {
    mapping_parameters.mismatch_penalty = result["mismatch-penalty"].as<int>();
  }
  if (result.count("O")) {
    mapping_parameters.gap_open_penalties =
        result["gap-open-penalties"].as<std::vector<int>>();
  }
  if (result.count("E")) {
    mapping_parameters.gap_extension_penalties =
        result["gap-extension-penalties"].as<std::vector<int>>();
  }
  if (result.count("s")) {
    mapping_parameters.min_num_seeds_required_for_mapping =
        result["min-num-seeds"].as<int>();
  }
  if (result.count("f")) {
    mapping_parameters.max_seed_frequencies =
        result["max-seed-frequencies"].as<std::vector<int>>();
  }
  if (result.count("n")) {
    mapping_parameters.max_num_best_mappings =
        result["max-num-best-mappings"].as<int>();
  }
  if (result.count("l")) {
    mapping_parameters.max_insert_size = result["max-insert-size"].as<int>();
  }
  if (result.count("q")) {
    mapping_parameters.mapq_threshold = result["MAPQ-threshold"].as<uint8_t>();
  }
  if (result.count("t")) {
    mapping_parameters.num_threads = result["num-threads"].as<int>();
  }


  // check cache-related parameters
  if (result.count("cache-update-param")) {
    mapping_parameters.cache_update_param = result["cache-update-param"].as<double>();
    if (mapping_parameters.cache_update_param < 0.0 || mapping_parameters.cache_update_param > 1.0){
      chromap::ExitWithMessage("cache update param is not approriate, must be in this range (0, 1]");
    }
  } 
  if (result.count("cache-size")) {
    mapping_parameters.cache_size = result["cache-size"].as<int>();
    if (mapping_parameters.cache_size < 2000000 || mapping_parameters.cache_size > 15000000) {
        chromap::ExitWithMessage("cache size is not in appropriate range\n");
    }
  }
  if (result.count("use-all-reads")) {
    mapping_parameters.use_all_reads = true;
  }
  if (result.count("debug-cache")) {
    mapping_parameters.debug_cache = true;
  }



  if (result.count("min-read-length")) {
    mapping_parameters.min_read_length = result["min-read-length"].as<int>();
  }
  if (result.count("bc-error-threshold")) {
    mapping_parameters.barcode_correction_error_threshold =
        result["bc-error-threshold"].as<int>();
  }
  if (result.count("bc-probability-threshold")) {
    mapping_parameters.barcode_correction_probability_threshold =
        result["bc-probability-threshold"].as<double>();
  }
  if (result.count("multi-mapping-allocation-distance")) {
    mapping_parameters.multi_mapping_allocation_distance =
        result["multi-mapping-allocation-distance"].as<int>();
  }
  if (result.count("multi-mapping-allocation-seed")) {
    mapping_parameters.multi_mapping_allocation_seed =
        result["multi-mapping-allocation-seed"].as<int>();
  }
  if (result.count("drop-repetitive-reads")) {
    mapping_parameters.drop_repetitive_reads =
        result["drop-repetitive-reads"].as<int>();
  }
  if (result.count("trim-adapters")) {
    mapping_parameters.trim_adapters = true;
  }
  if (result.count("remove-pcr-duplicates")) {
    mapping_parameters.remove_pcr_duplicates = true;
  }
  if (result.count("remove-pcr-duplicates-at-bulk-level")) {
    mapping_parameters.remove_pcr_duplicates_at_bulk_level = true;
  }
  if (result.count("remove-pcr-duplicates-at-cell-level")) {
    mapping_parameters.remove_pcr_duplicates_at_bulk_level = false;
  }
  if (result.count("allocate-multi-mappings")) {
    mapping_parameters.allocate_multi_mappings = true;
    mapping_parameters.only_output_unique_mappings = false;
  }
  if (result.count("Tn5-shift")) {
    mapping_parameters.Tn5_shift = true;
  }
  if (result.count("split-alignment")) {
    mapping_parameters.split_alignment = true;
  }
  if (result.count("output-mappings-not-in-whitelist")) {
    mapping_parameters.output_mappings_not_in_whitelist = true;
  }
  if (result.count("BED")) {
    mapping_parameters.mapping_output_format = MAPPINGFORMAT_BED;
  }
  if (result.count("TagAlign")) {
    mapping_parameters.mapping_output_format = MAPPINGFORMAT_TAGALIGN;
  }
  if (result.count("PAF")) {
    mapping_parameters.mapping_output_format = MAPPINGFORMAT_PAF;
  }
  if (result.count("pairs")) {
    mapping_parameters.mapping_output_format = MAPPINGFORMAT_PAIRS;
  }
  if (result.count("SAM")) {
    mapping_parameters.mapping_output_format = MAPPINGFORMAT_SAM;
  }
  if (result.count("low-mem")) {
    mapping_parameters.low_memory_mode = true;
  }
  if (result.count("cell-by-bin")) {
    mapping_parameters.cell_by_bin = true;
  }
  if (result.count("bin-size")) {
    mapping_parameters.bin_size = result["bin-size"].as<int>();
  }
  if (result.count("depth-cutoff")) {
    mapping_parameters.depth_cutoff_to_call_peak =
        result["depth-cutoff"].as<uint16_t>();
  }
  if (result.count("peak-min-length")) {
    mapping_parameters.peak_min_length = result["peak-min-length"].as<int>();
  }
  if (result.count("peak-merge-max-length")) {
    mapping_parameters.peak_merge_max_length =
        result["peak-merge-max-length"].as<int>();
  }

  std::cerr << std::setprecision(2) << std::fixed;
  if (result.count("i")) {
    if (result.count("r")) {
      index_parameters.reference_file_path = result["ref"].as<std::string>();
    } else {
      chromap::ExitWithMessage("No reference specified!");
    }
    if (result.count("o")) {
      index_parameters.index_output_file_path =
          result["output"].as<std::string>();
    } else {
      chromap::ExitWithMessage("No output file specified!");
    }
    std::cerr << "Build index for the reference.\n";
    std::cerr << "Kmer length: " << index_parameters.kmer_size
              << ", window size: " << index_parameters.window_size << "\n";
    std::cerr << "Reference file: " << index_parameters.reference_file_path
              << "\n";
    std::cerr << "Output file: " << index_parameters.index_output_file_path
              << "\n";
    chromap::Chromap chromap_for_indexing(index_parameters);
    chromap_for_indexing.ConstructIndex();
  } else if (result.count("1")) {
    std::cerr << "Start to map reads.\n";
    if (result.count("r")) {
      mapping_parameters.reference_file_path = result["ref"].as<std::string>();
    } else {
      chromap::ExitWithMessage("No reference specified!");
    }
    if (result.count("o")) {
      mapping_parameters.mapping_output_file_path =
          result["output"].as<std::string>();
    } else {
      chromap::ExitWithMessage("No output file specified!");
    }
    if (result.count("x")) {
      mapping_parameters.index_file_path = result["index"].as<std::string>();
    } else {
      chromap::ExitWithMessage("No index file specified!");
    }
    if (result.count("1")) {
      mapping_parameters.read_file1_paths =
          GetMatchedFilePaths(result["read1"].as<std::vector<std::string>>());
    } else {
      chromap::ExitWithMessage("No read file specified!");
    }
    if (result.count("2")) {
      mapping_parameters.read_file2_paths =
          GetMatchedFilePaths(result["read2"].as<std::vector<std::string>>());
    }

    if (result.count("b")) {
      mapping_parameters.is_bulk_data = false;
      mapping_parameters.barcode_file_paths =
          GetMatchedFilePaths(result["barcode"].as<std::vector<std::string>>());
      if (result.count("barcode-whitelist") == 0) {
        std::cerr << "WARNING: there are input barcode files but a barcode "
                     "whitelist file is missing!\n";
      }
    }

    if (result.count("barcode-whitelist")) {
      if (mapping_parameters.is_bulk_data) {
        chromap::ExitWithMessage(
            "No barcode file specified but the barcode whitelist file is "
            "given!");
      }
      mapping_parameters.barcode_whitelist_file_path =
          result["barcode-whitelist"].as<std::string>();
    }

    if (result.count("p")) {
      mapping_parameters.matrix_output_prefix =
          result["matrix-output-prefix"].as<std::string>();
      if (mapping_parameters.is_bulk_data) {
        chromap::ExitWithMessage(
            "No barcode file specified but asked to output matrix files!");
      }
    }
    if (result.count("read-format")) {
      mapping_parameters.read_format = result["read-format"].as<std::string>();
    }

    if (result.count("chr-order")) {
      mapping_parameters.custom_rid_order_file_path =
          result["chr-order"].as<std::string>();
    }

    if (result.count("pairs-natural-chr-order")) {
      mapping_parameters.pairs_flipping_custom_rid_order_file_path =
          result["pairs-natural-chr-order"].as<std::string>();
    }

    if (result.count("barcode-translate")) {
      mapping_parameters.barcode_translate_table_file_path =
          result["barcode-translate"].as<std::string>();
    }

    if (result.count("summary")) {
      mapping_parameters.summary_metadata_file_path =
          result["summary"].as<std::string>();
    }

    if (result.count("skip-barcode-check")) {
      mapping_parameters.skip_barcode_check = true;
    }

    // std::cerr << "Parameters: error threshold: " << error_threshold << ",
    // match score: " << match_score << ", mismatch_penalty: " <<
    // mismatch_penalty << ", gap open penalties for deletions and insertions: "
    // << gap_open_penalties[0] << "," << gap_open_penalties[1] << ", gap
    // extension penalties for deletions and insertions: " <<
    // gap_extension_penalties[0] << "," << gap_extension_penalties[1] << ",
    // min-num-seeds: " << min_num_seeds_required_for_mapping << ",
    // max-seed-frequency: " << max_seed_frequencies[0] << "," <<
    // max_seed_frequencies[1] << ", max-num-best-mappings: " <<
    // max_num_best_mappings << ", max-insert-size: " << max_insert_size << ",
    // MAPQ-threshold: " << (int)mapq_threshold << ", min-read-length: " <<
    // min_read_length << ", multi-mapping-allocation-distance: " <<
    // multi_mapping_allocation_distance << ", multi-mapping-allocation-seed: "
    // << multi_mapping_allocation_seed << ", drop-repetitive-reads: " <<
    // drop_repetitive_reads << "\n";
    std::cerr << "Parameters: error threshold: "
              << mapping_parameters.error_threshold << ", min-num-seeds: "
              << mapping_parameters.min_num_seeds_required_for_mapping
              << ", max-seed-frequency: "
              << mapping_parameters.max_seed_frequencies[0] << ","
              << mapping_parameters.max_seed_frequencies[1]
              << ", max-num-best-mappings: "
              << mapping_parameters.max_num_best_mappings
              << ", max-insert-size: " << mapping_parameters.max_insert_size
              << ", MAPQ-threshold: " << (int)mapping_parameters.mapq_threshold
              << ", min-read-length: " << mapping_parameters.min_read_length
              << ", bc-error-threshold: "
              << mapping_parameters.barcode_correction_error_threshold
              << ", bc-probability-threshold: "
              << mapping_parameters.barcode_correction_probability_threshold
              << "\n";
    std::cerr << "Number of threads: " << mapping_parameters.num_threads
              << "\n";
    if (mapping_parameters.is_bulk_data) {
      std::cerr << "Analyze bulk data.\n";
    } else {
      std::cerr << "Analyze single-cell data.\n";
    }
    if (mapping_parameters.trim_adapters) {
      std::cerr << "Will try to remove adapters on 3'.\n";
    } else {
      std::cerr << "Won't try to remove adapters on 3'.\n";
    }
    if (mapping_parameters.remove_pcr_duplicates) {
      std::cerr << "Will remove PCR duplicates after mapping.\n";
    } else {
      std::cerr << "Won't remove PCR duplicates after mapping.\n";
    }
    if (mapping_parameters.remove_pcr_duplicates_at_bulk_level) {
      std::cerr << "Will remove PCR duplicates at bulk level.\n";
    } else {
      std::cerr << "Will remove PCR duplicates at cell level.\n";
    }
    if (mapping_parameters.allocate_multi_mappings) {
      std::cerr << "Will allocate multi-mappings after mapping.\n";
    } else {
      std::cerr << "Won't allocate multi-mappings after mapping.\n";
    }
    if (mapping_parameters.only_output_unique_mappings) {
      std::cerr << "Only output unique mappings after mapping.\n";
    }
    if (!mapping_parameters.output_mappings_not_in_whitelist) {
      std::cerr << "Only output mappings of which barcodes are in whitelist.\n";
    } else {
      std::cerr << "No filtering of mappings based on whether their barcodes "
                   "are in whitelist.\n";
    }
    // if (allocate_multi_mappings && only_output_unique_mappings) {
    //  std::cerr << "WARNING: you want to output unique mappings only but you
    //  ask to allocate multi-mappings! In this case, it won't allocate
    //  multi-mappings and will only output unique mappings.\n";
    //  allocate_multi_mappings = false;
    //}
    if (mapping_parameters.max_num_best_mappings >
        mapping_parameters.drop_repetitive_reads) {
      std::cerr << "WARNING: you want to drop mapped reads with more than "
                << mapping_parameters.drop_repetitive_reads
                << " mappings. But you want to output top "
                << mapping_parameters.max_num_best_mappings
                << " best mappings. In this case, only reads with <="
                << mapping_parameters.drop_repetitive_reads
                << " best mappings will be output.\n";
      mapping_parameters.max_num_best_mappings =
          mapping_parameters.drop_repetitive_reads;
    }
    if (mapping_parameters.Tn5_shift) {
      std::cerr << "Perform Tn5 shift.\n";
    }
    if (mapping_parameters.split_alignment) {
      std::cerr << "Allow split alignment.\n";
    }

    switch (mapping_parameters.mapping_output_format) {
      case MAPPINGFORMAT_BED:
        std::cerr << "Output mappings in BED/BEDPE format.\n";
        break;
      case MAPPINGFORMAT_TAGALIGN:
        std::cerr << "Output mappings in TagAlign/PairedTagAlign format.\n";
        break;
      case MAPPINGFORMAT_PAF:
        std::cerr << "Output mappings in PAF format.\n";
        break;
      case MAPPINGFORMAT_SAM:
        std::cerr << "Output mappings in SAM format.\n";
        break;
      case MAPPINGFORMAT_PAIRS:
        std::cerr << "Output mappings in pairs format.\n";
        break;
      default:
        chromap::ExitWithMessage("Unknown mapping output format!");
        break;
    }

    std::cerr << "Reference file: " << mapping_parameters.reference_file_path
              << "\n";
    std::cerr << "Index file: " << mapping_parameters.index_file_path << "\n";
    for (size_t i = 0; i < mapping_parameters.read_file1_paths.size(); ++i) {
      std::cerr << i + 1
                << "th read 1 file: " << mapping_parameters.read_file1_paths[i]
                << "\n";
    }
    if (result.count("2") != 0) {
      for (size_t i = 0; i < mapping_parameters.read_file2_paths.size(); ++i) {
        std::cerr << i + 1 << "th read 2 file: "
                  << mapping_parameters.read_file2_paths[i] << "\n";
      }
    }
    if (result.count("b") != 0) {
      for (size_t i = 0; i < mapping_parameters.barcode_file_paths.size();
           ++i) {
        std::cerr << i + 1 << "th cell barcode file: "
                  << mapping_parameters.barcode_file_paths[i] << "\n";
      }
    }
    if (result.count("barcode-whitelist") != 0) {
      std::cerr << "Cell barcode whitelist file: "
                << mapping_parameters.barcode_whitelist_file_path << "\n";
    }
    std::cerr << "Output file: " << mapping_parameters.mapping_output_file_path
              << "\n";
    if (result.count("matrix-output-prefix") != 0) {
      std::cerr << "Matrix output prefix: "
                << mapping_parameters.matrix_output_prefix << "\n";
    }

    chromap::Chromap chromap_for_mapping(mapping_parameters);

    if (result.count("2") == 0) {
      // Single-end reads.
      switch (mapping_parameters.mapping_output_format) {
        case MAPPINGFORMAT_PAF: {
          chromap_for_mapping.MapSingleEndReads<chromap::PAFMapping>();
          break;
        }
        case MAPPINGFORMAT_SAM: {
          chromap_for_mapping.MapSingleEndReads<chromap::SAMMapping>();
          break;
        }
        case MAPPINGFORMAT_PAIRS:
          chromap::ExitWithMessage("No support for single-end HiC yet!");
          break;
        case MAPPINGFORMAT_BED:
        case MAPPINGFORMAT_TAGALIGN:
          if (result.count("b") != 0) {
            chromap_for_mapping
                .MapSingleEndReads<chromap::MappingWithBarcode>();
          } else {
            chromap_for_mapping
                .MapSingleEndReads<chromap::MappingWithoutBarcode>();
          }
          break;
        default:
          chromap::ExitWithMessage("Unknown mapping output format!");
          break;
      }
    } else {
      // Paired-end reads.
      switch (mapping_parameters.mapping_output_format) {
        case MAPPINGFORMAT_PAF: {
          chromap_for_mapping.MapPairedEndReads<chromap::PairedPAFMapping>();
          break;
        }
        case MAPPINGFORMAT_SAM: {
          chromap_for_mapping.MapPairedEndReads<chromap::SAMMapping>();
          break;
        }
        case MAPPINGFORMAT_PAIRS: {
          chromap_for_mapping.MapPairedEndReads<chromap::PairsMapping>();
          break;
        }
        case MAPPINGFORMAT_BED:
        case MAPPINGFORMAT_TAGALIGN:
          if (result.count("b") != 0) {
            chromap_for_mapping
                .MapPairedEndReads<chromap::PairedEndMappingWithBarcode>();
          } else {
            chromap_for_mapping
                .MapPairedEndReads<chromap::PairedEndMappingWithoutBarcode>();
          }
          break;
        default:
          chromap::ExitWithMessage("Unknown mapping output format!");
          break;
      }
    }
  } else {
    std::cerr << options.help(
        {"", "Indexing", "Mapping", "Peak", "Input", "Output"});
  }
}

}  // namespace chromap

int main(int argc, char *argv[]) {
  chromap::ChromapDriver chromap_driver;
  chromap_driver.ParseArgsAndRun(argc, argv);
  return 0;
}
