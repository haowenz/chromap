#include "chromap.h"

#include <assert.h>
#include <math.h>

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <random>
#include <sstream>
#include <system_error>
#include <thread>
#include <type_traits>
#include <unordered_map>

namespace chromap {

namespace {

// Loads sequences into slots 0, 1, ... of 'batch' until its file ends or
// 'max_num_sequences' are loaded, and returns how many it loaded.
uint32_t LoadUpToOneBatch(SequenceBatch &batch, uint32_t max_num_sequences) {
  uint32_t i = 0;
  for (i = 0; i < max_num_sequences; ++i) {
    if (batch.LoadOneSequenceAndSaveAt(i) == true) {  // true: no more read
      break;
    }
  }
  return i;
}

}  // namespace

// Reads the barcode files in order, in batches of up to batch_size records
// each, as SequenceBatch::LoadBatch would after InitializeLoading, but keeps
// only the whitelist key of each barcode without an N, in record order: the
// records go through one slot, so no batch of slots is built or filled.
// 'file_paths' must outlive it; each is read when its file is opened.
class BarcodeKeyBatchLoader {
 public:
  BarcodeKeyBatchLoader(const std::vector<std::string> &file_paths,
                        uint32_t batch_size,
                        const SequenceEffectiveRange &effective_range)
      : file_paths_(file_paths),
        batch_size_(batch_size),
        barcode_(/*max_num_sequences=*/1, effective_range) {}

  ~BarcodeKeyBatchLoader() {
    if (is_open_) {
      barcode_.FinalizeLoading();
    }
  }

  // Loads the next batch of the file being read, opening the next file first
  // when the last one ended, into 'keys', and writes the log line LoadBatch
  // writes. Returns the number of barcodes it read, N or not, 0 at the end of
  // a file. Exits (ExitWithMessage) where InitializeLoading and LoadBatch
  // would have exited, on a missing or truncated file.
  uint32_t LoadNextBatch(std::vector<uint64_t> &keys, std::ostream &log) {
    if (!is_open_) {
      barcode_.InitializeLoading(file_paths_[next_file_index_]);
      ++next_file_index_;
      is_open_ = true;
    }
    const double real_start_time = GetRealTime();
    keys.clear();
    uint32_t num_loaded_barcodes = 0;
    for (; num_loaded_barcodes < batch_size_; ++num_loaded_barcodes) {
      if (barcode_.LoadOneSequenceAndSaveAt(0)) {
        break;
      }
      std::vector<int> N_pos;  // position of Ns
      barcode_.GetSequenceNsAt(0, /*little_endian=*/true, N_pos);
      if (N_pos.size() > 0) continue;
      keys.push_back(barcode_.GenerateSeedFromSequenceAt(
          0, 0, barcode_.GetSequenceLengthAt(0)));
    }
    // As SequenceBatch::LoadBatch logs it.
    if (num_loaded_barcodes != 0) {
      log << "Loaded sequence batch successfully in "
          << GetRealTime() - real_start_time << "s, ";
      log << "number of sequences: " << num_loaded_barcodes << ".\n";
    } else {
      log << "No more sequences.\n";
      barcode_.FinalizeLoading();
      is_open_ = false;
    }
    return num_loaded_barcodes;
  }

 private:
  const std::vector<std::string> &file_paths_;
  const uint32_t batch_size_;
  SequenceBatch barcode_;
  size_t next_file_index_ = 0;
  bool is_open_ = false;
};

// Runs a BarcodeKeyBatchLoader over the first num_files of 'file_paths' on a
// thread of its own, up to two batches ahead of ComputeBarcodeAbundance, which
// takes them in order with TakeNextBatch, with their log lines and any exit
// the loading defers (a missing or truncated file) at the batch where the
// loader would have exited, so that the log and the first error are
// unchanged. Batches loaded ahead and never taken are dropped with their log
// lines and errors.
class BarcodeBatchReadAhead {
 public:
  BarcodeBatchReadAhead(const std::vector<std::string> &file_paths,
                        size_t num_files, uint32_t batch_size,
                        const SequenceEffectiveRange &effective_range,
                        const std::ostream &log_format)
      : num_files_(num_files),
        loader_(file_paths, batch_size, effective_range),
        log_flags_(log_format.flags()),
        log_precision_(log_format.precision()),
        log_locale_(log_format.getloc()) {}

  // Stops the loading after the batch it is loading and waits for it.
  ~BarcodeBatchReadAhead() {
    if (loading_thread_.joinable()) {
      {
        std::lock_guard<std::mutex> lock(mutex_);
        stop_ = true;
      }
      changed_.notify_all();
      loading_thread_.join();
    }
  }

  // Starts loading on a new thread. Returns false, having loaded nothing, when
  // the thread cannot be created.
  bool Start() {
    try {
      loading_thread_ = std::thread(&BarcodeBatchReadAhead::Load, this);
    } catch (const std::system_error &) {
      return false;
    }
    return true;
  }

  // Waits for the next batch and returns what BarcodeKeyBatchLoader::
  // LoadNextBatch(keys, log) returned for it, or exits as it did.
  uint32_t TakeNextBatch(std::vector<uint64_t> &keys, std::ostream &log) {
    std::string error;
    uint32_t num_loaded_barcodes = 0;
    {
      std::unique_lock<std::mutex> lock(mutex_);
      changed_.wait(lock, [this]() { return is_ready_; });
      log << log_;
      error.swap(error_);
      num_loaded_barcodes = num_loaded_barcodes_;
      keys.swap(keys_);
      is_ready_ = false;
    }
    changed_.notify_all();
    if (!error.empty()) {
      ExitWithMessage(error);
    }
    return num_loaded_barcodes;
  }

 private:
  // Hands a batch, or an error, to TakeNextBatch once it has taken the one
  // before. Returns false when the loading is to stop instead.
  bool Publish(uint32_t num_loaded_barcodes, const std::string &log,
               const std::string &error, std::vector<uint64_t> &keys) {
    std::unique_lock<std::mutex> lock(mutex_);
    changed_.wait(lock, [this]() { return !is_ready_ || stop_; });
    if (stop_) {
      return false;
    }
    num_loaded_barcodes_ = num_loaded_barcodes;
    log_ = log;
    error_ = error;
    keys_.swap(keys);
    is_ready_ = true;
    lock.unlock();
    changed_.notify_all();
    return true;
  }

  void Load() {
    ExitIsDeferredOnThisThread() = true;
    std::vector<uint64_t> keys;
    for (size_t file_index = 0; file_index < num_files_; ++file_index) {
      uint32_t num_loaded_barcodes = 0;
      do {
        std::ostringstream log;
        log.flags(log_flags_);
        log.precision(log_precision_);
        log.imbue(log_locale_);
        std::string error;
        num_loaded_barcodes = 0;
        try {
          num_loaded_barcodes = loader_.LoadNextBatch(keys, log);
        } catch (const DeferredExit &deferred_exit) {
          error = deferred_exit.message;
        }
        if (!Publish(num_loaded_barcodes, log.str(), error, keys) ||
            !error.empty()) {
          return;
        }
      } while (num_loaded_barcodes > 0);
    }
  }

  const size_t num_files_;
  BarcodeKeyBatchLoader loader_;
  const std::ios_base::fmtflags log_flags_;
  const std::streamsize log_precision_;
  const std::locale log_locale_;

  std::mutex mutex_;
  std::condition_variable changed_;
  // Set while a batch, or an error, waits for TakeNextBatch.
  bool is_ready_ = false;
  bool stop_ = false;
  uint32_t num_loaded_barcodes_ = 0;
  std::string log_;
  std::string error_;
  std::vector<uint64_t> keys_;
  std::thread loading_thread_;
};

void Chromap::ConstructIndex() {
  // TODO(Haowen): Need a faster algorithm
  // Load all sequences in the reference into one batch
  SequenceBatch reference;
  reference.InitializeLoading(index_parameters_.reference_file_path,
                              /*default_buffers=*/true);
  reference.LoadAllSequences();
  const uint32_t num_sequences = reference.GetNumSequences();
  Index index(index_parameters_);
  index.Construct(num_sequences, reference);
  index.Statistics(num_sequences, reference);
  index.Save();
  reference.FinalizeLoading();
}

uint32_t Chromap::LoadSingleEndReadsWithBarcodes(SequenceBatch &read_batch,
                                                 SequenceBatch &barcode_batch,
                                                 bool parallel_parsing) {
  //double real_start_time = GetRealTime();
  uint32_t num_loaded_reads = 0;

  if (!parallel_parsing || mapping_parameters_.is_bulk_data) {
    while (num_loaded_reads < read_batch_size_) {
      bool no_more_read = read_batch.LoadOneSequenceAndSaveAt(num_loaded_reads);
      bool no_more_barcode = no_more_read;
      if (!mapping_parameters_.is_bulk_data) {
        no_more_barcode =
          barcode_batch.LoadOneSequenceAndSaveAt(num_loaded_reads);
      }

      if (no_more_read && no_more_barcode) {
        break;
      } else if (no_more_read || no_more_barcode){
        ExitWithMessage("Numbers of reads and barcodes don't match!");
      }
      ++num_loaded_reads;
    }
  } else {
    uint32_t num_loaded_barcode = 0 ;
    
    if (omp_in_parallel()) {
#pragma omp task shared(num_loaded_reads, read_batch)
      {
        num_loaded_reads = LoadUpToOneBatch(read_batch, read_batch_size_);
      }

#pragma omp task shared(num_loaded_barcode, barcode_batch)
      { // bulk data will go to the other big branch
        num_loaded_barcode = LoadUpToOneBatch(barcode_batch, read_batch_size_);
      }

#pragma omp taskwait
    } else {
      // Outside a parallel region the tasks above run at once, one after the
      // other, on this thread. These are the same loads in the same order,
      // without the tasks, which an exit deferred on this thread cannot leave.
      num_loaded_reads = LoadUpToOneBatch(read_batch, read_batch_size_);
      num_loaded_barcode = LoadUpToOneBatch(barcode_batch, read_batch_size_);
    }

    if (num_loaded_reads != num_loaded_barcode) {
        ExitWithMessage("Numbers of reads and barcodes don't match!");
    }
  }
  /*if (num_loaded_reads > 0) {
    std::cerr << "Loaded " << num_loaded_reads << " reads in "
              << GetRealTime() - real_start_time << "s.\n";
  } else {
    std::cerr << "No more reads.\n";
  }*/
  return num_loaded_reads;
}

uint32_t Chromap::LoadPairedEndReadsWithBarcodes(SequenceBatch &read_batch1,
                                                 SequenceBatch &read_batch2,
                                                 SequenceBatch &barcode_batch,
                                                 bool parallel_parsing,
                                                 int num_loading_helpers) {
  // double real_start_time = Chromap<>::GetRealTime();
  uint32_t num_loaded_pairs = 0;
  
  if (!parallel_parsing) {
    while (num_loaded_pairs < read_batch_size_) {
      bool no_more_read1 = read_batch1.LoadOneSequenceAndSaveAt(num_loaded_pairs);
      bool no_more_read2 = read_batch2.LoadOneSequenceAndSaveAt(num_loaded_pairs);
      bool no_more_barcode = no_more_read2;
      if (!mapping_parameters_.is_bulk_data) {
        no_more_barcode =
          barcode_batch.LoadOneSequenceAndSaveAt(num_loaded_pairs);
      }
      
      if (no_more_read1 && no_more_read2 && no_more_barcode) {
        break;
      } else if (no_more_read1 || no_more_read2 || no_more_barcode){
        ExitWithMessage("Numbers of reads and barcodes don't match!");
      }
      ++num_loaded_pairs;
    }
  } else {
    uint32_t num_loaded_read1 = 0;
    uint32_t num_loaded_read2 = 0;
    uint32_t num_loaded_barcode = 0;
    
    if (omp_in_parallel()) {
#pragma omp task shared(num_loaded_read1, read_batch1)
      {
        num_loaded_read1 = LoadUpToOneBatch(read_batch1, read_batch_size_);
      }

#pragma omp task shared(num_loaded_read2, read_batch2)
      {
        num_loaded_read2 = LoadUpToOneBatch(read_batch2, read_batch_size_);
      }

#pragma omp task shared(num_loaded_barcode, barcode_batch)
      {
        if (!mapping_parameters_.is_bulk_data) {
          num_loaded_barcode =
              LoadUpToOneBatch(barcode_batch, read_batch_size_);
        }
      }

#pragma omp taskwait
    } else if (num_loading_helpers > 0) {
      // The same loads, each file by one thread: this one and up to
      // num_loading_helpers more, which take the files in order as they come
      // free. Each load defers its exit, and the first exit in that order is
      // taken once every load has returned: the one loading the files one
      // after the other, below, would take.
      SequenceBatch *const batches[3] = {&read_batch1, &read_batch2,
                                         &barcode_batch};
      uint32_t *const num_loaded[3] = {&num_loaded_read1, &num_loaded_read2,
                                       &num_loaded_barcode};
      const size_t num_files = mapping_parameters_.is_bulk_data ? 2 : 3;
      std::string errors[3];
      std::atomic<size_t> next_file(0);
      const auto load_files = [&]() {
        for (size_t i = next_file.fetch_add(1); i < num_files;
             i = next_file.fetch_add(1)) {
          errors[i] = RunWithExitDeferred([&]() {
            *num_loaded[i] = LoadUpToOneBatch(*batches[i], read_batch_size_);
          });
        }
      };
      std::vector<std::thread> loaders;
      for (int t = 0; t < num_loading_helpers && loaders.size() + 1 < num_files;
           ++t) {
        try {
          loaders.emplace_back(load_files);
        } catch (const std::system_error &) {
          // This thread loads what the missing ones would have.
          break;
        }
      }
      load_files();
      for (std::thread &loader : loaders) {
        loader.join();
      }
      for (size_t i = 0; i < num_files; ++i) {
        if (!errors[i].empty()) {
          ExitWithMessage(errors[i]);
        }
      }
    } else {
      // Outside a parallel region the tasks above run at once, one after the
      // other, on this thread. These are the same loads in the same order,
      // without the tasks, which an exit deferred on this thread cannot leave.
      num_loaded_read1 = LoadUpToOneBatch(read_batch1, read_batch_size_);
      num_loaded_read2 = LoadUpToOneBatch(read_batch2, read_batch_size_);
      if (!mapping_parameters_.is_bulk_data) {
        num_loaded_barcode = LoadUpToOneBatch(barcode_batch, read_batch_size_);
      }
    }
    if (mapping_parameters_.is_bulk_data) {
      num_loaded_barcode = num_loaded_read2;
    }
    if (num_loaded_read1 != num_loaded_read2 || num_loaded_read2 != num_loaded_barcode) {
        ExitWithMessage("Numbers of reads and barcodes don't match!");
    }

    num_loaded_pairs = num_loaded_read1 ;
  }
  // if (num_loaded_pairs > 0) {
  //  std::cerr << "Loaded " << num_loaded_pairs << " pairs in "<<
  //  Chromap<>::GetRealTime() - real_start_time << "s. ";
  //} else {
  //  std::cerr << "No more reads.\n";
  //}
  return num_loaded_pairs;
}

uint32_t Chromap::OpenAndLoadFirstSingleEndBatch(size_t read_file_index,
                                                 SequenceBatch &read_batch,
                                                 SequenceBatch &barcode_batch) {
  read_batch.InitializeLoading(
      mapping_parameters_.read_file1_paths[read_file_index]);
  if (!mapping_parameters_.is_bulk_data) {
    barcode_batch.InitializeLoading(
        mapping_parameters_.barcode_file_paths[read_file_index]);
  }
  return LoadSingleEndReadsWithBarcodes(
      read_batch, barcode_batch,
      mapping_parameters_.num_threads >= 3 ? true : false);
}

uint32_t Chromap::OpenAndLoadFirstPairedEndBatch(size_t read_file_index,
                                                 SequenceBatch &read_batch1,
                                                 SequenceBatch &read_batch2,
                                                 SequenceBatch &barcode_batch,
                                                 int num_loading_helpers) {
  read_batch1.InitializeLoading(
      mapping_parameters_.read_file1_paths[read_file_index]);
  read_batch2.InitializeLoading(
      mapping_parameters_.read_file2_paths[read_file_index]);
  if (!mapping_parameters_.is_bulk_data) {
    barcode_batch.InitializeLoading(
        mapping_parameters_.barcode_file_paths[read_file_index]);
  }
  return LoadPairedEndReadsWithBarcodes(
      read_batch1, read_batch2, barcode_batch,
      mapping_parameters_.num_threads >= 3 ? true : false,
      num_loading_helpers);
}

std::string Chromap::RunWithExitDeferred(const std::function<void()> &work) {
  std::string error;
  const bool was_deferred = ExitIsDeferredOnThisThread();
  ExitIsDeferredOnThisThread() = true;
  try {
    work();
  } catch (const DeferredExit &deferred_exit) {
    error = deferred_exit.message;
  }
  ExitIsDeferredOnThisThread() = was_deferred;
  return error;
}

void Chromap::TrimAdapterForPairedEndRead(uint32_t pair_index,
                                          SequenceBatch &read_batch1,
                                          SequenceBatch &read_batch2) {
  const uint32_t raw_read1_length = read_batch1.GetSequenceLengthAt(pair_index);
  const uint32_t raw_read2_length = read_batch2.GetSequenceLengthAt(pair_index);
  const char *raw_read1 = read_batch1.GetSequenceAt(pair_index);
  const char *raw_read2 = read_batch2.GetSequenceAt(pair_index);
  const std::string &raw_negative_read1 =
      read_batch1.GetNegativeSequenceAt(pair_index);
  const std::string &raw_negative_read2 =
      read_batch2.GetNegativeSequenceAt(pair_index);

  // In the actual adaptor trimming, we assuem length(read1)<=length(read2). So
  // we can have the case that read1 is a subset of read2.
  const char *read1 =
      raw_read1_length <= raw_read2_length ? raw_read1 : raw_read2;
  const std::string &negative_read2 = raw_read1_length <= raw_read2_length
                                          ? raw_negative_read2
                                          : raw_negative_read1;
  const uint32_t read1_length = raw_read1_length <= raw_read2_length
                                    ? raw_read1_length
                                    : raw_read2_length;
  const uint32_t read2_length = raw_read1_length <= raw_read2_length
                                    ? raw_read2_length
                                    : raw_read1_length;

  const int min_overlap_length = mapping_parameters_.min_read_length;
  const int seed_length = min_overlap_length / 2;
  const int error_threshold_for_merging = 1;
  bool is_merged = false;

  for (int si = 0; si < error_threshold_for_merging + 1; ++si) {
    size_t seed_start_position =
        negative_read2.find(read1 + si * seed_length, 0, seed_length);

    while (seed_start_position != std::string::npos) {
      const bool before_seed_is_enough_long =
          seed_start_position >= (size_t)(si * seed_length);
      const bool overlap_is_enough_long =
          (int)(read2_length - seed_start_position + seed_length * si) >=
          min_overlap_length;

      if (!before_seed_is_enough_long || !overlap_is_enough_long) {
        seed_start_position = negative_read2.find(
            read1 + si * seed_length, seed_start_position + 1, seed_length);
        continue;
      }

      bool can_merge = true;
      int num_errors = 0;

      // The bases before the seed.
      for (int i = 0; i < seed_length * si; ++i) {
        if (negative_read2[seed_start_position - si * seed_length + i] !=
            read1[i]) {
          ++num_errors;
        }
        if (num_errors > error_threshold_for_merging) {
          can_merge = false;
          break;
        }
      }

      // The bases after the seed.
      for (uint32_t i = seed_length; i + seed_start_position < read2_length &&
                                     si * seed_length + i < read1_length;
           ++i) {
        if (negative_read2[seed_start_position + i] !=
            read1[si * seed_length + i]) {
          ++num_errors;
        }
        if (num_errors > error_threshold_for_merging) {
          can_merge = false;
          break;
        }
      }

      if (can_merge) {
        // Trim adapters and TODO: fix sequencing errors
        int overlap_length =
            read2_length - seed_start_position + si * seed_length;
        int read2_offset = 0;
        // The case that read1 is strictly contained in read2. overlap_length is
        // inferred from the longer read2, which could be longer than read1. In
        // that case, we don't trim read1 (make overlap length equal to read1
        // length) and trim read2 as the original plan.
        if (overlap_length > (int)read1_length) {
          read2_offset = overlap_length - read1_length;
          overlap_length = read1_length;
        }

        if (raw_read1_length <= raw_read2_length) {
          read_batch1.TrimSequenceAt(pair_index, overlap_length);
          read_batch2.TrimSequenceAt(pair_index, overlap_length + read2_offset);
        } else {
          read_batch1.TrimSequenceAt(pair_index, overlap_length + read2_offset);
          read_batch2.TrimSequenceAt(pair_index, overlap_length);
        }

        is_merged = true;
        // std::cerr << "Trimed! overlap length: " << overlap_length << ", " <<
        // read1.GetLength() << " " << read2.GetLength() << "\n";
        break;
      }

      seed_start_position = negative_read2.find(
          read1 + si * seed_length, seed_start_position + 1, seed_length);
    }

    if (is_merged) {
      break;
    }
  }
}

bool Chromap::PairedEndReadWithBarcodeIsDuplicate(
    uint32_t pair_index, const SequenceBatch &barcode_batch,
    const SequenceBatch &read_batch1, const SequenceBatch &read_batch2) {
  int dedupe_seed_length = 16;
  uint32_t barcode_length = barcode_batch.GetSequenceLengthAt(pair_index);
  uint64_t barcode_key =
      barcode_batch.GenerateSeedFromSequenceAt(pair_index, 0, barcode_length);
  uint64_t read1_seed1 =
      read_batch1.GenerateSeedFromSequenceAt(pair_index, 0, dedupe_seed_length);
  uint64_t read2_seed1 =
      read_batch2.GenerateSeedFromSequenceAt(pair_index, 0, dedupe_seed_length);
  uint64_t read_seed_key =
      (read1_seed1 << (dedupe_seed_length * 2)) | read2_seed1;
  uint64_t read1_seed2 = read_batch1.GenerateSeedFromSequenceAt(
      pair_index, dedupe_seed_length, dedupe_seed_length * 2);
  uint64_t read2_seed2 = read_batch2.GenerateSeedFromSequenceAt(
      pair_index, dedupe_seed_length, dedupe_seed_length * 2);
  khiter_t barcode_table_iterator =
      kh_get(k64_seq, barcode_lookup_table_, barcode_key);
  if (barcode_table_iterator != kh_end(barcode_lookup_table_)) {
    uint32_t read_lookup_table_index =
        kh_value(barcode_lookup_table_, barcode_table_iterator);
    // std::cerr << "Have barcode, try to check read. " <<
    // read_lookup_table_index << "\n";
    khash_t(k128) *read_lookup_table =
        read_lookup_tables_[read_lookup_table_index];
    khiter_t read_lookup_table_iterator =
        kh_get(k128, read_lookup_table, read_seed_key);
    if (read_lookup_table_iterator != kh_end(read_lookup_table)) {
      // std::cerr << "Have barcode, have read, try whether match.\n";
      uint128_t read_seeds =
          kh_value(read_lookup_table, read_lookup_table_iterator);
      if (read_seeds.first == read1_seed2 && read_seeds.second == read2_seed2) {
        // std::cerr << "Have barcode, have read, and match.\n";
        return true;
      } else {
        // std::cerr << "Have barcode, have read, but don't match.\n";
        return false;
      }
    } else {
      // std::cerr << "Have barcode, no read.\n";
      uint128_t read_seeds = {.first = read1_seed2, .second = read2_seed2};
      int khash_return_code;
      khiter_t read_lookup_table_insert_iterator =
          kh_put(k128, read_lookup_table, read_seed_key, &khash_return_code);
      assert(khash_return_code != -1 && khash_return_code != 0);
      kh_value(read_lookup_table, read_lookup_table_insert_iterator) =
          read_seeds;
      // std::cerr << "Have barcode, no read.\n";
      return false;
    }
  } else {
    // insert the barcode and append a new read hash table to tables and then
    // insert the reads
    // std::cerr << "No barcode, no read.\n";
    int khash_return_code;
    khiter_t barcode_table_insert_iterator =
        kh_put(k64_seq, barcode_lookup_table_, barcode_key, &khash_return_code);
    assert(khash_return_code != -1 && khash_return_code != 0);
    kh_value(barcode_lookup_table_, barcode_table_insert_iterator) =
        read_lookup_tables_.size();
    khash_t(k128) *read_lookup_table = kh_init(k128);
    khiter_t read_lookup_table_iterator =
        kh_put(k128, read_lookup_table, read_seed_key, &khash_return_code);
    assert(khash_return_code != -1 && khash_return_code != 0);
    uint128_t read_seeds = {.first = read1_seed2, .second = read2_seed2};
    kh_value(read_lookup_table, read_lookup_table_iterator) = read_seeds;
    read_lookup_tables_.push_back(read_lookup_table);
    // std::cerr << "No barcode, no read.\n";
    return false;
  }
}

uint32_t Chromap::SampleInputBarcodesAndExamineLength(std::ostream &log) {
  if (mapping_parameters_.is_bulk_data) {
    return 0;
  }

  uint32_t sample_batch_size = 1000;
  SequenceBatch barcode_batch(sample_batch_size, barcode_effective_range_);

  barcode_batch.InitializeLoading(mapping_parameters_.barcode_file_paths[0]);

  uint32_t num_loaded_barcodes = barcode_batch.LoadBatch(log);

  uint32_t cell_barcode_length = barcode_batch.GetSequenceLengthAt(0);
  for (uint32_t i = 1; i < num_loaded_barcodes; ++i) {
    if (barcode_batch.GetSequenceLengthAt(i) != cell_barcode_length) {
      ExitWithMessage("ERROR: barcode lengths are not equal in the sample!");
    }
  }

  barcode_batch.FinalizeLoading();

  return cell_barcode_length;
}

void Chromap::LoadBarcodeWhitelist(std::ostream &log) {
  double real_start_time = GetRealTime();
  int num_barcodes = 0;

  if (1) {
    gzFile barcode_whitelist_file = 
      gzopen(mapping_parameters_.barcode_whitelist_file_path.c_str(), "r"); 
    const uint32_t barcode_buffer_size = 256;
    char barcode[barcode_buffer_size];
    while (gzgets(barcode_whitelist_file, barcode, barcode_buffer_size) != NULL) {
      size_t barcode_length = strlen(barcode);
      if (barcode[barcode_length - 1] == '\n') {
        barcode[barcode_length - 1] = '\0';
        --barcode_length;
      }
      if (barcode_length > 32) {
        ExitWithMessage("ERROR: barcode length is greater than 32!");
      }

      if (barcode_length != barcode_length_) {
        if (num_barcodes == 0) {
          ExitWithMessage(
              "ERROR: whitelist and input barcode lengths are not equal!");
        } else {
          ExitWithMessage(
              "ERROR: barcode lengths are not equal in the whitelist!");
        }
      }
      
      uint64_t barcode_key = GenerateSeedFromSequence(
          barcode, barcode_length, 0, barcode_length);
      
      int khash_return_code;
      khiter_t barcode_whitelist_lookup_table_iterator =
        kh_put(k64_seq, barcode_whitelist_lookup_table_, barcode_key,
            &khash_return_code);
      kh_value(barcode_whitelist_lookup_table_,
          barcode_whitelist_lookup_table_iterator) = 0;
      assert(khash_return_code != -1 && khash_return_code != 0);
      ++num_barcodes;
    }
    if (!gzeof(barcode_whitelist_file)) {
      ExitWithMessage("ERROR: barcode whitelist file does not exist or is truncated!");
    }
    gzclose(barcode_whitelist_file);
  } else { 
    std::ifstream barcode_whitelist_file_stream(
        mapping_parameters_.barcode_whitelist_file_path);
    std::string barcode_whitelist_file_line;
    // bool first_line = true;
    while (getline(barcode_whitelist_file_stream, barcode_whitelist_file_line)) {
      std::stringstream barcode_whitelist_file_line_string_stream(
          barcode_whitelist_file_line);
      //// skip the header
      // if (barcode_whitelist_file_line[0] == '#' ||
      // barcode_whitelist_file_line.find("kmer") == 0) {
      //  continue;
      //}
      std::string barcode;
      barcode_whitelist_file_line_string_stream >> barcode;
      size_t barcode_length = barcode.length();
      if (barcode_length > 32) {
        ExitWithMessage("ERROR: barcode length is greater than 32!");
      }

      if (barcode_length != barcode_length_) {
        if (num_barcodes == 0) {
          ExitWithMessage(
              "ERROR: whitelist and input barcode lengths are not equal!");
        } else {
          ExitWithMessage(
              "ERROR: barcode lengths are not equal in the whitelist!");
        }
      }

      // if (first_line) {
      //  //size_t barcode_length = kmer.length();
      //  // Allocate memory to save pore model parameters
      //  //size_t num_pore_models = 1 << (kmer_size_ * 2);
      //  //pore_models_.assign(num_pore_models, PoreModelParameters());
      //  //first_line = false;
      //}
      // assert(kmer.length() == (size_t)kmer_size_);
      uint64_t barcode_key = GenerateSeedFromSequence(
          barcode.data(), barcode_length, 0, barcode_length);
      // PoreModelParameters &pore_model_parameters =
      // pore_models_[kmer_hash_value]; barcode_whitelist_file_line_string_stream
      // >> pore_model_parameters.level_mean >> pore_model_parameters.level_stdv
      // >> pore_model_parameters.sd_mean >> pore_model_parameters.sd_stdv;
      int khash_return_code;
      khiter_t barcode_whitelist_lookup_table_iterator =
        kh_put(k64_seq, barcode_whitelist_lookup_table_, barcode_key,
            &khash_return_code);
      kh_value(barcode_whitelist_lookup_table_,
          barcode_whitelist_lookup_table_iterator) = 0;
      assert(khash_return_code != -1 && khash_return_code != 0);
      ++num_barcodes;
    }
    barcode_whitelist_file_stream.close();
  }
  log << "Loaded " << num_barcodes << " barcodes in "
      << GetRealTime() - real_start_time << "s.\n";
}

void Chromap::ComputeBarcodeAbundance(uint64_t max_num_sample_barcodes,
                                      std::ostream &log,
                                      BarcodeBatchReadAhead *read_ahead) {
  double real_start_time = GetRealTime();
  std::unique_ptr<BarcodeKeyBatchLoader> loader;
  if (read_ahead == nullptr) {
    loader.reset(new BarcodeKeyBatchLoader(
        mapping_parameters_.barcode_file_paths, read_batch_size_,
        barcode_effective_range_));
  }
  // The whitelist keys of a batch's barcodes without an N, in batch order.
  std::vector<uint64_t> barcode_keys;
  for (size_t read_file_index = 0;
       read_file_index < mapping_parameters_.read_file1_paths.size();
       ++read_file_index) {
    uint32_t num_loaded_barcodes =
        read_ahead != nullptr ? read_ahead->TakeNextBatch(barcode_keys, log)
                              : loader->LoadNextBatch(barcode_keys, log);
    while (num_loaded_barcodes > 0) {
      // The lookups, in barcode order, each with the first bucket kh_get
      // probes for the key kPrefetchDistance later prefetched (flags word, key
      // and value), so that their cache misses overlap.
      const size_t kPrefetchDistance = 16;
      const khash_t(k64_seq) *whitelist = barcode_whitelist_lookup_table_;
      for (size_t key_index = 0; key_index < barcode_keys.size();
           ++key_index) {
#if defined(__GNUC__) || defined(__clang__)
        if (key_index + kPrefetchDistance < barcode_keys.size() &&
            whitelist->n_buckets > 0) {
          const khint_t i =
              kh_int64_hash_func(barcode_keys[key_index + kPrefetchDistance]) &
              (whitelist->n_buckets - 1);
          __builtin_prefetch(&whitelist->flags[i >> 4]);
          __builtin_prefetch(&whitelist->keys[i]);
          __builtin_prefetch(&whitelist->vals[i]);
        }
#endif
        khiter_t barcode_whitelist_lookup_table_iterator =
            kh_get(k64_seq, barcode_whitelist_lookup_table_,
                   barcode_keys[key_index]);
        if (barcode_whitelist_lookup_table_iterator !=
            kh_end(barcode_whitelist_lookup_table_)) {
          // Correct barcode
          kh_value(barcode_whitelist_lookup_table_,
                   barcode_whitelist_lookup_table_iterator) += 1;
          ++num_sample_barcodes_;
        }
      }
      if (!mapping_parameters_.skip_barcode_check &&
          num_sample_barcodes_ * 20 < num_loaded_barcodes) {
        // Since num_loaded_pairs is a constant, this if is actuaclly only
        // effective in the first iteration
        ExitWithMessage(
            "Less than 5\% barcodes can be found or corrected based on the "
            "barcode whitelist.\nPlease check whether the barcode whitelist "
            "matches the data, e.g. length, reverse-complement. If this is a "
            "false warning, please run Chromap with the option "
            "--skip-barcode-check.");
      }

      if (num_sample_barcodes_ >= max_num_sample_barcodes) {
        break;
      }
      num_loaded_barcodes =
          read_ahead != nullptr ? read_ahead->TakeNextBatch(barcode_keys, log)
                                : loader->LoadNextBatch(barcode_keys, log);
    }
    if (num_sample_barcodes_ >= max_num_sample_barcodes) {
      break;
    }
  }

  log << "Compute barcode abundance using " << num_sample_barcodes_ << " in "
      << GetRealTime() - real_start_time << "s.\n";
}

std::string Chromap::PrepareBarcodes(std::ostream &log,
                                     int num_loading_helpers) {
  std::string error;
  const bool was_deferred = ExitIsDeferredOnThisThread();
  ExitIsDeferredOnThisThread() = true;
  try {
    barcode_length_ = SampleInputBarcodesAndExamineLength(log);
    if (!mapping_parameters_.barcode_whitelist_file_path.empty()) {
      // With a helper, the abundance pass's barcode batches load on it from
      // now on, while the whitelist loads; it is stopped and joined before
      // this returns or an exit leaves. The pass reads a barcode file per read
      // file, and with fewer barcode files it runs here, as it did.
      std::unique_ptr<BarcodeBatchReadAhead> read_ahead;
      const size_t num_files = mapping_parameters_.read_file1_paths.size();
      if (num_loading_helpers > 0 &&
          mapping_parameters_.barcode_file_paths.size() >= num_files) {
        read_ahead.reset(new BarcodeBatchReadAhead(
            mapping_parameters_.barcode_file_paths, num_files,
            read_batch_size_, barcode_effective_range_, log));
        if (!read_ahead->Start()) {
          // This thread loads them in the pass.
          read_ahead.reset();
        }
      }
      LoadBarcodeWhitelist(log);
      ComputeBarcodeAbundance(initial_num_sample_barcodes_, log,
                              read_ahead.get());
    }
  } catch (const DeferredExit &deferred_exit) {
    error = deferred_exit.message;
  }
  ExitIsDeferredOnThisThread() = was_deferred;
  return error;
}

void Chromap::UpdateBarcodeAbundance(uint32_t num_loaded_barcodes,
                                     const SequenceBatch &barcode_batch) {
  double real_start_time = GetRealTime();
  for (uint32_t barcode_index = 0; barcode_index < num_loaded_barcodes;
       ++barcode_index) {
    uint32_t barcode_length = barcode_batch.GetSequenceLengthAt(barcode_index);
    uint64_t barcode_key = barcode_batch.GenerateSeedFromSequenceAt(
        barcode_index, 0, barcode_length);
    khiter_t barcode_whitelist_lookup_table_iterator =
        kh_get(k64_seq, barcode_whitelist_lookup_table_, barcode_key);
    if (barcode_whitelist_lookup_table_iterator !=
        kh_end(barcode_whitelist_lookup_table_)) {
      // Correct barcode
      kh_value(barcode_whitelist_lookup_table_,
               barcode_whitelist_lookup_table_iterator) += 1;
      ++num_sample_barcodes_;
    }
  }
  std::cerr << "Update barcode abundance using " << num_sample_barcodes_
            << " in " << GetRealTime() - real_start_time << "s.\n";
}

bool Chromap::CorrectBarcodeAt(uint32_t barcode_index,
                               SequenceBatch &barcode_batch,
                               uint64_t &num_barcode_in_whitelist,
                               uint64_t &num_corrected_barcode) {
  const uint32_t barcode_length =
      barcode_batch.GetSequenceLengthAt(barcode_index);
  const uint64_t barcode_key = barcode_batch.GenerateSeedFromSequenceAt(
      barcode_index, 0, barcode_length);
  khiter_t barcode_whitelist_lookup_table_iterator =
      kh_get(k64_seq, barcode_whitelist_lookup_table_, barcode_key);
  std::vector<int> N_pos;  // position of Ns

  barcode_batch.GetSequenceNsAt(barcode_index, /*little_endian=*/true, N_pos);
  if (N_pos.size() >
      (uint32_t)mapping_parameters_.barcode_correction_error_threshold)
    return false;

  if (N_pos.size() == 0 && barcode_whitelist_lookup_table_iterator !=
                               kh_end(barcode_whitelist_lookup_table_)) {
    // Correct barcode
    ++num_barcode_in_whitelist;
    return true;
  } else if (mapping_parameters_.barcode_correction_error_threshold > 0) {
    // Need to correct this barcode
    // const char *barcode = barcode_batch->GetSequenceAt(barcode_index);
    // std::cerr << barcode_index << " barcode " << barcode << " needs
    // correction\n";
    const char *barcode_qual = barcode_batch.GetSequenceQualAt(barcode_index);
    std::vector<BarcodeWithQual> corrected_barcodes_with_quals;
    uint64_t mask = (uint64_t)3;
    uint32_t i_start = 0;
    uint32_t i_end = barcode_length;
    uint32_t ti_limit = 3;
    if (N_pos.size() > 0) {
      i_start = N_pos[0];
      i_end = N_pos[0] + 1;
      ti_limit = 4;
    }
    for (uint32_t i = i_start; i < i_end; ++i) {
      uint64_t barcode_key_to_change = mask << (2 * i);
      barcode_key_to_change = ~barcode_key_to_change;
      barcode_key_to_change &= barcode_key;
      uint64_t base_to_change1 = (barcode_key >> (2 * i)) & mask;
      for (uint32_t ti = 0; ti < ti_limit; ++ti) {
        // change the base
        base_to_change1 += 1;
        base_to_change1 &= mask;
        // generate the corrected key
        uint64_t corrected_barcode_key =
            barcode_key_to_change | (base_to_change1 << (2 * i));
        barcode_whitelist_lookup_table_iterator = kh_get(
            k64_seq, barcode_whitelist_lookup_table_, corrected_barcode_key);
        if (barcode_whitelist_lookup_table_iterator !=
            kh_end(barcode_whitelist_lookup_table_)) {
          // find one possible corrected barcode
          double barcode_abundance =
              kh_value(barcode_whitelist_lookup_table_,
                       barcode_whitelist_lookup_table_iterator) /
              (double)num_sample_barcodes_;
          int qual_offset = 33;
          int adjusted_qual =
              barcode_qual[barcode_length - 1 - i] - qual_offset;
          adjusted_qual = adjusted_qual > 40 ? 40 : adjusted_qual;
          adjusted_qual = adjusted_qual < 3 ? 3 : adjusted_qual;
          double score =
              pow(10.0, ((-adjusted_qual) / 10.0)) * barcode_abundance;
          corrected_barcodes_with_quals.emplace_back(
              BarcodeWithQual{barcode_length - 1 - i,
                              Uint8ToChar(base_to_change1), 0, 0, score});
          // std::cerr << "1score: " << score << " pos1: " << barcode_length - 1
          // - i << " b1: " << base_to_change1 << " pos2: " << 0 << " b2: " <<
          // (char)0 << "\n";
        }
        if (mapping_parameters_.barcode_correction_error_threshold == 2) {
          uint32_t j_start = i + 1;
          uint32_t j_end = barcode_length;
          uint32_t ti2_limit = 3;
          if (N_pos.size() == 2) {
            j_start = N_pos[1];
            j_end = N_pos[1] + 1;
            ti2_limit = 4;
          }
          for (uint32_t j = j_start; j < j_end; ++j) {
            uint64_t barcode_key_to_change2 = mask << (2 * i);
            barcode_key_to_change2 = mask << (2 * j);
            barcode_key_to_change2 = ~barcode_key_to_change2;
            barcode_key_to_change2 &= corrected_barcode_key;
            uint64_t base_to_change2 =
                (corrected_barcode_key >> (2 * j)) & mask;
            for (uint32_t ti2 = 0; ti2 < ti2_limit; ++ti2) {
              // change the base
              base_to_change2 += 1;
              base_to_change2 &= mask;
              // generate the corrected key
              uint64_t corrected_barcode_key2 =
                  barcode_key_to_change2 | (base_to_change2 << (2 * j));
              barcode_whitelist_lookup_table_iterator =
                  kh_get(k64_seq, barcode_whitelist_lookup_table_,
                         corrected_barcode_key2);
              if (barcode_whitelist_lookup_table_iterator !=
                  kh_end(barcode_whitelist_lookup_table_)) {
                // find one possible corrected barcode
                double barcode_abundance =
                    kh_value(barcode_whitelist_lookup_table_,
                             barcode_whitelist_lookup_table_iterator) /
                    (double)num_sample_barcodes_;
                int qual_offset = 33;
                int adjusted_qual =
                    barcode_qual[barcode_length - 1 - j] - qual_offset;
                adjusted_qual = adjusted_qual > 40 ? 40 : adjusted_qual;
                adjusted_qual = adjusted_qual < 3 ? 3 : adjusted_qual;
                int adjusted_qual1 =
                    barcode_qual[barcode_length - 1 - i] - qual_offset;
                adjusted_qual1 = adjusted_qual1 > 40 ? 40 : adjusted_qual1;
                adjusted_qual1 = adjusted_qual1 < 3 ? 3 : adjusted_qual1;
                adjusted_qual += adjusted_qual1;
                double score =
                    pow(10.0, ((-adjusted_qual) / 10.0)) * barcode_abundance;
                corrected_barcodes_with_quals.emplace_back(BarcodeWithQual{
                    barcode_length - 1 - i, Uint8ToChar(base_to_change1),
                    barcode_length - 1 - j, Uint8ToChar(base_to_change2),
                    score});
                // std::cerr << "2score: " << score << " pos1: " <<
                // barcode_length - 1 - i << " b1: " << base_to_change1 << "
                // pos2: " << barcode_length - 1 -j << " b2: " <<
                // base_to_change2
                // << "\n";
              }
            }
          }
        }
      }
    }
    size_t num_possible_corrected_barcodes =
        corrected_barcodes_with_quals.size();
    if (num_possible_corrected_barcodes == 0) {
      // Barcode cannot be corrected, leave it for downstream
      return false;
    } else if (num_possible_corrected_barcodes == 1) {
      // Just correct it
      // std::cerr << "Corrected the barcode from " << barcode << " to ";
      barcode_batch.CorrectBaseAt(
          barcode_index, corrected_barcodes_with_quals[0].corrected_base_index1,
          corrected_barcodes_with_quals[0].correct_base1);
      if (corrected_barcodes_with_quals[0].correct_base2 != 0) {
        barcode_batch.CorrectBaseAt(
            barcode_index,
            corrected_barcodes_with_quals[0].corrected_base_index2,
            corrected_barcodes_with_quals[0].correct_base2);
      }
      // std::cerr << barcode << "\n";
      // std::cerr << "score: " << corrected_barcodes_with_quals[0].score <<
      // "\n"; std::cerr << "score: " << corrected_barcodes_with_quals[0].score
      // << " pos1: " << corrected_barcodes_with_quals[0].corrected_base_index1
      // << " b1: " << corrected_barcodes_with_quals[0].correct_base1 << " pos2:
      // " << corrected_barcodes_with_quals[0].corrected_base_index2 << " b2: "
      // << corrected_barcodes_with_quals[0].correct_base2 << "\n";
      ++num_corrected_barcode;
      return true;
    } else {
      // Select the best correction
      std::sort(corrected_barcodes_with_quals.begin(),
                corrected_barcodes_with_quals.end(),
                std::greater<BarcodeWithQual>());
      // int num_ties = 0;
      double sum_score = 0;
      for (size_t ci = 0; ci < num_possible_corrected_barcodes; ++ci) {
        sum_score += corrected_barcodes_with_quals[ci].score;
        // std::cerr << ci << " score: " <<
        // corrected_barcodes_with_quals[ci].score << " pos1: " <<
        // corrected_barcodes_with_quals[ci].corrected_base_index1 << " b1: " <<
        // corrected_barcodes_with_quals[ci].correct_base1 << " pos2: " <<
        // corrected_barcodes_with_quals[ci].corrected_base_index2 << " b2: " <<
        // corrected_barcodes_with_quals[ci].correct_base2 << "\n"; if
        // (corrected_barcodes_with_quals[ci].qual ==
        // corrected_barcodes_with_quals[0].qual) {
        //  ++num_ties;
        //}
      }
      int best_corrected_barcode_index = 0;
      // if (num_ties > 0) {
      //  std::mt19937 tmp_generator(11);
      //  std::uniform_int_distribution<int> distribution(0, num_ties); //
      //  important: inclusive range best_corrected_barcode_index =
      //  distribution(tmp_generator);
      //}
      // std::cerr << "Corrected the barcode from " << barcode << " to ";
      double confidence_threshold =
          mapping_parameters_.barcode_correction_probability_threshold;
      if (corrected_barcodes_with_quals[best_corrected_barcode_index].score /
              sum_score >
          confidence_threshold) {
        barcode_batch.CorrectBaseAt(
            barcode_index,
            corrected_barcodes_with_quals[best_corrected_barcode_index]
                .corrected_base_index1,
            corrected_barcodes_with_quals[best_corrected_barcode_index]
                .correct_base1);
        if (corrected_barcodes_with_quals[best_corrected_barcode_index]
                .correct_base2 != 0) {
          barcode_batch.CorrectBaseAt(
              barcode_index,
              corrected_barcodes_with_quals[best_corrected_barcode_index]
                  .corrected_base_index2,
              corrected_barcodes_with_quals[best_corrected_barcode_index]
                  .correct_base2);
        }
        // std::cerr << barcode << "\n";
        // std::cerr << "score: " <<
        // corrected_barcodes_with_quals[best_corrected_barcode_index].score <<
        // "\n"; std::cerr << "best score: " <<
        // corrected_barcodes_with_quals[best_corrected_barcode_index].score <<
        // " sum score: " << sum_score << "\n";
        ++num_corrected_barcode;
        return true;
      } else {
        // std::cerr << "Didnt pass filter: " <<
        // corrected_barcodes_with_quals[best_corrected_barcode_index].score /
        // sum_score << "\n"; std::cerr << "best score: " <<
        // corrected_barcodes_with_quals[best_corrected_barcode_index].score <<
        // " sum score: " << sum_score << "\n";
        return false;
      }
    }
  } else {
    return false;
  }
}

void Chromap::OutputBarcodeStatistics() {
  std::cerr << "Number of barcodes in whitelist: " << num_barcode_in_whitelist_
            << ".\n";
  std::cerr << "Number of corrected barcodes: " << num_corrected_barcode_
            << ".\n";
}

void Chromap::OutputMappingStatistics() {
  std::cerr << "Number of reads: " << num_reads_ << ".\n";
  // std::cerr << "Number of duplicated reads: " << num_duplicated_reads_ <<
  // ".\n";
  std::cerr << "Number of mapped reads: " << num_mapped_reads_ << ".\n";
  std::cerr << "Number of uniquely mapped reads: " << num_uniquely_mapped_reads_
            << ".\n";
  std::cerr << "Number of reads have multi-mappings: "
            << num_mapped_reads_ - num_uniquely_mapped_reads_ << ".\n";
  std::cerr << "Number of candidates: " << num_candidates_ << ".\n";
  std::cerr << "Number of mappings: " << num_mappings_ << ".\n";
  std::cerr << "Number of uni-mappings: " << num_uniquely_mapped_reads_
            << ".\n";
  std::cerr << "Number of multi-mappings: "
            << num_mappings_ - num_uniquely_mapped_reads_ << ".\n";
}

void Chromap::ParseReadFormat(const std::string &read_format) {
  if (read_format.empty()) {
    return;
  }

  read1_effective_range_.InitializeParsing();
  read2_effective_range_.InitializeParsing();
  barcode_effective_range_.InitializeParsing();

  uint32_t i, j;
  for (i = 0; i < read_format.size();) {
    for (j = i + 1; j < read_format.size() && read_format[j] != ','; ++j)
      ;
    bool parse_success = true;
    if (read_format[i] == 'r' && read_format[i + 1] == '1') {
      parse_success =
          read1_effective_range_.ParseFormatStringAndAppendEffectiveRange(
              read_format.c_str() + i, j - i);
    } else if (read_format[i] == 'r' && read_format[i + 1] == '2') {
      parse_success =
          read2_effective_range_.ParseFormatStringAndAppendEffectiveRange(
              read_format.c_str() + i, j - i);
    } else if (read_format[i] == 'b' && read_format[i + 1] == 'c') {
      parse_success =
          barcode_effective_range_.ParseFormatStringAndAppendEffectiveRange(
              read_format.c_str() + i, j - i);
    } else {
      parse_success = false;
    }

    if (!parse_success) {
      ExitWithMessage("Unknown read format: " + read_format + "\n");
    }

    i = j + 1;
  }

  read1_effective_range_.FinalizeParsing();
  read2_effective_range_.FinalizeParsing();
  barcode_effective_range_.FinalizeParsing();
}

void Chromap::GenerateCustomRidRanks(
    const std::string &custom_rid_order_file_path,
    uint32_t num_reference_sequences, const SequenceBatch &reference,
    std::vector<int> &rid_ranks) {
  for (uint32_t i = 0; i < num_reference_sequences; ++i) {
    rid_ranks.emplace_back(i);
  }

  if (custom_rid_order_file_path.empty()) {
    return;
  }

  std::unordered_map<std::string, int> ref_name_to_rank;
  std::ifstream custom_rid_order_file_stream(custom_rid_order_file_path);
  std::string ref_name;
  uint32_t ref_rank = 0;
  while (getline(custom_rid_order_file_stream, ref_name)) {
    ref_name_to_rank[ref_name] = ref_rank;
    ref_rank += 1;
  }
  custom_rid_order_file_stream.close();

  // First, rank the chromosomes in the custom order provided by users.
  for (uint32_t i = 0; i < num_reference_sequences; ++i) {
    std::string ref_name(reference.GetSequenceNameAt(i));
    if (ref_name_to_rank.find(ref_name) != ref_name_to_rank.end()) {
      rid_ranks[i] = ref_name_to_rank[ref_name];
    } else {
      rid_ranks[i] = -1;
    }
  }

  // There might be some rids without any custom order. We just order them based
  // on their original order in the reference file.
  uint32_t k = ref_name_to_rank.size();
  // Rank the remaining chromosomes.
  for (uint32_t i = 0; i < num_reference_sequences; ++i) {
    if (rid_ranks[i] == -1) {
      rid_ranks[i] = k;
      ++k;
    }
  }

  if (k > num_reference_sequences) {
    ExitWithMessage(
        "ERROR: unknown chromsome names found in chromosome order file.");
  }
}

std::string Chromap::LoadReferenceAndIndex(
    SequenceBatch &reference, Index &index,
    const std::function<void(int)> &prepare_mapping,
    const std::function<void(int)> &prepare_barcodes,
    int max_num_preparing_helpers, int max_num_barcode_helpers) {
  // The log lines are formatted as std::cerr formats them.
  std::ostringstream index_log;
  index_log.flags(std::cerr.flags());
  index_log.precision(std::cerr.precision());
  index_log.imbue(std::cerr.getloc());

  if (mapping_parameters_.num_threads < 2) {
    reference.LoadAllSequences();
    index.Load(index_log);
    prepare_mapping(/*num_preparing_helpers=*/0);
    if (prepare_barcodes) {
      prepare_barcodes(/*num_barcode_helpers=*/0);
    }
    return index_log.str();
  }

  std::ostringstream reference_log;
  reference_log.flags(std::cerr.flags());
  reference_log.precision(std::cerr.precision());
  reference_log.imbue(std::cerr.getloc());

  // The loader thread must not call exit(), which would run the exit handlers
  // while this thread reads the index, so a truncated reference file is
  // reported after the join, as LoadAllSequences() reports it.
  bool is_reference_complete = false;
  // The threads besides the loader's own that may load the reference (a plain
  // FASTA file): none while the index loads, then those the index load used
  // but this one and those prepare_mapping runs, and all but the loader's
  // once this thread waits for it.
  std::atomic<int> num_reference_helpers(0);
  std::thread reference_loader([&reference, &reference_log,
                                &is_reference_complete,
                                &num_reference_helpers]() {
    is_reference_complete =
        reference.LoadAllSequences(reference_log, &num_reference_helpers);
  });
  // The index is read on the other mapping threads. With three threads or
  // more, one of them prepares the barcodes instead (with a whitelist, a pass
  // over the whole barcode file, longer than the index load), and with four or
  // more up to max_num_barcode_helpers more help it; with two, this thread
  // prepares them after the mapping's set-up.
  int num_index_threads = mapping_parameters_.num_threads - 1;
  std::thread barcode_preparer;
  if (prepare_barcodes && num_index_threads >= 2) {
    --num_index_threads;
    const int num_barcode_helpers =
        std::max(0, std::min(max_num_barcode_helpers, num_index_threads - 1));
    num_index_threads -= num_barcode_helpers;
    barcode_preparer = std::thread(prepare_barcodes, num_barcode_helpers);
  }
  index.Load(index_log, num_index_threads);
  // The index loads faster than the reference, so this thread, idle until the
  // join, does the mapping's set-up meanwhile, on up to
  // max_num_preparing_helpers of the threads the index load used but this
  // one; the reference loader may take the others.
  const int num_preparing_helpers =
      std::max(0, std::min(max_num_preparing_helpers, num_index_threads - 1));
  num_reference_helpers.store(num_index_threads - 1 - num_preparing_helpers);
  prepare_mapping(num_preparing_helpers);
  num_reference_helpers.store(num_index_threads - 1);
  if (barcode_preparer.joinable()) {
    barcode_preparer.join();
  } else if (prepare_barcodes) {
    prepare_barcodes(/*num_barcode_helpers=*/0);
  }
  num_reference_helpers.store(mapping_parameters_.num_threads - 1);
  reference_loader.join();

  if (!is_reference_complete) {
    ExitWithMessage(SequenceBatch::IncompleteFileMessage());
  }
  std::cerr << reference_log.str();
  return index_log.str();
}

void Chromap::RerankCandidatesRid(std::vector<Candidate> &candidates) {
  for (size_t i = 0; i < candidates.size(); ++i) {
    uint64_t rid = (uint32_t)(candidates[i].position >> 32);
    rid = custom_rid_rank_[rid];
    candidates[i].position =
        (candidates[i].position & (uint64_t)0xffffffff) | (rid << 32);
  }
}

}  // namespace chromap
