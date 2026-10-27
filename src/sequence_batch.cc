#include "sequence_batch.h"

#if defined(__linux__)
#include <errno.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#endif

#if defined(__SSSE3__)
#include <tmmintrin.h>
#elif defined(__SSE2__)
#include <emmintrin.h>
#endif

#include <algorithm>
#include <condition_variable>
#include <functional>
#include <memory>
#include <mutex>
#include <thread>
#include <tuple>

#include "utils.h"

namespace chromap {

namespace {

// The capacity the loader's buffers are given in the arena: more than a
// sequence can have under chromap's 32-bit lengths, so kseq_read never grows
// (reallocs) a buffer there. A sequence is read into the arena only while this
// much of each half is free; later ones are allocated one by one.
constexpr size_t kArenaBufferBytes = ((size_t)1 << 32) + 64;
// Each sequence and quality string in the arena starts on a cache line.
constexpr size_t kArenaAlignment = 64;

size_t ArenaBytes(size_t length) {
  // The string and its terminating null.
  return (length + kArenaAlignment) / kArenaAlignment * kArenaAlignment;
}

void PointAt(kstring_t &buffer, char *start) {
  buffer.s = start;
  buffer.m = kArenaBufferBytes;
  buffer.l = 0;
}

void Clear(kstring_t &buffer) {
  buffer.s = NULL;
  buffer.m = 0;
  buffer.l = 0;
}

#if defined(__linux__) && defined(MADV_HUGEPAGE)

// The parallel FASTA loader reads the file in pieces of this many bytes, each
// into a buffer of its thread that stays in the core's cache.
constexpr size_t kPieceBytes = (size_t)1 << 20;

// The longest header line it reads; a longer one sends the load to kseq.
constexpr size_t kMaxHeaderLineBytes = (size_t)1 << 26;

// Reads 'num_bytes' bytes at 'offset' of the file into 'buffer'. Returns
// false on an error or when the file ends first.
bool ReadAt(int fd, char *buffer, size_t num_bytes, uint64_t offset) {
  while (num_bytes > 0) {
    const ssize_t num_read = pread(fd, buffer, num_bytes, (off_t)offset);
    if (num_read < 0 && errno == EINTR) {
      continue;
    }
    if (num_read <= 0) {
      return false;
    }
    buffer += num_read;
    num_bytes -= num_read;
    offset += num_read;
  }
  return true;
}

// What the first pass finds in one piece of the file.
struct PieceScan {
  uint64_t num_newlines = 0;
  unsigned char first_byte = 0;
  unsigned char last_byte = 0;
  // A byte kseq treats in a way this loader does not: '\r', '\v' or '\f'
  // anywhere, or a '+' or '@' that starts a line (but the piece's first byte,
  // which the caller checks against the previous piece's last byte).
  bool needs_kseq = false;
  // Each '>' that starts a line, but the piece's first byte: its offset in the
  // piece and the number of '\n's before it in the piece.
  std::vector<std::pair<size_t, uint64_t>> headers;
};

// Scans bytes [begin, end) of a piece one by one, counting its '\n's in
// 'num_newlines'.
void ScanBytes(const unsigned char *piece, size_t begin, size_t end,
               uint64_t &num_newlines, PieceScan &scan) {
  for (size_t i = begin; i < end; ++i) {
    const unsigned char c = piece[i];
    if (c == '\n') {
      ++num_newlines;
    } else if (c == '\r' || c == '\v' || c == '\f') {
      scan.needs_kseq = true;
    } else if ((c == '>' || c == '+' || c == '@') && i > 0 &&
               piece[i - 1] == '\n') {
      if (c == '>') {
        scan.headers.emplace_back(i, num_newlines);
      } else {
        scan.needs_kseq = true;
      }
    }
  }
}

#if defined(__SSE2__)
uint64_t SumBytes(__m128i counts) {
  const __m128i sums = _mm_sad_epu8(counts, _mm_setzero_si128());
  return (uint64_t)_mm_cvtsi128_si64(sums) +
         (uint64_t)_mm_cvtsi128_si64(_mm_unpackhi_epi64(sums, sums));
}
#endif

// Counts the '\n's of a piece and finds its headers and the bytes that send
// the load to kseq. 16-byte blocks without '>', '+', '@', '\r', '\v' or '\f'
// only add to the count; the others are scanned byte by byte.
void ScanPiece(const unsigned char *piece, size_t num_bytes,
               PieceScan &scan) {
  scan.first_byte = piece[0];
  scan.last_byte = piece[num_bytes - 1];
  uint64_t num_newlines = 0;
  size_t i = 0;
#if defined(__SSE2__)
  const __m128i newline = _mm_set1_epi8('\n');
  const __m128i greater = _mm_set1_epi8('>');
  const __m128i plus = _mm_set1_epi8('+');
  const __m128i at = _mm_set1_epi8('@');
  const __m128i vertical_tab = _mm_set1_epi8('\v');
  const __m128i two = _mm_set1_epi8(2);
  // Per-byte counts of '\n', summed before any can pass 255.
  __m128i counts = _mm_setzero_si128();
  int num_blocks = 0;
  for (; i + 16 <= num_bytes; i += 16) {
    const __m128i x = _mm_loadu_si128((const __m128i *)(piece + i));
    // '\v', '\f' and '\r' are the three bytes from '\v' on.
    const __m128i from_vertical_tab = _mm_sub_epi8(x, vertical_tab);
    const __m128i is_carriage =
        _mm_cmpeq_epi8(_mm_min_epu8(from_vertical_tab, two), from_vertical_tab);
    const __m128i is_mark = _mm_or_si128(
        _mm_or_si128(_mm_cmpeq_epi8(x, greater), _mm_cmpeq_epi8(x, plus)),
        _mm_cmpeq_epi8(x, at));
    if (_mm_movemask_epi8(_mm_or_si128(is_mark, is_carriage)) != 0) {
      num_newlines += SumBytes(counts);
      counts = _mm_setzero_si128();
      num_blocks = 0;
      ScanBytes(piece, i, i + 16, num_newlines, scan);
      continue;
    }
    counts = _mm_sub_epi8(counts, _mm_cmpeq_epi8(x, newline));
    if (++num_blocks == 255) {
      num_newlines += SumBytes(counts);
      counts = _mm_setzero_si128();
      num_blocks = 0;
    }
  }
  num_newlines += SumBytes(counts);
#endif
  ScanBytes(piece, i, num_bytes, num_newlines, scan);
  scan.num_newlines = num_newlines;
}

// Copies the 'num_bytes' bytes at 'source' but the '\n's to 'destination',
// writing at most 'capacity' bytes there. Returns the number of bytes that are
// not '\n', which is more than 'capacity' when the source holds more.
size_t CopyWithoutNewlines(const char *source, size_t num_bytes,
                           char *destination, size_t capacity) {
  size_t i = 0;
  size_t num_copied = 0;
#if defined(__SSSE3__)
  const __m128i newline = _mm_set1_epi8('\n');
  const __m128i indexes =
      _mm_setr_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);
  // Each block stores 16 bytes, so only while 16 fit.
  for (; i + 16 <= num_bytes && num_copied + 16 <= capacity; i += 16) {
    const __m128i x = _mm_loadu_si128((const __m128i *)(source + i));
    const unsigned newlines = _mm_movemask_epi8(_mm_cmpeq_epi8(x, newline));
    if (newlines == 0) {
      _mm_storeu_si128((__m128i *)(destination + num_copied), x);
      num_copied += 16;
    } else if ((newlines & (newlines - 1)) == 0) {
      // One '\n', at k: byte j of the result is byte j + (j >= k) of the
      // block. Its 16th byte is overwritten by the next block or lies past
      // the copy's end, which is within 'capacity'.
      const int k = __builtin_ctz(newlines);
      const __m128i shuffle = _mm_sub_epi8(
          indexes, _mm_cmpgt_epi8(indexes, _mm_set1_epi8((char)(k - 1))));
      _mm_storeu_si128((__m128i *)(destination + num_copied),
                       _mm_shuffle_epi8(x, shuffle));
      num_copied += 15;
    } else {
      for (size_t j = i; j < i + 16; ++j) {
        if (source[j] != '\n') {
          destination[num_copied++] = source[j];
        }
      }
    }
  }
#endif
  for (; i < num_bytes; ++i) {
    if (source[i] != '\n') {
      if (num_copied < capacity) {
        destination[num_copied] = source[i];
      }
      ++num_copied;
    }
  }
  return num_copied;
}

// Runs jobs over numbered pieces on the calling thread and on helper threads,
// which it starts while a job runs, as many as *num_helper_threads allows at
// that moment (none when it is null). Each thread has its own piece buffer.
class PieceWorkers {
 public:
  typedef std::function<void(size_t, char *)> Work;

  explicit PieceWorkers(const std::atomic<int> *num_helper_threads)
      : num_helper_threads_(num_helper_threads),
        buffer_(new char[kPieceBytes]) {}

  ~PieceWorkers() {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stop_ = true;
    }
    has_job_.notify_all();
    for (std::thread &helper : helpers_) {
      helper.join();
    }
  }

  // Calls work(piece, buffer) once for each piece in [0, num_pieces) and
  // returns when every call has returned.
  void Run(size_t num_pieces, const Work &work) {
    jobs_.emplace_back(new Job(work, num_pieces));
    Job &job = *jobs_.back();
    {
      std::lock_guard<std::mutex> lock(mutex_);
      current_job_ = &job;
      ++job_number_;
    }
    has_job_.notify_all();
    for (;;) {
      StartHelpers();
      const size_t piece = job.next_piece.fetch_add(1);
      if (piece >= num_pieces) {
        break;
      }
      work(piece, buffer_.get());
      job.num_done.fetch_add(1, std::memory_order_release);
    }
    while (job.num_done.load(std::memory_order_acquire) < num_pieces) {
      std::this_thread::yield();
    }
    std::lock_guard<std::mutex> lock(mutex_);
    current_job_ = nullptr;
  }

 private:
  struct Job {
    Job(const Work &work, size_t num_pieces)
        : work(work), num_pieces(num_pieces) {}
    const Work &work;
    const size_t num_pieces;
    std::atomic<size_t> next_piece{0};
    std::atomic<size_t> num_done{0};
  };

  void StartHelpers() {
    if (num_helper_threads_ == nullptr) {
      return;
    }
    const int num_allowed = num_helper_threads_->load(std::memory_order_relaxed);
    while ((int)helpers_.size() < num_allowed) {
      try {
        helpers_.emplace_back(&PieceWorkers::Help, this);
      } catch (...) {
        num_helper_threads_ = nullptr;
        return;
      }
    }
  }

  void Help() {
    std::unique_ptr<char[]> buffer(new char[kPieceBytes]);
    uint64_t last_job_number = 0;
    for (;;) {
      Job *job = nullptr;
      {
        std::unique_lock<std::mutex> lock(mutex_);
        has_job_.wait(lock, [this, last_job_number]() {
          return stop_ ||
                 (current_job_ != nullptr && job_number_ != last_job_number);
        });
        if (stop_) {
          return;
        }
        job = current_job_;
        last_job_number = job_number_;
      }
      // A job stays allocated until the workers stop, so a helper that comes
      // to it after its last piece only finds no piece left.
      for (;;) {
        const size_t piece = job->next_piece.fetch_add(1);
        if (piece >= job->num_pieces) {
          break;
        }
        job->work(piece, buffer.get());
        job->num_done.fetch_add(1, std::memory_order_release);
      }
    }
  }

  const std::atomic<int> *num_helper_threads_;
  std::unique_ptr<char[]> buffer_;
  std::vector<std::unique_ptr<Job>> jobs_;
  std::vector<std::thread> helpers_;
  std::mutex mutex_;
  std::condition_variable has_job_;
  Job *current_job_ = nullptr;
  uint64_t job_number_ = 0;
  bool stop_ = false;
};

// Reads the line at 'offset' of the file, without its '\n', into 'line', and
// whether a '\n' ends it rather than the file. Returns false on a read error
// or a line longer than kMaxHeaderLineBytes.
bool ReadLineAt(int fd, uint64_t offset, uint64_t file_size, std::string &line,
                bool &has_newline) {
  line.clear();
  has_newline = false;
  char chunk[4096];
  while (offset < file_size) {
    const size_t num_bytes =
        (size_t)std::min<uint64_t>(sizeof(chunk), file_size - offset);
    if (!ReadAt(fd, chunk, num_bytes, offset)) {
      return false;
    }
    const char *newline = (const char *)memchr(chunk, '\n', num_bytes);
    if (newline != nullptr) {
      line.append(chunk, newline - chunk);
      has_newline = true;
      return true;
    }
    line.append(chunk, num_bytes);
    if (line.size() > kMaxHeaderLineBytes) {
      return false;
    }
    offset += num_bytes;
  }
  return true;
}

// Closes a file descriptor when it leaves scope.
struct FileCloser {
  int fd;
  ~FileCloser() {
    if (fd >= 0) {
      close(fd);
    }
  }
};

// A FASTA record as kseq_read reads it: its header line at 'header' and its
// sequence, the bytes from 'begin' to 'end' but the '\n's, of which there
// are 'newlines_before_begin' and 'newlines_before_end' before those offsets.
struct FastaRecord {
  uint64_t header = 0;
  uint64_t newlines_before_header = 0;
  uint64_t begin = 0;
  uint64_t end = 0;
  uint64_t newlines_before_begin = 0;
  uint64_t newlines_before_end = 0;
  uint64_t length = 0;
  // Where the sequence starts in the arena's first half.
  size_t arena_offset = 0;
  std::string name;
  std::string comment;
  bool has_comment = false;
};

void CopyToKstring(const std::string &text, kstring_t &kstring) {
  size_t capacity = text.size() + 1;
  kroundup32(capacity);
  kstring.s = (char *)malloc(capacity);
  kstring.m = capacity;
  kstring.l = text.size();
  memcpy(kstring.s, text.data(), text.size());
  kstring.s[text.size()] = '\0';
}

#endif  // defined(__linux__) && defined(MADV_HUGEPAGE)

}  // namespace

void SequenceBatch::InitializeLoading(const std::string &sequence_file_path,
                                      bool default_buffers) {
  sequence_file_ = gzopen(sequence_file_path.c_str(), "r");
  if (sequence_file_ == NULL) {
    ExitWithMessage("Cannot find sequence file " + sequence_file_path);
  }
  use_arena_ = !default_buffers;
  if (default_buffers) {
    sequence_kseq_ = kseq_init(sequence_file_);
    return;
  }
  // Only here: a copy of a long path is one more heap allocation, and index
  // construction keeps the heap's layout (see default_buffers).
  sequence_file_path_ = sequence_file_path;
  // gzbuffer takes effect only before the first read allocates zlib's buffers.
  gzbuffer(sequence_file_, kGzBufferSize);
  sequence_kseq_ = kseq_init2(sequence_file_, kKseqBufferSize);
}

void SequenceBatch::FinalizeLoading() {
  kseq_destroy(sequence_kseq_);
  gzclose(sequence_file_);
}

bool SequenceBatch::LoadOneSequenceAndSaveAt(uint32_t sequence_index) {
  if (sequence_index == 0) {
    num_loaded_sequences_ = 0;
  }

  int length = kseq_read(sequence_kseq_);
  while (length == 0) {
    length = kseq_read(sequence_kseq_);
  }

  if (length > 0) {
    kseq_t *sequence = sequence_batch_[sequence_index];
    std::swap(sequence_kseq_->seq, sequence->seq);
    ReplaceByEffectiveRange(sequence->seq, /*is_seq=*/true);
    std::swap(sequence_kseq_->name, sequence->name);
    std::swap(sequence_kseq_->comment, sequence->comment);
    sequence->id = total_num_loaded_sequences_;
    ++total_num_loaded_sequences_;

    if (sequence_index >= num_loaded_sequences_) {
      ++num_loaded_sequences_;
    } else if (sequence_index + 1 != num_loaded_sequences_) {
      std::cerr << sequence_index << " " << num_loaded_sequences_ << "\n";
      ExitWithMessage(
          "Shouldn't override other sequences rather than the last!");
    }

    if (sequence_kseq_->qual.l != 0) {  // fastq file
      std::swap(sequence_kseq_->qual, sequence->qual);
      ReplaceByEffectiveRange(sequence->qual, /*is_seq=*/false);
    }
    return false;
  }

  // Make sure to reach the end of the file rather than meet an error.
  if (length != -1) {
    ExitWithMessage(
        "Didn't reach the end of sequence file, which might be corrupted!");
  }
  return true;
}

uint32_t SequenceBatch::LoadBatch() { return LoadBatch(std::cerr); }

uint32_t SequenceBatch::LoadBatch(std::ostream &log) {
  double real_start_time = GetRealTime();
  num_loaded_sequences_ = 0;
  for (uint32_t sequence_index = 0; sequence_index < max_num_sequences_;
       ++sequence_index) {
    if (LoadOneSequenceAndSaveAt(sequence_index)) {
      break;
    }
  }

  if (num_loaded_sequences_ != 0) {
    log << "Loaded sequence batch successfully in "
        << GetRealTime() - real_start_time << "s, ";
    log << "number of sequences: " << num_loaded_sequences_ << ".\n";
  } else {
    log << "No more sequences.\n";
  }
  return num_loaded_sequences_;
}

void SequenceBatch::LoadAllSequences() {
  if (!LoadAllSequences(std::cerr)) {
    ExitWithMessage(IncompleteFileMessage());
  }
}

bool SequenceBatch::LoadAllSequences(
    std::ostream &log, const std::atomic<int> *num_helper_threads) {
  double real_start_time = GetRealTime();
  sequence_batch_.reserve(200);
  num_loaded_sequences_ = 0;
  num_bases_ = 0;
  if (use_arena_) {
    ReserveArena();
  }
  if (!LoadFastaInParallel(num_helper_threads)) {
    PointLoaderAtArena();
    int length = kseq_read(sequence_kseq_);
    while (length >= 0) {
      if (length > 0) {
        sequence_batch_.emplace_back((kseq_t *)calloc(1, sizeof(kseq_t)));
        kseq_t *sequence = sequence_batch_.back();
        if (IsInArena(sequence_kseq_->seq.s)) {
          arena_bases_used_ += ArenaBytes(sequence_kseq_->seq.l);
          if (sequence_kseq_->qual.l != 0) {
            arena_quals_used_ += ArenaBytes(sequence_kseq_->qual.l);
          }
        }
        std::swap(sequence_kseq_->seq, sequence->seq);
        ReplaceByEffectiveRange(sequence->seq, /*is_seq=*/true);
        std::swap(sequence_kseq_->name, sequence->name);
        std::swap(sequence_kseq_->comment, sequence->comment);
        if (sequence_kseq_->qual.l != 0) {  // fastq file
          std::swap(sequence_kseq_->qual, sequence->qual);
          ReplaceByEffectiveRange(sequence->qual, /*is_seq=*/false);
        }
        sequence->id = total_num_loaded_sequences_;
        ++total_num_loaded_sequences_;
        ++num_loaded_sequences_;
        num_bases_ += length;
        PointLoaderAtArena();
      }
      length = kseq_read(sequence_kseq_);
    }
    // FinalizeLoading frees the loader's buffers (kseq_destroy), so they
    // leave the arena.
    ReleaseLoaderFromArena();

    // Make sure to reach the end of the file rather than meet an error.
    if (length != -1) {
      return false;
    }
  }

  log << "Loaded all sequences successfully in "
      << GetRealTime() - real_start_time << "s, ";
  log << "number of sequences: " << num_loaded_sequences_ << ", ";
  log << "number of bases: " << num_bases_ << ".\n";
  return true;
}

void SequenceBatch::ReplaceByEffectiveRange(kstring_t &seq, bool is_seq) {
  seq.l = effective_range_.Replace(seq.s, seq.l, is_seq);
}

bool SequenceBatch::LoadFastaInParallel(
    const std::atomic<int> *num_helper_threads) {
#if defined(__linux__) && defined(MADV_HUGEPAGE)
  // Only with helper threads to come (on one thread the two passes take
  // longer than kseq_read: 0.98 s against 0.82-0.87 s for hg38), into the
  // arena of a new batch, where kseq_read would put the sequences, and before
  // anything has been read from the file.
  if (num_helper_threads == nullptr || arena_ == nullptr ||
      arena_bases_used_ != 0 || arena_quals_used_ != 0 ||
      !sequence_batch_.empty() || sequence_kseq_ == nullptr ||
      sequence_kseq_->last_char != 0 || sequence_kseq_->f->begin != 0 ||
      sequence_kseq_->f->end != 0 || sequence_kseq_->f->is_eof != 0) {
    return false;
  }
  FileCloser file{open(sequence_file_path_.c_str(), O_RDONLY | O_CLOEXEC)};
  struct stat file_stat;
  if (file.fd < 0 || fstat(file.fd, &file_stat) != 0 ||
      !S_ISREG(file_stat.st_mode) || file_stat.st_size <= 0) {
    return false;
  }
  const uint64_t file_size = (uint64_t)file_stat.st_size;
  // kseq_read skips the bytes before the first '>' or '@'; a gzip file starts
  // with 0x1f 0x8b, which zlib inflates.
  char first_byte = 0;
  if (!ReadAt(file.fd, &first_byte, 1, 0) || first_byte != '>') {
    return false;
  }

  const size_t num_pieces = (file_size + kPieceBytes - 1) / kPieceBytes;
  const auto piece_size = [file_size](size_t piece) {
    return (size_t)std::min<uint64_t>(kPieceBytes,
                                      file_size - (uint64_t)piece * kPieceBytes);
  };
  PieceWorkers workers(num_helper_threads);

  // The first pass: '\n' counts, line-start '>'s and the bytes kseq is needed
  // for, piece by piece.
  std::vector<PieceScan> scans(num_pieces);
  std::atomic<bool> needs_kseq(false);
  workers.Run(num_pieces, [&](size_t piece, char *buffer) {
    if (needs_kseq.load(std::memory_order_relaxed)) {
      return;
    }
    const size_t num_bytes = piece_size(piece);
    if (!ReadAt(file.fd, buffer, num_bytes, (uint64_t)piece * kPieceBytes)) {
      needs_kseq = true;
      return;
    }
    ScanPiece((const unsigned char *)buffer, num_bytes, scans[piece]);
    if (scans[piece].needs_kseq) {
      needs_kseq = true;
    }
  });
  if (needs_kseq) {
    return false;
  }

  // The records: each line-start '>' begins one, which kseq_read reads as the
  // header line and every byte up to the next one but the '\n's.
  std::vector<uint64_t> newlines_before_piece(num_pieces + 1, 0);
  for (size_t piece = 0; piece < num_pieces; ++piece) {
    newlines_before_piece[piece + 1] =
        newlines_before_piece[piece] + scans[piece].num_newlines;
  }
  std::vector<FastaRecord> records;
  for (size_t piece = 0; piece < num_pieces; ++piece) {
    const PieceScan &scan = scans[piece];
    const uint64_t piece_start = (uint64_t)piece * kPieceBytes;
    const bool starts_line = piece == 0 || scans[piece - 1].last_byte == '\n';
    if (starts_line && (scan.first_byte == '>' || scan.first_byte == '+' ||
                        scan.first_byte == '@')) {
      if (scan.first_byte != '>') {
        return false;
      }
      records.emplace_back();
      records.back().header = piece_start;
      records.back().newlines_before_header = newlines_before_piece[piece];
    }
    for (const std::pair<size_t, uint64_t> &header : scan.headers) {
      records.emplace_back();
      records.back().header = piece_start + header.first;
      records.back().newlines_before_header =
          newlines_before_piece[piece] + header.second;
    }
  }
  std::vector<PieceScan>().swap(scans);

  // Each record's name, comment and length, read as kseq_read reads them, and
  // its place in the arena, as LoadAllSequences gives it.
  const size_t arena_half = arena_size_ / 2;
  size_t arena_bytes = 0;
  std::vector<const FastaRecord *> kept_records;
  std::string line;
  for (size_t ri = 0; ri < records.size(); ++ri) {
    FastaRecord &record = records[ri];
    bool has_newline = false;
    if (!ReadLineAt(file.fd, record.header, file_size, line, has_newline)) {
      return false;
    }
    const bool is_last = ri + 1 == records.size();
    record.end = is_last ? file_size : records[ri + 1].header;
    record.newlines_before_end = is_last
                                     ? newlines_before_piece[num_pieces]
                                     : records[ri + 1].newlines_before_header;
    if (!has_newline) {
      // The header line ends the file: no bases.
      continue;
    }
    record.begin = record.header + line.size() + 1;
    record.newlines_before_begin = record.newlines_before_header + 1;
    if (record.begin > record.end ||
        record.newlines_before_begin > record.newlines_before_end ||
        record.newlines_before_end - record.newlines_before_begin >
            record.end - record.begin) {
      return false;
    }
    record.length = (record.end - record.begin) -
                    (record.newlines_before_end - record.newlines_before_begin);
    if (record.length == 0) {
      continue;
    }
    // kseq_read returns the length as an int, and LoadAllSequences puts a
    // sequence in the arena only while kArenaBufferBytes of it are free.
    if (record.length > (uint64_t)INT32_MAX ||
        arena_half - arena_bytes < kArenaBufferBytes) {
      return false;
    }
    record.arena_offset = arena_bytes;
    arena_bytes += ArenaBytes(record.length);
    // The name runs to the first white space; a comment follows unless that
    // is the line's end.
    size_t name_end = 1;
    while (name_end < line.size() &&
           !isspace((unsigned char)line[name_end])) {
      ++name_end;
    }
    record.name.assign(line, 1, name_end - 1);
    if (name_end < line.size()) {
      record.has_comment = true;
      record.comment.assign(line, name_end + 1, std::string::npos);
    }
    kept_records.push_back(&record);
  }

  // The second pass: each piece's bases copied to their place.
  char *const bases = arena_;
  std::atomic<bool> has_failed(false);
  workers.Run(num_pieces, [&](size_t piece, char *buffer) {
    if (has_failed.load(std::memory_order_relaxed)) {
      return;
    }
    const uint64_t piece_start = (uint64_t)piece * kPieceBytes;
    const size_t num_bytes = piece_size(piece);
    const uint64_t piece_end = piece_start + num_bytes;
    bool is_read = false;
    // From the first record that ends after the piece starts.
    for (auto it = std::upper_bound(kept_records.begin(), kept_records.end(),
                                    piece_start,
                                    [](uint64_t offset, const FastaRecord *r) {
                                      return offset < r->end;
                                    });
         it != kept_records.end() && (*it)->begin < piece_end; ++it) {
      const FastaRecord &record = **it;
      const uint64_t begin = std::max(record.begin, piece_start);
      const uint64_t end = std::min(record.end, piece_end);
      if (begin >= end) {
        continue;
      }
      if (!is_read) {
        if (!ReadAt(file.fd, buffer, num_bytes, piece_start)) {
          has_failed = true;
          return;
        }
        is_read = true;
      }
      const uint64_t newlines_before_begin =
          begin == record.begin ? record.newlines_before_begin
                                : newlines_before_piece[piece];
      const uint64_t newlines_before_end = end == record.end
                                               ? record.newlines_before_end
                                               : newlines_before_piece[piece + 1];
      const uint64_t position =
          (begin - record.begin) -
          (newlines_before_begin - record.newlines_before_begin);
      const uint64_t num_bases =
          (end - begin) - (newlines_before_end - newlines_before_begin);
      if (CopyWithoutNewlines(buffer + (begin - piece_start), end - begin,
                              bases + record.arena_offset + position,
                              num_bases) != num_bases) {
        has_failed = true;
        return;
      }
    }
  });
  if (has_failed) {
    return false;
  }

  for (const FastaRecord *record : kept_records) {
    sequence_batch_.emplace_back((kseq_t *)calloc(1, sizeof(kseq_t)));
    kseq_t *sequence = sequence_batch_.back();
    sequence->seq.s = bases + record->arena_offset;
    sequence->seq.s[record->length] = '\0';
    sequence->seq.l = record->length;
    sequence->seq.m = kArenaBufferBytes;
    ReplaceByEffectiveRange(sequence->seq, /*is_seq=*/true);
    CopyToKstring(record->name, sequence->name);
    if (record->has_comment) {
      CopyToKstring(record->comment, sequence->comment);
    }
    sequence->id = total_num_loaded_sequences_;
    ++total_num_loaded_sequences_;
    ++num_loaded_sequences_;
    num_bases_ += record->length;
  }
  arena_bases_used_ = arena_bytes;
  return true;
#else
  (void)num_helper_threads;
  return false;
#endif
}

void SequenceBatch::ReserveArena() {
#if defined(__linux__) && defined(MADV_HUGEPAGE)
  if (arena_ != nullptr) {
    return;
  }
  // Each half takes 64 GiB of address space, of which only the pages touched
  // are backed (MAP_NORESERVE).
  const size_t arena_size = (size_t)128 << 30;
  void *arena = mmap(NULL, arena_size, PROT_READ | PROT_WRITE,
                     MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
  if (arena == MAP_FAILED) {
    return;
  }
  // A hint: where it fails, the arena is backed by 4 KiB pages.
  madvise(arena, arena_size, MADV_HUGEPAGE);
  arena_ = static_cast<char *>(arena);
  arena_size_ = arena_size;
#endif
}

void SequenceBatch::PointLoaderAtArena() {
  if (arena_ == nullptr) {
    return;
  }
  kstring_t &seq = sequence_kseq_->seq;
  kstring_t &qual = sequence_kseq_->qual;
  // kseq_read grows the quality buffer to the sequence buffer's capacity, so
  // the two are in the arena together or not at all. A buffer kseq_read
  // allocated stays with the loader.
  const bool is_free = (seq.s == NULL || IsInArena(seq.s)) &&
                       (qual.s == NULL || IsInArena(qual.s));
  const size_t half = arena_size_ / 2;
  if (is_free && half - arena_bases_used_ >= kArenaBufferBytes &&
      half - arena_quals_used_ >= kArenaBufferBytes) {
    PointAt(seq, arena_ + arena_bases_used_);
    PointAt(qual, arena_ + half + arena_quals_used_);
  } else {
    ReleaseLoaderFromArena();
  }
}

void SequenceBatch::ReleaseLoaderFromArena() {
  if (sequence_kseq_ == nullptr) {
    return;
  }
  if (IsInArena(sequence_kseq_->seq.s)) {
    Clear(sequence_kseq_->seq);
  }
  if (IsInArena(sequence_kseq_->qual.s)) {
    Clear(sequence_kseq_->qual);
  }
}

void SequenceBatch::ReleaseArena() {
#if defined(__linux__) && defined(MADV_HUGEPAGE)
  if (arena_ == nullptr) {
    return;
  }
  // The loader's buffers left the arena when LoadAllSequences returned, and
  // FinalizeLoading may have freed the loader since.
  for (kseq_t *sequence : sequence_batch_) {
    if (IsInArena(sequence->seq.s)) {
      Clear(sequence->seq);
    }
    if (IsInArena(sequence->qual.s)) {
      Clear(sequence->qual);
    }
  }
  munmap(arena_, arena_size_);
  arena_ = nullptr;
  arena_size_ = 0;
  arena_bases_used_ = 0;
  arena_quals_used_ = 0;
#endif
}

}  // namespace chromap
