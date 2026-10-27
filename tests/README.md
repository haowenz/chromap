# chromap unit tests

A pytest suite for chromap in two layers. CLI tests run the chromap binary
on a small synthetic genome and simulated reads whose true origin is
recorded, and compare every output (BED, fragments, TagAlign, SAM, pairs,
PAF, summary CSV, stderr statistics) with values derived from that truth or
worked by hand. C++ tests reach chromap's public classes and functions
directly through small drivers in `tests/cpp/`, compiled with the
repository's own compiler and flags against the objects of `src/`. The suite
needs Python >= 3.8, pytest, `make` and the C++ toolchain that builds
chromap; samtools and pysam are optional and only cross-check the harness's
SAM parser.

## Run

```bash
python3 -m pytest -q tests                            # builds chromap from src/ into a cache first
CHROMAP_BIN=$PWD/chromap python3 -m pytest -q tests   # use a chromap already built (e.g. by `make`)
CHROMAP_TEST_THREADS=4 python3 -m pytest -q tests     # -t for every chromap run (default 1)
python3 -m pytest -q tests -m "not lowmem"            # skip the four 12 GB --low-mem runs on small machines
python3 -m pytest -q tests -m "not cpp"               # only the command-line tests (no C++ drivers built)
python3 -m pytest -q tests/test_c_*.py                # one area
```

To test a modified copy of the sources (a mutant, another branch's `src/`)
with this checkout's tests, build that copy and point both layers at it:

```bash
make -f Makefile src_dir=/path/to/copy/src objs_dir=/tmp/m/objs exec=/tmp/m/chromap -j8
CHROMAP_BIN=/tmp/m/chromap CHROMAP_TEST_SRC_DIR=/path/to/copy/src python3 -m pytest -q tests
```

### How the binary and the C++ drivers are built

- `CHROMAP_BIN` set: the CLI tests use that file as is (it must be an
  executable). Use it when chromap is already built, for example by `make`
  in a CI job, so that the suite does not build it again.
- Otherwise the suite builds chromap from this checkout with the Makefile
  itself (`make -f Makefile src_dir=src objs_dir=<cache>/objs
  exec=<cache>/chromap`), into a cache directory keyed by a hash of `src/`,
  the Makefile and the toolchain: a source change always rebuilds, an
  unchanged tree never does, and nothing is written into the repository.
- The C++ drivers are compiled from this checkout's `src/`, or from
  `$CHROMAP_TEST_SRC_DIR` when it is set: every `.cc` there except
  `chromap_driver.cc` (which holds `main` and takes ~15 s) is compiled once
  into the cache, and each `tests/cpp/<name>.cc` is linked against those
  objects with the Makefile's link line plus `-I <src> -I tests/cpp`.
  `CHROMAP_TEST_SRC_DIR` also sets the sources of the cached chromap build
  (when `CHROMAP_BIN` is unset) and of the version `test_d_cli.py` expects;
  the compiler flags still come from this checkout's Makefile.
- The compiler and flags are the Makefile's: `make` prints its `CXX`,
  `CXXFLAGS` and `LDFLAGS`, with `CHROMAP_TEST_CXX`, `CHROMAP_TEST_CXXFLAGS`
  and `CHROMAP_TEST_LDFLAGS` passed on make's command line when set, exactly
  like `make CXX=...` (the environment's `CXX` is ignored, as the Makefile
  ignores it).

| Variable | Default | Meaning |
|---|---|---|
| `CHROMAP_BIN` | unset | chromap binary for the CLI tests |
| `CHROMAP_TEST_SRC_DIR` | `<checkout>/src` | sources the C++ drivers (and a cached chromap build) are compiled from |
| `CHROMAP_TEST_THREADS` | 1 | `-t` (and `OMP_NUM_THREADS`) of every run that does not set `threads=`; 1 to 8 |
| `CHROMAP_TEST_BUILD_DIR` | `<tempdir>/chromap-tests-<user>` | build cache |
| `CHROMAP_TEST_BUILD_JOBS` | 4 | parallel compile jobs |
| `CHROMAP_TEST_CXX`, `CHROMAP_TEST_CXXFLAGS`, `CHROMAP_TEST_LDFLAGS` | the Makefile's | make overrides, e.g. `CHROMAP_TEST_CXX=clang++` |

The repository's CI builds with `make CXX=g++|clang++` on Ubuntu and with
`make CXX=clang++ CXXFLAGS="... -Xclang -fopenmp ..." LDFLAGS="... -lomp"` on
macOS. To run the suite there after `make`, pass the same values:

```bash
make CXX=clang++
CHROMAP_BIN=$PWD/chromap CHROMAP_TEST_CXX=clang++ python3 -m pytest -q tests
# macOS: also CHROMAP_TEST_CXXFLAGS="<the CXXFLAGS given to make>" CHROMAP_TEST_LDFLAGS="<the LDFLAGS>"
```

A clang++ that is not the system compiler is passed the same way. The checks
below used clang++ 18 from a conda-forge environment (`clangxx`, `llvm-openmp`)
against the system's GCC 13 libstdc++, with the environment's libomp on the
link line and in the binary's rpath (the C++ drivers and the cached chromap
build both take these values):

```bash
CLANG=<env>/bin/clang++ LIB=<env>/lib
CHROMAP_TEST_CXX="$CLANG --no-default-config --gcc-install-dir=/usr/lib/gcc/x86_64-linux-gnu/13" \
CHROMAP_TEST_LDFLAGS="-lm -lz -L$LIB -Wl,-rpath,$LIB" python3 -m pytest -q tests
```

`--no-default-config` keeps the environment's clang configuration file from
adding its own sysroot and flags; `--gcc-install-dir` picks the C++ standard
library. clang/libomp builds undercount chromap's per-thread statistics on
stderr at `-t` > 1 (known bug a-mapping-statistics-clang-threads), so tests
assert those counters only on runs pinned to `threads=1`.

### Checked configurations

The whole suite at chromap 949043c on Linux x86-64, all with `LC_ALL=C`:
g++ 13.3 (the Makefile's flags) at `CHROMAP_TEST_THREADS` 1 and 4 (three
runs at 4, one of them under CPU contention: 8 cores shared with 16 busy
processes); g++ 11.5 at 4; clang++ 18.1 with LLVM's libomp (conda-forge,
against the system libstdc++, see above) at 1 and 4; Python 3.8.20 (pytest
8.3.5, without pysam: its test skips) and 3.13 at 4; `CHROMAP_BIN` pointing at
a binary built once with `make`, at 4. The
harness and pilot tests also ran on Python 3.14 and without samtools and
pysam. macOS was not run.

## Runtime and memory

The whole suite (432 tests) took 85 to 96 s at `CHROMAP_TEST_THREADS=1`
and at 4 with a warm build cache on a shared 256-core server under other
load; about 13 s of it are the four `lowmem` tests
and about 5 s the three cache runs of `test_b_cache.py`. A cold cache adds
about 25 s of compiling at the default 4 build jobs. Areas C and D alone
(232 tests) take about 22 s.
Costs at 949043c that shape how tests are written:

| Step | Wall time | Peak memory |
|---|---|---|
| index of the 100 kb synthetic genome | ~0.01 s | small |
| a single-end mapping run | ~0.3 s | ~1.1 GB |
| a paired-end mapping run | ~0.5 s | ~2.0 GB |
| `--low-mem` with SAM, pairs or PAF output (also `--preset hic`) | ~0.8 s | ~3 GB |
| `--low-mem` with BED or TagAlign output (also `--preset atac`, `--preset chip`) | ~3 s | ~12 GB |
| 510,000 reads or pairs (two batches, so the candidate cache is used; `test_b_cache.py`) | 1-3 s | up to ~2 GB |
| one C++ driver call | milliseconds | small |

Most of the fixed cost of a mapping run is chromap allocating its candidate
cache (4,000,003 entries for paired-end, 2,000,003 for single-end) whatever
the input size: `--cache-size 2000000` took a paired-end run from 0.52 s and
1.96 GB to 0.33 s and 1.29 GB. `--low-mem` reads its temp files back through
buffers totalling 10 GB (BED, TagAlign) or 1 GB (SAM, pairs, PAF), sized with
`resize` and so constructed record by record
(`ProcessAndOutputMappingsInLowMemory`). Tests that use BED or
TagAlign with `--low-mem` carry `@pytest.mark.lowmem`, so a small machine can
deselect them; keep them few (four at this writing, two in area B and two
in area C). The whole suite should stay about a minute at
`CHROMAP_TEST_THREADS=4`, so that it can run on every change: share runs
(module-scoped fixtures, `shared=True`) and give each run many designed
reads.

Keep pytest's base temporary directory on a local disk (the default,
`/tmp/pytest-of-$USER`); chromap runs in its run directory.

## Coverage

Line coverage of `src/` (the vendored `cxxopts.hpp`, `khash.h` and `kseq.h`
excluded) with the whole suite at chromap 949043c, g++ 13 `--coverage` at
`-O0` for the binary and the drivers, a line counted once over template
instances (gcovr `--merge-lines`): 87.9% (5788 of 6585 lines; 93.0% without
`ksw.cc`, mostly klib functions chromap never calls); functions 88.0%,
branches 59.2%. The command-line tests alone cover 74.8% of those lines,
the C++ drivers alone 53.9%. Public names no test reaches:
`GetLongestMatchLength` (split alignment with more than 200 candidates),
`AdjustGapBeginning` (an alignment starting with a gap), and code chromap
never calls: `BandedTracebackToEnd`, `mm_cache::DirectUpdateWeight`,
`mm_cache::GetMemoryBytes`, `GetByteSize` of `SAMMapping`, `PAFMapping` and
`PairedPAFMapping`, `Chromap::PairedEndReadWithBarcodeIsDuplicate`,
`Chromap::UpdateBarcodeAbundance` and `AddPeakOptions`. To measure it, build
chromap with `make CXXFLAGS="-std=c++11 -O0 -fopenmp -msse4.1 --coverage"
LDFLAGS="-lm -lz --coverage"`, run the suite with `CHROMAP_BIN` at that
binary and the same values in `CHROMAP_TEST_CXXFLAGS` and
`CHROMAP_TEST_LDFLAGS`, and set `GCOV_PREFIX` so the drivers' data, compiled
in temporary directories, can be paired with their `.gcno` files afterwards.

## Known bugs are kept out of the suite

Every test asserts behaviour that follows from the simulated truth, from
chromap's documentation (README.md, chromap.1, `chromap -h`), from the
intent its code states (comments naming a rule) or from a value worked by
hand. Where chromap 949043c gets a case wrong, the test is not in the suite;
the case is left out, to be reported as a GitHub issue with its
reproduction. The neighbouring correct cases stay in the suite, and
TARGETS.md names each kept-out case in its row by the short names below.
At this writing:

| Area | Notes |
|---|---|
| A | a-sam-fasta-reads, a-truncated-gzip-input, a-lowercase-read-sam-nm-md, a-mapping-statistics-clang-threads |
| B | b-contig-end-reads, b-error-threshold-16, b-error-threshold-shortcut, r1-single-end-drop-repetitive-reads |
| C | paired-paf-coordinates, low-mem-dedup-representative, c-bulk-level-dedup-in-memory, c-low-mem-bulk-dedup-summary, c-low-mem-bulk-dedup-without-whitelist, c-allocate-without-multimappers, c-single-end-tn5-dedup-order, c-summary-duplicate-cap, c-sam-single-end-bulk-dedup, c-paf-residue-matches, r2-bulk-dedup-barcode-count |
| D | single-end-chr-order, d-read-format-strand-on-reads, d-bc-threshold2-n, d-option-file-handling, d-mapq-threshold-300, d-version-stdout, d-feature-matrix |

A test that pins output nothing documents says "pins the current output" in
its docstring.

## Files

| File | Covers |
|---|---|
| `chromap_helpers.py` | the harness: synthetic genome, read simulator, `run_chromap`/`build_index`/`run_cli`, output parsers, Python reference functions |
| `chromap_build.py` | toolchain discovery, cached chromap build, C++ unit objects and drivers |
| `conftest.py` | session fixtures (binary, genome, index, simulator, whitelist, C++ drivers) |
| `pytest.ini` | markers `lowmem` and `cpp`; disables the cache provider so runs leave nothing in the tree |
| `cpp/chromap_test.h` | the C++ test framework: `CT_TEST`, `CT_PROBE`, `CT_CHECK`, `CT_CHECK_EQ`, `CT_MAIN` |
| `test_harness.py` | the harness itself: genome layout, simulator truth, reference functions, parsers (and against samtools/pysam), the runner |
| `test_pilot_cli.py` | pilot of the CLI layer: single-end truth, paired-end fragments, Tn5 shift, SAM fields, thread equality |
| `test_pilot_unit.py`, `cpp/pilot.cc` | pilot of the C++ layer: hand-worked CT_TESTs, the minimizer generator against its definition through a probe |
| `TARGETS.md` | every public name and option, its area and the tests covering it |
| `test_a_*.py`, `cpp/a_*.cc`, `helpers_a.py` | area A: index build and load, minimizers, seed encoding, `SequenceBatch`, read input files |
| `test_b_*.py`, `cpp/b_*.cc`, `helpers_b.py` | area B: candidates, draft mappings, alignment, MAPQ, pairing, split alignment, the candidate cache, presets, thread equivalence |
| `test_c_*.py`, `cpp/c_*.cc`, `helpers_c.py` | area C: output records and writers, sorting, deduplication, Tn5 shift, `--low-mem`, multi-mapping allocation, the summary |
| `test_d_*.py`, `cpp/d_*.cc`, `helpers_d.py` | area D: barcodes and whitelist correction, `--read-format`, `--barcode-translate`, pairs output, `--chr-order`, the feature matrix, the command-line front end |

## The synthetic genome

`make_synthetic_genome()` builds the same genome every time (seed
`GENOME_SEED`), 0-based half-open coordinates:

| Contig | Length | Planted |
|---|---|---|
| chr1 | 50,000 | dup2 copy 1 at 20000-21000 (+); rep8 at 30000, 32000, 34000 (-); dupmm copy 1 at 40000-40600 |
| chr2 | 30,000 | N at 5000-5100; dup2 copy 2 at 10000-11000 (reverse complement); rep8 at 20000, 22000 (-) |
| chr3 | 15,000 | (CA)100 at 1000-1200; rep8 at 3000, 6000; dupmm copy 2 at 10000-10600 (12 substitutions, offsets 25, 75, ..., 575) |
| chrX | 5,000 | rep8 at 1000; poly-A at 4000-4100 |

rep8 copies are 400 bp, dup2 1 kb, dupmm 600 bp. The FASTA headers carry a
comment (`>chr1 synthetic`) and 60-column lines. Everything else is uniform
random sequence, so a 17-mer occurs once.

## Harness API

```python
import collections
import chromap_helpers as ch

def test_example(index, sim):
	"""Single-end --Tn5-shift moves a + strand read right by 4 and pulls a
	- strand read's end (its 5' end) left by 5."""
	reads = sim.se_reads(50)
	run = ch.run_chromap(index, reads, "--Tn5-shift")
	got = collections.Counter((b.chrom, b.start, b.end) for b in run.bed())
	want = collections.Counter((r.contig, r.start + 4, r.end + 4) if r.strand == "+"
		else (r.contig, r.start, r.end - 5) for r in reads)
	assert got == want
```

### Fixtures (conftest.py)

| Fixture | Scope | Gives |
|---|---|---|
| `chromap_bin` | session | path of the binary (builds it on first use) |
| `chromap_threads` | session | the default `-t` of the session |
| `genome` | session | the default `SyntheticGenome`, FASTA written |
| `index` | session | `ChromapIndex` of `genome`, chromap's default k 17, w 7 |
| `index_factory` | session | `make(k=None, w=None, genome=None, extra_args=())` -> `ChromapIndex`, cached per argument set |
| `sim` | function | `ReadSimulator` seeded from the test's module and name |
| `whitelist` | session | `(barcodes, path)`: 32 barcodes of 16 bp, pairwise >= 3 mismatches apart |
| `cpp` | session | `cpp(name)` -> `CppDriver` for `tests/cpp/<name>.cc`: `tests()`, `probes()`, `check(test)`, `probe(name, *args, input=None)` (stdout; fails on a non-zero exit), `run(argv, input_text=None)` (the `CompletedProcess`) |

### Running chromap

- `run_chromap(index, reads=None, extra_args=(), threads=None, name=None, outdir=None, output="out.txt", gz=False, fasta=False, name_style="plain", check=True, timeout=None, shared=False)` -> `ChromapRun`
  - `reads`: a `ReadSet` (written to the run directory and passed as `-1`/`-2`/`-b`), a dict `{"r1": path, "r2": path, "bc": path}` (values may be lists, joined with commas), or `None` when `extra_args` carries `-1`.
  - `extra_args`: list or one shell-like string. chromap runs in the run directory, so `"--summary summary.csv"` lands there.
  - `threads`: `-t` and `OMP_NUM_THREADS`; `None` means `$CHROMAP_TEST_THREADS` (else 1). `-t` in `extra_args` is refused, so every run's thread count is pinned. Tests that compare thread counts pass `threads=1` and `threads=4`; all others leave it unset so the suite runs at 1 and at 4.
  - `check=True` raises `ChromapError` (an `AssertionError` with argv and stderr) on a non-zero exit; use `check=False` and `run.returncode` for error paths.
  - `shared=True` returns the earlier run of an identical request (same index, reads, arguments, threads, writing options) from this session; treat it as read-only.
- `build_index(genome, outdir=None, k=None, w=None, extra_args=(), fasta=None, name="index", check=True)` -> `ChromapIndex(path, genome, fasta, k, w, args, returncode, stderr)`; prefer `index_factory`.
- `run_cli(args, cwd=None, threads=None, timeout=None, no_core_dump=False)` -> `ChromapRun`: chromap with exactly `args`, never raises; for help, version and error paths. `no_core_dump=True` runs the child with `RLIMIT_CORE` 1 on Linux (0 elsewhere, a no-op without the `resource` module), so a command line that makes chromap abort (cxxopts throwing) does not wait about a second for the system's crash handler.
- `ChromapRun`: `outdir`, `argv`, `returncode`, `stdout`, `stderr`, `output` (Path of `-o`), `reads`, `truth` (name -> `SimRead`/`SimPair`), `path(name)`, `exists(name)`, `text()`, `lines()`, `bed()`, `fragments()`, `tagalign()`, `sam()`, `pairs()`, `paf()`, `summary(name="summary.csv")`, `stats()` (integer statistics from stderr).

### Simulating reads

`ReadSimulator(genome, seed=1, read_len=50, qual_char="I")`; use the `sim`
fixture in tests and `ReadSimulator(genome, seed=ch.stable_seed("label"))`
in module fixtures, never one simulator shared between tests. Mismatch
offsets count along the read as sequenced (5' to 3').

- `se_read(contig, start, length=None, strand="+", mismatches=(), name=None, qual=None, tail="", barcode=None, barcode_qual=None)` -> `SimRead`: the read covers `contig[start:start+length]`; `tail` (adapter) is appended at the 3' end and is not genomic.
- `pair(contig, start, insert, read_len=None, strand="+", mismatches1=(), mismatches2=(), name=None, barcode=None, barcode_qual=None, qual1=None, qual2=None)` -> `SimPair`: fragment `contig[start:start+insert]`; `strand` is read 1's; inserts shorter than the read run into Nextera adapters (`NEXTERA_R1`, `NEXTERA_R2`).
- `segment_read(segment, offset, length=None, strand="+", copy=0, mismatches=())` -> `SimRead` with `extra["loci"]` (every copy's interval); `random_read(length=None)` (should not map; `contig` is None).
- batches returning a `ReadSet`: `se_reads(n, length=None, strands="+-", margin=200, barcodes=None)`, `pairs(n, insert=(150, 400), read_len=None, strands="+-", margin=200, barcodes=None)` at unique loci; `shuffled(readset)`.
- `unique_locus(length, margin=200, contigs=None, pad=50)` -> `(contig, start)` away from contig ends and planted features.
- `SimRead`: `name seq qual contig start end strand mismatches tail barcode barcode_qual kind extra`, `interval()`, `five_prime`. `SimPair`: `name r1 r2 contig start end strand barcode barcode_qual insert`, `fragment()` -> `(contig, start, end, strand of read 1)`.
- `ReadSet(records)`: iterable, `by_name`, `paired`, `has_barcodes`, `mates(i)`, `+`, `digest()`, `write(directory, prefix="reads", gz=False, fasta=False, name_style="plain"|"comment"|"slash")` -> `{"r1", "r2", "bc"}` paths.
- module level: `make_whitelist(n, length=16, seed=7, min_dist=3)`, `write_lines(path, lines)`, `stable_seed(*parts)`.

### Parsers

`parse_bed(path, allow_six=False)` -> `[BedRecord(chrom start end name mapq strand dups)]` (7-column BED/TagAlign; `allow_six` for TagAlign lines without a count), `parse_fragments` -> `[FragmentRecord(chrom start end barcode dups)]`, `parse_sam` -> `(header, [SamRecord])` (`tags` with ints for `i`, `is_reverse`, `is_read1`, `is_read2`, `is_secondary`), `parse_pairs` -> `(header, [PairsRecord(read_id chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2)])`, `parse_paf` -> `[PafRecord]`, `parse_summary` -> `Summary` (`header`, `rows` barcode -> dict; the bulk row's barcode is `""`), `parse_stats(text)`, `first_line_diff(a, b)`.

### Reference functions

`revcomp`, `hamming`, `substitute`, `encode_seed(seq, start=0, length=None)` (chromap's 2-bit seeds: A0 C1 G2 T3, N and padding as A), `decode_seed`, `hash64(key, mask)`, `kmer_minimizer_key(kmer)` -> `(hash, strand)`, `reference_minimizers(seq, k, w, seq_index=0)` (minimizers by definition, N-free sequence and odd k), constants `DEFAULT_KMER`, `DEFAULT_WINDOW`, `DEFAULT_MAPQ_THRESHOLD`, `MAX_THREADS`.

## The C++ layer

A driver `tests/cpp/<area>_<topic>.cc` includes `chromap_test.h` first (some
`src/` headers rely on standard headers included before them), then the
`src/` headers it needs, and ends with `CT_MAIN`:

```cpp
#include "chromap_test.h"
#include "bed_mapping.h"
using namespace chromap;

// + strand: start moves right by 4 (100 -> 104).
CT_TEST(tn5_shift_plus) {
  MappingWithoutBarcode m(0, 100, 50, 60, 1, 1, 1);
  m.Tn5Shift();
  CT_CHECK_EQ(m.GetStartPosition(), 104u);
}

CT_PROBE(encode) {  // driver --probe encode SEQ -> the key on stdout
  std::cout << GenerateSeedFromSequence(args[0].c_str(), args[0].size(), 0, args[0].size()) << "\n";
  return 0;
}

CT_MAIN
```

The matching Python file lists every `CT_TEST` and checks the list against
`driver.tests()`, so none is skipped silently, then runs each as its own
pytest item (`cpp("name").check(test)`); probes print values that the Python
test compares with expectations computed in Python (see `test_pilot_unit.py`).
Mark the module `pytestmark = pytest.mark.cpp`.

- Use only public members. Protected members (`MappingMetadata`'s vectors)
  are reached through a test subclass with `using` declarations; no
  `#define private public`.
- Every signature a driver calls becomes part of what a change to `src/`
  must keep compiling: such a change may edit `tests/cpp/*.cc` only to follow
  a signature it changed, never an input, an expected value or an assertion.
  Call the stable public entry points listed in TARGETS.md, avoid incidental
  helpers, and check what the caller sees rather than how the function gets
  there (no order a later sort erases, no wording of a diagnostic line).
- `ExitWithMessage` ends the process with status 255; test such paths
  through a probe and its exit code: `cpp(name).run(["--probe", "x", ...])`
  returns the `CompletedProcess` (`returncode`, `stdout`, `stderr` as bytes)
  without failing.
- Only the C++11 standard library; no downloads.

## Writing tests

- Expected values come from the simulation truth, a reference implementation
  in the test or `helpers_<area>.py`, or arithmetic in the docstring. A test
  that only records chromap's current output says "pins the current output";
  keep those rare. Checking only that chromap exits 0 is not coverage.
- When chromap's behaviour is wrong (it contradicts its documentation, the
  truth, or itself across modes), leave the test out, report the case as a
  GitHub issue with its reproduction, and keep the neighbouring correct cases
  in the suite.
- Use `sim` for reads, the session `index`, and `index_factory` for other k,
  w or genomes. Keep each test to a few chromap runs; share a module-scoped
  run between tests that only read it; respect the per-area run budget in
  TARGETS.md.
- Choose reads whose truth is unambiguous: `unique_locus`/`se_reads`/`pairs`
  for unique reads, planted segments for multimappers, substitutions placed
  so enough error-free 17-mers remain to seed the read.
- stderr: assert the statistics chromap reports ("Number of reads",
  "Number of mapped reads", "Number of output mappings (passed filters)",
  "Number of barcodes in whitelist", ...), error and warning messages, and
  the parameter report the driver prints before mapping; never timing
  lines, progress lines, or the counts chromap prints between pipeline
  steps ("After removing PCR duplications, ...", "Allocated N
  multi-mappings in ...", "# shifted mappings"), which a faster pipeline may
  print differently or not at all. Per-thread statistics only at
  `threads=1` (see the clang note above).
- Temp files: check that the run directory holds only what the test expects
  afterwards, never a temp file's name or how many there were.
- Python 3.8 syntax, standard library only; indent with tabs.
- Tick the rows you cover in `TARGETS.md`.

## chromap behaviour to know when writing tests

Observed at 949043c; the line numbers move, the behaviour should not.

- **Coordinates.** BED, TagAlign and PAF are 0-based half-open; SAM POS and
  pairs positions are 1-based. Pairs positions are each read's 5' end (the
  last base of a reverse-strand read). Single-end BED intervals are the read's
  aligned interval; paired-end BED intervals are the fragment, with read 1's
  strand.
- **Order.** Output is sorted per reference (FASTA order, or the
  `--chr-order` order) by each record class's `operator<`, which ends in the
  read id, so the order is total and outputs are byte-identical across
  thread counts (checked with 12,000 pairs and multimappers at `-q 0`). The
  summary CSV's rows come in hash-table order: compare them as a set.
- **Threads.** Single-end mapping is split into `t * t` OpenMP tasks, so any
  input exercises several threads. Paired-end mapping makes one task per
  5,000 to 9,999 pairs (`taskloop grainsize(5000)`), so fewer than 10,000
  pairs are mapped by one thread whatever `-t` is; thread tests need more.
  Multimappers' locations come from a generator seeded with 11 per read
  (single-end) or per task copy (paired-end), so they do not depend on `-t`.
- **Batches and the cache.** Reads are mapped in batches of 500,000; the
  candidate cache (`mm_cache`) is updated only after a batch, so a run of
  fewer reads never hits it and its summary's `cachehit`, `fric`, `estfrip`
  and `numcacheslots` are 0. Only `test_b_cache.py` maps more than one batch
  (510,000 reads or pairs), to check that reads answered from the cache map
  as they do without it.
- **Seeds.** A read maps only with at least `-s` (2) matching minimizers
  (k 17, w 7: about 2 per 8 bases). A 50 bp read has 34 17-mers; a
  substitution destroys the 17 that cover it. Place substitutions near the
  ends or keep one in the middle.
- **MAPQ.** Unique exact reads get 60; reads at several equally good loci
  get 0 and are dropped at the default `-q 30`.
- **Deduplication.** Single-end duplicates share a start position whatever
  their strand or length; paired-end duplicates share start and length;
  barcoded records also need the same barcode (cell level), or are collapsed
  across barcodes (bulk level, the default for single-cell data). The count
  column is the group size, capped at 255. MAPQ filtering applies to the
  record kept for each group, after deduplication. `--low-mem` and the
  default mode form the same groups but may keep a different record (known
  bug note); assert only on identical duplicates or on counts. Bulk-level
  removal (single-cell data with `--low-mem`) writes each group with the
  barcode that has the most pairs in it, ties broken by the barcode's
  abundance (exact whitelist matches in the input); chromap counts at most
  two pairs per barcode there (known bug r2-bulk-dedup-barcode-count).
- **Presets.** `atac` and `chip` turn on `--low-mem` (12 GB per run, see
  Runtime); `hic` writes pairs and turns on split alignment.
- **Names.** Read names are cut at the first whitespace; a `/1` or `/2`
  suffix is kept (SAM QNAMEs then differ between mates). Pairs use read 1's
  name; SAM and paired-end PAF write each mate's own.
- **Messages.** Help, version, progress and statistics go to stderr.
  `ExitWithMessage` exits with status 255. An option chromap does not
  register (the peak options, `-p`) or a malformed value makes cxxopts throw
  and the process abort: assert a non-zero exit and the message only.
- **Single-end only / paired-end only.** `--pairs` with single-end reads
  exits ("No support for single-end HiC yet!"); `--frip-est-params` is parsed
  only for paired-end runs.
