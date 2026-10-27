"""Area A: building, saving and loading the index.

The index file is parsed in Python (helpers_a.parse_index) and compared with
the index computed independently from the reference: the minimizers of
helpers_a.sketch, grouped by hash into singletons and occurrence lists
(helpers_a.expected_index). Index::Load, Save, CheckIndex, the index_utils.h
layouts and the heap merge are reached through tests/cpp/a_index.cc."""

import hashlib
import random
import re
import shutil

import pytest

import chromap_helpers as ch
import helpers_a as ha

A_INDEX_TESTS = [
	"lookup_table_layout_hand_worked",
	"candidate_position_layout_hand_worked",
	"khash_functions_ignore_singleton_bit",
	"heap_entry_order_and_seed_stats",
	"index_parameters_defaults",
]

A_INDEX_PROBES = ["heap_merge", "load_save", "check"]


# ------------------------------------------------------------ fixtures

def _rand(rng, n, alphabet="ACGT"):
	return "".join(rng.choice(alphabet) for _ in range(n))


@pytest.fixture(scope="module")
def edge_reference(tmp_path_factory):
	"""A reference FASTA with what genomes carry: a contig entirely in lower
	case; one with a mixed-case stretch, IUPAC codes, N runs of 1, 5 and 100
	bases, a leading run of 20 N and a trailing run of 30 N; a contig of 10
	bases (shorter than k, no minimizers but its own index); an empty record
	(skipped, no index); a contig carrying 300 bases of the first contig and
	their reverse complement (minimizers seen several times on both strands);
	and 250 contigs of 40 to 79 bases (more than the 200 LoadAllSequences
	reserves). Headers carry space and tab comments; lines are 61 bases.
	Returns (path, upper-case copy's path, [sequences of the kept records])."""
	rng = random.Random(ch.stable_seed("a_index", "edge"))
	lower = _rand(rng, 3000).lower()
	mixed = ("N" * 20 + _rand(rng, 300) + _rand(rng, 200).lower() + "N" + _rand(rng, 250) + "RYKM" + _rand(rng, 200)
		+ "N" * 5 + _rand(rng, 400) + "N" * 100 + _rand(rng, 600) + "N" * 30)
	copy = lower[1000:1300].upper()
	rep = _rand(rng, 500) + copy + _rand(rng, 200) + ch.revcomp(copy) + _rand(rng, 300) + copy + _rand(rng, 100)
	records = [("edge_lower lower case", lower), ("mixed\tN runs", mixed), ("short", _rand(rng, 10)), ("empty", ""),
		("rep copies", rep)]
	records += [("tiny%03d" % i, _rand(rng, 40 + i % 40)) for i in range(250)]
	d = tmp_path_factory.mktemp("edge_reference")
	path = ha.write_fasta(d / "edge.fa", records, line_width=61)
	upper = ha.write_fasta(d / "edge_upper.fa", [(h, s.upper()) for h, s in records], line_width=61)
	return path, upper, [s for _, s in records if s]


@pytest.fixture(scope="module")
def even_reference(tmp_path_factory):
	"""An N-free reference for even k: random contigs carrying palindromes
	(runs of AT, of ACGT and GAATTC sites), so symmetric k-mers occur."""
	rng = random.Random(ch.stable_seed("a_index", "even"))
	seqs = []
	for i in range(6):
		seqs.append(_rand(rng, 400) + "AT" * 20 + _rand(rng, 300) + "ACGT" * 10 + _rand(rng, 200) + "GAATTC"
			+ _rand(rng, 500))
	d = tmp_path_factory.mktemp("even_reference")
	path = ha.write_fasta(d / "even.fa", [("e%d" % i, s) for i, s in enumerate(seqs)])
	return path, seqs


def _build(fasta, outdir, k=None, w=None, extra_args=(), name="index"):
	return ch.build_index(None, outdir, k=k, w=w, extra_args=extra_args, fasta=fasta, name=name)


def _content(parsed):
	return {key: v for key, (_, v) in parsed.entries.items()}, list(parsed.occurrences)


@pytest.fixture(scope="module")
def default_expected(genome):
	"""The default genome's index as helpers_a.expected_index computes it."""
	return ha.expected_index(list(genome.contigs.values()), ch.DEFAULT_KMER, ch.DEFAULT_WINDOW)


# ------------------------------------------------------------ C++ unit tests

@pytest.mark.cpp
def test_a_index_driver_lists_every_test(cpp):
	"""Every CT_TEST in cpp/a_index.cc is run by test_a_index_cpp, and every
	probe used below exists."""
	driver = cpp("a_index")
	assert driver.tests() == A_INDEX_TESTS
	assert driver.probes() == A_INDEX_PROBES


@pytest.mark.cpp
@pytest.mark.parametrize("name", A_INDEX_TESTS)
def test_a_index_cpp(cpp, name):
	cpp("a_index").check(name)


@pytest.mark.cpp
def test_heap_merge_equals_sorted_concatenation(cpp):
	"""HeapMergeCandidatePositionLists merges sorted lists into their sorted
	concatenation, appended to what the output vector already holds: random
	lists with duplicates across and within lists, empty lists, one list, no
	list, and 64-bit positions."""
	rng = random.Random(ch.stable_seed("a_index", "heap"))
	driver = cpp("a_index")
	cases = [[], [[]], [[5]], [[1, 1, 2], [], [1, 3]], [[2 ** 63 + 5, 2 ** 64 - 1], [0, 2 ** 40]]]
	for _ in range(40):
		lists = []
		for _ in range(rng.randint(1, 8)):
			n = rng.choice([0, 1, 3, 10, 40])
			lists.append(sorted(rng.randint(0, rng.choice([5, 1000, 2 ** 40])) for _ in range(n)))
		cases.append(lists)
	for lists in cases:
		text = "".join(" ".join(map(str, lst)) + "\n" for lst in lists)
		got = [int(x) for x in driver.probe("heap_merge", input=text).split()]
		assert got == sorted(x for lst in lists for x in lst), lists
		got = [int(x) for x in driver.probe("heap_merge", "append", input=text).split()]
		assert got == [7] + sorted(x for lst in lists for x in lst), lists


# ------------------------------------------------------------ index content

def test_default_index_matches_reference(index, genome, default_expected):
	"""chromap -i on the default genome (k 17, w 7) writes k and w, then a
	lookup table holding every minimizer hash of the reference once: key
	hash << 1 with bit 0 set and the hit as value for hashes seen once, else
	offset << 32 | count, with the hits sorted ascending at that offset of
	the occurrence table (offsets in ascending hash order). The planted
	repeats make non-singletons on both strands and the N run on chr2 drops
	the minimizer pending before it (helpers_a.sketch)."""
	parsed = ha.parse_index(index.path)
	entries, occurrences = default_expected
	assert (parsed.k, parsed.w) == (17, 7)
	assert parsed.table_size == len(entries) == parsed.size == parsed.n_occupied
	got_entries, got_occurrences = _content(parsed)
	assert got_entries == entries
	assert got_occurrences == occurrences
	# The data exercise both kinds of key, hits on both strands of one hash,
	# and the minimizer dropped before the N run.
	assert sum(1 for key in entries if key & 1) < len(entries)
	assert any(len({hit & 1 for hit in occurrences[(v >> 32):(v >> 32) + (v & 0xFFFFFFFF)]}) == 2
		for key, v in entries.items() if not key & 1)
	chr2 = genome.contigs["chr2"]
	plain = ha.plain_segment_minimizers(chr2, 17, 7, genome.rid("chr2"))
	kept = ha.sketch(chr2, 17, 7, genome.rid("chr2"))
	assert len(plain) == len(kept) + 1 and set(kept) < set(plain)


def test_default_index_bytes(index):
	"""Pins the current output: the md5 of the default genome's index file
	(k 17, w 7) as chromap 949043c writes it. The other tests check what the
	file holds; this one also fixes the khash bucket layout and every other
	byte, so that indexes written before and after a change are the same
	file. (The genome is the same under every supported Python version.)"""
	assert hashlib.md5(index.path.read_bytes()).hexdigest() == "f15ce52dc146b9597ed2414041ae4af4"


def test_default_index_khash_layout(index, default_expected):
	"""The khash dump is a valid table: a power-of-two bucket count,
	size == n_occupied (nothing deleted) == the header's table size,
	upper_bound = n_buckets * 0.77 + 0.5 truncated, no deleted flags, and
	each key found where khash's probe (hash = key >> 1, quadratic steps)
	looks for it; nothing follows the occurrence table."""
	parsed = ha.parse_index(index.path)
	n = parsed.n_buckets
	assert n >= 4 and n & (n - 1) == 0
	assert parsed.size == parsed.n_occupied == parsed.table_size == len(parsed.entries)
	assert parsed.upper_bound == int(n * 0.77 + 0.5)
	assert parsed.size <= parsed.upper_bound
	deleted = [i for i in range(n) if (parsed.flags[i >> 4] >> ((i & 15) << 1)) & 3 == 1]
	assert deleted == []
	for key, (bucket, _) in parsed.entries.items():
		assert ha.khash_find(parsed, key) == bucket
	entries, _ = default_expected
	missing = [h for h in (5, 77, 1 << 20) if (h << 1) not in entries and ((h << 1) | 1) not in entries]
	assert all(ha.khash_find(parsed, h << 1) is None for h in missing)


def test_index_statistics(index, genome, default_expected):
	"""Index::Statistics, whose only output is its report on stderr when the
	index is built, gives k, w and the number of reference sequences, the
	number of distinct minimizers, the percentage that are singletons,
	the average occurrences per distinct minimizer and the average spacing
	(reference bases per minimizer occurrence), all computed here from the
	expected index. Other log lines of index building are not asserted."""
	entries, occurrences = default_expected
	err = index.stderr
	n = len(entries)
	n1 = sum(1 for key in entries if key & 1)
	total = n1 + len(occurrences)
	assert re.search(r"kmer size: 17; skip: 7; #seq: 4\b", err)
	stats = re.search(r"distinct minimizers: (\d+) \(([\d.]+)% are singletons\); average occurrences: ([\d.]+); "
		r"average spacing: ([\d.]+)", err)
	assert stats is not None
	assert int(stats.group(1)) == n
	assert stats.group(2) == "%.2f" % (100.0 * n1 / n)
	assert stats.group(3) == "%.3f" % (float(total) / n)
	assert stats.group(4) == "%.3f" % (float(genome.length) / total)


@pytest.mark.parametrize("k,w", [(17, 7), (21, 11), (11, 3), (27, 2), (15, 1)])
def test_index_of_edge_reference(edge_reference, tmp_path, k, w):
	"""chromap -i -k K -w W on a reference with lower case, IUPAC codes, N
	runs (leading, inner, trailing), a contig shorter than k, an empty record
	and 250 small contigs stores exactly helpers_a.expected_index: no k-mer
	spanning an ambiguous base, hits numbered by the non-empty records in
	FASTA order, lower case read as upper case. The index of the same FASTA
	in upper case is identical in content."""
	path, upper, seqs = edge_reference
	expected = ha.expected_index(seqs, k, w)
	parsed = ha.parse_index(_build(path, tmp_path / "lower", k, w).path)
	assert (parsed.k, parsed.w, parsed.table_size) == (k, w, len(expected[0]))
	assert _content(parsed) == expected
	parsed_upper = ha.parse_index(_build(upper, tmp_path / "upper", k, w).path)
	assert _content(parsed_upper) == expected
	hits = [v for key, v in expected[0].items() if key & 1] + expected[1]
	indices = {hit >> 33 for hit in hits}
	assert 2 not in indices and len(seqs) - 1 in indices


@pytest.mark.parametrize("k,w", [(16, 5), (18, 7), (4, 3)])
def test_index_with_even_k(even_reference, tmp_path, k, w):
	"""With even k the index skips k-mers equal to their reverse complement
	(helpers_a.expected_index), and k and w are written to the header."""
	path, seqs = even_reference
	parsed = ha.parse_index(_build(path, tmp_path, k, w).path)
	assert (parsed.k, parsed.w) == (k, w)
	assert _content(parsed) == ha.expected_index(seqs, k, w)


@pytest.mark.parametrize("args,k,w", [
	([], 17, 7),
	(["--min-frag-length", "30"], 17, 7),
	(["--min-frag-length", "60"], 17, 7),
	(["--min-frag-length", "61"], 19, 10),
	(["--min-frag-length", "80"], 19, 10),
	(["--min-frag-length", "81"], 23, 11),
	(["--min-frag-length", "150"], 23, 11),
	(["--min-frag-length", "81", "-k", "15"], 15, 11),
	(["--min-frag-length", "81", "-w", "5"], 23, 5),
	(["-k", "21", "-w", "9", "--min-frag-length", "61"], 21, 9),
	(["--kmer", "13", "--window", "4"], 13, 4),
])
def test_min_frag_length_k_and_w(even_reference, tmp_path, args, k, w):
	"""--min-frag-length picks k and w (<= 60: 17 and 7; <= 80: 19 and 10;
	above: 23 and 11), and -k/--kmer and -w/--window given as well win
	wherever they stand on the command line. The chosen values are written to
	the index header, and the index holds the minimizers of that k and w."""
	path, seqs = even_reference
	built = _build(path, tmp_path, extra_args=args)
	parsed = ha.parse_index(built.path)
	assert (parsed.k, parsed.w) == (k, w)
	assert _content(parsed) == ha.expected_index(seqs, k, w)


def test_build_index_needs_ref_and_output(genome, tmp_path):
	"""-i without -r exits 255 with "No reference specified!", without -o
	with "No output file specified!", and writes no index."""
	run = ch.run_cli(["-i", "-o", str(tmp_path / "idx")], cwd=tmp_path)
	assert run.returncode == 255
	assert "No reference specified!" in run.stderr
	run = ch.run_cli(["-i", "-r", str(genome.fasta_path)], cwd=tmp_path)
	assert run.returncode == 255
	assert "No output file specified!" in run.stderr
	assert not (tmp_path / "idx").exists()


# ------------------------------------------------------------ load and save

@pytest.mark.cpp
def test_load_then_save_is_byte_identical(cpp, index, edge_reference, tmp_path):
	"""Index::Load reads k, w and the lookup table (GetKmerSize,
	GetWindowSize, GetLookupTableSize = number of distinct minimizer
	hashes), and Index::Save writes the same bytes back, for the default
	index and a k 21, w 11 index; Destroy may be called repeatedly."""
	parsed = ha.parse_index(index.path)
	edge = _build(edge_reference[0], tmp_path / "edge", 21, 11)
	edge_parsed = ha.parse_index(edge.path)
	for src, p in ((index.path, parsed), (edge.path, edge_parsed)):
		copy = tmp_path / ("copy_" + src.name + str(p.k))
		shutil.copyfile(str(src), str(copy))
		out = cpp("a_index").probe("load_save", copy)
		assert out.split() == [str(p.k), str(p.w), str(len(p.entries))]
		assert copy.read_bytes() == src.read_bytes()


@pytest.mark.cpp
def test_check_index_accepts_built_index(cpp, index, genome):
	"""Index::CheckIndex regenerates the reference's minimizers and asserts
	(assert, so only in builds without NDEBUG, like the Makefile's) that each
	is stored with its hit at its place in the occurrence table; it returns
	for the index chromap built from the same reference."""
	assert cpp("a_index").probe("check", index.path, genome.fasta_path).strip() == "checked"
