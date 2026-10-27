"""Area D, C++ layer: pairs records, the barcode-correction candidate order,
--read-format ranges and the feature/barcode matrix classes, through
tests/cpp/d_units.cc.

The feature/barcode matrix (FeatureBarcodeMatrix, FeatureBarcodeMatrixWriter,
Peak) is not reachable from the command line at 949043c (the peak options
are never registered and -p is commented out), so it is tested here only.
Expected values come from hand-worked cases in the driver's comments and from
the Python references below."""

import collections
import gzip
import random

import pytest

import chromap_helpers as ch
import helpers_d as hd

pytestmark = pytest.mark.cpp

D_UNITS_TESTS = [
	"pairs_mapping_positions_and_strands",
	"pairs_mapping_barcode_and_tn5_shift",
	"pairs_mapping_order",
	"pairs_mapping_duplicate_keys",
	"pairs_mapping_write_load_round_trip",
	"pairs_mapping_byte_size_counts_the_name",
	"barcode_with_qual_order",
]


def test_d_units_driver_lists_every_test(cpp):
	"""Every CT_TEST in cpp/d_units.cc is run below, and the probes the
	tests use exist."""
	driver = cpp("d_units")
	assert driver.tests() == D_UNITS_TESTS
	assert driver.probes() == ["effective_range", "feature_matrix", "barcode_output"]


@pytest.mark.parametrize("name", D_UNITS_TESTS)
def test_d_units_cpp(cpp, name):
	"""The hand-worked C++ cases (see the comments in cpp/d_units.cc;
	pairs_mapping_barcode_and_tn5_shift pins the current output)."""
	cpp("d_units").check(name)


# ------------------------------------------------------------ read format

def _write_fastq(path, records):
	with open(str(path), "w") as f:
		for i, (s, q) in enumerate(records):
			f.write("@q%d\n%s\n+\n%s\n" % (i, s, q))
	return path


def _effective_range(cpp, tmp_path, records, fields):
	fq = _write_fastq(tmp_path / ("reads_%d.fq" % len(list(tmp_path.iterdir()))), records)
	out = cpp("d_units").probe("effective_range", fq, *fields)
	if out.strip() == "rejected":
		return "rejected"
	return [tuple(line.split("\t")) for line in out.splitlines()]


def test_read_format_ranges_hand_worked(cpp, tmp_path):
	"""Inclusive ranges, -1 as the read's last base, several ranges
	concatenated in order, and "-" reverse-complementing the sequence (N stays
	N) while reversing the quality string:
	ACGTNACGTT with bc:0:2 + bc:7:-1 is ACG + GTT = ACGGTT; with bc:2:5:-
	the range GTNA becomes TNAC and its qualities 2345 become 5432; with
	bc:0:-1:- the whole read becomes AACGTNACGT, qualities 9876543210."""
	recs = [("ACGTNACGTT", "0123456789")]
	assert _effective_range(cpp, tmp_path, recs, ["bc:0:2", "bc:7:-1"]) == [("ACGGTT", "012789")]
	assert _effective_range(cpp, tmp_path, recs, ["bc:2:5:-"]) == [("TNAC", "5432")]
	assert _effective_range(cpp, tmp_path, recs, ["bc:0:-1:-"]) == [("AACGTNACGT", "9876543210")]
	assert _effective_range(cpp, tmp_path, recs, ["r1:0:-1:+"]) == [recs[0]]
	assert _effective_range(cpp, tmp_path, recs, []) == [recs[0]]


def test_read_format_ranges_match_reference(cpp, tmp_path):
	"""For random reads (A, C, G, T and N) and random range lists (one to
	three ordered, non-overlapping ranges, ends given or -1, strand omitted,
	"+" or "-" on any range), the sequence and quality chromap extracts equal
	helpers_d.apply_read_format, the README's grammar applied by hand."""
	rng = random.Random(ch.stable_seed("test_d_unit", "read_format"))
	for case in range(40):
		n = rng.randint(1, 3)
		length = rng.randint(40, 90)
		cuts = sorted(rng.sample(range(length), 2 * n))
		fields = []
		for k in range(n):
			start, end = cuts[2 * k], cuts[2 * k + 1]
			end_txt = str(end)
			if k == n - 1 and rng.random() < 0.5:
				end_txt = "-1"
			f = "bc:%d:%s" % (start, end_txt)
			r = rng.random()
			if r < 0.25:
				f += ":-"
			elif r < 0.4:
				f += ":+"
			fields.append(f)
		reads = []
		for _ in range(3):
			s = "".join(rng.choice("ACGTACGTACGTN") for _ in range(length))
			q = "".join(chr(33 + rng.randint(0, 41)) for _ in range(length))
			reads.append((s, q))
		got = _effective_range(cpp, tmp_path, reads, fields)
		want = [hd.apply_read_format(s, q, fields) for s, q in reads]
		assert got == want, (fields, reads)


@pytest.mark.parametrize("fields", [
	["r1:0"],
	["r1:0:5:+:x"],
	["r1:3:9", "r1:12"],
])
def test_read_format_malformed_range_rejected(cpp, tmp_path, fields):
	"""A range without its end, or with a fifth field, is refused by the
	parser (ParseFormatStringAndAppendEffectiveRange returns false)."""
	assert _effective_range(cpp, tmp_path, [("ACGTACGTACGTACGT", "I" * 16)], fields) == "rejected"


# ------------------------------------------------------------ feature matrix

MATRIX_CONTIGS = (("cA", 2000), ("cB", 1500), ("cC", 700))


@pytest.fixture(scope="module")
def matrix_fasta(tmp_path_factory):
	rng = random.Random(ch.stable_seed("test_d_unit", "matrix_fasta"))
	path = tmp_path_factory.mktemp("matrix") / "ref.fa"
	lines = []
	for name, n in MATRIX_CONTIGS:
		lines += [">%s" % name, "".join(rng.choice("ACGT") for _ in range(n))]
	return ch.write_lines(path, lines)


def reference_feature_matrix(fragments, cell_by_bin, bin_size, overlap, depth):
	"""(peak lines, barcode lines, matrix lines) of the feature/barcode
	matrix, computed from the fragments (rid, start, length, barcode).

	Features are bins of bin_size per contig (1-based start, start + bin_size)
	when cell_by_bin, else peaks: maximal runs of positions covered by at
	least `depth` fragments, numbered in reference order. A fragment counts
	for a feature it overlaps once both are extended by `overlap` bases
	(peaks) or for every bin it covers. Barcodes are numbered in order of
	first appearance (reference by reference, input order within one) and
	written as SEQ-1; matrix lines are "feature barcode count", by barcode
	then feature, under a "features barcodes lines" header."""
	features = []
	if cell_by_bin:
		offsets, total = [], 0
		for name, n in MATRIX_CONTIGS:
			offsets.append(total)
			for pos in range(0, n, bin_size):
				features.append("%s\t%d\t%d" % (name, pos + 1, pos + bin_size))
			total += (n + bin_size - 1) // bin_size
		num_features = total

		def overlapped(rid, start, length):
			return [offsets[rid] + b for b in range(start // bin_size, (start + length - 1) // bin_size + 1)]
	else:
		peaks = []
		for rid, (name, n) in enumerate(MATRIX_CONTIGS):
			cov = [0] * n
			for r, s, l, _ in fragments:
				if r == rid:
					for p in range(s, s + l):
						cov[p] += 1
			p = 0
			while p < n:
				if cov[p] >= depth:
					q = p
					while q < n and cov[q] >= depth:
						q += 1
					peaks.append((rid, p, q, len(peaks)))
					features.append("%s\t%d\t%d" % (name, p + 1, q))
					p = q
				else:
					p += 1
		num_features = len(peaks)

		def overlapped(rid, start, length):
			lo, hi = max(0, start - overlap), start + length + overlap
			return [i for r, ps, pe, i in peaks if r == rid and ps < hi and lo < pe]
	order = []
	for rid in range(len(MATRIX_CONTIGS)):
		for r, s, l, bc in fragments:
			if r == rid and bc not in order:
				order.append(bc)
	counts = collections.Counter()
	for rid in range(len(MATRIX_CONTIGS)):
		for r, s, l, bc in fragments:
			if r == rid:
				for f in overlapped(r, s, l):
					counts[(order.index(bc), f)] += 1
	matrix = ["%d\t%d\t%d" % (num_features, len(order), len(counts))]
	matrix += ["%d\t%d\t%d" % (f, b, c) for (b, f), c in sorted(counts.items())]
	return features, ["%s-1" % bc for bc in order], matrix


def _feature_matrix(cpp, fasta, tmp_path, fragments, cell_by_bin, bin_size, overlap, depth):
	prefix = tmp_path / ("m_%d_%d_%d_%d" % (cell_by_bin, bin_size, overlap, depth))
	stdin = "".join("%d %d %d %s\n" % f for f in fragments)
	cpp("d_units").probe("feature_matrix", fasta, prefix, int(cell_by_bin), bin_size, overlap, depth, input=stdin)
	return tuple(open(str(prefix) + suffix).read().splitlines() for suffix in ("_peaks.bed", "_barcode.tsv",
		"_matrix.mtx"))


BC = ["AAAACCCCGGGGTTTT", "ACGTACGTACGTACGT", "TTTTGGGGCCCCAAAA", "CATGCATGCATGCATG"]

# Hand-made fragments (rid, start, length, barcode), none touching a contig
# end. cA: [100, 300), [150, 320) and [200, 400) pile up to depth 3 on
# [200, 300); [330, 410) starts 30 bp after that peak (depth 2 with
# [200, 400)); two fragments of depth 2 at 1000. cB: depth 4 on [500, 600)
# and two depth-3 runs, [800, 900) and [950, 1000).
HAND_FRAGMENTS = [
	(0, 100, 200, BC[0]), (0, 150, 170, BC[1]), (0, 200, 200, BC[0]),
	(0, 330, 80, BC[2]),
	(0, 1000, 100, BC[1]), (0, 1000, 100, BC[3]),
	(1, 500, 100, BC[1]), (1, 500, 100, BC[1]), (1, 500, 100, BC[2]), (1, 500, 100, BC[2]),
	(1, 800, 100, BC[3]), (1, 800, 100, BC[3]), (1, 800, 200, BC[0]),
	(1, 950, 50, BC[0]), (1, 950, 50, BC[2]),
]


def test_feature_matrix_peaks_hand_worked(cpp, matrix_fasta, tmp_path):
	"""Peak calling at depth 3 on the hand-made fragments: peaks cA [200,
	300), cB [500, 600), [800, 900) and [950, 1000), written 1-based (201-300,
	...); the depth-2 pile at cA 1000 is no peak, and [330, 410) misses the
	cA peak. Barcode indices follow first appearance (BC[0] 0, BC[1] 1, BC[2]
	2, BC[3] 3); matrix lines are "peak barcode count"."""
	peaks, barcodes, matrix = _feature_matrix(cpp, matrix_fasta, tmp_path, HAND_FRAGMENTS, False, 5000, 0, 3)
	assert peaks == ["cA\t201\t300", "cB\t501\t600", "cB\t801\t900", "cB\t951\t1000"]
	assert barcodes == ["%s-1" % b for b in BC]
	# Peak 0 (cA): BC[0] x2, BC[1] x1. Peak 1: BC[1] x2, BC[2] x2. Peak 2:
	# BC[3] x2, BC[0] x1 (its 200 bp fragment also reaches peak 3). Peak 3:
	# BC[0] x2 (the 800-1000 fragment and its own 950-1000), BC[2] x1.
	assert matrix == ["4\t4\t8",
		"0\t0\t2", "2\t0\t1", "3\t0\t2",
		"0\t1\t1", "1\t1\t2",
		"1\t2\t2", "3\t2\t1",
		"2\t3\t2"]


@pytest.mark.parametrize("cell_by_bin,bin_size,overlap,depth", [
	(False, 5000, 0, 3),
	(False, 5000, 50, 3),
	(False, 5000, 0, 2),
	(False, 5000, 25, 1),
	(True, 250, 0, 3),
	(True, 1000, 0, 3),
	(True, 600, 0, 3),
])
def test_feature_matrix_matches_reference(cpp, matrix_fasta, tmp_path, cell_by_bin, bin_size, overlap, depth):
	"""FeatureBarcodeMatrix::OutputFeatureMatrix writes the peaks (or bins),
	barcodes and matrix reference_feature_matrix computes, for the hand-made
	fragments plus 120 random ones (lengths 30-300, ending at least 10 bases
	before their contig's end, so no coverage run reaches it). With bins,
	fragments on the first reference cover up to three bins; fragments on
	the others stay within one bin. (Both restrictions avoid the known bug
	d-feature-matrix.)"""
	rng = random.Random(ch.stable_seed("test_d_unit", "matrix", cell_by_bin, bin_size))
	frags = list(HAND_FRAGMENTS)
	for i in range(120):
		rid = rng.randrange(len(MATRIX_CONTIGS))
		n = MATRIX_CONTIGS[rid][1]
		length = rng.randint(30, 300)
		start = rng.randint(5, n - length - 10)
		if cell_by_bin and rid > 0:
			b = rng.randrange(n // bin_size) if n // bin_size else 0
			length = min(length, bin_size - 1)
			start = b * bin_size + rng.randint(0, bin_size - length)
			if start + length > n - 10:
				continue
		frags.append((rid, start, length, rng.choice(BC)))
	got = _feature_matrix(cpp, matrix_fasta, tmp_path, frags, cell_by_bin, bin_size, overlap, depth)
	want = reference_feature_matrix(frags, cell_by_bin, bin_size, overlap, depth)
	assert got[0] == want[0]
	assert got[1] == want[1]
	assert got[2] == want[2]


# ------------------------------------------------------------ barcode output

def _barcode_output(cpp, tmp_path, length, table, barcodes):
	prefix = tmp_path / ("bo_%d" % len(list(tmp_path.iterdir())))
	cpp("d_units").probe("barcode_output", prefix, length, table, *barcodes)
	return open(str(prefix) + "_barcode.tsv").read().splitlines()


def test_matrix_barcodes_without_table(cpp, tmp_path):
	"""Without a translation table each barcode is written as its own
	sequence (of SetBarcodeLength bases) followed by -1."""
	assert _barcode_output(cpp, tmp_path, 8, "-", ["ACGTACGT", "TTTTGGGG"]) == ["ACGTACGT-1", "TTTTGGGG-1"]
	bcs = ch.make_whitelist(5, length=16, seed=3)
	assert _barcode_output(cpp, tmp_path, 16, "-", bcs) == ["%s-1" % b for b in bcs]


def test_matrix_barcodes_translated(cpp, tmp_path):
	"""--barcode-translate tables (README: "two-column tsv/csv", translated
	barcode first): a tab-separated table of whole barcodes, and a gzipped
	comma-separated table of 8 bp segments, where a 16 bp barcode is
	translated segment by segment and the IDs joined with "-"."""
	bcs = ch.make_whitelist(6, length=16, seed=5)
	tsv = ch.write_lines(tmp_path / "whole.tsv", ["CELL%d\t%s" % (i, b) for i, b in enumerate(bcs)])
	assert _barcode_output(cpp, tmp_path, 16, tsv, [bcs[3], bcs[0]]) == ["CELL3-1", "CELL0-1"]
	halves = sorted(set(b[:8] for b in bcs) | set(b[8:] for b in bcs))
	ids = dict((h, "S%02d" % i) for i, h in enumerate(halves))
	gz = tmp_path / "halves.csv.gz"
	with gzip.open(str(gz), "wt") as f:
		for h in halves:
			f.write("%s,%s\n" % (ids[h], h))
	want = ["%s-%s-1" % (ids[b[:8]], ids[b[8:]]) for b in bcs]
	assert _barcode_output(cpp, tmp_path, 16, gz, bcs) == want


def test_translation_of_unknown_barcode_exits(cpp, tmp_path):
	"""A barcode missing from the translation table ends the program with a
	non-zero status and says so on stderr."""
	bcs = ch.make_whitelist(3, length=16, seed=9)
	tsv = ch.write_lines(tmp_path / "t.tsv", ["CELL%d\t%s" % (i, b) for i, b in enumerate(bcs[:2])])
	proc = cpp("d_units").run(["--probe", "barcode_output", tmp_path / "x", 16, tsv, bcs[0], bcs[2]])
	assert proc.returncode != 0
	assert "Barcode does not exist in the translation table." in proc.stderr.decode()
