"""Area D: Hi-C output in the 4DN pairs format (--pairs, --preset hic,
--pairs-natural-chr-order) and --chr-order (paired-end; single-end is kept
out, see known bug single-end-chr-order).

Read pairs are simulated as Hi-C produces them (helpers_d.hic_reads): plain
fragments, mates on two contigs, and mates far apart on one contig, with
every strand combination and mates in both orders. Expected pairs records
follow the 4DN pairs specification: readID, chrom1, pos1, chrom2, pos2,
strand1, strand2, positions 1-based at each read's 5' end, the two sides in
"upper triangle" order (lower chromosome rank first, then lower position)."""

import collections

import pytest

import chromap_helpers as ch
import helpers_d as hd

PAIRS_COLUMNS = "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2"
PARTIAL_ORDER = ["chrX", "chr3"]
NATURAL_ORDER = ["chr2", "chrX", "chr1", "chr3"]


@pytest.fixture(scope="module")
def hic(genome):
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_d_hic", "reads"))
	return hd.hic_reads(sim, n_fragment=40, n_trans=40, n_far=40)


@pytest.fixture(scope="module")
def order_files(tmp_path_factory):
	d = tmp_path_factory.mktemp("d_orders")
	return dict(partial=ch.write_lines(d / "partial.txt", PARTIAL_ORDER),
		natural=ch.write_lines(d / "natural.txt", NATURAL_ORDER))


@pytest.fixture(scope="module")
def pairs_run(index, hic):
	"""The hic preset's mapping without --low-mem (-e 4 -q 1
	--split-alignment), --pairs followed by --BED: the in-memory pairs path."""
	return ch.run_chromap(index, hic, "--split-alignment -e 4 -q 1 --pairs --BED --cache-size 2000000",
		output="out.pairs", shared=True)


@pytest.fixture(scope="module")
def hic_sam_run(index, hic, order_files):
	"""--preset hic --SAM on the same reads (README: Hi-C reads and SAM
	output), with --chr-order listing chrX and chr3."""
	return ch.run_chromap(index, hic, ["--preset", "hic", "--SAM", "--chr-order", str(order_files["partial"]),
		"--cache-size", "2000000"], output="out.sam", shared=True)


@pytest.fixture(scope="module")
def ordered_run(index, hic, order_files):
	return ch.run_chromap(index, hic, ["--preset", "hic", "--chr-order", str(order_files["partial"]),
		"--cache-size", "2000000"], output="out.pairs", shared=True)


@pytest.fixture(scope="module")
def natural_run(index, hic, order_files):
	"""The hic preset's mapping in memory with both order files; read names
	end in /1 and /2."""
	return ch.run_chromap(index, hic, ["--split-alignment", "-e", "4", "-q", "1", "--pairs",
		"--chr-order", str(order_files["partial"]), "--pairs-natural-chr-order", str(order_files["natural"]),
		"--cache-size", "2000000"], output="out.pairs", name_style="slash", shared=True)


def _header_and_records(run):
	lines = run.lines()
	n_header = 0
	while n_header < len(lines) and lines[n_header].startswith("#"):
		n_header += 1
	assert not any(x.startswith("#") for x in lines[n_header:]), "header lines after records"
	header, recs = run.pairs()
	assert header == lines[:n_header]
	return header, recs


def _expected_header(genome, flip_order):
	return (["## pairs format v1.0.0", "#shape: upper triangle"]
		+ ["#chromsize: %s %d" % (c, len(genome.contigs[c])) for c in flip_order] + [PAIRS_COLUMNS])


def _sam_mapq(run):
	_, sam = run.sam()
	return dict(((s.qname, s.is_read1), s.mapq) for s in sam if not s.is_secondary)


def _check_records(recs, reads, flip_order, sort_order, suffix=""):
	rank = dict((c, i) for i, c in enumerate(flip_order))
	want = collections.Counter()
	for p in reads:
		r = hd.expected_pairs_record(p, rank)
		want[(r[0] + suffix,) + r[1:]] += 1
	got = collections.Counter(tuple(r[:7]) for r in recs)
	assert got == want
	assert set(r.pair_type for r in recs) == {"UU"}
	si = dict((c, i) for i, c in enumerate(sort_order))
	keys = [(si[r.chrom1], si[r.chrom2], r.pos1, r.pos2) for r in recs]
	assert keys == sorted(keys)


def test_pairs_header(genome, pairs_run):
	"""The header is "## pairs format v1.0.0", "#shape: upper triangle", one
	"#chromsize: NAME LENGTH" per contig in FASTA order, then the #columns
	line naming the ten columns; nothing else and nothing after the
	records start."""
	header, _ = _header_and_records(pairs_run)
	assert header == _expected_header(genome, genome.names)


def test_pairs_records_from_truth(genome, hic, pairs_run):
	"""One record per read pair with readID (read 1's name), both sides at
	their reads' 1-based 5' ends and strands, flipped into upper-triangle
	order by FASTA rank and position, pair type UU; records sorted by
	chrom1, chrom2 (FASTA order), pos1, pos2. --pairs given before --BED
	still writes pairs (pairs outranks BED whatever the order)."""
	_, recs = _header_and_records(pairs_run)
	assert len(recs) == len(hic)
	_check_records(recs, hic, genome.names, genome.names)
	assert sum(1 for r in recs if r.chrom1 != r.chrom2) == 40


def test_pairs_mapq_columns(hic, pairs_run, hic_sam_run):
	"""mapq1 and mapq2 carry the pair's MAPQ (chromap writes the pair's
	MAPQ in both columns). Every simulated mate maps uniquely, so here the
	pair's MAPQ equals each mate's: for every pair, mapq1 and mapq2 equal the
	SAM MAPQ of read 1 and read 2 from --preset hic --SAM on the same reads
	(the reference order does not change a MAPQ), and every pair passes the
	preset's -q 1. Pairs whose mates differ in MAPQ are not tested: the 4DN
	format defines mapq1 and mapq2 per side, chromap writes the pair's
	value."""
	_, recs = _header_and_records(pairs_run)
	sam_mapq = _sam_mapq(hic_sam_run)
	assert len(sam_mapq) == 2 * len(hic)
	_check_mapq(recs, sam_mapq)
	assert all(r.mapq1 >= 1 for r in recs)


def test_hic_sam_with_chr_order(genome, hic, hic_sam_run):
	"""--preset hic --SAM writes SAM, not pairs (explicit output options
	after a preset override it). With --chr-order listing chrX and chr3 the
	@SQ lines come in the order chrX, chr3, then the unlisted chr1, chr2 in
	FASTA order; every mate of every pair is at its true contig and 1-based
	start (mates on different contigs included), and records are sorted by
	that order, then position."""
	header, sam = hic_sam_run.sam()
	order = hd.full_order(PARTIAL_ORDER, genome.names)
	assert [h for h in header if h.startswith("@SQ")] == [
		"@SQ\tSN:%s\tLN:%d" % (c, len(genome.contigs[c])) for c in order]
	got = collections.Counter((s.qname, s.is_read1, s.rname, s.pos) for s in sam)
	want = collections.Counter()
	for p in hic:
		want[(p.name, True, p.r1.contig, p.r1.start + 1)] += 1
		want[(p.name, False, p.r2.contig, p.r2.start + 1)] += 1
	assert got == want
	oi = dict((c, i) for i, c in enumerate(order))
	keys = [(oi[s.rname], s.pos) for s in sam]
	assert keys == sorted(keys)


def _check_mapq(recs, sam_mapq, suffix=""):
	for r in recs:
		name = r.read_id[:len(r.read_id) - len(suffix)]
		assert (r.mapq1, r.mapq2) == (sam_mapq[(name, True)], sam_mapq[(name, False)]), r


def test_pairs_with_chr_order(genome, hic, ordered_run, hic_sam_run):
	"""--preset hic (pairs written through --low-mem temporary files) with
	--chr-order listing chrX and chr3: the order is chrX, chr3, then the
	unlisted chr1, chr2 in FASTA order. Without --pairs-natural-chr-order
	this order also flips the pairs (chromap.1): #chromsize lines in it,
	records flipped and sorted by it, MAPQ columns as in the SAM run."""
	order = hd.full_order(PARTIAL_ORDER, genome.names)
	assert order == ["chrX", "chr3", "chr1", "chr2"]
	header, recs = _header_and_records(ordered_run)
	assert header == _expected_header(genome, order)
	_check_records(recs, hic, order, order)
	_check_mapq(recs, _sam_mapq(hic_sam_run))


def test_pairs_natural_chr_order(genome, hic, natural_run, hic_sam_run):
	"""--pairs-natural-chr-order chr2, chrX, chr1, chr3 flips the pairs and
	orders the #chromsize lines; the records stay sorted by the --chr-order
	order (chrX, chr3, chr1, chr2) of chrom1, then chrom2, pos1, pos2. With
	mate names ending in /1 and /2 the readID is read 1's name, /1 kept."""
	header, recs = _header_and_records(natural_run)
	assert header == _expected_header(genome, NATURAL_ORDER)
	_check_records(recs, hic, NATURAL_ORDER, hd.full_order(PARTIAL_ORDER, genome.names), suffix="/1")
	_check_mapq(recs, _sam_mapq(hic_sam_run), suffix="/1")


@pytest.mark.parametrize("args", ["--pairs", "--preset hic"])
def test_single_end_pairs_exits(index, args, tmp_path):
	"""Pairs output needs read pairs: single-end input with --pairs (or the
	hic preset) stops with status 255 and "No support for single-end HiC
	yet!"."""
	files = ch.ReadSimulator(index.genome, seed=5).se_reads(3).write(tmp_path)
	run = ch.run_cli(["-x", index.path, "-r", index.fasta, "-1", files["r1"], "-o", tmp_path / "out.pairs"]
		+ args.split(), cwd=tmp_path)
	assert run.returncode == 255
	assert "No support for single-end HiC yet!" in run.stderr


# ------------------------------------------------------------ --chr-order

def test_chr_order_unknown_name_exits(index, tmp_path):
	"""A name in the --chr-order file that is not a reference sequence stops
	chromap with status 255 and its message."""
	files = ch.ReadSimulator(index.genome, seed=7).pairs(3).write(tmp_path)
	order = ch.write_lines(tmp_path / "order.txt", ["chr2", "chrUn", "chr1"])
	run = ch.run_cli(["-x", index.path, "-r", index.fasta, "-1", files["r1"], "-2", files["r2"],
		"-o", tmp_path / "out.bed", "--chr-order", order], cwd=tmp_path)
	assert run.returncode == 255
	assert "ERROR: unknown chromsome names found in chromosome order file." in run.stderr
