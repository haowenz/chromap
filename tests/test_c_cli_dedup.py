"""Area C, CLI layer: single-end deduplication (--remove-pcr-duplicates),
the duplicate count column and its cap, the MAPQ filter after
deduplication, --Tn5-shift on single-end reads, single-end TagAlign, the
bulk summary CSV (--summary) and --low-mem on BED output.

Reads come from helpers_c designs whose duplicate groups are known in
advance; expected lines and counts follow from the simulation truth and
chromap's documented rules (tests/README.md): single-end duplicates share a
start whatever their strand or length, the count column is the group size
capped at 255, MAPQ filtering applies to the record kept for each group."""

import pytest

import chromap_helpers as ch
import helpers_c as hc


@pytest.fixture(scope="module")
def se_design(genome):
	"""single_end_dedup_design: singletons, identical, mixed-strand,
	mixed-length and multimapper groups, random and too-short reads."""
	return hc.single_end_dedup_design(genome, ch.stable_seed("test_c_cli_dedup", "se_design"))


@pytest.fixture(scope="module")
def se_dedup_run(index, se_design):
	"""In-memory BED with --remove-pcr-duplicates and --summary. The
	malformed --frip-est-params is read only by paired-end runs."""
	reads, _, _ = se_design
	return ch.run_chromap(index, reads, "--remove-pcr-duplicates --summary summary.csv --frip-est-params 1;2",
		shared=True)


@pytest.fixture(scope="module")
def se_dedup_lowmem_run(index, se_design):
	"""The same reads with --low-mem (BED: about 3 s and 12 GB)."""
	reads, _, _ = se_design
	return ch.run_chromap(index, reads, "--remove-pcr-duplicates --low-mem --summary summary.csv", shared=True)


def _passing_groups(groups):
	return [g for g in groups if g.kind != "multimapper"]


def _check_se_groups(run, groups):
	"""Every BED line is one designed group (MAPQ 60, count = group size,
	any member's interval and strand); every group with MAPQ 60 is written
	once; multimapper groups (MAPQ 0 < 30), random and too-short reads are
	not written; lines are sorted by reference (FASTA order) and start."""
	lines = run.lines()
	passing = _passing_groups(groups)
	seen = hc.match_groups(lines, passing, lambda l: (l.split("\t")[0], int(l.split("\t")[1])),
		lambda g: hc.single_end_group_lines(g, g.size))
	assert len(seen) == len(passing)
	order = [(run.index.genome.rid(l.split("\t")[0]), int(l.split("\t")[1])) for l in lines]
	assert order == sorted(order)
	return seen


def _expected_bulk_summary(reads, groups, unmapped):
	"""The bulk row (README.md, *Summary File*): total = every read;
	duplicate = reads beyond the first of each written group; unmapped =
	reads without a mapping (random and too short); lowmapq = reads in
	groups below the MAPQ threshold; no cache hits in a first batch
	(tests/README.md), so cachehit, fric, estfrip and numcacheslots are 0."""
	dup = sum(g.size - 1 for g in _passing_groups(groups))
	low = sum(g.size for g in groups if g.kind == "multimapper")
	return {"": dict(total=len(reads), duplicate=dup, unmapped=len(unmapped), lowmapq=low, cachehit=0, fric=0.0,
		estfrip=0.0, numcacheslots=0)}


def test_single_end_duplicate_groups_and_counts(se_dedup_run, se_design):
	"""--remove-pcr-duplicates on single-end reads writes one line per
	duplicate group: identical reads, + and - reads sharing their leftmost
	base, and reads of 50, 60 and 70 bp sharing a start each collapse to one
	line whose count is the group size; singletons keep count 1."""
	_, groups, _ = se_design
	seen = _check_se_groups(se_dedup_run, groups)
	counts = sorted(int(l.split("\t")[6]) for l in seen.values())
	assert counts.count(1) == 40 and max(counts) == 12
	# Identical duplicates print the one possible line.
	for i, g in enumerate(_passing_groups(groups)):
		if g.kind == "identical":
			r = g.reads[0]
			assert seen[i] == "%s\t%d\t%d\tN\t60\t%s\t%d" % (r.contig, r.start, r.end, r.strand, g.size)


def test_single_end_dedup_statistics(se_dedup_run, se_design):
	"""chromap reports the number of mappings it wrote ("Number of output
	mappings (passed filters)"): one per duplicate group passing -q 30, the
	number of lines in the output. (The intermediate counts it prints after
	deduplication are progress diagnostics and are not checked.)"""
	_, groups, _ = se_design
	uni = len(_passing_groups(groups))
	assert se_dedup_run.stats()["Number of output mappings (passed filters)"] == uni == len(se_dedup_run.lines())


def test_single_end_summary(se_dedup_run, se_design):
	"""--summary for bulk single-end data: one row with an empty barcode and
	the counts of _expected_bulk_summary; the header has the numcacheslots
	column."""
	reads, groups, unmapped = se_design
	s = se_dedup_run.summary()
	assert s.header == hc.SUMMARY_HEADER + ["numcacheslots"]
	assert s.rows == _expected_bulk_summary(reads, groups, unmapped)


def test_single_end_ignores_frip_params(se_dedup_run):
	"""--frip-est-params is parsed only by paired-end mapping
	(Chromap::MapPairedEndReads): the single-end run above with a value of
	two numbers succeeds and writes estfrip 0 (pins the current behaviour)."""
	assert se_dedup_run.returncode == 0
	assert se_dedup_run.summary()[""]["estfrip"] == 0.0


@pytest.mark.lowmem
def test_low_mem_single_end_dedup_matches_truth(se_dedup_lowmem_run, se_dedup_run, se_design):
	"""--low-mem writes its mappings to a temp file and merges them back:
	the same duplicate groups and counts as the truth, the same lines as the
	default mode for groups of identical reads (other groups may keep a
	different member, see known-bugs/low-mem-dedup-representative.md), the
	same summary, and no temp file left next to the output."""
	reads, groups, unmapped = se_design
	seen = _check_se_groups(se_dedup_lowmem_run, groups)
	default = _check_se_groups(se_dedup_run, groups)
	for i, g in enumerate(_passing_groups(groups)):
		if g.kind in ("identical", "singleton"):
			assert seen[i] == default[i]
	assert se_dedup_lowmem_run.summary().rows == _expected_bulk_summary(reads, groups, unmapped)
	assert sorted(p.name for p in se_dedup_lowmem_run.outdir.iterdir()) == ["input", "out.txt", "summary.csv"]


@pytest.fixture(scope="module")
def cap_design(genome):
	return hc.single_end_cap_design(genome, ch.stable_seed("test_c_cli_dedup", "cap_design"))


def test_tagalign_tn5_and_count_cap(index, cap_design):
	"""Single-end --TagAlign has the 7 BED columns. With --Tn5-shift a +
	read's start moves right by 4 and a - read's end left by 5 (chromap.1);
	a + read keeps its length, so its end moves by 4 too (pins the current
	behaviour). With --remove-pcr-duplicates groups of 300 and 256
	identical reads are written once with the count capped at 255. No + and
	- read share a start before or after shifting here, so the groups are
	the same whichever comes first."""
	reads, groups = cap_design
	run = ch.run_chromap(index, reads, "--TagAlign --Tn5-shift --remove-pcr-duplicates")

	def shifted(r):
		return (r.contig, r.start + 4) if r.strand == "+" else (r.contig, r.start)

	groups = [hc.Group(g.kind, g.reads, key=shifted(g.reads[0])) for g in groups]
	starts = [g.reads[0].start for g in groups]
	assert len(set(starts)) == len(groups) and len(set(g.key for g in groups)) == len(groups)
	seen = hc.match_groups(run.lines(), groups, lambda l: (l.split("\t")[0], int(l.split("\t")[1])),
		lambda g: hc.single_end_group_lines(g, min(255, g.size), shift=True))
	assert len(seen) == len(groups)
	assert sorted(int(l.split("\t")[6]) for l in seen.values())[-2:] == [255, 255]
	order = [(index.genome.rid(l.split("\t")[0]), int(l.split("\t")[1])) for l in run.lines()]
	assert order == sorted(order)
