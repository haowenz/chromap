"""Area C, CLI layer: barcoded (single-cell) paired-end output: the fragment
file (5-column BED), duplicate removal at cell level
(--remove-pcr-duplicates-at-cell-level) and at bulk level
(--remove-pcr-duplicates-at-bulk-level, the default for single-cell data,
applied by --low-mem), --Tn5-shift on fragments, and the per-barcode
summary CSV with its non-whitelist row.

The design (helpers_c.barcoded_pair_design) reads each of four fragments
six times with three whitelisted barcodes (3, 1 and 2 pairs) and adds
singletons, a pair with a barcode far from the whitelist, a repeat pair and
a random pair, so every expected line and count follows from the design."""

import collections

import pytest

import chromap_helpers as ch
import helpers_c as hc


@pytest.fixture(scope="module")
def bc_design(genome, whitelist):
	return hc.barcoded_pair_design(genome, ch.stable_seed("test_c_cli_barcodes", "design"), whitelist[0],
		tie_breaks=True)


@pytest.fixture(scope="module")
def cell_level_run(index, bc_design, whitelist):
	reads = bc_design[0]
	return ch.run_chromap(index, reads, ["--remove-pcr-duplicates", "--remove-pcr-duplicates-at-cell-level",
		"--summary", "summary.csv", "--barcode-whitelist", str(whitelist[1])], shared=True)


def _fragment_line(contig, start, end, barcode, count):
	return "%s\t%d\t%d\t%s\t%d" % (contig, start, end, barcode, count)


def test_cell_level_dedup_fragments(cell_level_run, bc_design):
	"""Barcoded BED output is the fragment file (chrom, start, end, barcode,
	count; README.md). At cell level pairs are duplicates only with the same
	fragment and barcode: each designed fragment gives three lines (counts 3,
	1 and 2 for its three barcodes), sorted by start, length and then
	barcode; singletons give count 1; the non-whitelisted, repeat and random
	pairs are not written."""
	reads, fragments, singles, bad, rep, unm = bc_design
	want = []
	for g in fragments + singles:
		counts = {}
		for p in g.reads:
			counts[p.barcode] = counts.get(p.barcode, 0) + 1
		for bc, n in counts.items():
			want.append((g.key, bc, n))
	genome = cell_level_run.index.genome
	want.sort(key=lambda x: (genome.rid(x[0][0]), x[0][1], x[0][2] - x[0][1], ch.encode_seed(x[1])))
	assert cell_level_run.lines() == [_fragment_line(k[0], k[1], k[2], bc, n) for k, bc, n in want]


def test_cell_level_summary_per_barcode(cell_level_run, bc_design):
	"""With a whitelist the summary has one row per barcode seen (rows in no
	particular order) plus 'non-whitelist': total = pairs with that barcode,
	duplicate = pairs beyond the first per fragment, unmapped = pairs without
	a mapping, lowmapq = pairs with MAPQ below 30; the non-whitelisted pair
	counts as total and unmapped of the 'non-whitelist' row."""
	reads, fragments, singles, bad, rep, unm = bc_design
	rows = {}

	def row(bc):
		return rows.setdefault(bc, dict(total=0, duplicate=0, unmapped=0, lowmapq=0, cachehit=0, fric=0.0,
			estfrip=0.0, numcacheslots=0))

	for g in fragments + singles:
		seen = set()
		for p in g.reads:
			r = row(p.barcode)
			r["total"] += 1
			if p.barcode in seen:
				r["duplicate"] += 1
			seen.add(p.barcode)
	row(rep.barcode).update(total=1, lowmapq=1)
	row(unm.barcode).update(total=1, unmapped=1)
	row("non-whitelist").update(total=1, unmapped=1)
	s = cell_level_run.summary()
	assert s.header == hc.SUMMARY_HEADER + ["numcacheslots"]
	assert s.rows == rows
	assert s.lines[-1].startswith("non-whitelist,")


@pytest.mark.lowmem
def test_low_mem_bulk_level_dedup_with_tn5(index, bc_design, whitelist):
	"""--low-mem with --remove-pcr-duplicates on single-cell data removes
	duplicates at bulk level (the default): the pairs of a fragment collapse
	into one line with the whole group's count, carrying the barcode with the
	most pairs in the group, ties broken by the barcode's abundance (its
	exact matches among the input barcodes; the rule in
	MappingWriter::FindBestMappingIndexFromDuplicates). The designed
	fragments cover the most pairs winning over abundance (2 pairs against
	1) and abundance breaking a tie. chromap counts at most two pairs per
	barcode here (known bug r2-bulk-dedup-barcode-count), so wherever a
	barcode has three pairs it is also the most abundant one. --Tn5-shift
	moves every fragment's start right by 4 and its end left by 5."""
	reads, fragments, singles, bad, rep, unm = bc_design
	run = ch.run_chromap(index, reads, ["--remove-pcr-duplicates", "--low-mem", "--Tn5-shift",
		"--barcode-whitelist", str(whitelist[1])])
	listed = set(whitelist[0])
	abundance = collections.Counter(p.barcode for p in reads if p.barcode in listed)
	want = []
	for g in fragments:
		pairs = collections.Counter(p.barcode for p in g.reads)
		rank = dict((b, (pairs[b], abundance[b])) for b in pairs)
		best = max(pairs, key=rank.get)
		assert list(rank.values()).count(rank[best]) == 1, g
		assert best == max(pairs, key=lambda b: (min(pairs[b], 2), abundance[b])), g
		want.append((g.key, best, g.size))
	assert set(g.kind for g in fragments) == {"fragment", "most_pairs", "abundance_tie"}
	for g in singles:
		want.append((g.key, g.reads[0].barcode, 1))
	genome = index.genome
	want.sort(key=lambda x: (genome.rid(x[0][0]), x[0][1], x[0][2] - x[0][1]))
	assert run.lines() == [_fragment_line(k[0], k[1] + 4, k[2] - 5, bc, n) for k, bc, n in want]
	assert sorted(p.name for p in run.outdir.iterdir()) == ["input", "out.txt"]
