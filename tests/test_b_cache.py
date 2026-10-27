"""Area B, CLI layer: the candidate cache in the mapping loop.

chromap updates its candidate cache (src/mmcache.hpp) after each batch of
500,000 reads (single-end) or pairs, and queries it for every read of the
next batch; a hit replaces candidate generation for that read. Only an input
of more than one batch reaches a filled cache, so these tests map 510 copies
of a small designed read set: the first batch fills the cache with reads
seen 500 times each and every read of the second batch is a copy of one of
them, also on the other strand (the cache answers for a read's reverse
complement too). Whatever the cache holds, each copy of a read must map
exactly as the read does: the expected lines come from the simulation truth
times 510, and a multimapper's 510 copies all land at the one copy chosen
for that read. 3 chromap processes of about 1-3 s each (the paired-end ones
take about 2 GB)."""

import collections

import pytest

import chromap_helpers as ch

COPIES = 510


def _write_copies(path, reads):
	"""FASTQ of COPIES back-to-back copies of `reads` (in order)."""
	block = "".join("@%s\n%s\n+\n%s\n" % (r.name, r.seq, r.qual) for r in reads)
	with open(str(path), "w") as f:
		for _ in range(COPIES):
			f.write(block)
	return path


def _bed_counts(run):
	"""Counter of (chrom, start, end, strand, MAPQ) over the BED lines of a
	run (half a million lines: distinct lines are counted first, then
	parsed)."""
	out = collections.Counter()
	for line, n in collections.Counter(run.lines()).items():
		f = line.split("\t")
		out[(f[0], int(f[1]), int(f[2]), f[5], int(f[4]))] += n
	return out


def _segment_loci(genome, r):
	"""(contig, start, end, strand) of every copy of the planted segment a
	segment read comes from, with the read's strand at each copy."""
	seg = genome.segments[r.extra["segment"]]
	used = [c for c in seg.copies if c.contig == r.contig and c.start <= r.start < c.end][0]
	flip = {"+": "-", "-": "+"}
	rel = r.strand if used.strand == "+" else flip[r.strand]
	return [(contig, start, end, rel if cc.strand == "+" else flip[rel])
		for (contig, start, end), cc in zip(r.extra["loci"], seg.copies)]


def test_single_end_cache_hits_map_like_the_reads(index, genome, tmp_path):
	"""510,000 single-end reads: 1020 distinct reads (420 unique loci read on
	both strands, so each "-" read is the reverse complement of a "+" read;
	90 unique reads with one or two substitutions; 48 dup2 and 42 rep8 reads
	at distinct offsets, both strands) repeated 510 times. With -q 0 every
	unique read gives 510 lines at its true interval with MAPQ 60, and each
	multimapper 510 lines at one of its copies with MAPQ 0, its copies in the
	cached second batch included."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_b_cache", "se"))
	reads = []
	for i in range(420):
		c, s = sim.unique_locus(50)
		reads.append(sim.se_read(c, s, 50, "+"))
		reads.append(sim.se_read(c, s, 50, "-"))
	for i in range(90):
		c, s = sim.unique_locus(50)
		reads.append(sim.se_read(c, s, 50, "+-"[i % 2], mismatches=[(5,), (2, 47)][i % 2]))
	multi = []
	for seg, offsets in (("dup2", range(0, 950, 40)), ("rep8", range(0, 350, 17))):
		for off in offsets:
			for strand in "+-":
				multi.append(sim.segment_read(seg, off, 50, strand=strand, copy=len(multi) % 2))
	reads += multi
	reads = list(sim.shuffled(ch.ReadSet(reads)))
	path = _write_copies(tmp_path / "reads.fq", reads)
	run = ch.run_chromap(index, {"r1": path}, "-q 0", name="se_cache")
	got = _bed_counts(run)
	assert run.stats()["Number of reads"] == COPIES * len(reads)
	unique = collections.Counter()
	for r in reads:
		if not r.extra.get("segment"):
			unique[r.interval() + (60,)] += COPIES
	assert {k: got.pop(k, 0) for k in unique} == dict(unique)
	for r in multi:
		counts = [got.pop(locus + (0,), 0) for locus in _segment_loci(genome, r)]
		assert sorted(counts)[-2:] == [0, COPIES], (r, counts)
	assert not got


def _cache_pairs(genome):
	"""1000 distinct pairs: 300 fragments sequenced on both strands (read 1
	of one is read 2 of the other) and 200 couples of fragments placed so
	that read 2 of one is the reverse complement of read 1 of the other."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_b_cache", "pe"))
	pairs = []
	for i in range(300):
		ins = sim.rng.randint(150, 400)
		c, s = sim.unique_locus(ins)
		pairs.append(sim.pair(c, s, ins, strand="+"))
		pairs.append(sim.pair(c, s, ins, strand="-"))
	for i in range(200):
		ins = sim.rng.randint(150, 300)
		c, s = sim.unique_locus(ins + 200)
		# a: read 1 forward at [s + 200, s + 250); b: the fragment [s, s + 250),
		# whose read 2 is the reverse complement of a's read 1.
		a = sim.pair(c, s + 200, ins, strand="+")
		b = sim.pair(c, s, 250, strand="+")
		assert b.r2.seq == ch.revcomp(a.r1.seq)
		pairs += [a, b]
	return list(sim.shuffled(ch.ReadSet(pairs)))


@pytest.fixture(scope="module")
def cache_pair_files(genome, tmp_path_factory):
	pairs = _cache_pairs(genome)
	d = tmp_path_factory.mktemp("cache_pairs")
	paths = {"r1": _write_copies(d / "r1.fq", [p.r1 for p in pairs]),
		"r2": _write_copies(d / "r2.fq", [p.r2 for p in pairs])}
	return pairs, paths


def test_paired_end_cache_hits_map_like_the_pairs(index, cache_pair_files):
	"""510,000 pairs, 1000 distinct pairs repeated 510 times: every pair
	gives 510 fragment lines with MAPQ 60, the cached second batch included;
	the summary counts all pairs and reports cache hits."""
	pairs, paths = cache_pair_files
	run = ch.run_chromap(index, paths, "-q 0 --summary summary.csv", name="pe_cache")
	got = _bed_counts(run)
	want = collections.Counter()
	for p in pairs:
		want[p.fragment() + (60,)] += COPIES
	assert got == want
	row = run.summary().rows[""]
	assert row["total"] == COPIES * len(pairs)
	assert 0 < row["cachehit"] <= row["total"]


def test_paired_end_cache_hits_with_split_alignment(index, cache_pair_files):
	"""The same 510,000 pairs with --split-alignment, where a mate without
	candidates of its own is not rescued from its partner's, and a smaller
	cache (--cache-size 2000000, so other slots): each pair's 510 copies all
	give its true fragment, with one MAPQ, the cached second batch
	included."""
	pairs, paths = cache_pair_files
	run = ch.run_chromap(index, paths, "-q 0 --split-alignment --cache-size 2000000", name="pe_cache_split")
	got = collections.defaultdict(collections.Counter)
	for (chrom, start, end, strand, mapq), n in _bed_counts(run).items():
		got[(chrom, start, end, strand)][mapq] += n
	want = collections.Counter()
	for p in pairs:
		want[p.fragment()] += COPIES
	assert {k: sum(v.values()) for k, v in got.items()} == dict(want)
	assert all(len(v) == 1 for v in got.values()), [v for v in got.values() if len(v) > 1][:3]
