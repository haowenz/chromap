"""Area B, CLI layer, single-end reads: positions, CIGARs, MAPQ, multimappers
and the mapping options -e, -s, -f, -q, -n, --drop-repetitive-reads, -A, -B,
-O, -E and -t.

Every read is designed with a known origin. Expected positions come from the
simulation truth; expected CIGARs from the affine-gap reference of
helpers_b.semi_global over the window chromap searches; expected MAPQ values
from the MAPQ rule written out in helpers_b.mapq_single_end; seed counts from
the read's minimizers by definition (helpers_b.seed_support). Runs are shared
through module fixtures: 11 chromap processes in all."""

import collections

import pytest

import chromap_helpers as ch
import helpers_b as hb

E = hb.DEFAULT_E


def _flip(strand):
	return "-" if strand == "+" else "+"


def _segment_loci(genome, r):
	"""[(contig, start, end, strand)] of every copy of the planted segment a
	segment read comes from, with the read's strand at each copy."""
	seg = genome.segments[r.extra["segment"]]
	used = [c for c in seg.copies if c.contig == r.contig and c.start <= r.start < c.end][0]
	rel = r.strand if used.strand == "+" else _flip(r.strand)
	out = []
	for (contig, start, end), cc in zip(r.extra["loci"], seg.copies):
		out.append((contig, start, end, rel if cc.strand == "+" else _flip(rel)))
	return out


def _other_copy_view(genome, r):
	"""The dupmm read r seen from the other copy: same sequence, placed at the
	other copy, with the copies' substitution inside it as a mismatch."""
	seg = genome.segments["dupmm"]
	i = 0 if seg.copies[0].contig == r.contig else 1
	here, there = seg.copies[i], seg.copies[1 - i]
	off = r.start - here.start
	subs = [there.start + o for o in seg.substitutions if off <= o < off + len(r.seq)]
	return ch.SimRead(r.name, r.seq, r.qual, there.contig, there.start + off, there.start + off + len(r.seq), r.strand,
		mismatches=subs)


def _unique_sub_read(sim, genome, length, offsets, strand, min_hits=2):
	"""A read at a unique locus with substitutions at `offsets` that keeps at
	least `min_hits` seeds at its locus."""
	while True:
		c, s = sim.unique_locus(length)
		r = sim.se_read(c, s, length, strand, mismatches=offsets)
		if hb.seed_support(genome, r)[1] >= min_hits:
			return r


def _indel_reads(sim, genome, n, name):
	"""Reads with one 1-3 bp insertion or deletion at read offsets 15-35 at
	unique loci, kept when the draft mapping's end is the true end for every
	candidate chromap may pick (the start shifted by the indel or not) and the
	affine-gap optimum in chromap's SAM window is unique and equals the true
	CIGAR at the true start."""
	rng = sim.rng
	out = []
	while len(out) < n:
		kind = "DI"[len(out) % 2]
		size = rng.randint(1, 3)
		off = rng.randint(15, 35)
		strand = "+-"[(len(out) // 2) % 2]
		c, s = sim.unique_locus(60)
		ref = genome.seq(c, s, s + 60)
		if kind == "D" and not hb.indel_placement_unique(ref, off, kind, size):
			continue
		r = hb.indel_read(genome, c, s, 50, strand, kind, off, size, "%s%02d" % (name, len(out)))
		base = s if strand == "+" else r.end - 1
		drafts = {hb.draft_mapping(genome.contigs[c], r.seq, strand, base + d, E) for d in (-size, 0, size)}
		if drafts != {(size, r.end - 1)}:
			continue
		k = hb.expected_ksw(genome, c, r.end - 1, r.extra["fwd"], E)
		if k["count"] != 1 or k["alignments"][0][2] != r.extra["cigar"] or k["window_start"] + k["alignments"][0][0] != s:
			continue
		r.kind = "indel"
		out.append(r)
	return out


@pytest.fixture(scope="module")
def se_reads(genome):
	"""Single-end reads, 50 bp, both strands, labelled by r.kind:
	exact (16 unique exact), subs (16 unique with 1-3 substitutions at read
	offsets 8-41 keeping >= 2 seeds), indel (12, see _indel_reads), dupmm (8
	from either dupmm copy covering one of the copies' 12 differences at read
	offset 10-40, with >= 2 seeds at the other copy), dup2 (4), rep8 (4) and
	random (6, not from the genome)."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_b_single_end", "se_reads"), name_prefix="a")
	rng = sim.rng
	reads = []
	for i in range(16):
		c, s = sim.unique_locus(50)
		r = sim.se_read(c, s, 50, "+-"[i % 2])
		r.kind = "exact"
		reads.append(r)
	for i in range(16):
		r = _unique_sub_read(sim, genome, 50, sorted(rng.sample(range(8, 42), 1 + i % 3)), "+-"[i % 2])
		r.kind = "subs"
		reads.append(r)
	reads += _indel_reads(sim, genome, 12, "ind")
	seg = genome.segments["dupmm"]
	while sum(1 for r in reads if r.kind == "dupmm") < 8:
		i = sum(1 for r in reads if r.kind == "dupmm")
		sub = rng.choice(seg.substitutions)
		off = sub - rng.randint(10, 40)
		if off < 0 or off + 50 > seg.length:
			continue
		r = sim.segment_read("dupmm", off, 50, strand="+-"[i % 2], copy=i % 2)
		if hb.seed_support(genome, _other_copy_view(genome, r))[1] < 2:
			continue
		r.kind = "dupmm"
		reads.append(r)
	for i in range(4):
		r = sim.segment_read("dup2", rng.randint(0, 950), 50, strand="+-"[i % 2], copy=i % 2)
		r.kind = "dup2"
		reads.append(r)
	for i in range(4):
		r = sim.segment_read("rep8", rng.randint(0, 350), 50, strand="+-"[i % 2], copy=(3 * i) % 8)
		r.kind = "rep8"
		reads.append(r)
	for i in range(6):
		reads.append(sim.random_read())
	return sim.shuffled(ch.ReadSet(reads))


UNIQUE = ("exact", "subs", "indel")


def _errors(r):
	if r.kind == "indel":
		return hb.parse_cigar(r.extra["cigar"])[1][0]
	return len(r.mismatches)


def _expected_mapq(r, e=E):
	"""MAPQ of a read's best mapping (helpers_b.mapq_single_end): unique reads
	have no second-best mapping (second-best errors start at e + 1); a dupmm
	read has its other copy one error worse; repeats of equal copies get 0."""
	if r.kind in UNIQUE:
		span = r.end - r.start
		return hb.mapq_single_end(_errors(r), span, len(r.seq), e, e + 1, 1, 0, 0)
	if r.kind == "dupmm":
		return hb.mapq_single_end(0, 50, 50, e, 1, 1, 1, 0)
	return 0


@pytest.fixture(scope="module")
def se_sam(index, se_reads):
	return ch.run_chromap(index, se_reads, "--SAM -q 0 -n 8", name="se_sam", shared=True)


@pytest.fixture(scope="module")
def se_bed(index, se_reads):
	return ch.run_chromap(index, se_reads, name="se_bed", shared=True)


def test_single_end_sam_records_match_truth(genome, se_reads, se_sam):
	"""Unique reads (exact, substitutions, indels) give one SAM record each:
	FLAG 0 or 16 by strand, POS the true 1-based start, CIGAR 50M or the
	true indel CIGAR, SEQ in reference orientation, NM and MD computed from
	the reference and that CIGAR (SAM specification), MAPQ 60 (the rule of
	helpers_b.mapq_single_end with no second-best mapping caps at 60).
	Random reads give no record."""
	header, recs = se_sam.sam()
	by_name = collections.defaultdict(list)
	for rec in recs:
		by_name[rec.qname].append(rec)
	checked = 0
	for r in se_reads:
		got = by_name.get(r.name, [])
		if r.kind == "random":
			assert got == [], r
			continue
		if r.kind not in UNIQUE:
			continue
		assert len(got) == 1, r
		rec = got[0]
		cigar = r.extra["cigar"] if r.kind == "indel" else "50M"
		fwd = r.extra["fwd"] if r.kind == "indel" else (r.seq if r.strand == "+" else ch.revcomp(r.seq))
		nm, md = hb.nm_md(genome.seq(r.contig, r.start, r.end), fwd, cigar)
		assert (rec.rname, rec.pos, rec.cigar, rec.flag) == (r.contig, r.start + 1, cigar,
			0 if r.strand == "+" else 16), r
		assert rec.seq == fwd, r
		assert (rec.tags["NM"], rec.tags["MD"]) == (nm, md), r
		assert nm == _errors(r), r
		assert rec.mapq == _expected_mapq(r) == 60, r
		checked += 1
	assert checked == 44


def test_single_end_second_best_lowers_mapq(genome, se_reads, se_sam):
	"""A dupmm read matches its own copy exactly and the other copy with one
	mismatch: one record at its own copy with MAPQ
	int(30.1 * 1 * (3 / ln 50)^2 + 0.499) - int(4.343 ln 2 + 0.499) = 18 - 3 =
	15, the second-best penalty included."""
	header, recs = se_sam.sam()
	by_name = collections.defaultdict(list)
	for rec in recs:
		by_name[rec.qname].append(rec)
	reads = [r for r in se_reads if r.kind == "dupmm"]
	assert len(reads) == 8
	for r in reads:
		got = by_name[r.name]
		assert [(x.rname, x.pos, x.flag, x.mapq, x.cigar) for x in got] == [
			(r.contig, r.start + 1, 0 if r.strand == "+" else 16, 15, "50M")], r
		assert _expected_mapq(r) == 15


def test_multimappers_reported_at_every_locus(genome, se_reads, se_sam):
	"""With -n 8 and -q 0, a read from dup2 (2 identical copies) or rep8 (8)
	is reported at every copy, on the strand the copy implies, with MAPQ 0
	(several best mappings); exactly one record per read is primary, the
	others carry FLAG 256 (secondary)."""
	header, recs = se_sam.sam()
	by_name = collections.defaultdict(list)
	for rec in recs:
		by_name[rec.qname].append(rec)
	n = 0
	for r in se_reads:
		if r.kind not in ("dup2", "rep8"):
			continue
		got = by_name[r.name]
		want = sorted((c, s + 1, st) for c, s, e, st in _segment_loci(genome, r))
		assert sorted((x.rname, x.pos, "-" if x.is_reverse else "+") for x in got) == want, r
		assert {x.mapq for x in got} == {0}, r
		assert sum(1 for x in got if not x.is_secondary) == 1, r
		n += 1
	assert n == 8


def test_single_end_mapping_statistics(se_reads, se_sam, chromap_threads):
	"""The statistics on stderr follow from the design: every read but the
	random ones maps; unique reads and dupmm reads have one best mapping;
	dup2 and rep8 reads have 2 and 8, all reported with -n 8. The per-thread
	counters (mapped, unique, mappings) are checked at -t 1 only: clang/libomp
	builds at -t > 1 undercount them (known bug
	a-mapping-statistics-clang-threads); the read count and the number of
	output mappings are checked at any -t."""
	kinds = collections.Counter(r.kind for r in se_reads)
	unique = kinds["exact"] + kinds["subs"] + kinds["indel"] + kinds["dupmm"]
	multi = kinds["dup2"] + kinds["rep8"]
	mappings = unique + 2 * kinds["dup2"] + 8 * kinds["rep8"]
	stats = se_sam.stats()
	assert stats["Number of reads"] == len(se_reads)
	assert stats["Number of output mappings (passed filters)"] == mappings
	if chromap_threads != 1:
		return
	assert stats["Number of mapped reads"] == unique + multi
	assert stats["Number of uniquely mapped reads"] == unique
	assert stats["Number of reads have multi-mappings"] == multi
	assert stats["Number of mappings"] == mappings
	assert stats["Number of uni-mappings"] == unique
	assert stats["Number of multi-mappings"] == mappings - unique


def test_single_end_bed_default_threshold(se_reads, se_bed):
	"""Default BED output (-q 30, -n 1): each unique read once at its true
	aligned interval [start, end) (an indel read spans 50 + deleted - inserted
	reference bases) with its strand, MAPQ 60 and count 1. dupmm reads (MAPQ
	15), dup2 and rep8 reads (MAPQ 0) fall below -q 30 and random reads do not
	map."""
	want = collections.Counter((r.contig, r.start, r.end, r.strand, 60) for r in se_reads if r.kind in UNIQUE)
	got = collections.Counter((b.chrom, b.start, b.end, b.strand, b.mapq) for b in se_bed.bed())
	assert got == want
	assert {b.dups for b in se_bed.bed()} == {1}


CAP_ARGS = "--SAM -n 5 --drop-repetitive-reads 3 -q 0"


def test_max_num_best_mappings_capped_by_drop_repetitive_reads(index, genome, se_reads):
	"""-n 5 with --drop-repetitive-reads 3: chromap warns that only reads
	with at most 3 best mappings will be output. Reads with at most 3 are
	unaffected: a dup2 read is reported at both copies (one record primary),
	unique and dupmm reads once at their true position, with MAPQ as in the
	default run. (Single-end reads with more than 3 best mappings, the rep8
	reads, are still reported at 949043c: known bug
	r1-single-end-drop-repetitive-reads, not asserted here.)"""
	run = ch.run_chromap(index, se_reads, CAP_ARGS, name="se_cap", shared=True)
	assert ("WARNING: you want to drop mapped reads with more than 3 mappings. But you want to output top 5 "
		"best mappings. In this case, only reads with <=3 best mappings will be output.") in run.stderr
	header, recs = run.sam()
	by_name = collections.defaultdict(list)
	for rec in recs:
		by_name[rec.qname].append(rec)
	n = 0
	for r in se_reads:
		got = by_name.get(r.name, [])
		if r.kind == "random":
			assert got == [], r
			continue
		if r.kind == "rep8":
			continue
		if r.kind == "dup2":
			loci = sorted((c, s + 1, st) for c, s, e, st in _segment_loci(genome, r))
		else:
			loci = [(r.contig, r.start + 1, r.strand)]
		assert sorted((x.rname, x.pos, "-" if x.is_reverse else "+") for x in got) == loci, r
		assert {x.mapq for x in got} == {_expected_mapq(r)}, r
		assert sum(1 for x in got if not x.is_secondary) == 1, r
		n += 1
	assert n == len(se_reads) - 6 - 4


# --------------------------------------------- -e, -s, -q and MAPQ values

SUB_OFFSETS = (20, 45, 70, 95, 120)


@pytest.fixture(scope="module")
def threshold_reads(genome):
	"""150 bp reads with k = 0..5 substitutions at read offsets 20, 45, ...
	(4 per k, both strands; every gap keeps intact 17-mers); 50 bp reads with
	two substitutions placed so that exactly 1, 2, 3 or 4 of their minimizers
	still hit their locus (2 each, kind "seeds"); 4 dupmm reads."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_b_single_end", "threshold_reads"), name_prefix="t")
	rng = sim.rng
	reads = []
	for k in range(6):
		for i in range(4):
			r = _unique_sub_read(sim, genome, 150, SUB_OFFSETS[:k], "+-"[i % 2], min_hits=4)
			r.kind = "e%d" % k
			reads.append(r)
	for target in (1, 1, 2, 2, 3, 3, 4, 4):
		while True:
			c, s = sim.unique_locus(50)
			r = sim.se_read(c, s, 50, rng.choice("+-"), mismatches=(rng.randint(9, 14), rng.randint(34, 41)))
			nmin, hits = hb.seed_support(genome, r)
			if hits == target and nmin > 4:
				r.kind = "seeds"
				r.extra["hits"] = hits
				reads.append(r)
				break
	seg = genome.segments["dupmm"]
	while sum(1 for r in reads if r.kind == "dupmm") < 4:
		i = sum(1 for r in reads if r.kind == "dupmm")
		off = rng.choice(seg.substitutions) - rng.randint(10, 40)
		if off < 0 or off + 50 > seg.length:
			continue
		r = sim.segment_read("dupmm", off, 50, strand="+-"[i % 2], copy=i % 2)
		if hb.seed_support(genome, _other_copy_view(genome, r))[1] >= 2:
			r.kind = "dupmm"
			reads.append(r)
	return sim.shuffled(ch.ReadSet(reads))


@pytest.mark.parametrize("e,s,q", [(3, 1, 10), (4, 3, 16)])
def test_error_threshold_seeds_and_mapq_threshold(index, threshold_reads, chromap_threads, e, s, q):
	"""A read maps when its best alignment has at most -e errors and its
	locus is supported by at least min(-s, #minimizers) seeds; it is written
	when its MAPQ is at least -q. With -e 3 the reads with 4 and 5
	substitutions do not map, with -e 4 those with 5; with -s 1 every
	"seeds" read maps, with -s 3 only those keeping 3 or 4 seeds. MAPQ comes
	from helpers_b.mapq_single_end, which depends on -e: with e errors
	allowed and k made, the second-best error count is taken as min(e + 1,
	k + e). With -e 3 -q 10 the reads with 3 substitutions have MAPQ exactly
	10 and are written; with -q 16, dupmm reads (MAPQ 15) and reads with 4
	substitutions at -e 4 (MAPQ 10) map but are not written. "Number of
	mapped reads" (before the MAPQ filter) is checked at -t 1 only (known bug
	a-mapping-statistics-clang-threads)."""
	run = ch.run_chromap(index, threshold_reads, "-e %d -s %d -q %d" % (e, s, q), name="se_e%d" % e)
	want = collections.Counter()
	mapped = 0
	for r in threshold_reads:
		if r.kind == "dupmm":
			ok, mapq = True, hb.mapq_single_end(0, 50, 50, e, 1, 1, 1, 0)
		else:
			k = len(r.mismatches)
			seeds_ok = r.kind != "seeds" or r.extra["hits"] >= s
			ok = k <= e and seeds_ok
			mapq = hb.mapq_single_end(k, len(r.seq), len(r.seq), e, e + 1, 1, 0, 0)
		if ok:
			mapped += 1
			if mapq >= q:
				want[(r.contig, r.start, r.end, r.strand, mapq)] += 1
	got = collections.Counter((b.chrom, b.start, b.end, b.strand, b.mapq) for b in run.bed())
	assert got == want
	if chromap_threads == 1:
		assert run.stats()["Number of mapped reads"] == mapped
	# The design covers both sides of each threshold: reads with e + 1
	# substitutions, a MAPQ equal to -q (written), and mapped reads below it.
	assert any(r.kind == "e%d" % (e + 1) for r in threshold_reads)
	if q == 10:
		assert hb.mapq_single_end(e, 150, 150, e, e + 1, 1, 0, 0) == q
		assert any(m == q for (_, _, _, _, m) in want)
	else:
		assert mapped > sum(want.values())


# ------------------------------------------------- -f and repetitive seeds

@pytest.fixture(scope="module")
def repeat_reads(genome):
	"""Reads inside rep8 (4) and reads that start in unique sequence and end
	inside a forward rep8 copy, or the reverse, with 0-3 substitutions in the
	rep8 part (24 of 50 bp, kind "boundary"); 8 such reads of 100 bp with 3
	or 5 substitutions, whose identities 0.97 and 0.95 are the upper ends of
	two MAPQ identity bands; 4 unique exact reads."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_b_single_end", "repeat_reads"), name_prefix="p")
	rng = sim.rng
	reads = []
	for i in range(4):
		r = sim.segment_read("rep8", rng.randint(0, 350), 50, strand="+-"[i % 2], copy=(5 * i) % 8)
		r.kind = "rep8"
		reads.append(r)
	copies = [c for c in genome.segments["rep8"].copies if c.strand == "+"]
	counts = hb.segment_minimizer_counts(genome, "rep8")
	i = -1
	while sum(1 for r in reads if r.kind == "boundary") < 24:
		i += 1
		c = copies[i % len(copies)]
		d = rng.randint(12, 20)
		k = i % 4
		# Substitutions sit next to each other at the far end of the rep8
		# part, so a stretch of intact rep8 seeds remains.
		if i % 2 == 0:
			# [start, start + 50): d unique bases, then 50 - d in the copy
			start = c.start - d
			offs = [41 - j for j in range(k)]
		else:
			start = c.end - 50 + d
			offs = [8 + j for j in range(k)]
		strand = rng.choice("+-")
		read_offs = offs if strand == "+" else [49 - o for o in offs]
		r = sim.se_read(c.contig, start, 50, strand, mismatches=read_offs)
		if _unique_seed_hits(genome, r, counts) < 2:
			continue
		r.kind = "boundary"
		reads.append(r)
	i = -1
	while sum(1 for r in reads if len(r.seq) == 100) < 8:
		i += 1
		c = copies[i % len(copies)]
		d = rng.randint(12, 20)
		k = (3, 5)[(i // 2) % 2]
		if i % 2 == 0:
			start = c.start - d
			offs = [91 - j for j in range(k)]
		else:
			start = c.end - 100 + d
			offs = [8 + j for j in range(k)]
		strand = rng.choice("+-")
		read_offs = offs if strand == "+" else [99 - o for o in offs]
		r = sim.se_read(c.contig, start, 100, strand, mismatches=read_offs)
		if _unique_seed_hits(genome, r, counts) < 2:
			continue
		# Adjacent substitutions can be cheaper as an indel pair; keep reads
		# whose draft (edit distance) error count is the substitution count.
		cand = r.start if strand == "+" else r.end - 1
		if hb.draft_mapping(genome.contigs[c.contig], r.seq, strand, cand, E) != (k, r.end - 1):
			continue
		r.kind = "boundary"
		reads.append(r)
	for i in range(4):
		ctg, s = sim.unique_locus(50)
		r = sim.se_read(ctg, s, 50, "+-"[i % 2])
		r.kind = "exact"
		reads.append(r)
	return ch.ReadSet(reads)


def _unique_seed_hits(genome, r, counts):
	"""Seeds of an indel-free read that hit its locus through minimizers
	occurring fewer than 5 times (helpers_b.seed_support restricted to them)."""
	mins = hb.read_minimizers(r.seq)
	ref = hb.locus_minimizers(genome, r.contig, r.start, r.end)
	last = len(r.seq) - 1
	hits = 0
	for h, p, _ in mins:
		if counts.get(h, 0) >= 5:
			continue
		g0 = r.start + p - hb.K + 1 if r.strand == "+" else r.start + last - p
		if any(g0 <= m < g0 + hb.K for m in r.mismatches):
			continue
		ref_end = r.start + p if r.strand == "+" else r.start + last - p + hb.K - 1
		hits += (h, ref_end) in ref
	return hits


def _boundary_mapq(genome, r, counts, f0):
	"""MAPQ of a read whose minimizers inside rep8 occur at least f0 times:
	they are skipped and counted as repetitive, and the repetitive seed length
	(helpers_b.repetitive_seed_length) scales the MAPQ of its unique
	mapping."""
	rep = [p for h, p, _ in hb.read_minimizers(r.seq) if counts.get(h, 0) >= f0]
	n = len(r.seq)
	return hb.mapq_single_end(len(r.mismatches), n, n, E, E + 1, 1, 0, hb.repetitive_seed_length(rep))


def test_seed_frequency_second_round_and_repetitive_mapq(index, genome, repeat_reads):
	"""-f 5,1000: rep8 minimizers occur 8 times, at or above the first value,
	so a rep8 read has no seed in the first round; the second round (limit
	1000) finds all 8 copies and the read is reported at each with -n 8. A
	read half in unique sequence and half in a rep8 copy maps from its unique
	seeds; its rep8 seeds are counted as repetitive and lower its MAPQ by the
	repetitive seed length (helpers_b._boundary_mapq; with identity >= 0.999
	the factor is 1 - f^4, then 1 - f^2, 1 - f and 1 - sqrt(f) as
	substitutions lower the identity, the bands closed at 0.97 and 0.95:
	identity 0.97 still gets 1 - f and 0.95 1 - sqrt(f))."""
	run = ch.run_chromap(index, repeat_reads, "--SAM -f 5,1000 -n 8 -q 0", name="se_f5_1000")
	header, recs = run.sam()
	by_name = collections.defaultdict(list)
	for rec in recs:
		by_name[rec.qname].append(rec)
	counts = hb.segment_minimizer_counts(genome, "rep8")
	lowered = set()
	for r in repeat_reads:
		got = by_name[r.name]
		if r.kind == "rep8":
			want = sorted((c, s + 1) for c, s, e, st in _segment_loci(genome, r))
			assert sorted((x.rname, x.pos) for x in got) == want, r
			assert {x.mapq for x in got} == {0}, r
			continue
		mapq = _boundary_mapq(genome, r, counts, 5) if r.kind == "boundary" else 60
		assert [(x.rname, x.pos, x.mapq) for x in got] == [(r.contig, r.start + 1, mapq)], r
		if mapq < 60:
			lowered.add((len(r.seq), len(r.mismatches)))
	assert lowered == {(50, 0), (50, 1), (50, 2), (50, 3), (100, 3), (100, 5)}


def test_seed_frequency_both_rounds_skip(index, genome, repeat_reads, chromap_threads):
	"""-f 5,6: rep8 minimizers (8 occurrences) are too frequent for both
	rounds, so rep8 reads do not map at all; the unique and half-unique reads
	are unaffected (same MAPQ as with -f 5,1000)."""
	run = ch.run_chromap(index, repeat_reads, "-f 5,6 -n 8 -q 0", name="se_f5_6")
	counts = hb.segment_minimizer_counts(genome, "rep8")
	want = collections.Counter()
	for r in repeat_reads:
		if r.kind == "rep8":
			continue
		mapq = _boundary_mapq(genome, r, counts, 5) if r.kind == "boundary" else 60
		want[(r.contig, r.start, r.end, r.strand, mapq)] += 1
	got = collections.Counter((b.chrom, b.start, b.end, b.strand, b.mapq) for b in run.bed())
	assert got == want
	if chromap_threads == 1:
		assert run.stats()["Number of mapped reads"] == sum(1 for r in repeat_reads if r.kind != "rep8")


# ------------------------------------------------------ -A, -B, -O, -E

SCORING = hb.ScoreParams(2, 2, 12, 3, 12, 3)
SCORING_ARGS = "--SAM -A 2 -B 2 -O 12,12 -E 3,3"


@pytest.fixture(scope="module")
def scoring_reads(genome):
	"""Reads with a 1-4 bp indel 8-14 bases from either end whose best
	alignment under -A 2 -B 2 -O 12,12 -E 3,3 (helpers_b.semi_global in
	chromap's SAM window) is unique and changes when any one of the four
	options alone goes back to its default: at least two reads per
	option."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_b_single_end", "scoring_reads"), name_prefix="s")
	rng = sim.rng
	variants = {"A": SCORING._replace(match=1), "B": SCORING._replace(mismatch=4),
		"O": SCORING._replace(o_ins=6, o_del=6), "E": SCORING._replace(e_ins=1, e_del=1)}
	picked, cover = [], collections.Counter()
	for _ in range(5000):
		if len(cover) == 4 and min(cover.values()) >= 2:
			break
		kind = rng.choice("DI")
		size = rng.randint(1, 4)
		off = rng.choice([rng.randint(8, 14), rng.randint(36, 42)])
		strand = rng.choice("+-")
		c, s = sim.unique_locus(60)
		if kind == "D" and not hb.indel_placement_unique(genome.seq(c, s, s + 60), off, kind, size):
			continue
		r = hb.indel_read(genome, c, s, 50, strand, kind, off, size, "sc%02d" % len(picked))
		base = s if strand == "+" else r.end - 1
		drafts = {hb.draft_mapping(genome.contigs[c], r.seq, strand, base + d, E) for d in (-size, 0, size)}
		if len(drafts) != 1:
			continue
		errors, end = drafts.pop()
		if errors > E:
			continue
		res = hb.expected_ksw(genome, c, end, r.extra["fwd"], E, SCORING)
		if res["count"] != 1:
			continue
		best = res["alignments"][0]
		flips = [k for k, p in sorted(variants.items())
			if hb.expected_ksw(genome, c, end, r.extra["fwd"], E, p)["alignments"][:1] not in ([best], [])]
		if not flips or all(cover[f] >= 2 for f in flips):
			continue
		cover.update(flips)
		r.kind = "scoring"
		r.extra["expected"] = (res["window_start"] + best[0], best[2])
		picked.append(r)
	assert len(cover) == 4 and min(cover.values()) >= 2, cover
	return ch.ReadSet(picked)


def test_scoring_options_choose_the_alignment(index, genome, scoring_reads):
	"""-A, -B, -O and -E set the match score, mismatch penalty and gap
	open/extension penalties of the SAM alignment: each read's POS and CIGAR
	are the unique best alignment under those values (an indel near a read
	end becomes mismatches when gaps are dear enough), and NM/MD follow the
	CIGAR. Every option changes at least two of the expected alignments when
	set back alone, so each must reach the aligner."""
	run = ch.run_chromap(index, scoring_reads, SCORING_ARGS, name="se_scoring")
	header, recs = run.sam()
	assert sorted(rec.qname for rec in recs) == sorted(r.name for r in scoring_reads)
	for rec in recs:
		r = run.truth[rec.qname]
		start, cigar = r.extra["expected"]
		fwd = r.extra["fwd"]
		span = hb.cigar_spans(cigar)[1]
		nm, md = hb.nm_md(genome.seq(r.contig, start, start + span), fwd, cigar)
		assert (rec.rname, rec.pos, rec.cigar) == (r.contig, start + 1, cigar), r
		assert (rec.tags["NM"], rec.tags["MD"]) == (nm, md), r


# ----------------------------------------------------------- threads

def test_single_end_thread_count_does_not_change_output(index, genome, se_reads):
	"""Single-end mapping runs t * t OpenMP tasks, so -t 4 spreads even a
	hundred reads over threads. SAM output with -n 3 -q 0 (rep8 reads, 8 best
	mappings, reported at 3 distinct copies drawn by reservoir sampling from
	a generator seeded per read; dup2 reads at both copies), indels and
	second-best MAPQ is byte-identical at -t 1 and -t 4."""
	one = ch.run_chromap(index, se_reads, "--SAM -n 3 -q 0", threads=1, name="se_n3", shared=True)
	four = ch.run_chromap(index, se_reads, "--SAM -n 3 -q 0", threads=4, name="se_n3", shared=True)
	a, b = one.lines(), four.lines()
	assert a == b, ch.first_line_diff(a, b)
	header, recs = one.sam()
	by_name = collections.defaultdict(list)
	for rec in recs:
		by_name[rec.qname].append(rec)
	n = 0
	for r in se_reads:
		if r.kind not in ("dup2", "rep8"):
			continue
		got = by_name[r.name]
		keys = [(x.rname, x.pos, "-" if x.is_reverse else "+") for x in got]
		loci = {(c, s + 1, st) for c, s, e, st in _segment_loci(genome, r)}
		assert len(keys) == {"dup2": 2, "rep8": 3}[r.kind], r
		assert len(set(keys)) == len(keys) and set(keys) <= loci, r
		assert sum(1 for x in got if not x.is_secondary) == 1, r
		n += 1
	assert n == 8
