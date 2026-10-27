"""Area B, CLI layer, paired-end reads: fragments, pairing limits, paired-end
MAPQ, multimappers, the paired-end candidate filter and mate rescue, adapter
trimming, split alignment, the presets and -t.

Expected fragments come from the simulation truth, expected MAPQ values from
the rule written out in helpers_b.mapq_paired_end (fed with each mate's
error counts and best/second-best multiplicities, which follow from the
design), seed counts from the reads' minimizers by definition. 8 chromap
processes plus 2 --low-mem runs (marked lowmem)."""

import collections

import pytest

import chromap_helpers as ch
import helpers_b as hb

E = hb.DEFAULT_E
FAST = "--cache-size 2000000"


def _flip(strand):
	return "-" if strand == "+" else "+"


def _fragment_loci(genome, segment, pair):
	"""[(contig, start, end, read-1 strand)] of a fragment lying inside one
	copy of a planted segment, at every copy of that segment."""
	seg = genome.segments[segment]
	used = [c for c in seg.copies if c.contig == pair.contig and c.start <= pair.start and pair.end <= c.end][0]
	if used.strand == "+":
		off, rel = pair.start - used.start, pair.strand
	else:
		off, rel = used.end - pair.end, _flip(pair.strand)
	n = pair.end - pair.start
	out = []
	for cc in seg.copies:
		if cc.strand == "+":
			out.append((cc.contig, cc.start + off, cc.start + off + n, rel))
		else:
			out.append((cc.contig, cc.end - off - n, cc.end - off, _flip(rel)))
	return out


def _mate(errors, length, min_errors=None, second=E + 1, best=1, second_n=0):
	return hb.MateStats(errors, length, length, errors if min_errors is None else min_errors, second, best,
		second_n, 0)


def _unique_pair_mapq(p, lengths=None):
	"""MAPQ of a uniquely placed pair whose mates have only their
	substitutions as errors and no second-best placement."""
	l1, l2 = lengths or (len(p.r1.seq), len(p.r2.seq))
	k1, k2 = len(p.r1.mismatches), len(p.r2.mismatches)
	return hb.mapq_paired_end(_mate(k1, l1), _mate(k2, l2), k1 + k2, 2 * E + 1, 1, 0)


def _mates_seeded(genome, p, lengths=None):
	"""Both mates (cut to `lengths`, default their genomic part) keep the
	seeds chromap needs: min(2, #minimizers) hits at their locus."""
	for r, n in zip((p.r1, p.r2), lengths or (None, None)):
		n = n or (r.end - r.start)
		cut = ch.SimRead(r.name, r.seq[:n], r.qual[:n], r.contig, r.start, r.end, r.strand, r.mismatches)
		nmin, hits = hb.seed_support(genome, cut)
		if hits < min(2, max(1, nmin)):
			return False
	return True


def _pair(sim, name, contig, s1, l1, st1, s2, l2, st2, mm1=(), mm2=()):
	r1 = sim.se_read(contig, s1, l1, st1, mismatches=mm1, name=name)
	r2 = sim.se_read(contig, s2, l2, st2, mismatches=mm2, name=name)
	lo, hi = min(s1, s2), max(s1 + l1, s2 + l2)
	return ch.SimPair(name, r1, r2, contig, lo, hi, st1)


@pytest.fixture(scope="module")
def pe_reads(genome):
	"""Read pairs labelled by p.kind: unique (40, inserts 60-900, mates with
	0-2 interior substitutions keeping their seeds), ins1000 and ins1001 (3
	each, the default -l 1000 and one past it), dup2_anchor (4: one mate in a
	dup2 copy near its edge, the other in unique sequence next to it),
	dup2 (2, both mates inside one dup2 copy), rep8 (2, inside one rep8
	copy), dupmm (4: 40 bp mates inside a dupmm copy, one covering one of the
	copies' differences, both seeded at the other copy) and rescue (3: one
	mate unique, the other inside a rep8 copy)."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_b_paired_end", "pe_reads"), name_prefix="u")
	rng = sim.rng
	pairs = []

	def add(p, kind):
		p.kind = kind
		pairs.append(p)

	while sum(1 for p in pairs if p.kind == "unique") < 40:
		i = len(pairs)
		ins = rng.randint(60, 900)
		c, s = sim.unique_locus(ins)
		mm1 = sorted(rng.sample(range(8, 42), i % 3))
		mm2 = sorted(rng.sample(range(8, 42), (i // 3) % 3))
		p = sim.pair(c, s, ins, strand="+-"[i % 2], mismatches1=mm1, mismatches2=mm2)
		if _mates_seeded(genome, p):
			add(p, "unique")
	for i, ins in enumerate((1000, 1001) * 3):
		c, s = sim.unique_locus(ins)
		add(sim.pair(c, s, ins, strand="+-"[i % 2]), "ins%d" % ins)
	# dup2 copy 1 is chr1:20000-21000 (+), copy 2 chr2:10000-11000 (-).
	add(sim.pair("chr1", 20900, 250, strand="+"), "dup2_anchor")
	add(sim.pair("chr1", 20900, 250, strand="-"), "dup2_anchor")
	add(sim.pair("chr2", 9850, 250, strand="+"), "dup2_anchor")
	add(sim.pair("chr2", 9850, 250, strand="-"), "dup2_anchor")
	add(sim.pair("chr1", 20200, 300, strand="+"), "dup2")
	add(sim.pair("chr2", 10400, 250, strand="-"), "dup2")
	add(sim.pair("chr1", 30050, 300, strand="+"), "rep8")
	add(sim.pair("chr3", 6000, 150, strand="-"), "rep8")
	seg = genome.segments["dupmm"]
	while sum(1 for p in pairs if p.kind == "dupmm") < 4:
		j = sum(1 for p in pairs if p.kind == "dupmm")
		copy = seg.copies[j % 2]
		ins = rng.randint(80, 300)
		off = rng.randint(0, seg.length - ins)
		a, b = (off, off + 40), (off + ins - 40, off + ins)
		subs = [[o for o in seg.substitutions if lo <= o < hi] for lo, hi in (a, b)]
		if sorted(len(x) for x in subs) != [0, 1]:
			continue
		o = (subs[0] + subs[1])[0]
		lo = a[0] if subs[0] else b[0]
		if not 8 <= o - lo <= 31:
			continue
		# The forward mate is the left one; read 1 is the left mate for "+".
		name = "mm%d" % j
		if j // 2 % 2 == 0:
			p = _pair(sim, name, copy.contig, copy.start + a[0], 40, "+", copy.start + b[0], 40, "-")
			left, right = p.r1, p.r2
		else:
			p = _pair(sim, name, copy.contig, copy.start + b[0], 40, "-", copy.start + a[0], 40, "+")
			left, right = p.r2, p.r1
		other = seg.copies[1 - j % 2]
		mate = left if subs[0] else right
		view = ch.SimRead(mate.name, mate.seq, mate.qual, other.contig, other.start + lo, other.start + lo + 40,
			mate.strand, mismatches=[other.start + o])
		if hb.seed_support(genome, view)[1] < 2:
			continue
		add(p, "dupmm")
	# rep8 copies at chr1:32000 (+) and chr2:20000 (+); the other mate is in
	# unique sequence before or after the copy.
	add(sim.pair("chr1", 31800, 500, strand="+"), "rescue")
	add(sim.pair("chr1", 32300, 500, strand="-"), "rescue")
	add(sim.pair("chr2", 19700, 450, strand="+"), "rescue")
	return sim.shuffled(ch.ReadSet(pairs))


@pytest.fixture(scope="module")
def pe_run(index, pe_reads):
	return ch.run_chromap(index, pe_reads, "-q 0", name="pe_q0", shared=True)


def _dupmm_mapq_from_design(genome, p):
	seg = genome.segments["dupmm"]
	copy = [c for c in seg.copies if c.contig == p.contig][0]
	subs = {copy.start + o for o in seg.substitutions}
	k1 = sum(1 for g in range(p.r1.start, p.r1.end) if g in subs)
	k2 = sum(1 for g in range(p.r2.start, p.r2.end) if g in subs)
	assert sorted((k1, k2)) == [0, 1]
	# At the other copy the mate covering a difference has one error: the
	# second-best pair has error sum 1. The other mate maps equally well to
	# both copies (two best single-end placements).
	m_clean = _mate(0, 40, best=2)
	m_diff = _mate(0, 40, second=1, second_n=1)
	m1, m2 = (m_diff, m_clean) if k1 else (m_clean, m_diff)
	return hb.mapq_paired_end(m1, m2, 0, 1, 1, 1)


def test_paired_end_fragments_and_mapq(genome, pe_reads, pe_run):
	"""Each pair is one BED line: the fragment from read 1's and read 2's
	outer ends with read 1's strand. Unique pairs (also inserts of exactly
	-l 1000) have MAPQ 60 from the paired-end rule; a pair with one mate in a
	dup2 copy and the other next to it is unique as a pair (MAPQ 60: the
	paired-end filter discards the dup2 mate's other copy); a dupmm pair whose
	placement at the other copy costs one more error gets
	helpers_b.mapq_paired_end = 20; pairs inside dup2 or rep8 are reported once
	at one of the copies with MAPQ 0; a rescue pair (one mate in rep8) at its
	true fragment. Inserts of 1001 exceed -l 1000 and are not reported."""
	lines = collections.Counter((b.chrom, b.start, b.end, b.strand, b.mapq) for b in pe_run.bed())
	exact = collections.Counter()
	for p in pe_reads:
		if p.kind in ("unique", "ins1000", "dup2_anchor"):
			exact[p.fragment() + (_unique_pair_mapq(p),)] += 1
		elif p.kind == "dupmm":
			mapq = _dupmm_mapq_from_design(genome, p)
			assert mapq == 20
			exact[p.fragment() + (mapq,)] += 1
	assert {k: lines[k] for k in exact} == dict(exact)
	assert {m for (_, _, _, _, m) in exact} == {60, 20}
	rest = lines - exact
	for p in pe_reads:
		if p.kind in ("dup2", "rep8"):
			hits = [k for k in rest if k[:4] in _fragment_loci(genome, p.kind, p)]
			assert len(hits) == 1 and hits[0][4] == 0, p
			rest[hits[0]] -= 1
		elif p.kind == "rescue":
			hits = [k for k in rest if k[:4] == p.fragment()]
			assert len(hits) == 1, p
			rest[hits[0]] -= 1
	assert sum(rest.values()) == 0, rest
	assert {b.dups for b in pe_run.bed()} == {1}


def test_paired_end_mapping_statistics(pe_reads, pe_run, chromap_threads):
	"""Paired-end statistics count reads, two per pair: every pair but the
	1001 bp inserts maps; pairs with one best placement (unique, 1000 bp,
	dup2-anchored, dupmm, rescue) are unique, dup2 and rep8 pairs have
	several; with -n 1 each mapped pair gives one mapping per mate. The
	per-thread counters are checked at -t 1 only (known bug
	a-mapping-statistics-clang-threads)."""
	kinds = collections.Counter(p.kind for p in pe_reads)
	multi = kinds["dup2"] + kinds["rep8"]
	mapped = len(pe_reads) - kinds["ins1001"]
	stats = pe_run.stats()
	assert stats["Number of reads"] == 2 * len(pe_reads)
	assert stats["Number of output mappings (passed filters)"] == mapped
	if chromap_threads != 1:
		return
	assert stats["Number of mapped reads"] == 2 * mapped
	assert stats["Number of uniquely mapped reads"] == 2 * (mapped - multi)
	assert stats["Number of reads have multi-mappings"] == 2 * multi
	assert stats["Number of mappings"] == 2 * mapped
	assert stats["Number of uni-mappings"] == 2 * (mapped - multi)


def test_drop_repetitive_pairs_and_cache_options(index, genome, pe_reads, pe_run):
	"""--drop-repetitive-reads 2 drops pairs with more than 2 best placements
	(the rep8 pairs, 8) and keeps everything else, the dup2 pairs with
	exactly 2 included; --cache-size and --cache-update-param within their
	ranges do not change any line. The output is the default run's minus the
	rep8 lines, byte for byte."""
	run = ch.run_chromap(index, pe_reads, "-q 0 --drop-repetitive-reads 2 --cache-size 2000000 "
		"--cache-update-param 1", name="pe_drop2")
	rep8 = set()
	for p in pe_reads:
		if p.kind == "rep8":
			rep8.update(_fragment_loci(genome, "rep8", p))
	base = pe_run.lines()
	kept = [l for l in base if tuple([l.split("\t")[0], int(l.split("\t")[1]), int(l.split("\t")[2]),
		l.split("\t")[5]]) not in rep8]
	assert len(base) - len(kept) == 2
	assert sum(1 for p in pe_reads if p.kind == "dup2") == 2
	assert run.lines() == kept, ch.first_line_diff(run.lines(), kept)


# ------------------------------------------------- adapters and -l

@pytest.fixture(scope="module")
def trim_reads(genome):
	"""Pairs for --trim-adapters -l 400: inserts 30-49 (read-through into the
	Nextera adapters), some with one substitution in read 2's genomic part;
	inserts 30-36 with two substitutions in read 2 (no trimming, and the
	untrimmed mates have more than 8 errors); inserts 25-29 (overlap below
	--min-read-length 30, no trimming, untrimmed errors above 8); inserts
	50-400 and 401-600."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_b_paired_end", "trim_reads"), name_prefix="a")
	rng = sim.rng
	pairs = []

	def untrimmed_fails(p):
		for r in (p.r1, p.r2):
			cand = r.start if r.strand == "+" else r.end - 1
			d = hb.draft_mapping(genome.contigs[r.contig], r.seq, r.strand, cand, E)
			if d[0] > E:
				return True
		return False

	for ins in range(30, 50):
		c, s = sim.unique_locus(ins)
		p = sim.pair(c, s, ins, strand="+-"[ins % 2])
		p.kind = "trim"
		pairs.append(p)
	while sum(1 for p in pairs if p.kind == "trim1") < 6:
		ins = rng.randint(32, 46)
		c, s = sim.unique_locus(ins)
		p = sim.pair(c, s, ins, strand=rng.choice("+-"), mismatches2=(rng.randint(2, ins - 3),))
		if _mates_seeded(genome, p, (ins, ins)):
			p.kind = "trim1"
			pairs.append(p)
	while sum(1 for p in pairs if p.kind == "notrim") < 6:
		ins = rng.randint(25, 36)
		c, s = sim.unique_locus(ins)
		mm2 = () if ins < 30 else (rng.randint(2, ins // 2 - 2), rng.randint(ins // 2 + 2, ins - 3))
		p = sim.pair(c, s, ins, strand=rng.choice("+-"), mismatches2=mm2)
		if untrimmed_fails(p):
			p.kind = "notrim"
			pairs.append(p)
	for ins in (50, 120, 399, 400, 401, 450, 600):
		c, s = sim.unique_locus(ins)
		p = sim.pair(c, s, ins, strand=rng.choice("+-"))
		p.kind = "long"
		pairs.append(p)
	return ch.ReadSet(pairs)


def test_trim_adapters_and_max_insert_size(index, genome, trim_reads):
	"""--trim-adapters finds read 1 inside the reverse complement of read 2
	(seeds of --min-read-length / 2 = 15 bases, at most one mismatch over an
	overlap of at least 30) and cuts both reads to the overlap, so a fragment
	shorter than the reads is reported as [start, start + insert) even though
	the untrimmed reads would carry too many adapter errors to map. Two
	mismatches in the overlap, or an overlap below 30, leave the reads
	untrimmed and unmapped. -l 400 reports inserts up to 400 and drops longer
	ones. MAPQ from helpers_b.mapq_paired_end with the trimmed lengths."""
	run = ch.run_chromap(index, trim_reads, "--trim-adapters -l 400 -q 0 " + FAST, name="pe_trim")
	want = collections.Counter()
	for p in trim_reads:
		ins = p.end - p.start
		if p.kind in ("trim", "trim1") or (p.kind == "long" and ins <= 400):
			want[p.fragment() + (_unique_pair_mapq(p, (min(ins, 50), min(ins, 50))),)] += 1
	got = collections.Counter((b.chrom, b.start, b.end, b.strand, b.mapq) for b in run.bed())
	assert got == want
	assert sum(1 for p in trim_reads if p.kind == "notrim") == 6


def test_mate_rescue_of_repetitive_mate(index, genome, pe_reads):
	"""With -f 3,5 every minimizer of a mate inside rep8 (8 occurrences) is
	too frequent in both rounds, so that mate has no candidate of its own.
	SupplementCandidates then searches its minimizers' occurrences within
	2 * -l of the other mate's candidate and finds the rep8 copy next to it:
	the rescue pairs are reported at their true fragments. Pairs with both
	mates inside rep8 have no anchor and are not reported; pairs without
	frequent minimizers are unaffected."""
	rescue = ch.ReadSet([p for p in pe_reads if p.kind in ("rescue", "rep8", "unique")])
	run = ch.run_chromap(index, rescue, "-q 0 -f 3,5 " + FAST, name="pe_rescue")
	got = collections.Counter((b.chrom, b.start, b.end, b.strand) for b in run.bed())
	want = collections.Counter(p.fragment() for p in rescue if p.kind != "rep8")
	assert got == want
	assert sum(1 for p in rescue if p.kind == "rescue") == 3


# ------------------------------------------------ split alignment, presets

def _chimera(name, genome, parts):
	"""A read joined from (contig, start, length, strand) parts, 5' to 3';
	its interval is the first part's."""
	seq = ""
	for c, s, n, st in parts:
		x = genome.seq(c, s, s + n)
		seq += x if st == "+" else ch.revcomp(x)
	c, s, n, st = parts[0]
	return ch.SimRead(name, seq, "I" * len(seq), c, s, s + n, st, kind="chimera", extra={"parts": parts})


@pytest.fixture(scope="module")
def split_reads(genome):
	"""Pairs for split alignment: proper unique pairs (4), pairs whose mates
	come from unrelated unique loci, often on different contigs (4), pairs
	whose read 2 joins 35 bases of one locus to 15 of another (a ligation
	junction, 4), and pairs inside rep8 (2)."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_b_paired_end", "split_reads"), name_prefix="s")
	rng = sim.rng
	pairs = []
	for i in range(4):
		c, s = sim.unique_locus(300)
		p = sim.pair(c, s, 300, strand="+-"[i % 2])
		p.kind = "proper"
		pairs.append(p)
	for i in range(4):
		name = "inter%d" % i
		c1, s1 = sim.unique_locus(50)
		c2, s2 = sim.unique_locus(50)
		r1 = sim.se_read(c1, s1, 50, rng.choice("+-"), name=name)
		r2 = sim.se_read(c2, s2, 50, rng.choice("+-"), name=name)
		p = ch.SimPair(name, r1, r2, None, None, None, r1.strand)
		p.kind = "inter"
		pairs.append(p)
	for i in range(4):
		name = "chim%d" % i
		c1, s1 = sim.unique_locus(50)
		c2, s2 = sim.unique_locus(40)
		c3, s3 = sim.unique_locus(20)
		r1 = sim.se_read(c1, s1, 50, rng.choice("+-"), name=name)
		r2 = _chimera(name, genome, [(c2, s2, 35, rng.choice("+-")), (c3, s3, 15, rng.choice("+-"))])
		p = ch.SimPair(name, r1, r2, None, None, None, r1.strand)
		p.kind = "chimera"
		pairs.append(p)
	for i in range(2):
		p = sim.pair("chr1", 30050 + 40 * i, 300, strand="+-"[i])
		p.kind = "rep8"
		pairs.append(p)
	return ch.ReadSet(pairs)


def _five_prime(r):
	"""(contig, 1-based 5' position, strand) of a read's 5'-most aligned
	segment."""
	return (r.contig, r.start + 1 if r.strand == "+" else r.end, r.strand)


def _expected_pairs_record(genome, p):
	"""(read_id, chrom1, pos1, chrom2, pos2, strand1, strand2): each read's 5'
	end, ordered so the first has the lower contig rank, or the lower position
	on one contig."""
	a, b = _five_prime(p.r1), _five_prime(p.r2)
	if not (genome.rid(a[0]) < genome.rid(b[0]) or (a[0] == b[0] and a[1] < b[1])):
		a, b = b, a
	return (p.name, a[0], a[1], b[0], b[1], a[2], b[2])


@pytest.fixture(scope="module")
def split_run(index, split_reads):
	return ch.run_chromap(index, split_reads, "--split-alignment --pairs -q 0 " + FAST, name="pe_split", shared=True)


def test_split_alignment_reports_chimeric_pairs(genome, split_reads, split_run):
	"""--split-alignment pairs the best placement of each mate whatever their
	distance, contig or orientation, and aligns a read whose 3' part comes
	from elsewhere by its 5' part. Proper, inter-contig and chimeric pairs are
	each one pairs record at the mates' 5' ends; a pair inside rep8 has 64
	equally good mate combinations and MAPQ 0."""
	header, recs = split_run.pairs()
	got = {}
	for x in recs:
		assert x.read_id not in got
		got[x.read_id] = x
	assert set(got) == {p.name for p in split_reads}
	for p in split_reads:
		x = got[p.name]
		if p.kind == "rep8":
			assert (x.mapq1, x.mapq2) == (0, 0), x
			continue
		assert (x.read_id, x.chrom1, x.pos1, x.chrom2, x.pos2, x.strand1, x.strand2) == _expected_pairs_record(
			genome, p), x
		assert x.mapq1 > 0, x


def test_preset_hic(index, genome, split_reads, split_run):
	"""--preset hic turns on split alignment, pairs output, -e 4, -q 1 and
	--low-mem: every non-repetitive pair of the split-alignment run comes
	back at the same 5' ends in a pairs file, and the rep8 pairs (MAPQ 0)
	fall below -q 1."""
	run = ch.run_chromap(index, split_reads, "--preset hic " + FAST, name="pe_hic")
	header, recs = run.pairs()
	assert header[0] == "## pairs format v1.0.0"
	got = sorted((x.read_id, x.chrom1, x.pos1, x.chrom2, x.pos2, x.strand1, x.strand2) for x in recs)
	want = sorted(_expected_pairs_record(genome, p) for p in split_reads if p.kind != "rep8")
	assert got == want


def _duplicated(sim, pairs, kind, n):
	out = []
	for p in pairs:
		for j in range(n):
			name = "%s_d%d" % (p.name, j)
			q = sim.pair(p.contig, p.start, p.end - p.start, strand=p.strand, name=name)
			q.kind = kind
			out.append(q)
	return out


@pytest.fixture(scope="module")
def atac_reads(genome):
	"""Pairs for --preset atac: short inserts 30-49 (adapter read-through),
	inserts 60-900, an insert of 1500 (above the default -l 1000, within
	the preset's 2000) and 2100 (above it); three fragments sequenced three
	times each."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_b_paired_end", "atac_reads"), name_prefix="t")
	rng = sim.rng
	pairs = []
	for ins in list(range(30, 50, 3)) + [rng.randint(60, 900) for _ in range(12)] + [1500, 2100]:
		c, s = sim.unique_locus(ins)
		p = sim.pair(c, s, ins, strand=rng.choice("+-"))
		p.kind = "single"
		pairs.append(p)
	dups = []
	for i in range(3):
		c, s = sim.unique_locus(200)
		dups.append(sim.pair(c, s, 200, strand="+-"[i % 2]))
	return ch.ReadSet(pairs + _duplicated(sim, dups, "dup", 3))


@pytest.mark.lowmem
def test_preset_atac(index, atac_reads):
	"""--preset atac: -l 2000, adapter trimming, duplicate removal, Tn5
	shift, BED, --low-mem. Each fragment of up to 2000 bp (short inserts
	included) is one line [start + 4, end - 5) with read 1's strand and MAPQ
	60; the 2100 bp insert is not reported; a fragment sequenced three times
	is one line with count 3."""
	run = ch.run_chromap(index, atac_reads, "--preset atac " + FAST, name="pe_atac")
	want = collections.Counter()
	for p in atac_reads:
		if p.end - p.start <= 2000:
			want[(p.contig, p.start + 4, p.end - 5, p.strand)] += 1
	got = {(b.chrom, b.start, b.end, b.strand): (b.mapq, b.dups) for b in run.bed()}
	assert got == {k: (60, n) for k, n in want.items()}
	assert sorted(n for n in want.values()) == [1] * (len(want) - 3) + [3, 3, 3]


@pytest.fixture(scope="module")
def chip_reads(genome):
	"""Pairs for --preset chip -l 1500: inserts 60-900, 1500 and 1800, and
	three fragments sequenced twice each."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_b_paired_end", "chip_reads"), name_prefix="c")
	rng = sim.rng
	pairs = []
	for ins in [rng.randint(60, 900) for _ in range(12)] + [1500, 1800]:
		c, s = sim.unique_locus(ins)
		p = sim.pair(c, s, ins, strand=rng.choice("+-"))
		pairs.append(p)
	dups = []
	for i in range(3):
		c, s = sim.unique_locus(300)
		dups.append(sim.pair(c, s, 300, strand="+-"[i % 2]))
	return ch.ReadSet(pairs + _duplicated(sim, dups, "dup", 2))


@pytest.mark.lowmem
def test_preset_chip_with_override(index, chip_reads):
	"""--preset chip (-l 2000, duplicate removal, BED, --low-mem) with -l 1500
	given as well: the explicit option wins, so the 1500 bp insert is
	reported and the 1800 bp one is not; no Tn5 shift, so fragments are their
	true [start, end); duplicated fragments come back once with count 2."""
	run = ch.run_chromap(index, chip_reads, "--preset chip -l 1500 " + FAST, name="pe_chip")
	want = collections.Counter(p.fragment() for p in chip_reads if p.end - p.start <= 1500)
	got = {(b.chrom, b.start, b.end, b.strand): (b.mapq, b.dups) for b in run.bed()}
	assert got == {k: (60, n) for k, n in want.items()}
	assert sorted(want.values()).count(2) == 3


# ----------------------------------------------------------- threads

def test_paired_end_thread_count_does_not_change_output(index, genome):
	"""Paired-end mapping makes one OpenMP task per 5000-9999 pairs, so 11000
	pairs give two tasks that -t 4 runs on two threads. SAM output with
	--trim-adapters, -n 2 and -q 0 (CIGARs, adapter-trimmed short inserts,
	dup2 and rep8 multimappers drawn by a per-task generator seeded with 11)
	is byte-identical at -t 1 and -t 4."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_b_paired_end", "threads"), name_prefix="h")
	recs = list(sim.pairs(8500)) + list(sim.pairs(500, insert=(30, 49)))
	for i in range(2000):
		seg = genome.segments[("dup2", "rep8")[i % 2]]
		copy = seg.copies[i % len(seg.copies)]
		offset = sim.rng.randint(0, seg.length - 300)
		recs.append(sim.pair(copy.contig, copy.start + offset, 250, strand="+-"[i % 2]))
	reads = sim.shuffled(ch.ReadSet(recs))
	args = "--SAM --trim-adapters -n 2 -q 0 " + FAST
	one = ch.run_chromap(index, reads, args, threads=1, name="pe_t1")
	four = ch.run_chromap(index, reads, args, threads=4, name="pe_t4")
	a, b = one.lines(), four.lines()
	assert len(a) > 2 * 11000
	assert a == b, ch.first_line_diff(a, b)
