"""Pilot of the CLI layer: the pattern for tests that run the chromap binary.

Reads are simulated from the synthetic genome with their origin recorded,
and every expected value below is derived from that truth (or worked by hand
in the docstring), never copied from chromap's output."""

import collections

import pytest

import chromap_helpers as ch


@pytest.fixture(scope="module")
def pe_reads(genome):
	"""120 read pairs (50 bp, inserts 150-400) at unique loci, read 1 on
	either strand. Module-scoped: seeded from its own name, never from a
	simulator shared with tests."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_pilot_cli", "pe_reads"))
	return sim.pairs(120)


@pytest.fixture(scope="module")
def pe_bed(index, pe_reads):
	"""One default BED run on pe_reads, shared by the tests that only read it
	(shared=True also lets other modules reuse the identical run)."""
	return ch.run_chromap(index, pe_reads, name="pe_bed", shared=True)


def test_single_end_reads_map_to_their_origin(index, sim):
	"""Exact 50 bp reads from unique loci on both strands come back as one BED
	line each at their true interval [start, start + 50) and strand. MAPQ is
	60: with no second-best hit the error difference is capped at -e 8, so
	MAPQ = 5 * 6.02 * 8 * (int(ln 50) / ln 50)^2 = 141.6, clipped to 60
	(src/mapping_generator.h, GetMAPQForSingleEndRead). Without
	--remove-pcr-duplicates every line's count is 1. Lines are sorted by
	reference (FASTA order) and start."""
	reads = sim.se_reads(80)
	run = ch.run_chromap(index, reads)
	recs = run.bed()
	expected = collections.Counter((r.contig, r.start, r.end, r.strand) for r in reads)
	assert collections.Counter((b.chrom, b.start, b.end, b.strand) for b in recs) == expected
	assert {(b.name, b.mapq, b.dups) for b in recs} == {("N", 60, 1)}
	order = [(index.genome.rid(b.chrom), b.start) for b in recs]
	assert order == sorted(order)


def test_paired_end_fragments(index, pe_reads, pe_bed):
	"""Each pair is reported once as its fragment [start, end) with read 1's
	strand, MAPQ 60 and count 1."""
	expected = collections.Counter(p.fragment() for p in pe_reads)
	recs = pe_bed.bed()
	assert collections.Counter((b.chrom, b.start, b.end, b.strand) for b in recs) == expected
	assert {(b.mapq, b.dups) for b in recs} == {(60, 1)}


def test_tn5_shift_on_fragments(index, pe_reads, pe_bed):
	"""--Tn5-shift moves a fragment's start right by 4 and its end left by 5
	(Tn5 inserts its adapters 9 bp apart), whatever read 1's strand."""
	run = ch.run_chromap(index, pe_reads, "--Tn5-shift")
	shifted = collections.Counter((b.chrom, b.start, b.end) for b in run.bed())
	expected = collections.Counter((c, s + 4, e - 5) for c, s, e, _ in (p.fragment() for p in pe_reads))
	assert shifted == expected
	assert len(run.bed()) == len(pe_bed.bed())


def test_sam_fields_from_truth(index, sim, genome):
	"""Single-end SAM: POS is the 1-based start, CIGAR 50M, FLAG 0 or 16 by
	strand, SEQ the read in reference orientation and QUAL reversed with it
	for reverse-strand reads, NM the number of substitutions and MD the
	reference bases at the substituted positions with the match runs between
	them (SAM specification). The @SQ lines list every contig with its
	length in FASTA order. Substitutions sit where enough error-free 17-mers
	remain to seed the read (one at offset 25 leaves 17 of 34 17-mers intact;
	offsets 6 and 43 leave 20)."""
	reads = []
	for i in range(24):
		contig, start = sim.unique_locus(50)
		r = sim.se_read(contig, start, 50, "+-"[i % 2], mismatches=[(), (25,), (6, 43)][i % 3],
			qual="".join(chr(35 + (j * 7 + i) % 38) for j in range(50)))
		reads.append(r)
	run = ch.run_chromap(index, ch.ReadSet(reads), "--SAM")
	header, recs = run.sam()
	assert [h for h in header if h.startswith("@SQ")] == [
		"@SQ\tSN:%s\tLN:%d" % (c, len(s)) for c, s in genome.contigs.items()]
	assert sorted(r.qname for r in recs) == sorted(r.name for r in reads)
	for rec in recs:
		r = run.truth[rec.qname]
		ref = genome.seq(r.contig, r.start, r.end)
		md, run_len = "", 0
		for i in range(50):
			if r.start + i in r.mismatches:
				md += "%d%s" % (run_len, ref[i])
				run_len = 0
			else:
				run_len += 1
		md += str(run_len)
		assert (rec.rname, rec.pos, rec.cigar) == (r.contig, r.start + 1, "50M"), r
		assert rec.flag == (0 if r.strand == "+" else 16), r
		assert rec.seq == (r.seq if r.strand == "+" else ch.revcomp(r.seq)), r
		assert rec.qual == (r.qual if r.strand == "+" else r.qual[::-1]), r
		assert (rec.tags["NM"], rec.tags["MD"]) == (len(r.mismatches), md), r


def test_thread_count_does_not_change_output(index, genome):
	"""The BED output is byte-identical at -t 1 and -t 4. Paired-end mapping
	runs one OpenMP task per 5000 to 9999 pairs (taskloop grainsize 5000 in
	MapPairedEndReads), so 12000 pairs make two tasks that four threads can
	run at once. A quarter of the pairs come from the two-copy repeat dup2 and
	the eight-copy rep8 and are written with -q 0, so the choice among
	equally good locations is part of the comparison."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_pilot_cli", "threads"))
	recs = list(sim.pairs(9000))
	for i in range(3000):
		seg = genome.segments[("dup2", "rep8")[i % 2]]
		copy = seg.copies[0]
		offset = sim.rng.randint(0, seg.length - 300)
		recs.append(sim.pair(copy.contig, copy.start + offset, 250, strand="+-"[i % 2]))
	reads = sim.shuffled(ch.ReadSet(recs))
	one = ch.run_chromap(index, reads, "-q 0", threads=1)
	four = ch.run_chromap(index, reads, "-q 0", threads=4)
	assert len(one.lines()) == 12000
	assert one.lines() == four.lines(), ch.first_line_diff(one.lines(), four.lines())
