"""Area C, CLI layer: the output formats chromap writes and the options that
shape them: paired-end --TagAlign (two lines per fragment), paired-end and
barcoded single-end --SAM (flags, POS, RNEXT, PNEXT, TLEN, CIGAR, NM, MD,
CB, secondary records), --PAF for single-end reads, --summary for
paired-end data with --turn-off-num-uniq-cache-slots, --k-for-minhash and
--frip-est-params, and --low-mem with SAM and PAF output, including
-o /dev/stdout.

Expected fields come from the simulation truth and the format definitions
(README.md; the SAM specification for SAM; the PAF definition of minimap2
for PAF)."""

import pytest

import chromap_helpers as ch
import helpers_c as hc


def _f(line):
	return line.split("\t")


# ------------------------------------------------------- paired TagAlign


@pytest.fixture(scope="module")
def pe_tag_design(genome):
	return hc.paired_dedup_design(genome, ch.stable_seed("test_c_cli_formats", "pe_tag"))


@pytest.fixture(scope="module")
def pe_tag_run(index, pe_tag_design):
	reads = pe_tag_design[0]
	return ch.run_chromap(index, reads,
		"--TagAlign --remove-pcr-duplicates --summary summary.csv --turn-off-num-uniq-cache-slots", shared=True)


def test_paired_tagalign_lines(pe_tag_run, pe_tag_design):
	"""Paired-end --TagAlign writes two lines per fragment (README.md): read
	1's interval and strand, then read 2's, MAPQ 60, the duplicate count only
	on the second line. With --remove-pcr-duplicates pairs from one fragment
	are written once with the group size as the count, whichever strand read
	1 is on. Fragments come in reference and start order; repeat (MAPQ below
	30) and random pairs are not written."""
	reads, groups, reps, unm = pe_tag_design
	lines = pe_tag_run.lines()
	assert len(lines) % 2 == 0
	records = [tuple(lines[i:i + 2]) for i in range(0, len(lines), 2)]
	assert all(len(_f(a)) == 6 and len(_f(b)) == 7 for a, b in records)

	def key(rec):
		a, b = _f(rec[0]), _f(rec[1])
		return (a[0], min(int(a[1]), int(b[1])), max(int(a[2]), int(b[2])))

	seen = hc.match_groups(records, groups, key,
		lambda g: set(hc.pair_tagalign_lines(p, g.size) for p in g.reads))
	assert len(seen) == len(groups)
	order = [(pe_tag_run.index.genome.rid(k[0]), k[1], k[2] - k[1]) for k in (key(r) for r in records)]
	assert order == sorted(order)


def test_paired_summary_without_cache_slot_column(pe_tag_run, pe_tag_design):
	"""--summary for bulk paired-end data counts pairs: total = every pair,
	duplicate = pairs beyond the first of each written fragment, unmapped =
	random pairs, lowmapq = repeat pairs; no cache hits in a first batch, so
	cachehit, fric and estfrip are 0. --turn-off-num-uniq-cache-slots drops
	the numcacheslots column."""
	reads, groups, reps, unm = pe_tag_design
	s = pe_tag_run.summary()
	assert s.header == hc.SUMMARY_HEADER
	assert s.rows == {"": dict(total=len(reads), duplicate=sum(g.size - 1 for g in groups), unmapped=len(unm),
		lowmapq=len(reps), cachehit=0, fric=0.0, estfrip=0.0)}


def test_paired_frip_params_must_be_five_numbers(index, genome):
	"""--frip-est-params takes five numbers separated by ';' (README.md,
	*Estimating FRiP*); a paired-end run given four stops before mapping with
	a non-zero exit and says so."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_c_cli_formats", "frip"))
	run = ch.run_chromap(index, sim.pairs(5), "--frip-est-params 1;2;3;4 --summary summary.csv", check=False)
	assert run.returncode != 0
	assert "expecting 5 parameters but found 4" in run.stderr
	assert not run.exists("summary.csv")


# ------------------------------------------------------- paired SAM


@pytest.fixture(scope="module")
def pe_sam_design(genome):
	return hc.paired_dedup_design(genome, ch.stable_seed("test_c_cli_formats", "pe_sam"), identical_sizes=(2, 3, 4),
		mixed_groups=0, uniques=40)


@pytest.fixture(scope="module")
def pe_sam_run(index, pe_sam_design):
	return ch.run_chromap(index, pe_sam_design[0],
		"--SAM --remove-pcr-duplicates --summary summary.csv --k-for-minhash 1", shared=True)


def _expected_pair_record(pair, read1):
	"""SAM fields of one mate of a properly paired, uniquely mapped pair
	(SAM specification): FLAG = paired (1) + proper pair (2) + reverse (16)
	if this mate is on - + mate reverse (32) if the other is on - + read 1
	(64) or read 2 (128); POS and PNEXT 1-based leftmost positions; RNEXT
	'='; TLEN the insert, positive for the leftmost (+) mate; SEQ and QUAL in
	reference orientation."""
	mate, other = (pair.r1, pair.r2) if read1 else (pair.r2, pair.r1)
	flag = 1 | 2 | (16 if mate.strand == "-" else 0) | (32 if other.strand == "-" else 0) | (64 if read1 else 128)
	tlen = pair.insert if mate.strand == "+" else -pair.insert
	seq = mate.seq if mate.strand == "+" else ch.revcomp(mate.seq)
	qual = mate.qual if mate.strand == "+" else mate.qual[::-1]
	return [str(flag), pair.contig, str(mate.start + 1), "60", "50M", "=", str(other.start + 1), str(tlen), seq, qual,
		"NM:i:0", "MD:Z:50"]


def test_paired_sam_records(pe_sam_run, pe_sam_design, genome):
	"""Paired-end --SAM writes both mates of every mapped pair with the
	fields of _expected_pair_record; with --remove-pcr-duplicates a fragment
	read several times keeps one read-1 and one read-2 record. Repeat (MAPQ
	below 30) and random pairs are not written; bulk data has no CB tag; the
	header lists every reference."""
	reads, groups, reps, unm = pe_sam_design
	header, recs = pe_sam_run.sam()
	assert header == ["@SQ\tSN:%s\tLN:%d" % (c, len(s)) for c, s in genome.contigs.items()]
	by_fragment = {}
	for line in pe_sam_run.lines()[len(header):]:
		f = _f(line)
		assert len(f) == 13, line
		left = min(int(f[3]), int(f[7])) - 1
		by_fragment.setdefault((f[2], left, left + abs(int(f[8]))), []).append(f)
	assert set(by_fragment) == set(g.key for g in groups)
	for g in groups:
		got = by_fragment[g.key]
		assert sorted(int(f[1]) & 192 for f in got) == [64, 128], got
		for f in got:
			assert f[0] in g.names
			pair = pe_sam_run.truth[f[0]]
			assert f[1:] == _expected_pair_record(pair, bool(int(f[1]) & 64)), f


def test_paired_sam_summary_counts_pairs(pe_sam_run, pe_sam_design):
	"""With SAM output each pair gives two records, and the summary still
	counts pairs: duplicate, unmapped and lowmapq are per pair (chromap
	halves its per-record counts). The numcacheslots column is written (0:
	no cache hits); --k-for-minhash 1 is accepted and echoed."""
	reads, groups, reps, unm = pe_sam_design
	s = pe_sam_run.summary()
	assert s.header == hc.SUMMARY_HEADER + ["numcacheslots"]
	assert s.rows == {"": dict(total=len(reads), duplicate=sum(g.size - 1 for g in groups), unmapped=len(unm),
		lowmapq=len(reps), cachehit=0, fric=0.0, estfrip=0.0, numcacheslots=0)}
	assert "K for MinHash: 1" in pe_sam_run.stderr


# ------------------------------------------------------- single-end SAM


def _indel_read(sim, genome, kind, p, strand, barcode, name):
	"""A 50 bp read with one deletion (kind "D": the reference base at
	offset p, counted along the + strand, is missing) or one insertion ("I":
	an extra base after offset p) at a place where the gap cannot be shifted
	(the neighbouring bases differ), so its alignment is unique:
	CIGAR pM1D(50-p)M with MD p^X(50-p), or pM1I(49-p)M with MD 49."""
	while True:
		contig, start = sim.unique_locus(60)
		g = genome.seq(contig, start, start + 60)
		if kind == "D" and g[p - 1] != g[p] and g[p] != g[p + 1]:
			fwd = g[:p] + g[p + 1:51]
			span, cigar, md = 51, "%dM1D%dM" % (p, 50 - p), "%d^%s%d" % (p, g[p], 50 - p)
			break
		if kind == "I":
			ins = [b for b in "ACGT" if b != g[p - 1] and b != g[p]][0]
			fwd = g[:p] + ins + g[p:49]
			span, cigar, md = 49, "%dM1I%dM" % (p, 49 - p), "49"
			break
	seq = fwd if strand == "+" else ch.revcomp(fwd)
	qual = "".join(chr(35 + (i * 3 + p) % 40) for i in range(50))
	r = ch.SimRead(name, seq, qual, contig, start, start + span, strand, (), kind="indel", barcode=barcode)
	r.extra.update(cigar=cigar, md=md, nm=1)
	return r


@pytest.fixture(scope="module")
def se_sam_design(genome, whitelist):
	"""Barcoded single-end reads: 18 with 0, 1 or 2 substitutions on both
	strands, six with one insertion or deletion, four from the two-copy
	repeat dup2 (two best locations each) and two random reads."""
	barcodes = whitelist[0]
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_c_cli_formats", "se_sam"))
	reads = []
	for i in range(18):
		contig, start = sim.unique_locus(50)
		r = sim.se_read(contig, start, 50, "+-"[i % 2], mismatches=[(), (25,), (6, 43)][i % 3],
			qual="".join(chr(35 + (j * 7 + i) % 38) for j in range(50)), barcode=barcodes[i % 5])
		ref = genome.seq(r.contig, r.start, r.end)
		md, run_len = "", 0
		for j in range(50):
			if r.start + j in r.mismatches:
				md += "%d%s" % (run_len, ref[j])
				run_len = 0
			else:
				run_len += 1
		r.extra.update(cigar="50M", md=md + str(run_len), nm=len(r.mismatches))
		reads.append(r)
	for i, (kind, p, strand) in enumerate([("D", 20, "+"), ("D", 25, "-"), ("I", 20, "+"), ("I", 30, "-"),
			("D", 12, "+"), ("I", 38, "-")]):
		reads.append(_indel_read(sim, genome, kind, p, strand, barcodes[5 + i], "indel%d" % i))
	for i, offset in enumerate((100, 300, 500, 700)):
		r = sim.segment_read("dup2", offset, 50, "+-"[i % 2])
		r.barcode = barcodes[11 + i]
		reads.append(r)
	for i in range(2):
		r = sim.random_read(50)
		r.barcode = barcodes[20 + i]
		reads.append(r)
	return ch.ReadSet(reads)


SE_SAM_ARGS = "--SAM --Tn5-shift -n 2 -q 0"


@pytest.fixture(scope="module")
def se_sam_run(index, se_sam_design):
	return ch.run_chromap(index, se_sam_design, SE_SAM_ARGS, shared=True)


def test_single_end_sam_fields(se_sam_run, se_sam_design, genome):
	"""Single-end --SAM with barcodes: RNEXT '*', PNEXT 0 and TLEN 0; CB:Z:
	the read's barcode; POS the 1-based leftmost base; CIGAR, NM and MD from
	the read's true substitutions and indels; SEQ and QUAL in reference
	orientation. --Tn5-shift does not apply to SAM output (chromap.1), so
	positions are unshifted. -n 2 writes both copies of a dup2 read, one of
	them flagged secondary (256), both with MAPQ 0; random reads are not
	written."""
	header, recs = se_sam_run.sam()
	lines = se_sam_run.lines()[len(header):]
	by_name = {}
	for line in lines:
		by_name.setdefault(_f(line)[0], []).append(_f(line))
	unique = [r for r in se_sam_design if r.kind in ("genomic", "indel")]
	repeats = [r for r in se_sam_design if r.kind == "segment"]
	assert set(by_name) == set(r.name for r in unique + repeats)
	for r in unique:
		(f,) = by_name[r.name]
		seq = r.seq if r.strand == "+" else ch.revcomp(r.seq)
		qual = r.qual if r.strand == "+" else r.qual[::-1]
		assert f[1:] == [str(0 if r.strand == "+" else 16), r.contig, str(r.start + 1), "60", r.extra["cigar"], "*",
			"0", "0", seq, qual, "NM:i:%d" % r.extra["nm"], "MD:Z:%s" % r.extra["md"], "CB:Z:%s" % r.barcode], r
	for r in repeats:
		got = by_name[r.name]
		assert len(got) == 2
		assert sorted(int(f[1]) & 256 for f in got) == [0, 256]
		loci = set((c, s + 1) for c, s, e in r.extra["loci"])
		assert set((f[2], int(f[3])) for f in got) == loci
		for f in got:
			assert (f[4], f[5], f[6:9], f[-1]) == ("0", "50M", ["*", "0", "0"], "CB:Z:" + r.barcode)
			copy_forward = (f[2], int(f[3]) - 1) == (r.contig, r.start)
			reverse = (r.strand == "-") == copy_forward
			assert int(f[1]) & 16 == (16 if reverse else 0)
			assert f[9] == (ch.revcomp(r.seq) if reverse else r.seq)
	order = [(genome.rid(f[2]), int(f[3])) for f in (_f(l) for l in lines)]
	assert order == sorted(order)


def test_low_mem_sam_to_stdout(index, se_sam_design, se_sam_run):
	"""--low-mem with SAM output and -o /dev/stdout: the records go to
	standard output, identical to the default mode's file (no deduplication
	here, so the modes have nothing to choose), and the temp file, written to
	the working directory when the output is /dev/stdout, is removed."""
	run = ch.run_chromap(index, se_sam_design, SE_SAM_ARGS + " --low-mem", output="/dev/stdout")
	assert run.stdout.splitlines() == se_sam_run.lines()
	assert sorted(p.name for p in run.outdir.iterdir()) == ["input"]


# ------------------------------------------------------- PAF


def _paf_fields(r, genome, mapq=60):
	"""PAF columns for an exactly matching single-end read: query name,
	length, start 0, end, strand, target name, length, start, end, residue
	matches (all bases), alignment block length, MAPQ."""
	n = r.end - r.start
	return [r.name, str(n), "0", str(n), r.strand, r.contig, str(len(genome.contigs[r.contig])), str(r.start),
		str(r.end), str(n), str(n), str(mapq)]


def test_single_end_paf(index, genome):
	"""--PAF with exact single-end reads of 40, 50 and 75 bp on both strands:
	one 12-column line per read, sorted by reference and start."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_c_cli_formats", "paf"))
	reads = ch.ReadSet([sim.se_read(*sim.unique_locus(75), length=(40, 50, 75)[i % 3], strand="+-"[i % 2])
		for i in range(30)])
	run = ch.run_chromap(index, reads, "--PAF")
	lines = [_f(l) for l in run.lines()]
	assert sorted(lines) == sorted(_paf_fields(r, genome) for r in reads)
	order = [(genome.rid(f[5]), int(f[7])) for f in lines]
	assert order == sorted(order)


def test_low_mem_paf_dedup(index, genome):
	"""--PAF --low-mem --remove-pcr-duplicates: one line per duplicate group
	of identical reads (the query name is one of the group's), the other
	columns those of the group's interval; the temp file is removed."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_c_cli_formats", "paf_lowmem"))
	groups = [hc.Group("singleton", [r]) for r in sim.se_reads(20)]
	for n in (2, 3, 4, 5):
		contig, start = sim.unique_locus(50)
		strand = "+-"[n % 2]
		groups.append(hc.Group("identical", [sim.se_read(contig, start, 50, strand) for _ in range(n)]))
	reads = sim.shuffled(ch.ReadSet([r for g in groups for r in g.reads]))
	run = ch.run_chromap(index, reads, "--PAF --low-mem --remove-pcr-duplicates")
	lines = [tuple(_f(l)) for l in run.lines()]
	seen = hc.match_groups(lines, groups, lambda f: (f[5], int(f[7])),
		lambda g: set(tuple(_paf_fields(r, genome)) for r in g.reads))
	assert len(seen) == len(groups)
	assert sorted(p.name for p in run.outdir.iterdir()) == ["input", "out.txt"]
