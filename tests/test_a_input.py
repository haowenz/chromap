"""Area A: read input files and the index on the mapping command line.

Reads come from the simulator with their origin recorded; each test writes
them to the files it needs (several files, comma lists and glob patterns,
FASTQ and FASTA, plain and gzipped) and checks chromap's output against that
truth. Runs: 12 mapping processes in all (see tests/TARGETS.md for the
budget).

Of chromap's stderr statistics only "Number of reads" is asserted: built with
clang++ and libomp, the mapped-read counts come out low at -t > 1 while the
output is right (known bug a-mapping-statistics-clang-threads)."""

import collections

import pytest

import chromap_helpers as ch
import helpers_a as ha

MISMATCH_MESSAGE = "Numbers of reads and barcodes don't match!"
# Lengths below, at and above the default --min-read-length (30) and 45.
SE_LENGTHS = (26, 27, 28, 29, 30, 31, 44, 45, 46, 50)
PE_SHORT_LENGTHS = (26, 28, 29, 30, 31)


def _bed_set(run):
	return collections.Counter((b.chrom, b.start, b.end, b.strand) for b in run.bed())


def _truth(reads):
	return collections.Counter(r.interval() for r in reads)


# ------------------------------------------------------------ single-end

@pytest.fixture(scope="module")
def se_reads(genome):
	"""80 single-end reads at unique loci on both strands: 40 of the lengths
	in SE_LENGTHS and 40 of 50 bp, 12 of them in lower case."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_a_input", "se_reads"))
	recs = []
	for i, length in enumerate(SE_LENGTHS * 4):
		contig, start = sim.unique_locus(length)
		recs.append(sim.se_read(contig, start, length, "+-"[i % 2]))
	for i in range(40):
		contig, start = sim.unique_locus(50)
		r = sim.se_read(contig, start, 50, "+-"[(i // 2) % 2])
		recs.append(ha.lower_read(r) if i % 3 == 0 and i < 36 else r)
	return sim.shuffled(ch.ReadSet(recs))


@pytest.fixture(scope="module")
def se_files(se_reads, tmp_path_factory):
	"""se_reads in four files, a quarter each: p1.fq (FASTQ, headers with a
	comment), p2.fq.gz (gzipped FASTQ), q_b.fa (FASTA, 17 bases per line)
	and q_a.fa.gz (gzipped FASTA, headers with a tab comment). Returns
	(directory, [the four paths in that order], [the four read lists])."""
	d = tmp_path_factory.mktemp("se_files")
	reads = list(se_reads)
	q = len(reads) // 4
	parts = [reads[:q], reads[q:2 * q], reads[2 * q:3 * q], reads[3 * q:]]
	paths = [
		ha.write_reads(d / "p1.fq", parts[0], header="{name} comment 1:N:0"),
		ha.write_reads(d / "p2.fq.gz", parts[1], gz=True),
		ha.write_reads(d / "q_b.fa", parts[2], fasta=True, line_width=17),
		ha.write_reads(d / "q_a.fa.gz", parts[3], fasta=True, gz=True, header="{name}\ttab comment"),
	]
	return d, paths, parts


@pytest.fixture(scope="module")
def se_multi_run(index, se_files):
	"""One run over the four files given as a comma-separated list whose
	last element is the glob pattern q_* (matching q_a.fa.gz and q_b.fa)."""
	d, paths, _ = se_files
	return ch.run_chromap(index, {"r1": [paths[0], paths[1], d / "q_*"]}, name="se_multi")


def test_multiple_files_fastq_fasta_gzip_map_to_truth(se_multi_run, se_reads):
	"""-1 takes a comma-separated list whose last element is the glob
	pattern q_*: reads from plain and gzipped FASTQ and FASTA files
	(multi-line FASTA, headers with comments, reads in lower case) of every
	listed and matched file are all loaded ("Number of reads" counts every
	read of every file) and those of at least 30 bases (the default
	--min-read-length) map to their true interval and strand with MAPQ 60
	and count 1; shorter ones are dropped before mapping."""
	kept = [r for r in se_reads if len(r.seq) >= 30]
	assert _bed_set(se_multi_run) == _truth(kept)
	assert {(b.mapq, b.dups) for b in se_multi_run.bed()} == {(60, 1)}
	assert se_multi_run.stats()["Number of reads"] == len(se_reads)
	assert any(r.seq.islower() for r in kept)


@pytest.mark.parametrize("min_length", [20, 45])
def test_min_read_length(index, se_reads, min_length):
	"""--min-read-length L drops reads shorter than L (45: the 26 to 44 bp
	reads; 20: none, so the 26 to 29 bp reads dropped by default map at their
	origin); every read is still counted in "Number of reads"."""
	run = ch.run_chromap(index, se_reads, ["--min-read-length", str(min_length)], name="se_min_%d" % min_length)
	kept = [r for r in se_reads if len(r.seq) >= min_length]
	assert _bed_set(run) == _truth(kept)
	assert run.stats()["Number of reads"] == len(se_reads)
	if min_length == 45:
		assert len(kept) < len(se_reads)
	else:
		assert len(kept) == len(se_reads) and min(len(r.seq) for r in kept) < 30


def test_read_ids_continue_across_files_and_names(index, genome, tmp_path):
	"""Read ids count on from one input file to the next, in the order the
	files are given, a glob pattern contributing its matches in sorted
	order. chromap sorts SAM records with equal position, strand and MAPQ by
	read id, so of three identical reads at one locus, one in each file, the
	one in the first file is written first; were ids to restart at 0 in each
	file (each tie read is placed first there), or the pattern's matches
	read in another order, the order would differ. -1 is "x.fq,y_*" and
	y_* matches y_a.fq.gz and y_b.fq (written in the other order). Names are
	cut at the first space or tab; a /1 suffix is kept."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_a_input", "read_ids"))
	first = [sim.se_read(*sim.unique_locus(50), strand="+-"[i % 2], name="x%02d/1" % i if i % 4 == 0 else None)
		for i in range(20)]
	second = [sim.se_read(*sim.unique_locus(50), strand="+-"[i % 2]) for i in range(20)]
	third = [sim.se_read(*sim.unique_locus(50), strand="+-"[i % 2]) for i in range(10)]
	ties = []
	for i in range(3):
		contig, start = sim.unique_locus(50)
		a = sim.se_read(contig, start, 50, "+-"[i % 2], name="tie%d_first" % i)
		b = ch.SimRead("tie%d_second" % i, a.seq, a.qual, a.contig, a.start, a.end, a.strand)
		c = ch.SimRead("tie%d_third" % i, a.seq, a.qual, a.contig, a.start, a.end, a.strand)
		ties.append((a, b, c))
	f1 = ha.write_reads(tmp_path / "x.fq", first + [t[0] for t in ties], header="{name} first file")
	ha.write_reads(tmp_path / "y_b.fq", [t[2] for t in ties] + third)
	ha.write_reads(tmp_path / "y_a.fq.gz", [t[1] for t in ties] + second, gz=True, header="{name}\tsecond")
	run = ch.run_chromap(index, {"r1": [f1, tmp_path / "y_*"]}, "--SAM", name="se_ids")
	header, recs = run.sam()
	truth = {r.name: r for r in first + second + third + [x for t in ties for x in t]}
	assert sorted(rec.qname for rec in recs) == sorted(truth)
	assert sum(1 for rec in recs if rec.qname.endswith("/1")) == 5
	for rec in recs:
		r = truth[rec.qname]
		assert (rec.rname, rec.pos, rec.flag) == (r.contig, r.start + 1, 0 if r.strand == "+" else 16)
	names = [rec.qname for rec in recs]
	for i in range(3):
		at = names.index("tie%d_first" % i)
		assert names[at:at + 3] == ["tie%d_%s" % (i, x) for x in ("first", "second", "third")]


def test_mapping_takes_k_and_w_from_index(index_factory, se_reads):
	"""Mapping computes read minimizers with the k and w stored in the index
	file (-x), not the defaults 17 and 7. With an index built with -k 21
	-w 11 and -s 1000, a read maps only when every one of its minimizers
	hits its locus (chromap requires min(-s, number of minimizers) seeds).
	Every minimizer of an exact read with k 21 and w 11 is a reference
	minimizer of its locus, so all 50 bp reads map to their origin; with k
	17 no read minimizer would be in the index, and with w 7 a read has
	minimizers that the w 11 index does not hold."""
	idx = index_factory(k=21, w=11)
	reads = ch.ReadSet([r for r in se_reads if len(r.seq) == 50])
	run = ch.run_chromap(idx, reads, "-s 1000", name="se_k21")
	assert _bed_set(run) == _truth(reads)


# ------------------------------------------------------------ paired-end

@pytest.fixture(scope="module")
def pe_reads(genome):
	"""160 pairs at unique loci (50 bp reads, inserts 150-400) and 20 more
	with one mate cut at its 3' end to a length in PE_SHORT_LENGTHS (read 1
	or read 2, both strands); a cut mate keeps its 5' end, so the fragment is
	unchanged."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_a_input", "pe_reads"))
	recs = list(sim.pairs(160))
	for i, length in enumerate(PE_SHORT_LENGTHS * 4):
		contig, start = sim.unique_locus(300)
		p = sim.pair(contig, start, 300, strand="+-"[(i // 2) % 2])
		recs.append(ha.shorten_mate(p, 1 + i % 2, length))
	return sim.shuffled(ch.ReadSet(recs))


@pytest.fixture(scope="module")
def pe_files(pe_reads, tmp_path_factory):
	"""pe_reads in two file pairs: a_1.fq/a_2.fq (plain FASTQ) and
	b_1.fq.gz/b_2.fq.gz (gzipped)."""
	d = tmp_path_factory.mktemp("pe_files")
	pairs = list(pe_reads)
	half = len(pairs) // 2
	out = {"r1": [], "r2": []}
	for prefix, chunk, gz in (("a", pairs[:half], False), ("b", pairs[half:], True)):
		ext = ".fq.gz" if gz else ".fq"
		out["r1"].append(ha.write_reads(d / (prefix + "_1" + ext), [p.r1 for p in chunk], gz=gz))
		out["r2"].append(ha.write_reads(d / (prefix + "_2" + ext), [p.r2 for p in chunk], gz=gz))
	return out


@pytest.fixture(scope="module")
def pe_runs(index, pe_files):
	"""The same paired-end run at -t 1 (reads parsed by one thread) and -t 4
	(read 1 and read 2 files parsed in parallel tasks)."""
	return {t: ch.run_chromap(index, pe_files, "--cache-size 2000000", threads=t, name="pe_multi_t%d" % t)
		for t in (1, 4)}


def test_paired_end_files_map_to_truth(pe_runs, pe_reads):
	"""Paired-end reads from a list of two file pairs (plain and gzipped)
	map as one run: each pair whose mates are both at least 30 bases is
	reported once as its fragment with read 1's strand; a pair with either
	mate shorter than --min-read-length (30) is dropped. "Number of reads"
	counts both mates of every pair."""
	kept = [p for p in pe_reads if min(len(p.r1.seq), len(p.r2.seq)) >= 30]
	assert len(kept) == len(pe_reads) - 12
	for run in pe_runs.values():
		assert _bed_set(run) == collections.Counter(p.fragment() for p in kept)
		assert run.stats()["Number of reads"] == 2 * len(pe_reads)


def test_parallel_parsing_gives_same_output(pe_runs):
	"""At -t 3 or more chromap parses the read 1 and read 2 files in parallel
	tasks; the output is byte-identical to the -t 1 run."""
	one, four = pe_runs[1].lines(), pe_runs[4].lines()
	assert one == four, ch.first_line_diff(one, four)


def test_paired_end_min_read_length_either_mate(index, pe_files, pe_reads):
	"""With --min-read-length 20 the pairs with a 26 to 29 bp mate are kept
	and map at their fragments, so the default drops them for their length:
	one short mate is enough."""
	run = ch.run_chromap(index, pe_files, "--cache-size 2000000 --min-read-length 20", name="pe_min20")
	assert _bed_set(run) == collections.Counter(p.fragment() for p in pe_reads)


# ------------------------------------------------------------ errors

@pytest.mark.parametrize("threads,longer", [(1, "r1"), (4, "r2")])
def test_unequal_read_counts_exit(index, genome, tmp_path, threads, longer):
	"""Read 1 and read 2 files with different numbers of reads end the run
	with "Numbers of reads and barcodes don't match!" and a non-zero exit,
	both when one thread reads them in turn (-t 1) and when they are parsed
	in parallel (-t 4)."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_a_input", "unequal", threads))
	pairs = list(sim.pairs(30))
	n1 = 30 if longer == "r1" else 29
	n2 = 29 if longer == "r1" else 30
	files = {"r1": ha.write_reads(tmp_path / "r1.fq", [p.r1 for p in pairs[:n1]]),
		"r2": ha.write_reads(tmp_path / "r2.fq", [p.r2 for p in pairs[:n2]])}
	run = ch.run_chromap(index, files, "--cache-size 2000000", threads=threads, check=False, name="unequal")
	assert run.returncode != 0
	assert MISMATCH_MESSAGE in run.stderr


@pytest.mark.parametrize("threads", [1, 4])
def test_fewer_barcodes_than_reads_exit(index, genome, tmp_path, threads):
	"""A barcode file (-b) with fewer records than the read file ends a
	single-end run with "Numbers of reads and barcodes don't match!" and a
	non-zero exit, read in turn (-t 1) or in parallel tasks (-t 4)."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_a_input", "barcodes", threads))
	barcodes = ch.make_whitelist(4, seed=3)
	reads = list(sim.se_reads(30, barcodes=barcodes))
	files = {"r1": ha.write_reads(tmp_path / "r.fq", reads),
		"bc": ha.write_fastq(tmp_path / "bc.fq", [(r.name, r.barcode, "I" * len(r.barcode)) for r in reads[:-1]])}
	run = ch.run_chromap(index, files, threads=threads, check=False, name="fewer_barcodes")
	assert run.returncode != 0
	assert MISMATCH_MESSAGE in run.stderr


@pytest.mark.parametrize("flag", ["-1", "-2", "-b"])
def test_pattern_matching_nothing_exits(index, se_files, tmp_path, flag):
	"""A read or barcode file pattern that matches no file (a glob pattern or
	a plain path that does not exist) exits non-zero with "glob() failed"
	(the value glob() returns differs between platforms)."""
	d, paths, _ = se_files
	for missing in (str(d / "nothing_*.fq"), str(d / "absent.fq")):
		files = {"-1": str(paths[0]), "-2": str(paths[1]), "-b": str(paths[1])}
		files[flag] = missing
		argv = ["-x", str(index.path), "-r", str(index.fasta), "-o", str(tmp_path / "out.bed"), "-1", files["-1"]]
		if flag != "-1":
			argv += [flag, files[flag]]
		run = ch.run_cli(argv, cwd=tmp_path)
		assert run.returncode != 0
		assert "glob() failed" in run.stderr
		assert not (tmp_path / "out.bed").exists()
