"""Checks of the harness itself (synthetic genome, simulator truth, the
Python reference functions, parsers and the runner), so a failure elsewhere
points at chromap rather than at chromap_helpers."""

import gzip
import shutil
import subprocess
import sys

import pytest

import chromap_helpers as ch


def test_synthetic_genome_layout(genome):
	"""The default genome is deterministic and 100 kb on four contigs; every
	copy of dup2 and rep8 carries the same bases in its own orientation,
	dupmm's two copies differ at exactly its 12 listed offsets, and the N run
	and low-complexity stretches are where the specs put them."""
	again = ch.make_synthetic_genome()
	assert again.contigs == genome.contigs
	assert [(c, len(s)) for c, s in genome.contigs.items()] == list(ch.CONTIG_SPECS)
	assert genome.length == 100000
	for name in ("dup2", "rep8"):
		seg = genome.segments[name]
		copies = {genome.segment_seq(name, i) for i in range(len(seg.copies))}
		assert len(copies) == 1, name
	a, b = genome.segment_seq("dupmm", 0), genome.segment_seq("dupmm", 1)
	assert [i for i in range(len(a)) if a[i] != b[i]] == list(range(25, 600, 50))
	assert genome.seq("chr2", 4999, 5101) == genome.seq("chr2", 4999, 5000) + "N" * 100 + genome.seq("chr2", 5100, 5101)
	assert genome.seq("chr3", 1000, 1200) == "CA" * 100
	assert genome.seq("chrX", 4000, 4100) == "A" * 100
	assert genome.seq("chr1", 20000, 21000) == ch.revcomp(genome.seq("chr2", 10000, 11000))
	other = ch.make_synthetic_genome(seed=1)
	assert other.contigs["chr1"] != genome.contigs["chr1"]


def test_fasta_written_and_read_back(genome):
	"""genome.fa has one '>name comment' header per contig, 60-column lines,
	and the sequence reads back unchanged."""
	text = genome.fasta_path.read_text()
	seqs, name = {}, None
	for line in text.splitlines():
		if line.startswith(">"):
			name = line[1:].split()[0]
			assert line == ">%s synthetic" % name
			seqs[name] = []
		else:
			assert len(line) <= 60
			seqs[name].append(line)
	assert {k: "".join(v) for k, v in seqs.items()} == dict(genome.contigs)


def _aligned_part(read):
	return read.seq[:len(read.seq) - len(read.tail)]


def test_single_end_truth_is_consistent(genome, sim):
	"""A read's genomic part equals the genome over [start, end) (reverse-
	complemented for '-') except at the recorded mismatch positions, and its
	tail is not genomic; segment reads list every copy as a locus."""
	reads = list(sim.se_reads(40))
	contig, start = sim.unique_locus(60)
	reads.append(sim.se_read(contig, start, 60, "-", mismatches=(0, 7, 59), tail="ACGTACGT"))
	reads.append(sim.se_read(contig, start, 60, "+", mismatches=(30,)))
	seg = sim.segment_read("dup2", 100, 50, strand="+", copy=1)
	reads.append(seg)
	for r in reads:
		g = genome.seq(r.contig, r.start, r.end)
		s = _aligned_part(r)
		s = s if r.strand == "+" else ch.revcomp(s)
		assert len(s) == r.end - r.start
		diffs = [r.start + i for i in range(len(g)) if g[i] != s[i]]
		assert diffs == list(r.mismatches), r
		assert len(r.qual) == len(r.seq)
		assert r.five_prime == (r.start if r.strand == "+" else r.end - 1)
	assert reads[-3].mismatches == (start, start + 52, start + 59)
	loci = seg.extra["loci"]
	assert len(loci) == 2
	for c, s, e in loci:
		x = genome.seq(c, s, e)
		assert seg.seq in (x, ch.revcomp(x))


def test_pair_truth_is_consistent(genome, sim):
	"""Read 1 on '+' starts at the fragment's left end and read 2 is the
	reverse complement of its right end (the other way round for '-');
	inserts shorter than the read continue into the Nextera adapters."""
	for strand in "+-":
		contig, start = sim.unique_locus(300)
		p = sim.pair(contig, start, 300, read_len=50, strand=strand, mismatches1=(5,))
		left = genome.seq(contig, start, start + 50)
		right = ch.revcomp(genome.seq(contig, start + 250, start + 300))
		want1, want2 = (left, right) if strand == "+" else (right, left)
		assert p.r2.seq == want2
		assert [i for i in range(50) if p.r1.seq[i] != want1[i]] == [5]
		assert p.fragment() == (contig, start, start + 300, strand)
		assert (p.r1.start, p.r1.end) == ((start, start + 50) if strand == "+" else (start + 250, start + 300))
		assert len(p.r1.mismatches) == 1
	contig, start = sim.unique_locus(30)
	short = sim.pair(contig, start, 30, read_len=50)
	frag = genome.seq(contig, start, start + 30)
	assert short.r1.seq == frag + ch.NEXTERA_R1[:20]
	assert short.r2.seq == ch.revcomp(frag) + ch.NEXTERA_R2[:20]
	assert short.r1.tail == ch.NEXTERA_R1[:20] and short.insert == 30


def test_readset_write_round_trip(tmp_path, sim, whitelist):
	"""ReadSet.write gives FASTQ (or gzipped FASTQ, or FASTA) that reads back
	to the same names, sequences and qualities, plus a barcode FASTQ in the
	same order when records carry barcodes."""
	barcodes, _ = whitelist
	reads = sim.pairs(5, barcodes=barcodes)
	files = reads.write(tmp_path / "plain")
	assert sorted(files) == ["bc", "r1", "r2"]

	def fastq(path, opener=open):
		with opener(str(path), "rt") as f:
			lines = f.read().splitlines()
		return [tuple(lines[i:i + 4]) for i in range(0, len(lines), 4)]

	assert fastq(files["r1"]) == [("@" + p.name, p.r1.seq, "+", p.r1.qual) for p in reads]
	assert fastq(files["r2"]) == [("@" + p.name, p.r2.seq, "+", p.r2.qual) for p in reads]
	assert [rec[1] for rec in fastq(files["bc"])] == [p.barcode for p in reads]
	gz = reads.write(tmp_path / "gz", gz=True, name_style="comment")
	assert fastq(gz["r2"], gzip.open)[0][0] == "@%s 2:N:0:1" % reads.records[0].name
	fa = ch.ReadSet(sim.se_reads(3)).write(tmp_path / "fa", fasta=True)
	assert fa["r1"].read_text().splitlines()[0].startswith(">")
	with pytest.raises(ValueError):
		ch.ReadSet(list(sim.se_reads(2)) + list(sim.pairs(2)))


def test_seed_encoding_and_hash():
	"""encode_seed follows chromap's 2-bit layout (A0 C1 G2 T3, first base
	highest, N and padding as A); hash64 is a bijection on 2k-bit keys
	(checked exhaustively for k = 5), which is what makes minimizer ties come
	only from identical k-mers."""
	assert ch.encode_seed("ACGT") == 27
	assert ch.encode_seed("ANGT") == 11
	assert ch.encode_seed("AC", length=4) == 16
	assert ch.decode_seed(ch.encode_seed("GATTACA"), 7) == "GATTACA"
	mask = (1 << 10) - 1
	assert len({ch.hash64(x, mask) for x in range(1 << 10)}) == 1 << 10
	assert ch.revcomp("AACGTN") == "NACGTT"


def test_reference_minimizers_definition():
	"""reference_minimizers: every window of w k-mers contributes its
	smallest-hash k-mer, w = 1 selects every k-mer, and positions are the
	0-based end of each k-mer with the strand of the smaller hash."""
	seq = "ACGTTGCAAGGCTTAGCATCGATCGGATACGATTACGCGAT"
	k, w = 5, 4
	mins = ch.reference_minimizers(seq, k, w)
	keys = [ch.kmer_minimizer_key(seq[e - k + 1:e + 1]) for e in range(k - 1, len(seq))]
	ends = [(hit >> 1) & 0xffffffff for _, hit in mins]
	for s in range(len(keys) - w + 1):
		best = min(range(s, s + w), key=lambda i: keys[i][0])
		assert best + k - 1 in ends
	assert len(ch.reference_minimizers(seq, k, 1)) == len(seq) - k + 1
	for h, hit in mins:
		e = (hit >> 1) & 0xffffffff
		assert (h, hit & 1) == ch.kmer_minimizer_key(seq[e - k + 1:e + 1])
	with pytest.raises(ValueError):
		ch.reference_minimizers("ACGNACGT", 3, 2)


def test_parsers_on_handwritten_files(tmp_path):
	"""BED, fragments, SAM, pairs, PAF, summary and stderr-statistics parsers
	on small handwritten files."""
	(tmp_path / "a.bed").write_text("chr1\t10\t60\tN\t60\t+\t1\nchr2\t5\t55\tN\t3\t-\t2\n")
	assert ch.parse_bed(tmp_path / "a.bed") == [ch.BedRecord("chr1", 10, 60, "N", 60, "+", 1),
		ch.BedRecord("chr2", 5, 55, "N", 3, "-", 2)]
	(tmp_path / "t.bed").write_text("chr1\t10\t60\tN\t60\t+\nchr1\t200\t250\tN\t60\t-\t1\n")
	with pytest.raises(ValueError):
		ch.parse_bed(tmp_path / "t.bed")
	assert ch.parse_bed(tmp_path / "t.bed", allow_six=True)[0].dups is None
	(tmp_path / "f.tsv").write_text("chrX\t100\t300\tACGTACGTACGTACGT\t3\n")
	assert ch.parse_fragments(tmp_path / "f.tsv") == [ch.FragmentRecord("chrX", 100, 300, "ACGTACGTACGTACGT", 3)]
	(tmp_path / "x.sam").write_text("@SQ\tSN:chr1\tLN:100\n"
		"r1\t99\tchr1\t11\t60\t50M\t=\t61\t100\tACGT\tIIII\tNM:i:1\tMD:Z:2A1\tCB:Z:AAAA\n")
	header, recs = ch.parse_sam(tmp_path / "x.sam")
	assert header == ["@SQ\tSN:chr1\tLN:100"]
	r = recs[0]
	assert (r.qname, r.flag, r.pos, r.tlen, r.tags["NM"], r.tags["MD"], r.tags["CB"]) == (
		"r1", 99, 11, 100, 1, "2A1", "AAAA")
	assert r.is_read1 and not r.is_reverse and not r.is_secondary
	(tmp_path / "x.pairs").write_text("## pairs format v1.0.0\n#columns: readID chrom1\n"
		"p1\tchr1\t11\tchr2\t40\t+\t-\tUU\t60\t50\n")
	header, recs = ch.parse_pairs(tmp_path / "x.pairs")
	assert len(header) == 2 and recs == [ch.PairsRecord("p1", "chr1", 11, "chr2", 40, "+", "-", "UU", 60, 50)]
	(tmp_path / "x.paf").write_text("q\t50\t0\t50\t+\tchr1\t1000\t10\t60\t50\t50\t60\n")
	assert ch.parse_paf(tmp_path / "x.paf")[0] == ch.PafRecord("q", 50, 0, 50, "+", "chr1", 1000, 10, 60, 50, 50, 60)
	(tmp_path / "s.csv").write_text("barcode,total,duplicate,unmapped,lowmapq,cachehit,fric,estfrip\n"
		",10,1,2,3,0,0.00000,0.00000\nnon-whitelist,4,0,4,0,0,0.00000,0.00000\n")
	s = ch.parse_summary(tmp_path / "s.csv")
	assert s[""] == {"total": 10, "duplicate": 1, "unmapped": 2, "lowmapq": 3, "cachehit": 0, "fric": 0.0,
		"estfrip": 0.0}
	assert "non-whitelist" in s and s.header[0] == "barcode"
	stats = ch.parse_stats("Number of reads: 400.\nMapped all reads in 0.01s.\n"
		"Number of output mappings (passed filters): 200\n")
	assert stats == {"Number of reads": 400, "Number of output mappings (passed filters)": 200}


def test_runner_pins_threads(index, sim):
	"""run_chromap passes -t itself and refuses a -t in extra_args, so the
	thread count of every run is the one the session asked for."""
	with pytest.raises(ValueError):
		ch.run_chromap(index, sim.se_reads(2), "-t 2")
	with pytest.raises(ValueError):
		ch.run_chromap(index, sim.se_reads(2), threads=ch.MAX_THREADS + 1)
	run = ch.run_chromap(index, sim.se_reads(5), threads=2)
	i = run.argv.index("-t")
	assert run.argv[i + 1] == "2"
	assert "Number of threads: 2" in run.stderr


def test_no_core_dump_lowers_the_child_core_limit(tmp_path):
	"""run_cli(no_core_dump=True) runs the child with RLIMIT_CORE 1 on Linux
	(0 elsewhere, never above the hard limit); without it the child keeps
	the parent's limit. Skipped where the resource module does not exist."""
	resource = pytest.importorskip("resource")
	argv = [sys.executable, "-c", "import resource; print(resource.getrlimit(resource.RLIMIT_CORE)[0])"]
	soft, hard = resource.getrlimit(resource.RLIMIT_CORE)
	want = 1 if sys.platform.startswith("linux") else 0
	if hard != resource.RLIM_INFINITY:
		want = min(want, hard)
	assert int(ch._execute(argv, tmp_path, 1, None, no_core_dump=True).stdout) == want
	assert int(ch._execute(argv, tmp_path, 1, None).stdout) == soft


def test_shared_runs_are_reused(index, genome):
	"""shared=True returns the same run for an identical request and a new one
	when anything differs."""
	reads = ch.ReadSimulator(genome, seed=ch.stable_seed("harness", "shared")).se_reads(4)
	a = ch.run_chromap(index, reads, shared=True)
	b = ch.run_chromap(index, ch.ReadSet(list(reads)), shared=True)
	c = ch.run_chromap(index, reads, "--TagAlign", shared=True)
	assert a is b and a is not c


def _sam_run(index, genome):
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("harness", "sam"))
	return ch.run_chromap(index, sim.pairs(30), "--SAM", shared=True)


def test_sam_parser_matches_samtools(index, genome):
	"""parse_sam reads chromap's SAM the way `samtools view` does (skipped
	without samtools on PATH)."""
	samtools = shutil.which("samtools")
	if samtools is None:
		pytest.skip("samtools not on PATH")
	run = _sam_run(index, genome)
	text = subprocess.run([samtools, "view", str(run.output)], stdout=subprocess.PIPE,
		stderr=subprocess.PIPE, check=True).stdout.decode()
	theirs = [line.split("\t") for line in text.splitlines()]
	ours = run.sam()[1]
	assert [(r.qname, r.flag, r.rname, r.pos, r.cigar, r.tlen) for r in ours] == [
		(x[0], int(x[1]), x[2], int(x[3]), x[5], int(x[8])) for x in theirs]


def test_sam_parser_matches_pysam(index, genome):
	"""parse_sam agrees with pysam on every field (skipped without pysam)."""
	pysam = pytest.importorskip("pysam")
	run = _sam_run(index, genome)
	with pysam.AlignmentFile(str(run.output), "r", check_sq=False) as f:
		theirs = [(a.query_name, a.flag, a.reference_name, a.reference_start + 1, a.mapping_quality,
			a.cigarstring, a.template_length, a.get_tag("NM"), a.get_tag("MD")) for a in f]
	ours = [(r.qname, r.flag, r.rname, r.pos, r.mapq, r.cigar, r.tlen, r.tags["NM"], r.tags["MD"])
		for r in run.sam()[1]]
	assert ours == theirs
