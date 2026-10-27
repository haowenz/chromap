"""Area D: --read-format through the command line. Barcodes and genomic
sequence are cut out of the read files by the README's grammar
[r1|r2|bc]:start:end:strand (inclusive ranges, -1 the read's end, several
ranges concatenated, "-" reverse-complementing the barcode), and the
fragments must come back at the genomic part's true position with the
barcode helpers_d.correct_barcode gives for the extracted barcode and its
(reversed, for "-") qualities."""

import collections

import pytest

import chromap_helpers as ch
import helpers_d as hd

BC_LEN = 16

# Single-end layout (71 bp): barcode part 1 (8) | genomic 0-24 | 5 bp spacer
# | genomic 25-49 | barcode part 2 (8). The barcode is the reverse
# complement of part 1 + part 2.
SE_FORMAT = "bc:0:7,bc:63:70:-,r1:8:32,r1:38:62"


def _raw_index(k):
	"""Raw read position of base k of the extracted single-end barcode:
	extracted = revcomp(raw[0:8] + raw[63:71]), so base k is concatenated
	base 15 - k."""
	m = BC_LEN - 1 - k
	return m if m < 8 else 63 + (m - 8)


@pytest.fixture(scope="module")
def se_layout(genome, tmp_path_factory):
	"""58 reads: 37 with whitelisted barcodes, 10 one substitution away, 10
	exact reads of two designed whitelist barcodes a and b (5 each), and one
	read one substitution from both a and b whose only low quality ('$')
	sits under a's mismatch in extracted orientation.
	With the quality reversed as documented, a's score is 0.501 * 5 against
	1e-4 * 5 for b, so it is corrected to a; read unreversed both
	qualities would be 'I', a tie, and the read would be dropped."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_d_read_format", "se"))
	rng = hd.stable_rng("test_d_read_format", "se")
	picker = hd.LocusPicker(sim)
	base = ch.make_whitelist(12, length=BC_LEN, seed=67, min_dist=4)
	x = hd.isolated_barcodes(1, BC_LEN, rng, base, 5)[0]
	pa, pb = sorted(rng.sample(range(BC_LEN), 2))
	wa, wb = hd.mutate(x, [pa]), hd.mutate(x, [pb])
	whitelist = base + [wa, wb]
	plan = [(w, "I" * BC_LEN) for w in base for _ in range(2)] + [(base[i], "I" * BC_LEN) for i in range(6)]
	plan += [(hd.mutate(base[i], [rng.randrange(BC_LEN)]), hd.random_quality(rng, BC_LEN)) for i in range(10)]
	plan += [(wa, "I" * BC_LEN)] * 5 + [(wb, "I" * BC_LEN)] * 5
	q = ["I"] * BC_LEN
	q[pa] = "$"
	plan.append((x, "".join(q)))
	plan += [(base[i], hd.random_quality(rng, BC_LEN)) for i in range(7)]
	records, fastq = [], []
	for i, (bc, bq) in enumerate(plan):
		c, s = picker.take(50)
		g = sim.se_read(c, s, 50, rng.choice("+-"), name="f%03d" % i, barcode=bc, barcode_qual=bq)
		rc = ch.revcomp(bc)
		raw_q = ["I"] * 71
		for k in range(BC_LEN):
			raw_q[_raw_index(k)] = bq[k]
		gq = hd.random_quality(rng, 50, "5?I")
		seq = rc[:8] + g.seq[:25] + "".join(rng.choice("ACGT") for _ in range(5)) + g.seq[25:] + rc[8:]
		qual = "".join(raw_q[:8]) + gq[:25] + "IIIII" + gq[25:] + "".join(raw_q[63:])
		assert hd.apply_read_format(seq, qual, SE_FORMAT.split(",")[:2]) == (bc, bq)
		assert hd.apply_read_format(seq, qual, SE_FORMAT.split(",")[2:])[0] == g.seq
		records.append(g)
		fastq += ["@%s" % g.name, seq, "+", qual]
	d = tmp_path_factory.mktemp("d_read_format_se")
	path = ch.write_lines(d / "reads.fq", fastq)
	wl_path = ch.write_lines(d / "whitelist.txt", whitelist)
	return dict(records=records, whitelist=whitelist, path=path, wl_path=wl_path, design=(x, wa, wb))


def test_single_end_barcode_and_read_cut_from_one_file(index, se_layout):
	"""-1 and -b name the same file; --read-format bc:0:7,bc:63:70:-,
	r1:8:32,r1:38:62 takes the barcode from both ends, reverse-complemented,
	and the read from two ranges around a spacer. Every fragment is the
	genomic part's true interval with the corrected barcode, including the
	designed read that only the reversed qualities resolve (to a)."""
	run = ch.run_chromap(index, {"r1": se_layout["path"], "bc": se_layout["path"]},
		["--read-format", SE_FORMAT, "--barcode-whitelist", str(se_layout["wl_path"])])
	want = hd.expected_fragments(se_layout["records"], se_layout["whitelist"])
	got = collections.Counter((f.chrom, f.start, f.end, f.barcode) for f in run.fragments())
	assert got == want
	x, wa, wb = se_layout["design"]
	assert sum(n for k, n in want.items() if k[3] == wa) == 6
	assert len(run.fragments()) == len(se_layout["records"])


def test_paired_end_barcode_in_read_one(index, genome, tmp_path):
	"""The README's 10x layout: read 1 is a 16 bp barcode followed by
	genomic sequence, read 2 is genomic sequence followed by 12 bases of
	junk; -b names the read 1 file and --read-format bc:0:15,r1:16:-1,r2:0:49
	cuts both. Without --barcode-whitelist chromap warns and corrects
	nothing: fragments are the true fragments with the barcodes as read
	(10 of the 40 are one substitution from another barcode of the set)."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_d_read_format", "pe"))
	rng = hd.stable_rng("test_d_read_format", "pe")
	picker = hd.LocusPicker(sim)
	wl = ch.make_whitelist(10, length=BC_LEN, seed=71, min_dist=4)
	recs, r1, r2 = [], [], []
	for i in range(40):
		ins = rng.randint(150, 400)
		c, s = picker.take(ins)
		w = wl[i % len(wl)]
		bc = w if i < 30 else hd.mutate(w, [rng.randrange(BC_LEN)])
		bq = hd.random_quality(rng, BC_LEN)
		p = sim.pair(c, s, ins, strand="+-"[i % 2], barcode=bc, barcode_qual=bq)
		recs.append(p)
		junk = "".join(rng.choice("ACGT") for _ in range(12))
		r1 += ["@%s" % p.name, bc + p.r1.seq, "+", bq + p.r1.qual]
		r2 += ["@%s" % p.name, p.r2.seq + junk, "+", p.r2.qual + "I" * 12]
	p1 = ch.write_lines(tmp_path / "r1.fq", r1)
	p2 = ch.write_lines(tmp_path / "r2.fq", r2)
	run = ch.run_chromap(index, {"r1": p1, "r2": p2, "bc": p1}, ["--read-format", "bc:0:15,r1:16:-1,r2:0:49",
		"--cache-size", "2000000"])
	want = collections.Counter((p.contig, p.start, p.end, p.barcode) for p in recs)
	got = collections.Counter((f.chrom, f.start, f.end, f.barcode) for f in run.fragments())
	assert got == want
	assert "WARNING: there are input barcode files but a barcode whitelist file is missing!" in run.stderr
	stats = run.stats()
	assert (stats["Number of barcodes in whitelist"], stats["Number of corrected barcodes"]) == (0, 0)


@pytest.mark.parametrize("fmt", [
	"r1:0",
	"bc:0:15:+:x",
	"x1:0:5",
	"bc:0:15,r3:0:-1",
	"r1:0:-1,,",
])
def test_malformed_read_format_exits(index, fmt, tmp_path):
	"""A field without an end, with a fifth part, naming something other
	than r1, r2 or bc, or empty, stops chromap with status 255 and "Unknown
	read format: " followed by the whole string."""
	reads = ch.ReadSimulator(index.genome, seed=3).se_reads(2)
	files = reads.write(tmp_path)
	run = ch.run_cli(["-x", index.path, "-r", index.fasta, "-1", files["r1"], "-o", tmp_path / "out.bed",
		"--read-format", fmt], cwd=tmp_path)
	assert run.returncode == 255
	assert "Unknown read format: %s" % fmt in run.stderr
	assert not (tmp_path / "out.bed").exists()
