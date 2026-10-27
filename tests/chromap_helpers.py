"""Harness for the chromap unit tests.

A deterministic synthetic genome with planted repeats, a read simulator that
records each read's true origin, a runner for the chromap binary (index
building and mapping) and parsers for every output format chromap writes
(BED, fragments, TagAlign, SAM, pairs, PAF, summary CSV, stderr statistics).
Only the standard library is used.

Coordinates are 0-based and half-open throughout ([start, end)), as in
chromap's BED output; SAM POS and pairs positions are 1-based and are
converted by the tests that read them.
"""

import gzip
import hashlib
import os
import random
import re
import shlex
import subprocess
import sys
import tempfile
import zlib
from collections import OrderedDict, namedtuple
from pathlib import Path

try:
	import resource as _resource  # POSIX only
except ImportError:
	_resource = None

import chromap_build

REPO_ROOT = chromap_build.REPO_ROOT
GENOME_SEED = 20261024
DEFAULT_READ_LEN = 50
MAX_THREADS = 8
DEFAULT_TIMEOUT = 300
# chromap's defaults (src/mapping_parameters.h, src/index_parameters.h).
DEFAULT_KMER = 17
DEFAULT_WINDOW = 7
DEFAULT_MAPQ_THRESHOLD = 30
# Nextera (Tn5) adapter read-through, as in ATAC-seq libraries.
NEXTERA_R1 = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"
NEXTERA_R2 = "CTGTCTCTTATACACATCTGACGCTGCCGACGA"


# ---------------------------------------------------------------- sequences

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
NT4 = {"A": 0, "C": 1, "G": 2, "T": 3, "a": 0, "c": 1, "g": 2, "t": 3}


def revcomp(seq):
	return seq.translate(_COMP)[::-1]


def hamming(a, b):
	if len(a) != len(b):
		raise ValueError("hamming() needs equal lengths")
	return sum(1 for x, y in zip(a, b) if x != y)


def substitute(base):
	"""A fixed different base for `base` (A->C->G->T->A), so simulated
	mismatches are deterministic."""
	return {"A": "C", "C": "G", "G": "T", "T": "A"}.get(base.upper(), "A")


def encode_seed(seq, start=0, length=None):
	"""2-bit encoding of seq[start:start+length], first base in the highest
	bits: A=0, C=1, G=2, T=3. N (any non-ACGT) encodes as A, and positions past
	the end of `seq` pad with A, as chromap documents for its seeds and
	barcode keys (src/utils.h)."""
	if length is None:
		length = len(seq) - start
	key = 0
	for i in range(length):
		p = start + i
		key <<= 2
		if p < len(seq):
			key |= NT4.get(seq[p], 0)
	return key


def decode_seed(key, length):
	"""Inverse of encode_seed for an N-free sequence of `length` bases."""
	return "".join("ACGT"[(key >> (2 * (length - 1 - i))) & 3] for i in range(length))


def hash64(key, mask):
	"""Thomas Wang's invertible 64-bit integer hash restricted to `mask`, the
	k-mer hash of minimap2 and chromap."""
	m64 = (1 << 64) - 1
	key = ((~key) + (key << 21)) & mask
	key = key ^ (key >> 24)
	key = ((key + (key << 3)) + (key << 8)) & mask
	key = key ^ (key >> 14)
	key = ((key + (key << 2)) + (key << 4)) & mask
	key = key ^ (key >> 28)
	key = (key + (key << 31)) & mask
	return key & m64


def kmer_minimizer_key(kmer):
	"""(hash, strand) chromap assigns to an N-free k-mer: both strands are
	hashed, the strand with the smaller hash is chosen (0 forward, 1 reverse;
	ties go to reverse), and the stored hash is that hash hashed once more."""
	k = len(kmer)
	mask = (1 << (2 * k)) - 1
	fwd = encode_seed(kmer)
	rev = encode_seed(revcomp(kmer))
	h = (hash64(fwd, mask), hash64(rev, mask))
	strand = 0 if h[0] < h[1] else 1
	return hash64(h[strand], mask), strand


def reference_minimizers(seq, k, w, seq_index=0):
	"""Minimizers of `seq` by definition: every k-mer whose hash is the
	smallest (ties included) in at least one window of w consecutive k-mers.
	Returns a list of (hash, hit) sorted by position, where hit = ((seq_index
	<< 32 | end position) << 1) | strand and the end position is the 0-based
	index of the k-mer's last base (chromap's Minimizer encoding,
	src/minimizer.h). For a sequence shorter than w + k - 1 the one smallest
	k-mer is returned.

	Valid for N-free sequences with odd k (no k-mer equals its reverse
	complement). Around N bases minimap2-style winnowing differs from this
	definition; tests of that behaviour need their own expectations."""
	if any(c not in NT4 for c in seq):
		raise ValueError("reference_minimizers() needs an N-free sequence")
	if k % 2 == 0:
		raise ValueError("reference_minimizers() needs an odd k")
	keys = []
	for end in range(k - 1, len(seq)):
		h, strand = kmer_minimizer_key(seq[end - k + 1:end + 1])
		keys.append((h, ((seq_index << 32 | end) << 1) | strand))
	if not keys:
		return []
	chosen = set()
	if len(keys) < w:
		best = min(h for h, _ in keys)
		chosen.update(i for i, (h, _) in enumerate(keys) if h == best)
	else:
		for s in range(len(keys) - w + 1):
			window = keys[s:s + w]
			best = min(h for h, _ in window)
			chosen.update(s + i for i, (h, _) in enumerate(window) if h == best)
	return [keys[i] for i in sorted(chosen)]


def stable_seed(*parts):
	"""A seed that depends only on `parts` (not on Python's hash salt)."""
	return zlib.crc32("::".join(str(p) for p in parts).encode())


# ------------------------------------------------------- synthetic genome

CONTIG_SPECS = (("chr1", 50000), ("chr2", 30000), ("chr3", 15000), ("chrX", 5000))

# name, length, copies [(contig, 0-based start, strand)], substitution offsets
# (within the segment, in its own orientation) applied to copies 2+.
SEGMENT_SPECS = (
	("dup2", 1000, (("chr1", 20000, "+"), ("chr2", 10000, "-")), ()),
	("rep8", 400, (("chr1", 30000, "+"), ("chr1", 32000, "+"), ("chr1", 34000, "-"),
		("chr2", 20000, "+"), ("chr2", 22000, "-"), ("chr3", 3000, "+"), ("chr3", 6000, "+"),
		("chrX", 1000, "+")), ()),
	("dupmm", 600, (("chr1", 40000, "+"), ("chr3", 10000, "+")), tuple(range(25, 600, 50))),
)

# Runs of N: (contig, start, end).
N_RUN_SPECS = (("chr2", 5000, 5100),)

# Low-complexity sequence: (contig, start, repeat unit, length).
LOW_COMPLEXITY_SPECS = (("chr3", 1000, "CA", 200), ("chrX", 4000, "A", 100))

Segment = namedtuple("Segment", "name length copies substitutions")
SegmentCopy = namedtuple("SegmentCopy", "contig start end strand")
Feature = namedtuple("Feature", "contig start end kind name")


class SyntheticGenome(object):
	"""Uniform random chromosomes with planted repeats, N runs and
	low-complexity stretches.

	contigs: OrderedDict name -> sequence (upper case, may contain N);
	segments: OrderedDict name -> Segment, every copy listed; features: list
	of Feature(contig, start, end, kind, name) for everything planted
	(kind is "segment", "N" or "low_complexity"); fasta_path (set by write()).
	"""

	def __init__(self, contigs, segments, features, seed):
		self.contigs = contigs
		self.segments = segments
		self.features = features
		self.seed = seed
		self.fasta_path = None

	@property
	def names(self):
		return list(self.contigs)

	@property
	def length(self):
		return sum(len(s) for s in self.contigs.values())

	def rid(self, contig):
		"""Index of `contig` in the FASTA (chromap's reference id)."""
		return self.names.index(contig)

	def seq(self, contig, start, end):
		"""Forward-strand sequence of contig[start:end] (0-based, half-open)."""
		return self.contigs[contig][start:end]

	def segment_seq(self, name, copy=0):
		"""Sequence of one copy of a planted segment, in the segment's own
		orientation."""
		c = self.segments[name].copies[copy]
		s = self.seq(c.contig, c.start, c.end)
		return s if c.strand == "+" else revcomp(s)

	def overlaps_feature(self, contig, start, end, pad=0):
		return any(f.contig == contig and f.start - pad < end and start < f.end + pad for f in self.features)

	def write(self, directory, name="genome.fa", line_width=60):
		"""Write the FASTA into `directory` and set fasta_path; returns self.
		Headers carry a comment ("chr1 synthetic") so name parsing at the
		first whitespace is exercised."""
		directory = Path(directory).resolve()
		directory.mkdir(parents=True, exist_ok=True)
		self.fasta_path = directory / name
		with open(str(self.fasta_path), "w") as f:
			for contig, s in self.contigs.items():
				f.write(">%s synthetic\n" % contig)
				for i in range(0, len(s), line_width):
					f.write(s[i:i + line_width] + "\n")
		return self


def make_synthetic_genome(seed=GENOME_SEED, contig_specs=CONTIG_SPECS, segment_specs=SEGMENT_SPECS,
		n_run_specs=N_RUN_SPECS, low_complexity_specs=LOW_COMPLEXITY_SPECS):
	"""Build a SyntheticGenome deterministically from `seed`. The default is
	chr1 50 kb, chr2 30 kb, chr3 15 kb, chrX 5 kb (100 kb in all) with:
	dup2 (1 kb, two copies, the second reverse-complemented), rep8 (400 bp,
	eight copies on both strands), dupmm (600 bp, two copies differing at 12
	positions), 100 N at chr2:5000-5100, (CA)100 at chr3:1000-1200 and
	poly-A at chrX:4000-4100."""
	rng = random.Random(seed)
	seqs = OrderedDict((name, [rng.choice("ACGT") for _ in range(n)]) for name, n in contig_specs)
	painted = {}

	def paint(contig, start, bases, what):
		if start < 0 or start + len(bases) > len(seqs[contig]):
			raise ValueError("%s does not fit on %s" % (what, contig))
		for i, b in enumerate(bases):
			key = (contig, start + i)
			if key in painted:
				raise ValueError("layout conflict at %s:%d (%s)" % (contig, start + i, what))
			painted[key] = b
			seqs[contig][start + i] = b

	features = []
	segments = OrderedDict()
	for name, length, copies, subs in segment_specs:
		base = "".join(rng.choice("ACGT") for _ in range(length))
		placed = []
		for ci, (contig, start, strand) in enumerate(copies):
			s = list(base)
			if ci > 0:
				for off in subs:
					s[off] = substitute(s[off])
			s = "".join(s)
			paint(contig, start, s if strand == "+" else revcomp(s), name)
			placed.append(SegmentCopy(contig, start, start + length, strand))
			features.append(Feature(contig, start, start + length, "segment", name))
		segments[name] = Segment(name, length, tuple(placed), tuple(subs))
	for contig, start, end in n_run_specs:
		paint(contig, start, "N" * (end - start), "N run")
		features.append(Feature(contig, start, end, "N", "N run"))
	for contig, start, unit, length in low_complexity_specs:
		bases = (unit * (length // len(unit) + 1))[:length]
		paint(contig, start, bases, "low complexity")
		features.append(Feature(contig, start, start + length, "low_complexity", unit))
	contigs = OrderedDict((name, "".join(s)) for name, s in seqs.items())
	return SyntheticGenome(contigs, segments, features, seed)


# -------------------------------------------------------------- reads

class SimRead(object):
	"""A simulated read and where it came from.

	name, seq, qual: as written to FASTQ (seq 5' to 3' as sequenced).
	contig, start, end, strand: the genomic interval the read's genomic part
	covers (0-based, half-open) and its strand; contig is None for reads that
	should not map. tail: bases appended at the 3' end that are not from the
	genome (adapter); they are part of seq but not of [start, end).
	mismatches: genomic positions where the read differs from the genome.
	barcode: the cell barcode read as sequenced (None for bulk reads).
	kind, extra: free-form labels (e.g. extra["loci"] for repeat reads).
	"""

	def __init__(self, name, seq, qual, contig, start, end, strand, mismatches=(), tail="", kind="genomic",
			barcode=None, barcode_qual=None, extra=None):
		self.name = name
		self.seq = seq
		self.qual = qual
		self.contig = contig
		self.start = start
		self.end = end
		self.strand = strand
		self.mismatches = tuple(mismatches)
		self.tail = tail
		self.kind = kind
		self.barcode = barcode
		self.barcode_qual = barcode_qual
		self.extra = extra or {}

	@property
	def paired(self):
		return False

	def interval(self):
		"""(contig, start, end, strand) of the genomic part."""
		return (self.contig, self.start, self.end, self.strand)

	@property
	def five_prime(self):
		"""0-based genomic position of the read's 5'-most genomic base."""
		return self.start if self.strand == "+" else self.end - 1

	def __repr__(self):
		return "SimRead(%s %s:%s-%s%s mm=%s)" % (self.name, self.contig, self.start, self.end, self.strand,
			list(self.mismatches))


class SimPair(object):
	"""A simulated read pair: r1 and r2 (SimRead, with the same name), the
	fragment (contig, start, end) both came from, r1's strand, and the
	barcode read as sequenced (None for bulk pairs)."""

	def __init__(self, name, r1, r2, contig, start, end, strand, kind="genomic", barcode=None, barcode_qual=None,
			extra=None):
		self.name = name
		self.r1 = r1
		self.r2 = r2
		self.contig = contig
		self.start = start
		self.end = end
		self.strand = strand
		self.kind = kind
		self.barcode = barcode
		self.barcode_qual = barcode_qual
		self.extra = extra or {}

	@property
	def paired(self):
		return True

	@property
	def insert(self):
		return None if self.contig is None else self.end - self.start

	def fragment(self):
		"""(contig, start, end, strand of read 1)."""
		return (self.contig, self.start, self.end, self.strand)

	def __repr__(self):
		return "SimPair(%s %s:%s-%s r1%s)" % (self.name, self.contig, self.start, self.end, self.strand)


def _fastq_name(name, mate, style):
	if style == "plain":
		return name
	if style == "comment":
		return "%s %d:N:0:1" % (name, mate)
	if style == "slash":
		return "%s/%d" % (name, mate)
	raise ValueError("unknown name style %r" % style)


def _open_out(path, gz):
	return gzip.open(str(path), "wt") if gz else open(str(path), "w")


class ReadSet(object):
	"""An ordered collection of SimRead (single-end) or SimPair (paired-end)
	records. Iterable; `by_name` maps names to records; `+` concatenates."""

	def __init__(self, records):
		self.records = list(records)
		kinds = {r.paired for r in self.records}
		if len(kinds) > 1:
			raise ValueError("cannot mix single-end and paired records")
		self.paired = kinds.pop() if kinds else False
		names = [r.name for r in self.records]
		if len(set(names)) != len(names):
			raise ValueError("read names must be unique")
		self.by_name = OrderedDict((r.name, r) for r in self.records)

	def __iter__(self):
		return iter(self.records)

	def __len__(self):
		return len(self.records)

	def __add__(self, other):
		return ReadSet(self.records + list(other))

	@property
	def has_barcodes(self):
		with_bc = [r.barcode is not None for r in self.records]
		if any(with_bc) and not all(with_bc):
			raise ValueError("either every record or none carries a barcode")
		return bool(with_bc) and all(with_bc)

	def mates(self, i):
		"""SimRead list of mate i (1 or 2) for paired sets."""
		if not self.paired:
			raise ValueError("single-end read set")
		return [r.r1 if i == 1 else r.r2 for r in self.records]

	def digest(self):
		"""Content hash, used to share identical runs."""
		h = hashlib.sha1()
		for r in self.records:
			reads = [r.r1, r.r2] if r.paired else [r]
			for x in reads:
				h.update(("%s\t%s\t%s\n" % (x.name, x.seq, x.qual)).encode())
			h.update(("%s\t%s\n" % (r.barcode, r.barcode_qual)).encode())
		return h.hexdigest()

	def write(self, directory, prefix="reads", gz=False, fasta=False, name_style="plain"):
		"""Write FASTQ (or FASTA) files and return {"r1": path[, "r2": path]
		[, "bc": path]}. The barcode file is written when records carry
		barcodes; it is always FASTQ (chromap reads barcode qualities)."""
		directory = Path(directory)
		directory.mkdir(parents=True, exist_ok=True)
		ext = (".fa" if fasta else ".fq") + (".gz" if gz else "")
		mates = [1, 2] if self.paired else [1]
		paths = OrderedDict()
		for m in mates:
			p = directory / ("%s_%d%s" % (prefix, m, ext))
			with _open_out(p, gz) as f:
				for rec in self.records:
					r = (rec.r1 if m == 1 else rec.r2) if self.paired else rec
					name = _fastq_name(r.name, m, name_style)
					if fasta:
						f.write(">%s\n%s\n" % (name, r.seq))
					else:
						f.write("@%s\n%s\n+\n%s\n" % (name, r.seq, r.qual))
			paths["r%d" % m] = p
		if self.has_barcodes:
			p = directory / ("%s_bc.fq%s" % (prefix, ".gz" if gz else ""))
			with _open_out(p, gz) as f:
				for rec in self.records:
					qual = rec.barcode_qual or "I" * len(rec.barcode)
					f.write("@%s\n%s\n+\n%s\n" % (rec.name, rec.barcode, qual))
			paths["bc"] = p
		return paths


class ReadSimulator(object):
	"""Reads from a SyntheticGenome with their true origin recorded.

	Offsets given to `mismatches` count along the read as sequenced (5' to
	3'). Use the `sim` fixture (seeded per test) or
	ReadSimulator(genome, seed=stable_seed("label")) in module fixtures."""

	def __init__(self, genome, seed=1, read_len=DEFAULT_READ_LEN, qual_char="I", name_prefix="r"):
		self.g = genome
		self.rng = random.Random(seed)
		self.read_len = read_len
		self.qual_char = qual_char
		self.name_prefix = name_prefix
		self._count = 0
		self.adapter1 = NEXTERA_R1 + "".join(self.rng.choice("ACGT") for _ in range(300))
		self.adapter2 = NEXTERA_R2 + "".join(self.rng.choice("ACGT") for _ in range(300))

	def next_name(self, prefix=None):
		self._count += 1
		return "%s%06d" % (prefix or self.name_prefix, self._count)

	def unique_locus(self, length, margin=200, contigs=None, pad=50):
		"""(contig, start) of a random interval of `length` bases at least
		`margin` from the contig ends and `pad` from any planted feature."""
		names = list(contigs) if contigs else self.g.names
		weights = [len(self.g.contigs[c]) for c in names]
		for _ in range(10000):
			contig = self.rng.choices(names, weights)[0]
			hi = len(self.g.contigs[contig]) - margin - length
			if hi < margin:
				continue
			start = self.rng.randint(margin, hi)
			if not self.g.overlaps_feature(contig, start, start + length, pad=pad):
				return contig, start
		raise RuntimeError("no unique locus of length %d" % length)

	def _apply_mismatches(self, seq, offsets):
		s = list(seq)
		for off in offsets:
			if not 0 <= off < len(s):
				raise ValueError("mismatch offset %d outside read of length %d" % (off, len(s)))
			s[off] = substitute(s[off])
		return "".join(s)

	def se_read(self, contig, start, length=None, strand="+", mismatches=(), name=None, qual=None, tail="",
			barcode=None, barcode_qual=None, kind="genomic"):
		"""A read covering contig[start:start+length] on `strand`, with
		substitutions at read offsets `mismatches` and `tail` appended at its
		3' end."""
		length = length or self.read_len
		g = self.g.seq(contig, start, start + length)
		if len(g) != length:
			raise ValueError("read runs off %s" % contig)
		s = g if strand == "+" else revcomp(g)
		s = self._apply_mismatches(s, mismatches)
		genomic_mm = sorted(start + o if strand == "+" else start + length - 1 - o for o in mismatches)
		seq = s + tail
		return SimRead(name or self.next_name(), seq, qual or self.qual_char * len(seq), contig, start,
			start + length, strand, genomic_mm, tail=tail, kind=kind, barcode=barcode, barcode_qual=barcode_qual)

	def pair(self, contig, start, insert, read_len=None, strand="+", mismatches1=(), mismatches2=(), name=None,
			barcode=None, barcode_qual=None, qual1=None, qual2=None, kind="genomic"):
		"""A read pair from the fragment contig[start:start+insert]. strand
		"+": read 1 reads the fragment's left end forward and read 2 its right
		end reverse-complemented; "-": the other way round. When the insert is
		shorter than the read, each read continues into its Nextera adapter."""
		L = read_len or self.read_len
		end = start + insert
		frag = self.g.seq(contig, start, end)
		if len(frag) != insert:
			raise ValueError("fragment runs off %s" % contig)
		n = min(L, insert)
		left = frag[:n]
		right = revcomp(frag[insert - n:])
		left_iv = (start, start + n)
		right_iv = (end - n, end)
		if strand == "+":
			s1, iv1, st1, s2, iv2, st2 = left, left_iv, "+", right, right_iv, "-"
		else:
			s1, iv1, st1, s2, iv2, st2 = right, right_iv, "-", left, left_iv, "+"
		t1 = self.adapter1[:L - n]
		t2 = self.adapter2[:L - n]
		name = name or self.next_name("p")

		def mate(seq, iv, st, mm, tail, qual):
			s = self._apply_mismatches(seq, [o for o in mm if o < len(seq)])
			tail_mm = [o - len(seq) for o in mm if o >= len(seq)]
			tail = self._apply_mismatches(tail, tail_mm) if tail_mm else tail
			gm = sorted(iv[0] + o if st == "+" else iv[1] - 1 - o for o in mm if o < len(seq))
			full = s + tail
			return SimRead(name, full, qual or self.qual_char * len(full), contig, iv[0], iv[1], st, gm,
				tail=tail, kind=kind)

		r1 = mate(s1, iv1, st1, mismatches1, t1, qual1)
		r2 = mate(s2, iv2, st2, mismatches2, t2, qual2)
		return SimPair(name, r1, r2, contig, start, end, strand, kind=kind, barcode=barcode,
			barcode_qual=barcode_qual)

	def segment_read(self, segment, offset, length=None, strand="+", copy=0, mismatches=(), name=None):
		"""A read from one copy of a planted repeat; extra["loci"] lists the
		(contig, start, end) intervals of the same bases in every copy."""
		seg = self.g.segments[segment]
		length = length or self.read_len
		if offset + length > seg.length:
			raise ValueError("read runs off segment %s" % segment)
		c = seg.copies[copy]
		if c.strand == "+":
			gstart = c.start + offset
			gstrand = strand
		else:
			gstart = c.end - offset - length
			gstrand = "-" if strand == "+" else "+"
		r = self.se_read(c.contig, gstart, length, gstrand, mismatches, name=name, kind="segment")
		loci = []
		for cc in seg.copies:
			s = cc.start + offset if cc.strand == "+" else cc.end - offset - length
			loci.append((cc.contig, s, s + length))
		r.extra["loci"] = loci
		r.extra["segment"] = segment
		return r

	def random_read(self, length=None, name=None):
		"""A read of random sequence that should not map (contig None)."""
		length = length or self.read_len
		seq = "".join(self.rng.choice("ACGT") for _ in range(length))
		return SimRead(name or self.next_name(), seq, self.qual_char * length, None, None, None, None,
			kind="random")

	def se_reads(self, n, length=None, strands="+-", margin=200, barcodes=None):
		"""n reads at unique loci, strands drawn from `strands`; barcodes, if
		given, is a list of barcode sequences drawn uniformly per read."""
		length = length or self.read_len
		out = []
		for _ in range(n):
			contig, start = self.unique_locus(length, margin=margin)
			bc = self.rng.choice(barcodes) if barcodes else None
			out.append(self.se_read(contig, start, length, self.rng.choice(strands), barcode=bc))
		return ReadSet(out)

	def pairs(self, n, insert=(150, 400), read_len=None, strands="+-", margin=200, barcodes=None):
		"""n read pairs from fragments at unique loci with inserts drawn
		uniformly from `insert` (inclusive)."""
		out = []
		for _ in range(n):
			ins = self.rng.randint(insert[0], insert[1])
			contig, start = self.unique_locus(ins, margin=margin)
			bc = self.rng.choice(barcodes) if barcodes else None
			out.append(self.pair(contig, start, ins, read_len=read_len, strand=self.rng.choice(strands), barcode=bc))
		return ReadSet(out)

	def shuffled(self, readset):
		"""The same records in a random order (seeded)."""
		recs = list(readset)
		self.rng.shuffle(recs)
		return ReadSet(recs)


def make_whitelist(n, length=16, seed=7, min_dist=3):
	"""n random barcodes, pairwise at least `min_dist` mismatches apart."""
	rng = random.Random(seed)
	out = []
	while len(out) < n:
		b = "".join(rng.choice("ACGT") for _ in range(length))
		if all(hamming(b, x) >= min_dist for x in out):
			out.append(b)
	return out


def write_lines(path, lines):
	"""Write `lines` one per line (newline-terminated) and return the path."""
	path = Path(path)
	path.parent.mkdir(parents=True, exist_ok=True)
	path.write_text("".join(str(x) + "\n" for x in lines))
	return path


# ------------------------------------------------------------ running

def find_chromap():
	"""Path of the binary under test ($CHROMAP_BIN or a cached build)."""
	return chromap_build.find_chromap()


def default_threads():
	"""-t used when a test does not pass threads: $CHROMAP_TEST_THREADS or 1."""
	n = int(os.environ.get("CHROMAP_TEST_THREADS", "1"))
	if not 1 <= n <= MAX_THREADS:
		raise ValueError("CHROMAP_TEST_THREADS must be between 1 and %d" % MAX_THREADS)
	return n


def subprocess_env(threads):
	"""Environment for chromap and the C++ drivers: OpenMP pinned to
	`threads` so no pool sizes itself to the machine."""
	env = dict(os.environ)
	env["OMP_NUM_THREADS"] = str(threads)
	env.pop("OMP_DYNAMIC", None)
	return env


_RUN_ROOT = [None]
_SHARED_RUNS = {}


def set_run_root(path):
	"""Directory under which each run gets a fresh directory (conftest points
	it at a pytest tmp directory on local disk)."""
	_RUN_ROOT[0] = Path(path).resolve()
	_RUN_ROOT[0].mkdir(parents=True, exist_ok=True)


def fresh_dir(name=None):
	root = _RUN_ROOT[0]
	if root is None:
		root = Path(tempfile.mkdtemp(prefix="chromap_tests_"))
		_RUN_ROOT[0] = root
	safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name or "run")[:60]
	return Path(tempfile.mkdtemp(prefix=safe + "_", dir=str(root)))


def split_args(extra_args):
	if extra_args is None:
		return []
	if isinstance(extra_args, str):
		return shlex.split(extra_args)
	return [str(a) for a in extra_args]


class ChromapError(AssertionError):
	"""chromap exited non-zero where the test expected success."""


def _failure(argv, proc):
	return ChromapError("chromap exited %d\nargv: %s\nstderr (tail):\n%s\nstdout (tail):\n%s" % (
		proc.returncode, " ".join(shlex.quote(a) for a in argv),
		proc.stderr.decode(errors="replace")[-3000:], proc.stdout.decode(errors="replace")[-1000:]))


def _limit_core_dump():
	"""preexec_fn that keeps an aborting child from dumping core. On Linux
	RLIMIT_CORE 1 (not 0) also stops the kernel from handing the dump to a
	core_pattern pipe helper such as apport or systemd-coredump (core(5)),
	which can take about a second per abort; elsewhere the limit is 0. The
	hard limit is respected."""
	soft, hard = _resource.getrlimit(_resource.RLIMIT_CORE)
	limit = 1 if sys.platform.startswith("linux") else 0
	if hard != _resource.RLIM_INFINITY and hard < limit:
		limit = hard
	_resource.setrlimit(_resource.RLIMIT_CORE, (limit, hard))


def _core_limit_preexec():
	"""The preexec_fn for no_core_dump=True, or None where the resource
	module is not available (not POSIX): there the option does nothing."""
	if os.name != "posix" or _resource is None:
		return None
	return _limit_core_dump


def _execute(argv, cwd, threads, timeout, no_core_dump=False):
	return subprocess.run(argv, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
		env=subprocess_env(threads), timeout=timeout or DEFAULT_TIMEOUT,
		preexec_fn=_core_limit_preexec() if no_core_dump else None)


class ChromapIndex(object):
	"""A chromap index: path (the index file), genome (SyntheticGenome or
	None), fasta (reference FASTA passed as -r when mapping), k, w (as
	requested; None means chromap's default), args (extra index arguments),
	returncode/stderr of the build."""

	def __init__(self, path, genome, fasta, k=None, w=None, args=(), returncode=0, stderr=""):
		self.path = Path(path)
		self.genome = genome
		self.fasta = Path(fasta)
		self.k = k
		self.w = w
		self.args = list(args)
		self.returncode = returncode
		self.stderr = stderr

	def __repr__(self):
		return "ChromapIndex(%s, k=%s, w=%s)" % (self.path, self.k, self.w)


def build_index(genome, outdir=None, k=None, w=None, extra_args=(), fasta=None, name="index", check=True,
		timeout=None):
	"""Run `chromap -i -r FASTA -o INDEX [-k k] [-w w] extra_args`.

	genome: a SyntheticGenome that has been write()-n (its fasta_path is
	used), or None with `fasta` given. Returns a ChromapIndex; raises
	ChromapError on a non-zero exit when check is True."""
	exe = find_chromap()
	fasta = Path(fasta or genome.fasta_path).resolve()
	out = Path(outdir).resolve() if outdir is not None else fresh_dir("index")
	out.mkdir(parents=True, exist_ok=True)
	path = out / name
	argv = [str(exe), "-i", "-r", str(fasta), "-o", str(path)]
	if k is not None:
		argv += ["-k", str(k)]
	if w is not None:
		argv += ["-w", str(w)]
	args = split_args(extra_args)
	argv += args
	proc = _execute(argv, out, 1, timeout)
	if check and proc.returncode != 0:
		raise _failure(argv, proc)
	return ChromapIndex(path, genome, fasta, k, w, args, proc.returncode, proc.stderr.decode(errors="replace"))


def _check_threads(threads, args):
	for a in args:
		if a in ("-t", "--num-threads") or a.startswith("--num-threads="):
			raise ValueError("pass threads=, not -t/--num-threads, so the thread count stays pinned")
	threads = default_threads() if threads is None else int(threads)
	if not 1 <= threads <= MAX_THREADS:
		raise ValueError("threads must be between 1 and %d" % MAX_THREADS)
	return threads


def run_chromap(index, reads=None, extra_args=(), threads=None, name=None, outdir=None, output="out.txt",
		gz=False, fasta=False, name_style="plain", check=True, timeout=None, shared=False):
	"""Map reads with chromap in a fresh run directory and return a
	ChromapRun.

	index: ChromapIndex (gives -x and -r).
	reads: a ReadSet (written as FASTQ into the run directory and passed as
	-1/-2/-b); a dict {"r1": path, "r2": path, "bc": path} of existing files
	(values may be lists, joined with commas); or None when extra_args
	carries -1 itself.
	extra_args: list of strings or one shell-like string. chromap runs in the
	run directory, so relative output names such as "--summary summary.csv"
	land there.
	threads: -t; None means default_threads() ($CHROMAP_TEST_THREADS, else
	1). Never put -t in extra_args.
	output: the -o file name inside the run directory, an absolute path
	("/dev/stdout" sends the mappings to run.stdout), or None for no -o.
	gz/fasta/name_style: how a ReadSet is written (see ReadSet.write).
	check: raise ChromapError (an AssertionError) on a non-zero exit.
	shared: reuse the result of an identical earlier run in this session
	(same index, reads, arguments, threads and writing options). Shared runs
	must be treated as read-only.
	"""
	exe = find_chromap()
	args = split_args(extra_args)
	threads = _check_threads(threads, args)
	key = None
	if shared:
		if outdir is not None or not (reads is None or isinstance(reads, ReadSet)):
			raise ValueError("shared=True needs a ReadSet (or None) and no outdir")
		key = (str(exe), str(index.path), None if reads is None else reads.digest(), tuple(args), threads,
			output, gz, fasta, name_style)
		if key in _SHARED_RUNS:
			run = _SHARED_RUNS[key]
			if check and run.returncode != 0:
				raise _failure(run.argv, run._proc)
			return run
	out = Path(outdir).resolve() if outdir is not None else fresh_dir(name)
	out.mkdir(parents=True, exist_ok=True)
	argv = [str(exe), "-x", str(index.path), "-r", str(index.fasta)]
	readset = None
	if isinstance(reads, ReadSet):
		readset = reads
		files = reads.write(out / "input", gz=gz, fasta=fasta, name_style=name_style)
		argv += ["-1", str(files["r1"])]
		if "r2" in files:
			argv += ["-2", str(files["r2"])]
		if "bc" in files:
			argv += ["-b", str(files["bc"])]
	elif isinstance(reads, dict):
		for key_, flag in (("r1", "-1"), ("r2", "-2"), ("bc", "-b")):
			if reads.get(key_) is not None:
				v = reads[key_]
				v = ",".join(str(x) for x in v) if isinstance(v, (list, tuple)) else str(v)
				argv += [flag, v]
	elif reads is not None:
		raise TypeError("reads must be a ReadSet, a dict of paths or None")
	if output is not None:
		argv += ["-o", output if output.startswith("/") else str(out / output)]
	argv += ["-t", str(threads)] + args
	proc = _execute(argv, out, threads, timeout)
	run = ChromapRun(out, argv, proc, readset, index, output)
	if key is not None:
		_SHARED_RUNS[key] = run
	if check and proc.returncode != 0:
		raise _failure(argv, proc)
	return run


def run_cli(args, cwd=None, threads=None, timeout=None, no_core_dump=False):
	"""Run chromap with exactly `args` (list or string) in `cwd` (default a
	fresh directory) and return a ChromapRun; never raises on a non-zero exit.
	For help, version and error-path tests. OMP_NUM_THREADS is pinned to
	`threads` (default default_threads()). no_core_dump=True lowers the
	child's RLIMIT_CORE (see _limit_core_dump) for command lines expected to
	abort, so the system's crash handler does not add a second per abort; it
	is a no-op where the resource module is unavailable."""
	exe = find_chromap()
	args = split_args(args)
	threads = default_threads() if threads is None else int(threads)
	out = Path(cwd).resolve() if cwd is not None else fresh_dir("cli")
	argv = [str(exe)] + args
	proc = _execute(argv, out, threads, timeout, no_core_dump=no_core_dump)
	return ChromapRun(out, argv, proc, None, None, None)


_STAT_RE = re.compile(r"^(Number of [^:]+|# [^:]+): (\d+)\.?$")


class ChromapRun(object):
	"""Result of run_chromap or run_cli.

	outdir, argv, returncode, stdout, stderr (decoded text), stdout_bytes,
	reads (the ReadSet mapped, or None), index, output (Path of -o, or None),
	truth (read name -> SimRead/SimPair)."""

	def __init__(self, outdir, argv, proc, reads, index, output):
		self.outdir = Path(outdir)
		self.argv = argv
		self._proc = proc
		self.returncode = proc.returncode
		self.stdout_bytes = proc.stdout
		self.stdout = proc.stdout.decode(errors="replace")
		self.stderr = proc.stderr.decode(errors="replace")
		self.reads = reads
		self.index = index
		if output is None:
			self.output = None
		else:
			self.output = Path(output) if output.startswith("/") else self.outdir / output
		self.truth = reads.by_name if reads is not None else {}

	def path(self, name):
		return self.outdir / name

	def exists(self, name):
		return (self.outdir / name).exists()

	def text(self, name=None):
		p = self.output if name is None else self.path(name)
		return Path(p).read_text()

	def lines(self, name=None):
		return self.text(name).splitlines()

	def bed(self, name=None):
		return parse_bed(self.output if name is None else self.path(name))

	def fragments(self, name=None):
		return parse_fragments(self.output if name is None else self.path(name))

	def tagalign(self, name=None):
		return parse_bed(self.output if name is None else self.path(name), allow_six=True)

	def sam(self, name=None):
		return parse_sam(self.output if name is None else self.path(name))

	def pairs(self, name=None):
		return parse_pairs(self.output if name is None else self.path(name))

	def paf(self, name=None):
		return parse_paf(self.output if name is None else self.path(name))

	def summary(self, name="summary.csv"):
		return parse_summary(self.path(name))

	def stats(self):
		"""Integer statistics chromap prints on stderr, e.g.
		{"Number of reads": 400, "Number of mapped reads": 398, ...}."""
		return parse_stats(self.stderr)

	def __repr__(self):
		return "ChromapRun(%s, exit %s)" % (self.outdir, self.returncode)


# ------------------------------------------------------------ parsers

BedRecord = namedtuple("BedRecord", "chrom start end name mapq strand dups")
FragmentRecord = namedtuple("FragmentRecord", "chrom start end barcode dups")
PairsRecord = namedtuple("PairsRecord", "read_id chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2")
PafRecord = namedtuple("PafRecord", "qname qlen qstart qend strand tname tlen tstart tend nmatch alen mapq")


def _rows(path):
	with open(str(path)) as f:
		return [line.rstrip("\n").split("\t") for line in f if line.strip()]


def parse_bed(path, allow_six=False):
	"""Records of chromap's 7-column BED/TagAlign lines (chrom, start, end,
	"N", mapq, strand, duplicate count). With allow_six, 6-column lines
	(TagAlign lines without a count) are accepted with dups None."""
	out = []
	for f in _rows(path):
		if len(f) == 7:
			out.append(BedRecord(f[0], int(f[1]), int(f[2]), f[3], int(f[4]), f[5], int(f[6])))
		elif len(f) == 6 and allow_six:
			out.append(BedRecord(f[0], int(f[1]), int(f[2]), f[3], int(f[4]), f[5], None))
		else:
			raise ValueError("%s: unexpected BED line %r" % (path, "\t".join(f)))
	return out


def parse_fragments(path):
	"""Records of the 5-column barcoded BED (fragments) format: chrom, start,
	end, barcode, duplicate count."""
	out = []
	for f in _rows(path):
		if len(f) != 5:
			raise ValueError("%s: unexpected fragment line %r" % (path, "\t".join(f)))
		out.append(FragmentRecord(f[0], int(f[1]), int(f[2]), f[3], int(f[4])))
	return out


class SamRecord(object):
	"""One SAM line: qname flag rname pos mapq cigar rnext pnext tlen seq qual
	and tags (OrderedDict name -> value; i-typed values are ints)."""

	__slots__ = ("qname", "flag", "rname", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq", "qual",
		"tags", "tag_types")

	def __init__(self, fields):
		(self.qname, flag, self.rname, pos, mapq, self.cigar, self.rnext, pnext, tlen, self.seq,
			self.qual) = fields[:11]
		self.flag, self.pos, self.mapq, self.pnext, self.tlen = int(flag), int(pos), int(mapq), int(pnext), int(tlen)
		self.tags = OrderedDict()
		self.tag_types = OrderedDict()
		for t in fields[11:]:
			name, typ, val = t.split(":", 2)
			self.tags[name] = int(val) if typ == "i" else val
			self.tag_types[name] = typ

	@property
	def is_reverse(self):
		return bool(self.flag & 16)

	@property
	def is_read1(self):
		return bool(self.flag & 64)

	@property
	def is_read2(self):
		return bool(self.flag & 128)

	@property
	def is_secondary(self):
		return bool(self.flag & 256)

	def __repr__(self):
		return "SamRecord(%s %d %s:%d %s)" % (self.qname, self.flag, self.rname, self.pos, self.cigar)


def parse_sam(path):
	"""(header lines, [SamRecord]) of a SAM file."""
	header, recs = [], []
	with open(str(path)) as f:
		for line in f:
			line = line.rstrip("\n")
			if not line:
				continue
			if line.startswith("@"):
				header.append(line)
			else:
				recs.append(SamRecord(line.split("\t")))
	return header, recs


def parse_pairs(path):
	"""(header lines, [PairsRecord]) of a 4DN pairs file; positions are
	1-based as written."""
	header, recs = [], []
	with open(str(path)) as f:
		for line in f:
			line = line.rstrip("\n")
			if not line:
				continue
			if line.startswith("#"):
				header.append(line)
				continue
			x = line.split("\t")
			if len(x) != 10:
				raise ValueError("%s: unexpected pairs line %r" % (path, line))
			recs.append(PairsRecord(x[0], x[1], int(x[2]), x[3], int(x[4]), x[5], x[6], x[7], int(x[8]), int(x[9])))
	return header, recs


def parse_paf(path):
	"""[PafRecord] of chromap's 12-column PAF lines."""
	out = []
	for x in _rows(path):
		if len(x) != 12:
			raise ValueError("%s: unexpected PAF line %r" % (path, "\t".join(x)))
		out.append(PafRecord(x[0], int(x[1]), int(x[2]), int(x[3]), x[4], x[5], int(x[6]), int(x[7]), int(x[8]),
			int(x[9]), int(x[10]), int(x[11])))
	return out


class Summary(object):
	"""The --summary CSV: header (column names) and rows (OrderedDict barcode
	-> dict column -> value, ints for counts and floats for fric/estfrip). The
	bulk row's barcode is the empty string; "non-whitelist" is its own row."""

	def __init__(self, header, rows, lines):
		self.header = header
		self.rows = rows
		self.lines = lines

	def __getitem__(self, barcode):
		return self.rows[barcode]

	def __contains__(self, barcode):
		return barcode in self.rows


def parse_summary(path):
	with open(str(path)) as f:
		lines = [line.rstrip("\n") for line in f if line.strip()]
	header = lines[0].split(",")
	rows = OrderedDict()
	for line in lines[1:]:
		vals = line.split(",")
		if len(vals) != len(header):
			raise ValueError("%s: unexpected summary line %r" % (path, line))
		rec = {}
		for k, v in zip(header[1:], vals[1:]):
			rec[k] = float(v) if k in ("fric", "estfrip") else int(v)
		if vals[0] in rows:
			raise ValueError("%s: barcode %r twice" % (path, vals[0]))
		rows[vals[0]] = rec
	return Summary(header, rows, lines)


def parse_stats(text):
	"""{label: int} from chromap's stderr lines such as "Number of reads:
	400." and "Number of output mappings (passed filters): 200"."""
	out = OrderedDict()
	for line in text.splitlines():
		m = _STAT_RE.match(line.strip())
		if m:
			out[m.group(1)] = int(m.group(2))
	return out


def first_line_diff(a, b, limit=5):
	"""Readable description of where two lists of lines first differ."""
	for i, (x, y) in enumerate(zip(a, b)):
		if x != y:
			lo = max(0, i - 1)
			return "line %d differs:\n  %s\n  %s\n(context a: %s)" % (i + 1, x, y, a[lo:i + limit])
	if len(a) != len(b):
		return "lengths differ: %d vs %d" % (len(a), len(b))
	return "identical"
