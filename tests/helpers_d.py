"""Helpers for the area D tests: cell barcodes, whitelist correction,
--read-format, --barcode-translate, the pairs output, --chr-order and the
command-line front end.

Everything here is an independent reference: barcode correction is written
from the rule chromap documents (README "Map scATAC-seq reads", chromap.1
--bc-error-threshold and --bc-probability-threshold), the pairs records from
the 4DN pairs specification, the read-format extraction from the README's
[r1|r2|bc]:start:end:strand grammar. Only the standard library is used.
"""

import collections
import random

import chromap_helpers as ch

QUAL_OFFSET = 33
# chromap clamps each barcode base quality to [3, 40] before turning it into
# an error probability (src/chromap.cc, CorrectBarcodeAt).
MIN_BARCODE_QUAL = 3
MAX_BARCODE_QUAL = 40

UNCORRECTABLE = "uncorrectable"
EXACT = "exact"
CORRECTED = "corrected"


def error_probability(qual_char):
	"""Phred+33 error probability of one base, quality clamped to [3, 40]."""
	q = ord(qual_char) - QUAL_OFFSET
	q = min(MAX_BARCODE_QUAL, max(MIN_BARCODE_QUAL, q))
	return 10.0 ** (-q / 10.0)


def whitelist_abundance(barcodes, whitelist):
	"""Counter whitelist barcode -> number of input barcodes equal to it.
	Only N-free exact matches count: this is the abundance chromap estimates
	from the barcode files before mapping (README: "estimate barcode
	abundance and use this information to perform barcode correction")."""
	wl = set(whitelist)
	return collections.Counter(b for b in barcodes if "N" not in b and b in wl)


def correction_candidates(barcode, qual, whitelist, abundance, max_errors):
	"""Every whitelist barcode within Hamming distance 1..max_errors of
	`barcode` (an N matches no base), as (sort key, whitelist barcode). The
	score of a candidate is the product of the error probabilities of the
	mismatched bases times the candidate's abundance (its share of all exact
	whitelist matches). The sort key is chromap's BarcodeWithQual order
	(src/utils.h): score, then the larger and the smaller mismatched
	position with the whitelist bases there (0 and NUL for a single
	mismatch)."""
	total = float(sum(abundance.values()))
	out = []
	for w in whitelist:
		diff = [i for i in range(len(barcode)) if barcode[i] != w[i]]
		if not 1 <= len(diff) <= max_errors:
			continue
		p = 1.0
		for i in diff:
			p *= error_probability(qual[i])
		ab = abundance.get(w, 0) / total if total > 0 else float("nan")
		score = p * ab
		if len(diff) == 1:
			key = (score, diff[0], w[diff[0]], 0, "\0")
		else:
			hi, lo = max(diff), min(diff)
			key = (score, hi, w[hi], lo, w[lo])
		out.append((key, w))
	return out


def correct_barcode(barcode, qual, whitelist, abundance, max_errors=1, min_probability=0.9):
	"""(status, barcode) for one cell barcode under chromap's documented
	rule. status is EXACT (the barcode is in the whitelist), CORRECTED (one
	candidate within max_errors, or several and the best one's share of the
	summed scores is above min_probability) or UNCORRECTABLE (more Ns than
	max_errors, no candidate, or no candidate probable enough); the barcode
	returned is the whitelist barcode, or None when uncorrectable."""
	wl = set(whitelist)
	n_count = barcode.count("N")
	if n_count > max_errors:
		return UNCORRECTABLE, None
	if n_count == 0 and barcode in wl:
		return EXACT, barcode
	if max_errors == 0:
		return UNCORRECTABLE, None
	cands = correction_candidates(barcode, qual, whitelist, abundance, max_errors)
	if not cands:
		return UNCORRECTABLE, None
	if len(cands) == 1:
		return CORRECTED, cands[0][1]
	cands.sort(key=lambda kw: kw[0], reverse=True)
	total = 0.0
	for key, _ in cands:
		total += key[0]
	best = cands[0][0][0]
	if total > 0 and best / total > min_probability:
		return CORRECTED, cands[0][1]
	return UNCORRECTABLE, None


def seed_sequence(barcode):
	"""The barcode as chromap writes an uncorrected one: its 2-bit key
	decoded, so N becomes A (chromap_helpers.encode_seed)."""
	return ch.decode_seed(ch.encode_seed(barcode), len(barcode))


def expected_fragments(records, whitelist, max_errors=1, min_probability=0.9, keep_uncorrectable=False):
	"""Counter of (contig, start, end, barcode) -> number of barcoded records
	(SimRead: the read's interval; SimPair: the fragment) chromap should
	write with a whitelist: each line once per record without deduplication,
	or one line with that count as its duplicate count with cell-level
	deduplication."""
	recs = list(records)
	ab = whitelist_abundance([r.barcode for r in recs], whitelist)
	out = collections.Counter()
	for r in recs:
		status, bc = correct_barcode(r.barcode, r.barcode_qual or "I" * len(r.barcode), whitelist, ab,
			max_errors, min_probability)
		if status == UNCORRECTABLE:
			if not keep_uncorrectable:
				continue
			bc = seed_sequence(r.barcode)
		out[(r.contig, r.start, r.end, bc)] += 1
	return out


def correction_outcomes(records, whitelist, max_errors=1, min_probability=0.9):
	"""{read name: (status, barcode)} for barcoded records."""
	recs = list(records)
	ab = whitelist_abundance([r.barcode for r in recs], whitelist)
	return dict((r.name, correct_barcode(r.barcode, r.barcode_qual or "I" * len(r.barcode), whitelist, ab,
		max_errors, min_probability)) for r in recs)


def fragments_by_locus(run):
	"""{(chrom, start, end): [barcode, ...]} of a barcoded fragment output."""
	out = collections.defaultdict(list)
	for f in run.fragments():
		out[(f.chrom, f.start, f.end)].append(f.barcode)
	return out


# ------------------------------------------------------------ barcodes

def mutate(barcode, positions):
	"""`barcode` with a different base at each position (chromap_helpers
	.substitute: A->C->G->T->A)."""
	b = list(barcode)
	for p in positions:
		b[p] = ch.substitute(b[p])
	return "".join(b)


def with_n(barcode, positions):
	b = list(barcode)
	for p in positions:
		b[p] = "N"
	return "".join(b)


def isolated_barcodes(n, length, rng, avoid, min_dist):
	"""n random barcodes at least `min_dist` from each other and from every
	barcode in `avoid`."""
	out = []
	pool = list(avoid)
	while len(out) < n:
		b = "".join(rng.choice("ACGT") for _ in range(length))
		if all(ch.hamming(b, x) >= min_dist for x in pool):
			out.append(b)
			pool.append(b)
	return out


class LocusPicker(object):
	"""Draws unique loci with distinct (contig, start) and (contig, end), so
	every read or fragment of a run has its own output line key."""

	def __init__(self, sim):
		self.sim = sim
		self.starts = set()
		self.ends = set()

	def take(self, length, contigs=None):
		for _ in range(1000):
			c, s = self.sim.unique_locus(length, contigs=contigs)
			if (c, s) not in self.starts and (c, s + length) not in self.ends:
				self.starts.add((c, s))
				self.ends.add((c, s + length))
				return c, s
		raise RuntimeError("no free locus")


def barcoded_read(sim, picker, barcode, barcode_qual=None, strand=None):
	"""A 50 bp single-end read at a fresh unique locus carrying `barcode`."""
	c, s = picker.take(50)
	return sim.se_read(c, s, 50, strand or sim.rng.choice("+-"), barcode=barcode, barcode_qual=barcode_qual)


# ------------------------------------------------------------ pairs

def pairs_position(read):
	"""1-based position of a read's 5' end, as the 4DN pairs format records
	it: the first aligned base of a + strand read, the last of a - strand
	read."""
	return read.start + 1 if read.strand == "+" else read.end


def expected_pairs_record(pair, rank):
	"""(readID, chrom1, pos1, chrom2, pos2, strand1, strand2) of a pair in
	the 4DN "upper triangle" shape: the two sides are ordered by chromosome
	rank (`rank`: contig -> rank in the flipping order), then by position on
	the same chromosome. readID is read 1's name."""
	a = (pair.r1.contig, pairs_position(pair.r1), pair.r1.strand)
	b = (pair.r2.contig, pairs_position(pair.r2), pair.r2.strand)
	if (rank[a[0]], a[1]) > (rank[b[0]], b[1]):
		a, b = b, a
	return (pair.name, a[0], a[1], b[0], b[1], a[2], b[2])


def full_order(listed, fasta_order):
	"""The reference order --chr-order gives: the listed contigs first, the
	others after them in FASTA order (chromap.1, --chr-order; src/chromap.cc
	GenerateCustomRidRanks)."""
	return list(listed) + [c for c in fasta_order if c not in listed]


def hic_pair(sim, name, picker, kind, strands):
	"""A read pair as Hi-C produces it: read 1 and read 2 (50 bp each, exact)
	from two unique loci. kind "trans": two different contigs; "far": the
	same contig (one of 15 kb or more), more than 3 kb apart (beyond any
	insert size). strands:
	(read 1 strand, read 2 strand)."""
	long_contigs = [c for c in sim.g.names if len(sim.g.contigs[c]) >= 15000]
	ca, sa = picker.take(50, contigs=None if kind == "trans" else long_contigs)
	for _ in range(1000):
		if kind == "trans":
			cb, sb = picker.take(50)
			if cb != ca:
				break
		else:
			cb, sb = picker.take(50, contigs=[ca])
			if abs(sb - sa) > 3000:
				break
	else:
		raise RuntimeError("no partner locus")
	r1 = sim.se_read(ca, sa, 50, strands[0], name=name)
	r2 = sim.se_read(cb, sb, 50, strands[1], name=name)
	return ch.SimPair(name, r1, r2, ca, None, None, strands[0], kind="hic_" + kind)


def hic_reads(sim, n_fragment, n_trans, n_far):
	"""A ReadSet of Hi-C style pairs: ordinary fragments (insert 200-400)
	and pairs whose mates are on different contigs or far apart on one
	contig, every strand combination, mates in both orders."""
	picker = LocusPicker(sim)
	recs = []
	combos = [("+", "+"), ("+", "-"), ("-", "+"), ("-", "-")]
	for i in range(n_fragment):
		ins = sim.rng.randint(200, 400)
		c, s = picker.take(ins)
		recs.append(sim.pair(c, s, ins, strand="+-"[i % 2], name="frag%04d" % i))
	for i in range(n_trans):
		recs.append(hic_pair(sim, "trans%04d" % i, picker, "trans", combos[i % 4]))
	for i in range(n_far):
		recs.append(hic_pair(sim, "far%04d" % i, picker, "far", combos[i % 4]))
	return sim.shuffled(ch.ReadSet(recs))


# ------------------------------------------------------------ read format

def apply_read_format(seq, qual, fields):
	"""Extraction by the --read-format grammar [r1|r2|bc]:start:end[:strand]
	(README): the inclusive ranges of `fields` (each "xx:start:end" or
	"xx:start:end:strand", end -1 meaning the read's last base) are
	concatenated in order; with strand "-" the result is reverse-complemented
	(and the quality string reversed). The last strand given applies."""
	parts_s, parts_q, strand = [], [], "+"
	for f in fields:
		x = f.split(":")
		start, end = int(x[1]), int(x[2])
		if end == -1:
			end = len(seq) - 1
		parts_s.append(seq[start:end + 1])
		parts_q.append(qual[start:end + 1])
		if len(x) > 3:
			strand = x[3]
	s, q = "".join(parts_s), "".join(parts_q)
	if strand == "-":
		s, q = ch.revcomp(s), q[::-1]
	return s, q


def random_quality(rng, n, alphabet="#+5?I"):
	return "".join(rng.choice(alphabet) for _ in range(n))


def stable_rng(*parts):
	return random.Random(ch.stable_seed("helpers_d", *parts))


# ------------------------------------------------------------ aborting runs

def run_cli_expect_abort(args, cwd=None):
	"""chromap_helpers.run_cli for command lines the option parser rejects by
	throwing (chromap does not catch it and aborts): the same ChromapRun, run
	without a core dump (run_cli's no_core_dump)."""
	return ch.run_cli(args, cwd=cwd, no_core_dump=True)
