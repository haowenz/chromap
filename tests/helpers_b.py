"""Reference implementations and read designers for the area B tests
(mapping core: candidates, draft mappings, alignment, MAPQ, pairing).

Everything here is computed from definitions written out in the docstrings,
never from chromap's output: a banded edit distance (the draft mapping's
semi-global distance), an affine-gap semi-global aligner with chromap's
scoring matrix (the SAM CIGAR), the MAPQ rules of src/mapping_generator.h
written out as arithmetic, the repetitive-seed-length rule of
Index::UpdateRepetitiveSeedStats, and minimizer bookkeeping on top of
chromap_helpers.reference_minimizers. Python 3.8, standard library only.
"""

import collections
import math
import re

import chromap_helpers as ch

K = ch.DEFAULT_KMER
W = ch.DEFAULT_WINDOW
DEFAULT_E = 8
INF = float("inf")


def nt4(c):
	"""chromap's base code (CharToUint8): ACGT/acgt -> 0-3, anything else 4."""
	return ch.NT4.get(c, 4)


# ------------------------------------------------------------ edit distance

def banded_edit_distance(window, read, e):
	"""Semi-global edit distance of the whole `read` against `window`
	(len(read) + 2e bases) as chromap's draft mapping defines it: the read
	may start at any of the first 2e + 1 window positions for free, the
	alignment stays on window diagonals 0..2e (window index minus read index),
	and it may end at any window index L - 1 .. L - 1 + 2e. Bases compare by
	their 2-bit code, so N matches N and lower case matches upper case.

	Returns (errors, end) where end is the window index of the last aligned
	window base. Among equal minima the end is the first one scanning from
	L - 1 upwards, except that the centre L - 1 + e wins a tie with the
	minimum found before it (chromap prefers the unshifted diagonal)."""
	L = len(read)
	B = 2 * e + 1
	if len(window) != L + 2 * e:
		raise ValueError("window must be len(read) + 2e")
	wcodes = [nt4(c) for c in window]
	# prev[d]: distance of read[:i] against a window substring ending before
	# window index i + d (d = diagonal, window prefix length minus read prefix
	# length, 0..2e). Row 0 is free: the read may start anywhere in the band.
	prev = [0] * B
	for i in range(1, L + 1):
		rc = nt4(read[i - 1])
		cur = [INF] * B
		for d in range(B):
			# diagonal: D(i-1, i-1+d) + mismatch of read[i-1], window[i-1+d]
			best = prev[d] + (0 if wcodes[i - 1 + d] == rc else 1)
			# vertical: D(i-1, i+d), a read base without window base
			if d + 1 < B and prev[d + 1] + 1 < best:
				best = prev[d + 1] + 1
			# horizontal: D(i, i+d-1), a window base without read base
			if d > 0 and cur[d - 1] + 1 < best:
				best = cur[d - 1] + 1
			cur[d] = best
		prev = cur
	vals = prev
	best, pos = vals[0], 0
	for t in range(1, 2 * e + 1):
		if vals[t] < best or (vals[t] == best and t == e):
			best, pos = vals[t], t
	return best, L - 1 + pos


def draft_mapping(ref, read, strand, candidate, e):
	"""Draft mapping of `read` (as sequenced) at a candidate on contig
	sequence `ref`. `candidate` is the candidate position chromap derives
	from seeds: the start of the read's alignment for strand "+", the last
	aligned reference base for strand "-". Returns (errors, end) with `end`
	the reference position of the last aligned base, or None when the
	candidate's window does not fit on the contig."""
	L = len(read)
	text = read if strand == "+" else ch.revcomp(read)
	pos = candidate if strand == "+" else candidate - L + 1
	if pos < e or pos + L + e >= len(ref):
		return None
	ws = pos - e
	errors, end = banded_edit_distance(ref[ws:ws + L + 2 * e], text, e)
	return errors, ws + end


# --------------------------------------------------- affine-gap alignment

ScoreParams = collections.namedtuple("ScoreParams", "match mismatch o_ins e_ins o_del e_del")
# chromap's defaults (-A 1 -B 4 -O 6,6 -E 1,1); o_ins/e_ins price CIGAR I
# (read bases without reference), o_del/e_del CIGAR D.
DEFAULT_SCORES = ScoreParams(1, 4, 6, 1, 6, 1)


def substitution_score(a, b, p):
	x, y = nt4(a), nt4(b)
	if x == 4 or y == 4:
		return 0
	return p.match if x == y else -p.mismatch


def semi_global(window, read, p, w):
	"""Best affine-gap alignment of the whole `read` against `window` with
	chromap's matrix (match +A, mismatch -B, N 0): a gap of l read bases
	(CIGAR I) costs o_ins + l * e_ins, a gap of l window bases (CIGAR D) costs
	o_del + l * e_del, and gaps open only after an aligned (M) column, so I
	and D are never adjacent. The read's first base aligns to a window index
	<= w (free leading window bases), the alignment ends anywhere (free
	trailing window bases), and every aligned cell (read i, window j) keeps
	i <= j <= i + w.

	Returns dict(score, count, alignments): `count` optimal alignments, and
	when it is 1, alignments[0] = (start, end, cigar) with start/end the
	window interval [start, end) and cigar a string."""
	L, Q = len(read), len(window)
	NEG = -10 ** 9
	NONE = (NEG, 0)
	rc = [nt4(c) for c in read]
	wc = [nt4(c) for c in window]
	oi, ei = p.o_ins + p.e_ins, p.e_ins
	od, ed = p.o_del + p.e_del, p.e_del
	# Rows of (score, number of optimal paths) per diagonal d = j - i in
	# [0, w]: M (read[i] aligned to window[j]), I (read[i] unaligned after
	# window[j]), D (window[j] unaligned after read[i]).
	Ms, Is, Ds = [], [], []
	for i in range(L):
		M = [NONE] * (w + 1)
		I = [NONE] * (w + 1)
		D = [NONE] * (w + 1)
		for d in range(w + 1):
			j = i + d
			if j >= Q:
				break
			a, b = rc[i], wc[j]
			s = 0 if a == 4 or b == 4 else (p.match if a == b else -p.mismatch)
			if i == 0:
				M[d] = (s, 1)
			else:
				x, y, z = Ms[i - 1][d], Is[i - 1][d], Ds[i - 1][d]
				best = max(x[0], y[0], z[0])
				n = (x[1] if x[0] == best else 0) + (y[1] if y[0] == best else 0) + (z[1] if z[0] == best else 0)
				M[d] = (best + s, n) if n else NONE
				if d + 1 <= w:
					x, y = Ms[i - 1][d + 1], Is[i - 1][d + 1]
					o, c = x[0] - oi, y[0] - ei
					best = max(o, c)
					n = (x[1] if o == best else 0) + (y[1] if c == best else 0)
					I[d] = (best, n) if n else NONE
			if d > 0:
				x, y = M[d - 1], D[d - 1]
				o, c = x[0] - od, y[0] - ed
				best = max(o, c)
				n = (x[1] if o == best else 0) + (y[1] if c == best else 0)
				D[d] = (best, n) if n else NONE
		Ms.append(M)
		Is.append(I)
		Ds.append(D)
	ends = []
	for d in range(w + 1):
		if L - 1 + d >= Q:
			break
		for state, row in (("M", Ms[L - 1]), ("I", Is[L - 1])):
			v = row[d]
			if v[1]:
				ends.append((v[0], v[1], d, state))
	score = max(x[0] for x in ends)
	count = sum(x[1] for x in ends if x[0] == score)
	out = dict(score=score, count=count, alignments=[])
	if count == 1:
		_, _, d, state = [x for x in ends if x[0] == score][0]
		out["alignments"].append(_traceback(rc, wc, p, Ms, Is, Ds, L - 1, d, state))
	return out


def _traceback(rc, wc, p, Ms, Is, Ds, i, d, state):
	"""Walk back the single optimal path from (read i, diagonal d)."""
	ops = []
	end = i + d + 1
	tables = {"M": Ms, "I": Is, "D": Ds}
	while True:
		score = tables[state][i][d][0]
		j = i + d
		if state == "M":
			ops.append("M")
			if i == 0:
				break
			a, b = rc[i], wc[j]
			s = 0 if a == 4 or b == 4 else (p.match if a == b else -p.mismatch)
			state = [st for st in ("M", "I", "D") if tables[st][i - 1][d][1] and tables[st][i - 1][d][0] + s == score][0]
			i -= 1
		elif state == "I":
			ops.append("I")
			state = "M" if Ms[i - 1][d + 1][0] - p.o_ins - p.e_ins == score else "I"
			i, d = i - 1, d + 1
		else:
			ops.append("D")
			state = "M" if Ms[i][d - 1][0] - p.o_del - p.e_del == score else "D"
			d -= 1
	ops.reverse()
	return (i + d, end, compress_cigar(ops))


def compress_cigar(ops):
	out = []
	for op in ops:
		if out and out[-1][1] == op:
			out[-1][0] += 1
		else:
			out.append([1, op])
	return "".join("%d%s" % (n, op) for n, op in out)


def parse_cigar(cigar):
	return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]


def cigar_spans(cigar):
	"""(read bases, reference bases) a CIGAR consumes."""
	q = r = 0
	for n, op in parse_cigar(cigar):
		if op in "MI=XS":
			q += n
		if op in "MDN=X":
			r += n
	return q, r


def cigar_score(window, start, read, cigar, p):
	"""Score of the alignment `cigar` of `read` placed at window[start]."""
	i, j, score = 0, start, 0
	for n, op in parse_cigar(cigar):
		if op == "M":
			for _ in range(n):
				score += substitution_score(read[i], window[j], p)
				i += 1
				j += 1
		elif op == "I":
			score -= p.o_ins + n * p.e_ins
			i += n
		elif op == "D":
			score -= p.o_del + n * p.e_del
			j += n
		else:
			raise ValueError("unexpected CIGAR op %s" % op)
	return score


def nm_md(ref_seq, read_seq, cigar):
	"""NM and MD (SAM specification) of `read_seq` aligned by `cigar` to
	`ref_seq` starting at its first base."""
	nm, md, run = 0, "", 0
	i = j = 0
	for n, op in parse_cigar(cigar):
		if op == "M":
			for _ in range(n):
				if nt4(ref_seq[j]) == nt4(read_seq[i]):
					run += 1
				else:
					nm += 1
					md += "%d%s" % (run, ref_seq[j])
					run = 0
				i += 1
				j += 1
		elif op == "I":
			nm += n
			i += n
		elif op == "D":
			nm += n
			md += "%d^%s" % (run, ref_seq[j:j + n])
			run = 0
			j += n
	return nm, md + str(run)


# ------------------------------------------------------------------- MAPQ

def mapq_single_end(num_errors, alignment_length, read_length, max_diff, second_min_errors, num_best,
		num_second_best, repetitive_length):
	"""MAPQ of one read's best mapping (GetMAPQForSingleEndRead, non-split):

	length = max(alignment length, read length); identity = 1 - errors /
	length. With one best mapping, the error gap d = min(second-best errors,
	errors + max_diff) - errors and
	MAPQ = int(30.1 * d * t^2 + 0.499) with t = (int(ln 50) / ln length, or 1
	below 50 bp) * identity^2. Second-best mappings subtract
	int(4.343 * ln(n + 1) + 0.499); the value is clipped to [0, 60]. Several
	best mappings give 0. A repetitive seed length r scales MAPQ by
	(1 - sqrt(f)), (1 - f), (1 - f^2) or (1 - f^4) (f = r / read length,
	0.999 when r >= read length) for identity <= 0.95, <= 0.97, < 0.999 and
	>= 0.999, rounding with + 0.499. max_diff is -e for single-end reads and 2
	for each mate of a pair."""
	length = max(alignment_length, read_length)
	identity = 1 - float(num_errors) / length
	mapq = 0
	if num_best <= 1:
		second = min(second_min_errors, num_errors + max_diff)
		tmp = 1.0 if length < 50 else int(math.log(50)) / math.log(length)
		tmp = tmp * (identity * identity)
		mapq = int(5 * 6.02 * (second - num_errors) * tmp * tmp + 0.499)
	if num_second_best > 0:
		mapq -= int(4.343 * math.log(num_second_best + 1) + 0.499)
	mapq = max(0, min(60, mapq))
	if repetitive_length > 0:
		frac = 0.999 if repetitive_length >= read_length else repetitive_length / float(read_length)
		mapq = _scale_repetitive(mapq, frac, identity)
	return mapq


def _scale_repetitive(mapq, frac, identity):
	if identity <= 0.95:
		return int(mapq * (1 - math.sqrt(frac)) + 0.499)
	if identity <= 0.97:
		return int(mapq * (1 - frac) + 0.499)
	if identity >= 0.999:
		return int(mapq * (1 - frac * frac * frac * frac) + 0.499)
	return int(mapq * (1 - frac * frac) + 0.499)


MateStats = collections.namedtuple("MateStats",
	"errors alignment_length read_length min_errors second_min_errors num_best num_second_best repetitive_length")


def mapq_paired_end(m1, m2, min_sum, second_min_sum, num_best_pairs, num_second_best_pairs, e=DEFAULT_E,
		force_mapq=-1):
	"""MAPQ of a pair (GetMAPQForPairedEndRead, non-split). The pair score
	q = min(30.1 * d + 0.499, 60) with d = min(second-best pair error sum,
	min errors 1 + min errors 2 + 3) - best pair error sum (0 when the pair
	has several best placements), minus int(4.343 ln(n + 1) + 0.499) for n
	second-best pairs, scaled for the mates' summed repetitive seed length as
	in mapq_single_end. Each mate's single-end MAPQ s (max_diff 2) becomes s
	if s > q, else min(q, s + 0.65 q); then times 1.2, capped at 60. The pair
	gets the smaller mate value (all conversions truncate to an unsigned
	byte)."""
	mapq_pe = 0
	if num_best_pairs <= 1:
		adjusted = min(second_min_sum, m1.min_errors + m2.min_errors + 3)
		mapq_pe = int(5 * 6.02 * (adjusted - min_sum) / 1 + .499) % 256
		if num_second_best_pairs > 0:
			mapq_pe = (mapq_pe - int(4.343 * math.log(num_second_best_pairs + 1) + 0.499)) % 256
		mapq_pe = min(mapq_pe, 60)
		rep = m1.repetitive_length + m2.repetitive_length
		if rep > 0:
			total = float(m1.read_length + m2.read_length)
			frac = 0.999 if rep >= total else rep / total
			id1 = 1 - float(m1.errors) / max(m1.read_length, m1.alignment_length)
			id2 = 1 - float(m2.errors) / max(m2.read_length, m2.alignment_length)
			mapq_pe = _scale_repetitive(mapq_pe, frac, min(id1, id2)) % 256
	out = []
	for m in (m1, m2):
		s = mapq_single_end(m.errors, m.alignment_length, m.read_length, 2, m.second_min_errors, m.num_best,
			m.num_second_best, m.repetitive_length)
		if s > mapq_pe:
			v = s
		elif mapq_pe < s + mapq_pe * 0.65:
			v = mapq_pe
		else:
			v = int(s + mapq_pe * 0.65)
		v = int(v * 1.2) % 256
		out.append(min(v, 60))
	mapq = min(out)
	if mapq < 60 and 0 <= force_mapq < mapq:
		mapq = force_mapq
	return mapq


def repetitive_seed_length(positions, k=K, w=W):
	"""Index::UpdateRepetitiveSeedStats over the read end positions of the
	repetitive minimizers, in read order: the first adds k; each next one adds
	its distance to the previous one when that is below k + w - 1, else k."""
	total, prev = 0, None
	for pos in positions:
		if prev is None or prev > pos:
			total += k
		elif pos < prev + k + w - 1:
			total += pos - prev
		else:
			total += k
		prev = pos
	return total


# ------------------------------------------------------------- minimizers

def read_minimizers(seq):
	"""[(hash, end position, strand)] of a read as chromap computes them
	(N-free reads, k 17, w 7)."""
	return [(h, (hit >> 1) & 0xffffffff, hit & 1) for h, hit in ch.reference_minimizers(seq, K, W)]


def locus_minimizers(genome, contig, start, end, pad=40):
	"""{(hash, reference end position)} of the reference minimizers whose
	k-mer lies in contig[start:end]; computed on a padded stretch so windows
	reaching past [start, end) are included."""
	a = max(0, start - pad)
	b = min(len(genome.contigs[contig]), end + pad)
	out = set()
	for h, hit in ch.reference_minimizers(genome.seq(contig, a, b), K, W):
		pos = a + ((hit >> 1) & 0xffffffff)
		if start + K - 1 <= pos < end:
			out.add((h, pos))
	return out


def seed_support(genome, read):
	"""(number of read minimizers, number that hit the read's true locus) for
	an indel-free SimRead: a read minimizer hits when its k-mer has no
	substitution and the same hash is a reference minimizer at the matching
	reference position."""
	mins = read_minimizers(read.seq)
	ref = locus_minimizers(genome, read.contig, read.start, read.end)
	L = len(read.seq) - len(read.tail)
	bad = set()
	for g in read.mismatches:
		bad.add(g - read.start if read.strand == "+" else read.start + L - 1 - g)
	hits = 0
	for h, p, _ in mins:
		if p >= L or any(p - K + 1 <= b <= p for b in bad):
			continue
		ref_end = read.start + p if read.strand == "+" else read.start + L - 1 - p + K - 1
		if (h, ref_end) in ref:
			hits += 1
	return len(mins), hits


def segment_minimizer_counts(genome, segment, pad=40):
	"""Counter hash -> number of copies of a planted segment in which the
	hash is a reference minimizer."""
	counts = collections.Counter()
	for c in genome.segments[segment].copies:
		hashes = {h for h, _ in locus_minimizers(genome, c.contig, c.start, c.end, pad=pad)}
		counts.update(hashes)
	return counts


# ----------------------------------------------------------- read design

def indel_read(genome, contig, start, length, strand, kind, offset, size, name):
	"""A read with one insertion or deletion of `size` bases after read
	offset `offset` (counted on the forward strand of the locus). Deletion:
	the read skips reference bases [start + offset, start + offset + size)
	and covers length + size reference bases. Insertion: `size` bases not in
	the reference follow offset - 1; the read covers length - size reference
	bases. Returns a SimRead whose [start, end) is the reference span and
	extra["cigar"] the CIGAR in reference orientation."""
	if kind == "D":
		ref = genome.seq(contig, start, start + length + size)
		fwd = ref[:offset] + ref[offset + size:]
		end = start + length + size
		cigar = "%dM%dD%dM" % (offset, size, length - offset)
	else:
		ref = genome.seq(contig, start, start + length - size)
		ins = "".join(ch.substitute(ref[offset - 1 + (i % 2)]) for i in range(size))
		fwd = ref[:offset] + ins + ref[offset:]
		end = start + length - size
		cigar = "%dM%dI%dM" % (offset, size, length - offset - size)
	seq = fwd if strand == "+" else ch.revcomp(fwd)
	r = ch.SimRead(name, seq, "I" * len(seq), contig, start, end, strand, kind="indel",
		extra={"cigar": cigar, "fwd": fwd})
	return r


def indel_placement_unique(ref, offset, kind, size):
	"""True when an indel of `size` bases at `offset` of `ref` (the reference
	stretch for a deletion) cannot slide left or right to an equivalent
	position."""
	if kind == "D":
		return ref[offset - 1] != ref[offset + size - 1] and ref[offset] != ref[offset + size]
	return True


def expected_ksw(genome, contig, draft_end, read_fwd, e, params=DEFAULT_SCORES):
	"""The semi-global optimum chromap's SAM path searches: the window of
	L + 2e reference bases ending e bases after the draft mapping's last base
	(verification window start = end + 1 - L - e), band 2e + 1."""
	L = len(read_fwd)
	ws = draft_end + 1 - L - e
	window = genome.seq(contig, ws, ws + L + 2 * e)
	res = semi_global(window, read_fwd, params, 2 * e + 1)
	res["window_start"] = ws
	res["window"] = window
	return res


# ------------------------------------------------- candidate generation

SMALL_GENOME_SEED = 4242
# Three N-free contigs with a 300 bp segment in six copies on both strands
# (rep6) and a 200 bp segment in two (rep2); everything else is uniform
# random sequence, so a 17-mer outside the segments occurs once.
SMALL_CONTIGS = (("s1", 4000), ("s2", 3000), ("s3", 2000))
SMALL_SEGMENTS = (
	("rep6", 300, (("s1", 500, "+"), ("s1", 1500, "-"), ("s1", 2600, "+"), ("s2", 400, "+"), ("s2", 1800, "-"),
		("s3", 900, "+")), ()),
	("rep2", 200, (("s1", 3400, "+"), ("s2", 2500, "-")), ()),
)


def make_small_genome():
	return ch.make_synthetic_genome(seed=SMALL_GENOME_SEED, contig_specs=SMALL_CONTIGS,
		segment_specs=SMALL_SEGMENTS, n_run_specs=(), low_complexity_specs=())


def minimizer_table(genome):
	"""hash -> [(rid, end position, strand)] of every reference minimizer of
	an N-free genome, by definition (chromap_helpers.reference_minimizers)."""
	table = collections.defaultdict(list)
	for rid, (name, seq) in enumerate(genome.contigs.items()):
		for h, hit in ch.reference_minimizers(seq, K, W):
			table[h].append((rid, (hit >> 1) & 0xffffffff, hit & 1))
	for v in table.values():
		v.sort()
	return table


def candidate_position(ref_hit, read_pos, read_strand, k=K):
	"""Candidate position of a read minimizer hit (Index): the reference
	position minus the read position when both k-mers are on the same strand
	(where the read starts), else the reference position plus the read
	position - k + 1 (where a reverse-strand read ends). Returns (strand,
	rid << 32 | position)."""
	rid, pos, strand = ref_hit
	if strand == read_strand:
		return "+", (rid << 32) | ((pos - read_pos) & 0xffffffff)
	return "-", (rid << 32) | ((pos + read_pos - k + 1) & 0xffffffff)


def candidate_hits(table, read_mins, max_freq, rep_freq):
	"""(positive hits, negative hits, #repetitive seeds, repetitive length)
	of Index::GenerateCandidatePositions: every occurrence of a read minimizer
	whose hash occurs once in the reference or fewer than max_freq times;
	hashes with at least rep_freq occurrences (and more than one) are
	repetitive. Hits are sorted."""
	pos, neg, rep = [], [], []
	for h, p, s in read_mins:
		occ = table.get(h, [])
		if not occ:
			continue
		if len(occ) == 1 or len(occ) < max_freq:
			for o in occ:
				strand, c = candidate_position(o, p, s)
				(pos if strand == "+" else neg).append(c)
		if len(occ) > 1 and len(occ) >= rep_freq:
			rep.append(p)
	return sorted(pos), sorted(neg), len(rep), repetitive_seed_length(rep)


def cluster_hits(hits, e, required, num_minimizers):
	"""Candidates from sorted hits (CandidateProcessor): consecutive hits on
	one reference within e of each other form a group, which also ends once
	it holds num_minimizers hits and the next hit lies more than e past the
	group's best hit. The best hit is the most frequent exact value of the
	group (the first on ties); a group with at least `required` hits gives a
	candidate (best hit, its frequency)."""
	out = []
	if not hits:
		return out
	groups = []
	cur = [hits[0]]
	best, best_n = hits[0], 1
	for h in hits[1:]:
		prev = cur[-1]
		if (h >> 32) != (prev >> 32) or (h & 0xffffffff) > (prev & 0xffffffff) + e or (
				len(cur) >= num_minimizers and (h & 0xffffffff) > (best & 0xffffffff) + e):
			groups.append(cur)
			cur = [h]
			best, best_n = h, 1
		else:
			cur.append(h)
		n = cur.count(h)
		if n > best_n:
			best, best_n = h, n
	groups.append(cur)
	for g in groups:
		if len(g) >= required:
			counts = collections.Counter(g)
			top = max(counts.values())
			first = min(v for v in counts if counts[v] == top)
			out.append((first, top))
	return out


def generate_candidates(table, read_mins, e, s, f0, f1):
	"""CandidateProcessor::GenerateCandidates: a first round of
	candidate_hits with max frequency f0 (repetitive at f0); when it finds no
	hit at all, a second round with max frequency f1. Seeds required per
	candidate: min(s, max(1, #minimizers - #repetitive)), or s when the
	second round found hits on both strands. Returns (positive, negative,
	repetitive length, hits)."""
	pos, neg, nrep, rlen = candidate_hits(table, read_mins, f0, f0)
	high = False
	if not pos and not neg:
		pos, neg, nrep, rlen = candidate_hits(table, read_mins, f1, f0)
		high = bool(pos) and bool(neg)
	n = len(read_mins)
	required = min(s, max(1, n - nrep))
	if high:
		required = s
	return (cluster_hits(pos, e, required, n), cluster_hits(neg, e, required, n), rlen, (pos, neg))
