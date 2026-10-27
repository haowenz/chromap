"""Helpers for the area A tests (index build and load, minimizers, seed
encoding, SequenceBatch, read input files).

Independent references in Python: chromap's minimizer sketch around N bases
and for even k, the index chromap should write for a reference, a parser of
the index file and khash's lookup; plus writers for the FASTA/FASTQ files the
tests read in, and read edits (shortened mates, lower case). Only the
standard library is used."""

import gzip
import struct
from collections import OrderedDict, namedtuple
from pathlib import Path

import chromap_helpers as ch

ACGT = "ACGTacgt"


# ------------------------------------------------------------- minimizers

def _kmer_keys(seq, k):
	"""(end position, hash, strand) of every k-mer of an N-free `seq` that is
	not its own reverse complement, in position order (chromap skips such
	"symmetric" k-mers, which exist only for even k, because their strand is
	undefined)."""
	out = []
	for end in range(k - 1, len(seq)):
		kmer = seq[end - k + 1:end + 1].upper()
		if kmer == ch.revcomp(kmer):
			continue
		h, strand = ch.kmer_minimizer_key(kmer)
		out.append((end, h, strand))
	return out


def _window_minima(keys, w):
	"""Indices (into keys) of every k-mer that is the smallest hash, ties
	included, in at least one window of w consecutive keys; for fewer than w
	keys, the smallest ones."""
	if not keys:
		return []
	chosen = set()
	if len(keys) < w:
		best = min(h for _, h, _ in keys)
		chosen.update(i for i, (_, h, _) in enumerate(keys) if h == best)
	else:
		for s in range(len(keys) - w + 1):
			best = min(h for _, h, _ in keys[s:s + w])
			chosen.update(s + i for i in range(w) if keys[s + i][1] == best)
	return sorted(chosen)


def segments(seq):
	"""(start, end) of every maximal run of A, C, G, T (either case) in seq;
	every other byte breaks a run, as chromap's base table maps it to 4."""
	out = []
	start = None
	for i, c in enumerate(seq):
		if c in ACGT:
			if start is None:
				start = i
		elif start is not None:
			out.append((start, i))
			start = None
	if start is not None:
		out.append((start, len(seq)))
	return out


def sketch(seq, k, w, seq_index=0):
	"""The minimizers chromap's MinimizerGenerator should return for `seq`, as
	[(hash, hit)] in position order, from the definition rather than from
	chromap's code:

	- Within each maximal N-free segment (any byte other than ACGTacgt is an
	  N), a k-mer is a minimizer when its hash is the smallest, ties
	  included, in some window of w consecutive k-mers of that segment; a
	  segment with fewer than w k-mers keeps its smallest. Windows never span
	  an N, and lower case is the same base.
	- Symmetric k-mers (equal to their reverse complement; even k only) are
	  not k-mers for this purpose: windows count the remaining ones.
	- minimap2-style winnowing writes a window's minimizer only once it is
	  superseded or leaves the window, so the minimizer still pending when a
	  segment runs into an N, the smallest k-mer of the segment's last
	  window, is dropped. The last segment of the sequence keeps it, since
	  the pending minimizer is written at the end of the sequence.

	Exact (ValueError otherwise) when no segment has a tie for the smallest
	hash in its first window or, before an N, in its last window, when the
	last segment that runs into an N ends at least w bases before the end of
	the sequence (otherwise its pending minimizer, or another of its k-mers
	still in the window buffer, is written at the end of the sequence), and,
	for even k, when the sequence is N-free (after an N, chromap tests the
	k-mer that skips the N for symmetry)."""
	segs = segments(seq)
	if k % 2 == 0 and segs != [(0, len(seq))] and seq:
		raise ValueError("even k needs an N-free sequence")
	ends_before_n = [b for _, b in segs if b != len(seq)]
	if ends_before_n and len(seq) - ends_before_n[-1] < w:
		raise ValueError("a segment runs into an N fewer than w bases before the end")
	out = []
	for si, (a, b) in enumerate(segs):
		keys = _kmer_keys(seq[a:b], k)
		chosen = _window_minima(keys, w)
		first = [h for _, h, _ in keys[:w]]
		if first and first.count(min(first)) > 1:
			raise ValueError("tie for the smallest hash in a segment's first window")
		if b != len(seq) and keys:
			last = keys[-w:]
			best = min(h for _, h, _ in last)
			if sum(1 for _, h, _ in last if h == best) > 1:
				raise ValueError("tie for the smallest hash in the last window before an N")
			pending = len(keys) - len(last) + [h for _, h, _ in last].index(best)
			chosen = [i for i in chosen if i != pending]
		for i in chosen:
			end, h, strand = keys[i]
			out.append((h, ((seq_index << 32 | (a + end)) << 1) | strand))
	return out


def plain_segment_minimizers(seq, k, w, seq_index=0):
	"""Window minimizers of each N-free segment of `seq` by the plain
	definition (chromap_helpers.reference_minimizers, odd k), without the
	drop of the minimizer pending before an N; tests use it to show that
	their data exercise that drop."""
	out = []
	for a, b in segments(seq):
		if b - a >= k:
			out.extend((h, hit + (a << 1)) for h, hit in ch.reference_minimizers(seq[a:b].upper(), k, w, seq_index))
	return out


# ------------------------------------------------------------- sequence files

def write_fasta(path, records, line_width=60, gz=False):
	"""Write (header, sequence) pairs as FASTA with `line_width` bases per
	line (0: one line). The header is written as given (it may carry a
	comment). Returns the path."""
	path = Path(path)
	path.parent.mkdir(parents=True, exist_ok=True)
	opener = gzip.open if gz else open
	with opener(str(path), "wt") as f:
		for header, seq in records:
			f.write(">%s\n" % header)
			if not seq:
				continue
			step = line_width or len(seq)
			for i in range(0, len(seq), step):
				f.write(seq[i:i + step] + "\n")
	return path


def write_fastq(path, records, gz=False):
	"""Write (header, sequence, quality) triples as four-line FASTQ."""
	path = Path(path)
	path.parent.mkdir(parents=True, exist_ok=True)
	opener = gzip.open if gz else open
	with opener(str(path), "wt") as f:
		for header, seq, qual in records:
			f.write("@%s\n%s\n+\n%s\n" % (header, seq, qual))
	return path


# ------------------------------------------------------------- index files

IndexFile = namedtuple("IndexFile",
	"k w table_size n_buckets size n_occupied upper_bound flags keys vals occurrences entries size_bytes")


def parse_index(path):
	"""Parse a chromap index file (little-endian): int32 k, int32 w, uint32
	lookup table size, the khash dump (uint32 n_buckets, size, n_occupied,
	upper_bound; uint32 flags[max(1, n_buckets / 16)]; uint64 keys[n_buckets];
	uint64 vals[n_buckets]), uint32 occurrence count, uint64 occurrences.

	entries: OrderedDict key -> (bucket, value) for the occupied buckets
	(flag bits 2 = empty, 1 = deleted, two bits per bucket)."""
	data = Path(path).read_bytes()
	k, w, table_size, n_buckets, size, n_occupied, upper_bound = struct.unpack_from("<iiIIIII", data, 0)
	off = 28
	nflags = 1 if n_buckets < 16 else n_buckets >> 4
	if n_buckets == 0:
		flags, keys, vals = (), (), ()
	else:
		flags = struct.unpack_from("<%dI" % nflags, data, off)
		off += 4 * nflags
		keys = struct.unpack_from("<%dQ" % n_buckets, data, off)
		off += 8 * n_buckets
		vals = struct.unpack_from("<%dQ" % n_buckets, data, off)
		off += 8 * n_buckets
	(nocc,) = struct.unpack_from("<I", data, off)
	off += 4
	occurrences = struct.unpack_from("<%dQ" % nocc, data, off)
	off += 8 * nocc
	if off != len(data):
		raise ValueError("%s: %d trailing bytes" % (path, len(data) - off))
	entries = OrderedDict()
	for i in range(n_buckets):
		if (flags[i >> 4] >> ((i & 15) << 1)) & 3 == 0:
			entries[keys[i]] = (i, vals[i])
	return IndexFile(k, w, table_size, n_buckets, size, n_occupied, upper_bound, flags, keys, vals, occurrences,
		entries, len(data))


def khash_find(index_file, key):
	"""Bucket of `key` found the way khash's kh_get probes the table (hash =
	key >> 1 truncated to 32 bits, quadratic probing i += ++step, keys equal
	when equal after >> 1), or None."""
	n = index_file.n_buckets
	if n == 0:
		return None
	mask = n - 1
	flags = index_file.flags
	i = ((key >> 1) & 0xFFFFFFFF) & mask
	last = i
	step = 0
	while True:
		f = (flags[i >> 4] >> ((i & 15) << 1)) & 3
		if f & 2:
			return None
		if not (f & 1) and index_file.keys[i] >> 1 == key >> 1:
			return i
		step += 1
		i = (i + step) & mask
		if i == last:
			return None


def expected_index(sequences, k, w):
	"""What a chromap index of `sequences` (the list of sequences in FASTA
	order, empty records already removed) should hold: (entries, occurrences)
	with entries a dict lookup key -> value. Every minimizer hash h of
	sketch() becomes key h << 1; a hash seen once gets key bit 0 set and its
	hit as the value; the others get value offset << 32 | count, where
	offset counts the hits of all smaller non-singleton hashes, and their
	hits sorted ascending fill occurrences[offset:offset + count]."""
	by_hash = {}
	for i, s in enumerate(sequences):
		for h, hit in sketch(s, k, w, seq_index=i):
			by_hash.setdefault(h, []).append(hit)
	entries = {}
	occurrences = []
	for h in sorted(by_hash):
		hits = sorted(by_hash[h])
		if len(hits) == 1:
			entries[(h << 1) | 1] = hits[0]
		else:
			entries[h << 1] = (len(occurrences) << 32) | len(hits)
			occurrences.extend(hits)
	return entries, occurrences


# ------------------------------------------------------------- reads

def shorten_read(read, length):
	"""A copy of SimRead `read` cut to its first `length` bases (the 3' end is
	removed, as by quality trimming): the genomic interval keeps the read's 5'
	end. Needs a read without an adapter tail."""
	if read.tail:
		raise ValueError("shorten_read() needs a read without a tail")
	if read.strand == "+":
		start, end = read.start, read.start + length
	else:
		start, end = read.end - length, read.end
	mm = [p for p in read.mismatches if start <= p < end]
	return ch.SimRead(read.name, read.seq[:length], read.qual[:length], read.contig, start, end, read.strand, mm,
		kind=read.kind, barcode=read.barcode, barcode_qual=read.barcode_qual, extra=dict(read.extra))


def shorten_mate(pair, mate, length):
	"""A copy of SimPair `pair` with mate 1 or 2 cut to `length` bases at its
	3' end. The fragment is unchanged: it is bounded by the 5' ends."""
	r1 = shorten_read(pair.r1, length) if mate == 1 else pair.r1
	r2 = shorten_read(pair.r2, length) if mate == 2 else pair.r2
	return ch.SimPair(pair.name, r1, r2, pair.contig, pair.start, pair.end, pair.strand, kind=pair.kind,
		barcode=pair.barcode, barcode_qual=pair.barcode_qual, extra=dict(pair.extra))


def lower_read(read):
	"""A copy of SimRead `read` with its sequence in lower case."""
	return ch.SimRead(read.name, read.seq.lower(), read.qual, read.contig, read.start, read.end, read.strand,
		read.mismatches, tail=read.tail, kind=read.kind, barcode=read.barcode, barcode_qual=read.barcode_qual,
		extra=dict(read.extra))


def write_reads(path, reads, fasta=False, gz=False, header="{name}", line_width=0):
	"""Write SimRead records to one FASTQ (or FASTA, with `line_width` bases
	per line, 0 for one line) file; `header` formats each record's header
	from its name (to add comments). Returns the path."""
	if fasta:
		return write_fasta(path, [(header.format(name=r.name), r.seq) for r in reads], line_width=line_width, gz=gz)
	return write_fastq(path, [(header.format(name=r.name), r.seq, r.qual) for r in reads], gz=gz)
