"""Area A, C++ layer: base and seed encoding, Hash64, hit encoding,
Minimizer, MinimizerGenerator and SequenceBatch, through tests/cpp/a_sequence.cc.

Hand-worked values live in the driver's CT_TESTs; the probes print what
chromap computes and the expected values are computed here, in Python
(chromap_helpers.encode_seed and hash64, helpers_a.sketch for minimizers, and
the records each test writes for SequenceBatch)."""

import gzip
import random

import pytest

import chromap_helpers as ch
import helpers_a as ha

pytestmark = pytest.mark.cpp

A_SEQUENCE_TESTS = [
	"base_tables_hand_worked",
	"seed_lengths_one_and_32",
	"hit_encoding_extremes",
	"minimizer_class",
]

A_SEQUENCE_PROBES = ["seeds", "hash64", "minimizers", "load", "slots", "strands", "trim", "edit", "reorder", "swap",
	"batch_seeds"]

CORRUPT_MESSAGE = "Didn't reach the end of sequence file, which might be corrupted!"


def _driver(cpp):
	return cpp("a_sequence")


def test_a_sequence_driver_lists_every_test(cpp):
	"""Every CT_TEST in cpp/a_sequence.cc is run by test_a_sequence_cpp, and
	every probe used below exists."""
	driver = _driver(cpp)
	assert driver.tests() == A_SEQUENCE_TESTS
	assert driver.probes() == A_SEQUENCE_PROBES


@pytest.mark.parametrize("name", A_SEQUENCE_TESTS)
def test_a_sequence_cpp(cpp, name):
	_driver(cpp).check(name)


def _random_seq(rng, n, alphabet="ACGT"):
	return "".join(rng.choice(alphabet) for _ in range(n))


def _revcomp_chromap(seq):
	"""Reverse complement as chromap writes it: A, C, G, T in either case
	complement to upper case, anything else becomes N."""
	comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
	return "".join(comp.get(c.upper(), "N") for c in reversed(seq))


# ------------------------------------------------------------ encodings

def test_seed_encoding_matches_reference(cpp):
	"""GenerateSeedFromSequence equals chromap_helpers.encode_seed: two bits
	per base, first base highest, A0 C1 G2 T3 in either case, any other byte
	(N, n, IUPAC codes) encodes as A, positions past the end pad with A; for
	start offsets inside, at and past the end and seed lengths 1 to 32."""
	rng = random.Random(ch.stable_seed("a_sequence", "seeds"))
	cases = [("ACGT", 0, 4), ("ACGT", 2, 2), ("AC", 0, 4), (".", 0, 5), ("T", 0, 1)]
	for _ in range(300):
		seq = _random_seq(rng, rng.randint(1, 40), "ACGTACGTacgtNnRY")
		length = rng.randint(1, 32)
		start = rng.randint(0, len(seq) + 3)
		cases.append((seq, start, length))
	text = "".join("%s %d %d\n" % c for c in cases)
	got = [int(x) for x in _driver(cpp).probe("seeds", input=text).split()]
	want = [ch.encode_seed("" if s == "." else s, start, length) for s, start, length in cases]
	assert got == want
	assert got[:5] == [27, 11, 16, 0, 3]


def test_batch_seed_encoding_matches_reference(cpp, tmp_path):
	"""SequenceBatch::GenerateSeedFromSequenceAt encodes the stored sequence
	like GenerateSeedFromSequence: for several start offsets and lengths
	(including a range running past the end, padded with A)."""
	rng = random.Random(ch.stable_seed("a_sequence", "batch_seeds"))
	seqs = [_random_seq(rng, n, "ACGTacgtN") for n in (5, 17, 31, 32, 40)]
	fasta = ha.write_fasta(tmp_path / "s.fa", [("s%d" % i, s) for i, s in enumerate(seqs)])
	for start, length in ((0, 1), (0, 17), (3, 32), (20, 15), (38, 4)):
		got = [int(x) for x in _driver(cpp).probe("batch_seeds", fasta, start, length).split()]
		assert got == [ch.encode_seed(s, start, length) for s in seqs], (start, length)


def test_hash64_matches_reference(cpp):
	"""Hash64 equals Thomas Wang's 64-bit integer hash restricted to the mask
	(chromap_helpers.hash64) for keys 0, mask and 30 random keys under every
	k-mer mask of 2k bits, k = 1 to 28."""
	rng = random.Random(ch.stable_seed("a_sequence", "hash64"))
	cases = []
	for k in range(1, 29):
		mask = (1 << (2 * k)) - 1
		cases.extend((key, 2 * k) for key in [0, mask] + [rng.randint(0, mask) for _ in range(30)])
	text = "".join("%d %d\n" % c for c in cases)
	got = [int(x) for x in _driver(cpp).probe("hash64", input=text).split()]
	assert got == [ch.hash64(key, (1 << bits) - 1) for key, bits in cases]


def test_hash64_is_a_bijection_on_small_masks(cpp):
	"""Under a mask of 2k bits Hash64 permutes the 4^k keys (Wang's hash is
	invertible), so distinct k-mers never collide; checked for k = 1 to 6."""
	cases = [(key, 2 * k) for k in range(1, 7) for key in range(4 ** k)]
	text = "".join("%d %d\n" % c for c in cases)
	got = [int(x) for x in _driver(cpp).probe("hash64", input=text).split()]
	pos = 0
	for k in range(1, 7):
		n = 4 ** k
		assert sorted(got[pos:pos + n]) == list(range(n)), k
		pos += n


# ------------------------------------------------------------ minimizers

def _minimizers(cpp, tmp_path, seqs, k, w, append=False, name="m.fa"):
	fasta = ha.write_fasta(tmp_path / name, [("s%d" % i, s) for i, s in enumerate(seqs)], line_width=37)
	extra = ["append"] if append else []
	out = _driver(cpp).probe("minimizers", fasta, k, w, *extra)
	got = {}
	for line in out.splitlines():
		i, h, hit = (int(x) for x in line.split("\t"))
		got.setdefault(i, []).append((h, hit))
	return got


def _sequences_with_n(rng, n, k, w):
	"""Random sequences of N-free segments (upper and lower case) separated by
	runs of N, n and IUPAC codes, kept when helpers_a.sketch is exact for them."""
	out = []
	while len(out) < n:
		s = ""
		for j in range(rng.randint(2, 5)):
			s += _random_seq(rng, rng.randint(0, 3 * (w + k)), rng.choice(["ACGT", "acgt", "ACGTacgt"]))
			s += rng.choice(["N", "n", "R", "Y", "NNNNN"]) * rng.randint(1, 30)
		s += _random_seq(rng, rng.randint(w + k, 4 * (w + k)))
		try:
			ha.sketch(s, k, w)
		except ValueError:
			continue
		out.append(s)
	return out


@pytest.mark.parametrize("k,w", [(17, 7), (11, 5), (21, 11), (19, 10)])
def test_minimizers_restart_at_n_bases(cpp, tmp_path, k, w):
	"""Around ambiguous bases MinimizerGenerator returns helpers_a.sketch: the
	window minimizers of each N-free segment (lower case is the same base;
	N, n, R and Y all break the sequence), no k-mer that spans an ambiguous
	base, and, as minimap2-style winnowing does, without the minimizer still
	pending when a segment runs into an N. The sequences are chosen so that
	this dropped minimizer exists (the plain per-segment definition differs
	from the expectation for most of them)."""
	rng = random.Random(ch.stable_seed("a_sequence", "n_runs", k, w))
	seqs = _sequences_with_n(rng, 40, k, w)
	got = _minimizers(cpp, tmp_path, seqs, k, w)
	differs = 0
	for i, s in enumerate(seqs):
		want = ha.sketch(s, k, w, seq_index=i)
		assert got.get(i, []) == want, "sequence %d: %r" % (i, s)
		differs += ha.plain_segment_minimizers(s, k, w, i) != want
		for h, hit in want:
			end = (hit >> 1) & 0xFFFFFFFF
			assert all(c in ha.ACGT for c in s[end - k + 1:end + 1])
	assert differs > len(seqs) // 2


def test_minimizers_lower_case_equals_upper_case(cpp, tmp_path):
	"""A sequence in lower case (or mixed case) has exactly the minimizers of
	its upper-case version: same hashes, positions and strands."""
	rng = random.Random(ch.stable_seed("a_sequence", "lower"))
	upper = [_random_seq(rng, n) for n in (30, 200, 1000)]
	mixed = ["".join(c.lower() if rng.random() < 0.5 else c for c in s) for s in upper]
	lower = [s.lower() for s in upper]
	got = _minimizers(cpp, tmp_path, upper + mixed + lower, 17, 7)
	for i in range(3):
		want = ch.reference_minimizers(upper[i], 17, 7, seq_index=0)
		for j in (i, i + 3, i + 6):
			strip = [(h, hit - ((j << 32) << 1)) for h, hit in got[j]]
			assert strip == want, j


@pytest.mark.parametrize("k,w", [(2, 3), (4, 2), (8, 5), (16, 7), (28, 4)])
def test_minimizers_even_k_skip_symmetric_kmers(cpp, tmp_path, k, w):
	"""For even k a k-mer can equal its reverse complement (ACGT for k 4); its
	strand is undefined and chromap skips it: no such k-mer is returned and
	windows count only the other k-mers (helpers_a.sketch). The sequences
	carry palindromic stretches (repeats of ACGT and of AT) so the case
	occurs."""
	rng = random.Random(ch.stable_seed("a_sequence", "even", k, w))
	seqs = []
	while len(seqs) < 12:
		n = len(seqs)
		parts = [_random_seq(rng, rng.randint(0, 40)), ("ACGT" if n % 2 else "AT") * rng.randint(1, 12),
			_random_seq(rng, rng.randint(0, 40)), ("GAATTC" * 6)[:k], _random_seq(rng, rng.randint(k, 60))]
		try:
			ha.sketch("".join(parts), k, w)
		except ValueError:
			continue
		seqs.append("".join(parts))
	got = _minimizers(cpp, tmp_path, seqs, k, w)
	num_symmetric = 0
	for i, s in enumerate(seqs):
		assert got.get(i, []) == ha.sketch(s, k, w, seq_index=i), "sequence %d: %r" % (i, s)
		for h, hit in got.get(i, []):
			end = (hit >> 1) & 0xFFFFFFFF
			kmer = s[end - k + 1:end + 1]
			assert kmer != ch.revcomp(kmer)
		num_symmetric += sum(1 for e in range(k - 1, len(s)) if s[e - k + 1:e + 1] == ch.revcomp(s[e - k + 1:e + 1]))
	assert num_symmetric > 0


def test_minimizers_short_sequences(cpp, tmp_path):
	"""A sequence shorter than k has no minimizer; one of exactly k bases has
	one, its only k-mer, at end position k - 1; a sequence of fewer than w
	k-mers keeps its smallest k-mer."""
	rng = random.Random(ch.stable_seed("a_sequence", "short"))
	k, w = 17, 7
	seqs = [_random_seq(rng, n) for n in (1, 5, 16, 17, 18, 22)]
	got = _minimizers(cpp, tmp_path, seqs, k, w)
	assert not any(i in got for i in range(3))
	h, strand = ch.kmer_minimizer_key(seqs[3])
	assert got[3] == [(h, ((3 << 32 | 16) << 1) | strand)]
	for i in (4, 5):
		keys = [ch.kmer_minimizer_key(seqs[i][e - k + 1:e + 1]) + (e,) for e in range(k - 1, len(seqs[i]))]
		h, strand, e = min(keys)
		assert got[i] == [(h, ((i << 32 | e) << 1) | strand)]


@pytest.mark.parametrize("k,w", [(17, 7), (11, 5), (15, 10)])
def test_minimizers_first_window_ties(cpp, tmp_path, k, w):
	"""Identical k-mers tied for the smallest hash in a sequence's first
	window (from a tandem repeat of period below w, or a k-mer and its reverse
	complement in an inverted repeat, at the start) are all minimizers, as in
	every later window: MinimizerGenerator returns the window minimizers by
	definition, ties included (chromap_helpers.reference_minimizers). Kept to
	sequences whose w-th k-mer does not itself tie the minimum: there chromap
	writes only that last copy for the first window, a case the definition
	leaves open."""
	rng = random.Random(ch.stable_seed("a_sequence", "first_window_ties", k, w))
	seqs = []
	for attempt in range(20000):
		if len(seqs) == 30:
			break
		kind = attempt % 3
		lead = _random_seq(rng, rng.randint(0, 4))
		if kind == 0:
			unit = _random_seq(rng, rng.randint(2, w - 1))
			s = lead + (unit * (k + w))[:rng.randint(k + len(unit), k + w + 4)]
		elif kind == 1:
			s = lead + ("CA" * 30)[:rng.randint(k + 1, k + w + 3)]
		else:
			# Inverted repeat A + revcomp(A): the k-mers at offsets 0 and d are
			# each other's reverse complement.
			d = rng.choice([x for x in (1, 3, 5) if x < w])
			a = _random_seq(rng, (k + d) // 2)
			s = lead + a + ch.revcomp(a)
		s += _random_seq(rng, rng.randint(10, 80))
		keys = [ch.kmer_minimizer_key(s[e - k + 1:e + 1])[0] for e in range(k - 1, len(s))]
		first = keys[:w]
		if first.count(min(first)) > 1 and keys[w - 1] != min(first):
			seqs.append((s, min(first)))
	assert len(seqs) == 30
	got = _minimizers(cpp, tmp_path, [s for s, _ in seqs], k, w)
	for i, (s, smallest) in enumerate(seqs):
		want = ch.reference_minimizers(s, k, w, seq_index=i)
		assert got.get(i, []) == want, "sequence %d: %r" % (i, s)
		# Several tied minimizers end inside the first window.
		assert sum(1 for h, hit in want if h == smallest and ((hit >> 1) & 0xFFFFFFFF) < k - 1 + w) > 1


def test_minimizers_append_to_non_empty_vector(cpp, tmp_path):
	"""GenerateMinimizers appends to the vector it is given: an entry already
	there stays first and the sequence's minimizers follow unchanged."""
	rng = random.Random(ch.stable_seed("a_sequence", "append"))
	seqs = [_random_seq(rng, 300), _random_seq(rng, 60)]
	plain = _minimizers(cpp, tmp_path, seqs, 15, 6)
	appended = _minimizers(cpp, tmp_path, seqs, 15, 6, append=True, name="m2.fa")
	for i in range(2):
		assert appended[i] == [(12345, 678)] + plain[i]
		assert plain[i] == ch.reference_minimizers(seqs[i], 15, 6, seq_index=i)


# ------------------------------------------------------------ SequenceBatch

def _records(out):
	"""Parse the probe's rec lines into (slot, id, name, name_length, length,
	seq, qual) tuples and the other lines into a dict."""
	recs, info = [], {}
	for line in out.splitlines():
		f = line.split("\t")
		if f[0] == "rec":
			recs.append((int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5]), f[6], f[7]))
		elif f[0] in ("num_sequences", "num_bases", "max_batch"):
			info[f[0]] = int(f[1])
	return recs, info


def _fasta_text():
	"""A FASTA file exercising the parser and the records it should give:
	headers with a space and a tab comment, sequence lines of different
	widths, a blank line, lower case and N, an empty record in the middle and
	one at the end (both skipped without using an id)."""
	text = (">s1 first record\nACGTA\nCC\nggtt\n\n>empty1\n>s2\tTabbed comment\nNNNNacgtnRYACGT\n>s3\n"
		+ "A" * 70 + "\n" + "C" * 30 + "\n>empty2 no bases\n")
	want = [("s1", "ACGTACCggtt"), ("s2", "NNNNacgtnRYACGT"), ("s3", "A" * 70 + "C" * 30)]
	return text, want


def _fastq_text():
	"""FASTQ records: a comment, an empty record (skipped), a multi-line
	record, a quality line starting with @ and a + line repeating the name."""
	text = ("@q1 comment here\nACGTN\n+\nIIII#\n@q2\n\n+\n\n@q3\nACGT\nTTGG\n+\n!!!!\n####\n"
		+ "@q4 x\nacgt\n+q4 x\n@III\n")
	want = [("q1", "ACGTN", "IIII#"), ("q3", "ACGTTTGG", "!!!!####"), ("q4", "acgt", "@III")]
	return text, want


@pytest.mark.parametrize("gz", [False, True])
def test_load_all_sequences_fasta(cpp, tmp_path, gz):
	"""LoadAllSequences (the reference path) keeps every non-empty FASTA
	record in file order, plain or gzipped: the name up to the first
	whitespace (space or tab) and its length, the sequence lines joined with
	case and N kept, ids 0, 1, ... skipping empty records, GetNumSequences
	and GetNumBases over the kept records; there is no quality (null) and the
	default constructor has no batch limit (GetMaxBatchSize 0)."""
	text, want = _fasta_text()
	path = tmp_path / ("ref.fa" + (".gz" if gz else ""))
	if gz:
		with gzip.open(str(path), "wt") as f:
			f.write(text)
	else:
		path.write_text(text)
	recs, info = _records(_driver(cpp).probe("load", "all", path))
	assert recs == [(i, i, n, len(n), len(s), s, "*") for i, (n, s) in enumerate(want)]
	assert info == {"num_sequences": 3, "num_bases": sum(len(s) for _, s in want), "max_batch": 0}


@pytest.mark.parametrize("gz", [False, True])
def test_load_all_sequences_fastq(cpp, tmp_path, gz):
	"""FASTQ records load with their qualities (multi-line sequence and
	quality joined; a quality line may start with @ and the + line may repeat
	the name); an empty record is skipped without an id."""
	text, want = _fastq_text()
	path = tmp_path / ("reads.fq" + (".gz" if gz else ""))
	if gz:
		with gzip.open(str(path), "wt") as f:
			f.write(text)
	else:
		path.write_text(text)
	recs, info = _records(_driver(cpp).probe("load", "all", path))
	assert recs == [(i, i, n, len(n), len(s), s, q) for i, (n, s, q) in enumerate(want)]
	assert info["num_sequences"] == 3


def test_load_batch_limits_and_ids(cpp, tmp_path):
	"""LoadBatch fills at most GetMaxBatchSize records per call, returns 0 at
	the end of the file and keeps counting ids across calls (7 records in
	batches of 3: 3, 3, 1, 0 with ids 0-6); GetNumSequences is the size of
	the current batch. (GetNumBases is documented only for LoadAllSequences
	and the end-of-file message is a log line, so neither is asserted
	here.)"""
	rng = random.Random(ch.stable_seed("a_sequence", "batch"))
	reads = [("r%d extra" % i, _random_seq(rng, 10 + i), "".join(chr(33 + (i + j) % 40) for j in range(10 + i)))
		for i in range(7)]
	path = ha.write_fastq(tmp_path / "r.fq.gz", reads, gz=True)
	proc = _driver(cpp).run(["--probe", "load", "batch", path, 3])
	assert proc.returncode == 0
	out = proc.stdout.decode()
	batches = [(int(f.split("\t")[1]), int(f.split("\t")[2])) for f in out.splitlines() if f.startswith("batch")]
	assert batches == [(3, 3), (3, 3), (1, 1), (0, 0)]
	recs, info = _records(out)
	assert recs == [(i % 3, i, "r%d" % i, len("r%d" % i), 10 + i, reads[i][1], reads[i][2]) for i in range(7)]
	assert (info["num_sequences"], info["max_batch"]) == (0, 3)


def test_load_one_sequence_at_slots(cpp, tmp_path):
	"""LoadOneSequenceAndSaveAt returns false for each record it stores and
	true at the end of the file; GetNumSequences counts the slots filled."""
	text, want = _fasta_text()
	path = tmp_path / "ref.fa"
	path.write_text(text)
	out = _driver(cpp).probe("load", "one", path, 5)
	calls = [tuple(int(x) for x in line.split("\t")[1:]) for line in out.splitlines() if line.startswith("one")]
	assert calls == [(0, 0, 1), (1, 0, 2), (2, 0, 3), (3, 1, 3)]
	recs, info = _records(out)
	assert [(r[0], r[1], r[2], r[5]) for r in recs] == [(i, i, n, s) for i, (n, s) in enumerate(want)]
	assert info["max_batch"] == 5


def test_load_one_sequence_overwrite_rules(cpp, tmp_path):
	"""Slot 0 starts a new batch and the last loaded slot may be loaded again
	(ids still advance: slots 0 1 2 2 0 get ids 0-4); asking for an earlier
	slot (1 after 0 1 2) exits with chromap's message and status 255."""
	path = ha.write_fasta(tmp_path / "five.fa", [(c, c * 4) for c in "ACGTA"[:4]] + [("e", "ACGT")])
	out = _driver(cpp).probe("slots", path, 8, 0, 1, 2, 2, 0, 1)
	slots = [tuple(int(x) for x in line.split("\t")[1:]) for line in out.splitlines() if line.startswith("slot")]
	assert slots == [(0, 0, 1), (1, 0, 2), (2, 0, 3), (2, 0, 3), (0, 0, 1), (1, 1, 1)]
	recs, _ = _records(out)
	assert [(r[0], r[1], r[5]) for r in recs] == [(0, 0, "AAAA"), (1, 1, "CCCC"), (2, 2, "GGGG"), (2, 3, "TTTT"),
		(0, 4, "ACGT")]
	proc = _driver(cpp).run(["--probe", "slots", path, 8, 0, 1, 2, 1])
	assert proc.returncode == 255
	assert "Shouldn't override other sequences rather than the last!" in proc.stderr.decode()
	assert proc.stdout.decode().count("slot\t") == 3


@pytest.mark.parametrize("mode", ["all", "batch", "one"])
def test_missing_file_exits(cpp, tmp_path, mode):
	"""InitializeLoading on a path that does not exist exits with status 255
	and "Cannot find sequence file <path>"."""
	path = tmp_path / "missing.fq"
	proc = _driver(cpp).run(["--probe", "load", mode, path, 4])
	assert proc.returncode == 255
	assert ("Cannot find sequence file %s" % path) in proc.stderr.decode()
	assert proc.stdout == b""


@pytest.mark.parametrize("mode", ["all", "batch", "one"])
@pytest.mark.parametrize("tail", ["IIII\n", "\n", "", "IIIIIIII\n"])
def test_truncated_fastq_exits(cpp, tmp_path, mode, tail):
	"""A FASTQ record whose quality is shorter than its sequence, empty, or
	missing at the end of the file, or longer than the sequence (8 for 6
	bases), ends loading with "Didn't reach the end of sequence file, which
	might be corrupted!" and status 255, in every loading mode."""
	path = tmp_path / "trunc.fq"
	path.write_text("@ok\nACGTAC\n+\nIIIIII\n@bad\nACGTAC\n+\n" + tail)
	proc = _driver(cpp).run(["--probe", "load", mode, path, 4])
	assert proc.returncode == 255
	assert CORRUPT_MESSAGE in proc.stderr.decode()


def test_negative_sequence_and_n_positions(cpp, tmp_path):
	"""PrepareNegativeSequenceAt stores the reverse complement (A, C, G, T of
	either case complemented to upper case, anything else N), readable
	through GetNegativeSequenceAt and GetNegativeSequenceBatch. GetSequenceNsAt
	lists the N positions from the left, or (seed order) as distances from
	the right end, and IsNInSequenceAt says whether there is one."""
	seqs = ["ACGTNacgtRYK", "AAAAAN", "TTTT", "NGATTACAN", "gattaca"]
	path = ha.write_fasta(tmp_path / "s.fa", [("s%d" % i, s) for i, s in enumerate(seqs)])
	lines = _driver(cpp).probe("strands", path, 8).splitlines()
	assert len(lines) == len(seqs)
	for s, line in zip(seqs, lines):
		seq, neg, neg_batch, has_n, big, little = line.split("\t")
		assert (seq, neg, neg_batch) == (s, _revcomp_chromap(s), _revcomp_chromap(s))
		npos = [i for i, c in enumerate(s) if c == "N"]
		assert has_n == ("1" if npos else "0")
		assert big == (",".join(map(str, npos)) or "-")
		assert little == (",".join(str(len(s) - 1 - i) for i in reversed(npos)) or "-")


@pytest.mark.parametrize("length", [0, 3, 9, 10, 25])
def test_trim_sequence(cpp, tmp_path, length):
	"""TrimSequenceAt(L) keeps the first L bases and qualities and drops the
	same bases from the reverse complement (its start); an L at or above the
	length changes nothing."""
	reads = [("a", "ACGTNacgta", "0123456789"), ("b", "GGGATTTCCA", "ABCDEFGHIJ")]
	path = ha.write_fastq(tmp_path / "t.fq", reads)
	lines = _driver(cpp).probe("trim", path, 4, length).splitlines()
	for (_, s, q), line in zip(reads, lines):
		n = min(length, len(s))
		assert line.split("\t") == [str(n), s[:n], q[:n], _revcomp_chromap(s[:n])]


def test_correct_base_at(cpp, tmp_path):
	"""CorrectBaseAt replaces one base of one record and nothing else."""
	path = ha.write_fasta(tmp_path / "s.fa", [("a", "AAAA"), ("b", "CCCC"), ("c", "GGGG")])
	assert _driver(cpp).probe("edit", path, 4, 1, 2, "T").split() == ["AAAA", "CCTC", "GGGG"]
	assert _driver(cpp).probe("edit", path, 4, 2, 0, "n").split() == ["AAAA", "CCCC", "nGGG"]


@pytest.mark.parametrize("mode", ["all", "batch"])
def test_reorder_sequences(cpp, tmp_path, mode):
	"""ReorderSequences(rank) moves sequence i to slot rank[i] (with its name
	and, for a batch, its reverse complement)."""
	names = ["a", "b", "c", "d", "e"]
	seqs = ["AAAC", "CCCA", "GGGT", "TTTG", "ACGG"]
	path = ha.write_fasta(tmp_path / "s.fa", list(zip(names, seqs)))
	rank = [2, 0, 4, 1, 3]
	lines = _driver(cpp).probe("reorder", mode, path, *rank).splitlines()
	slot = [None] * 5
	for i, r in enumerate(rank):
		slot[r] = i
	assert [line.split("\t")[:2] for line in lines] == [[names[i], seqs[i]] for i in slot]
	if mode == "batch":
		assert [line.split("\t")[2] for line in lines] == [ch.revcomp(seqs[i]) for i in slot]


def test_swap_sequence_batch(cpp, tmp_path):
	"""SwapSequenceBatch exchanges the records and reverse complements of two
	batches."""
	one = ha.write_fasta(tmp_path / "one.fa", [("a", "AACG"), ("b", "TTGCA")])
	two = ha.write_fasta(tmp_path / "two.fa", [("c", "GATC"), ("d", "CCAT"), ("e", "acgtt")])
	lines = _driver(cpp).probe("swap", one, two, 3).splitlines()
	want = [("first", s) for s in ("GATC", "CCAT", "acgtt")] + [("second", s) for s in ("AACG", "TTGCA")]
	assert [tuple(line.split("\t")) for line in lines] == [(w, s, _revcomp_chromap(s)) for w, s in want]
