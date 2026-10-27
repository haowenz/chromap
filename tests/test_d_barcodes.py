"""Area D: cell barcodes through the command line: -b, whitelist correction
(--barcode-whitelist, --bc-error-threshold, --bc-probability-threshold,
--output-mappings-not-in-whitelist, --skip-barcode-check), the barcode
statistics, cell-level deduplication of corrected barcodes, the fragment
output and --barcode-translate.

Every expected barcode comes from helpers_d.correct_barcode, chromap's
documented rule written in Python (README "Map scATAC-seq reads", chromap.1):
a barcode in the whitelist is kept; otherwise the whitelist barcodes within
--bc-error-threshold mismatches (an N matches nothing) are candidates; one
candidate is taken; several are scored by the error probability of the
mismatched bases (Phred+33, clamped to [3, 40]) times the candidate's
abundance (its exact matches among all input barcodes), and the best is
taken when its share of the summed scores is above
--bc-probability-threshold. Reads whose barcode cannot be corrected are not
written unless --output-mappings-not-in-whitelist is given. Every read comes
from its own unique locus, so each output line is matched to its read by
position. Whitelists are synthetic."""

import collections
import gzip

import pytest

import chromap_helpers as ch
import helpers_d as hd

BC_LEN = 16


class Design(object):
	"""A barcoded single-end read set with its whitelist. `names[category]`
	lists the read names of each category; `source[name]` is the whitelist
	barcode a read's barcode was made from; `designs[label]` is (read name,
	[candidate whitelist barcodes in order a, b, c])."""

	def __init__(self):
		self.records = []
		self.whitelist = []
		self.names = collections.defaultdict(list)
		self.source = {}
		self.designs = {}

	@property
	def reads(self):
		return ch.ReadSet(self.records)


def build_design(genome):
	"""Categories (counts): exact (72: 3 per base whitelist barcode), one (40:
	one substitution), two (24), three (10), n1 (22: one N, quality '#' there),
	n2 (10: two Ns, half over an A of the whitelist barcode at the second),
	n1sub (10: an N and a substitution before it), far (8:
	random barcodes at least 4 from every whitelist barcode), and designed
	ambiguous reads, each one substitution away from two or three whitelist
	barcodes of its own (helpers_d.isolated_barcodes keeps the designs 5 or
	more mismatches from everything else), with exact reads of those barcodes
	("abundance") to set their abundances:

	- abundance: a:b = 19:1, equal qualities -> a (0.95 > 0.9)
	- quality: 5:5, '$' (Q3) at a's mismatch, 'I' (Q40) at b's -> a (0.9998)
	- tie: 5:5, equal qualities -> 0.5, uncorrectable at 0.9; two more
	  reads of a with an N over one of its A bases (its 2-bit key, N read
	  as A, equals a) are corrected to a but add nothing to a's abundance
	- low_clamp: 17:2, '!' (Q0, clamped to 3) and '$' (Q3) -> 17/19 = 0.895,
	  uncorrectable (unclamped it would be 0.944)
	- high_clamp: 17:2, 'I' (Q40) and 'K' (Q42, clamped to 40) -> 0.895,
	  uncorrectable (unclamped 0.931)
	- quality_over_abundance: 1:5, '$' at a's mismatch, 'I' at b's -> a
	  (0.501 * 1 against 1e-4 * 5: 0.999)
	- three_way: 17:1:1, equal qualities -> 17/19 = 0.895, uncorrectable
	  (the best against the second alone would be 0.944)
	- unseen: one candidate, never seen exactly (abundance 0) -> corrected
	"""
	d = Design()
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_d_barcodes", "design"))
	rng = hd.stable_rng("test_d_barcodes", "design")
	picker = hd.LocusPicker(sim)
	base = ch.make_whitelist(24, length=BC_LEN, seed=31, min_dist=4)
	used = list(base)
	wl = list(base)

	def add(category, barcode, qual, source=None):
		r = hd.barcoded_read(sim, picker, barcode, qual)
		d.records.append(r)
		d.names[category].append(r.name)
		if source is not None:
			d.source[r.name] = source
		return r

	def quals():
		return hd.random_quality(rng, BC_LEN)

	for label, abund, qchars in (
			("abundance", (19, 1), ("I", "I")),
			("quality", (5, 5), ("$", "I")),
			("tie", (5, 5), ("I", "I")),
			("low_clamp", (17, 2), ("!", "$")),
			("high_clamp", (17, 2), ("I", "K")),
			("quality_over_abundance", (1, 5), ("$", "I")),
			("three_way", (17, 1, 1), ("I", "I", "I")),
			("unseen", (0,), ("I",))):
		x = hd.isolated_barcodes(1, BC_LEN, rng, used, 5)[0]
		used.append(x)
		positions = sorted(rng.sample(range(BC_LEN), len(abund)))
		cands = [hd.mutate(x, [p]) for p in positions]
		used += cands
		wl += cands
		q = ["I"] * BC_LEN
		for p, c in zip(positions, qchars):
			q[p] = c
		r = add("design", x, "".join(q))
		d.designs[label] = (r.name, cands)
		for w, n in zip(cands, abund):
			for _ in range(n):
				add("abundance", w, quals(), w)
		if label == "tie":
			a_bases = [k for k in range(BC_LEN) if cands[0][k] == "A" and k not in positions]
			assert a_bases
			for k in a_bases[:2]:
				q = list(quals())
				q[k] = "#"
				add("n1", hd.with_n(cands[0], [k]), "".join(q), cands[0])
	for w in base:
		for _ in range(3):
			add("exact", w, quals(), w)
	for i in range(40):
		w = base[i % len(base)]
		add("one", hd.mutate(w, [rng.randrange(BC_LEN)]), quals(), w)
	for i in range(24):
		w = base[i % len(base)]
		add("two", hd.mutate(w, rng.sample(range(BC_LEN), 2)), quals(), w)
	for i in range(10):
		w = base[(i * 7) % len(base)]
		add("three", hd.mutate(w, rng.sample(range(BC_LEN), 3)), quals(), w)
	for i in range(20):
		w = base[(i * 5) % len(base)]
		p = rng.randrange(BC_LEN)
		q = list(quals())
		q[p] = "#"
		add("n1", hd.with_n(w, [p]), "".join(q), w)
	for i in range(10):
		# Two Ns where the whitelist base under the first N is not A (half
		# of them with an A under the second), and an N with a substitution
		# before it: the cases --bc-error-threshold 2 handles correctly (see
		# known bug d-bc-threshold2-n).
		w = base[(i * 3 + 1) % len(base)]
		pairs = [(a, b) for a in range(BC_LEN) for b in range(a + 1, BC_LEN)
			if w[a] != "A" and (w[b] == "A") == (i % 2 == 0)]
		ps = list(rng.choice(pairs))
		q = list(quals())
		for p in ps:
			q[p] = "#"
		add("n2", hd.with_n(w, ps), "".join(q), w)
		w = base[(i * 3 + 2) % len(base)]
		pn = rng.randrange(4, BC_LEN)
		ps = rng.randrange(0, pn)
		q = list(quals())
		q[pn] = "#"
		add("n1sub", hd.with_n(hd.mutate(w, [ps]), [pn]), "".join(q), w)
	for b in hd.isolated_barcodes(8, BC_LEN, rng, used, 4):
		add("far", b, quals())
	d.records = list(sim.shuffled(ch.ReadSet(d.records)))
	rng.shuffle(wl)
	d.whitelist = wl
	return d


@pytest.fixture(scope="module")
def design(genome, tmp_path_factory):
	d = build_design(genome)
	d.whitelist_path = ch.write_lines(tmp_path_factory.mktemp("d_whitelist") / "whitelist.txt", d.whitelist)
	return d


def _run(index, design, extra=()):
	args = ["--barcode-whitelist", str(design.whitelist_path)] + list(extra)
	return ch.run_chromap(index, design.reads, args, shared=True)


@pytest.fixture(scope="module")
def default_run(index, design):
	return _run(index, design)


@pytest.fixture(scope="module")
def keep_run(index, design):
	"""--output-mappings-not-in-whitelist, at -t 1. The barcode statistics on
	stderr (counted before the option decides what is written) are checked
	on this run: they are summed from per-thread counters, which clang/libomp
	builds undercount at -t > 1 (known bug a-mapping-statistics-clang-threads)."""
	args = ["--barcode-whitelist", str(design.whitelist_path), "--output-mappings-not-in-whitelist"]
	return ch.run_chromap(index, design.reads, args, threads=1, shared=True)


@pytest.fixture(scope="module")
def two_error_run(index, design):
	return _run(index, design, ["--bc-error-threshold", "2"])


@pytest.fixture(scope="module")
def low_probability_run(index, design):
	return _run(index, design, ["--bc-probability-threshold", "0.45"])


def _lines(run):
	"""Counter of (chrom, start, end, barcode) of a fragment output, after
	checking every count column is 1 (no deduplication requested)."""
	recs = run.fragments()
	assert set(f.dups for f in recs) <= {1}
	return collections.Counter((f.chrom, f.start, f.end, f.barcode) for f in recs)


def _barcode_at(run, read):
	got = hd.fragments_by_locus(run).get((read.contig, read.start, read.end), [])
	assert len(got) <= 1, (read, got)
	return got[0] if got else None


# ------------------------------------------------------------ default run

def test_fragment_lines_follow_correction_rule(design, default_run):
	"""Every fragment line (chrom, start, end, barcode) is the read's true
	interval with the barcode helpers_d.correct_barcode gives at the default
	thresholds (1 mismatch, probability 0.9); uncorrectable reads are
	absent."""
	want = hd.expected_fragments(design.records, design.whitelist)
	assert _lines(default_run) == want


def test_exact_barcodes_are_kept(design, default_run):
	"""Reads whose barcode is in the whitelist are written with it."""
	reads = default_run.truth
	for name in design.names["exact"] + design.names["abundance"]:
		assert _barcode_at(default_run, reads[name]) == reads[name].barcode, reads[name]


def test_one_mismatch_is_corrected(design, default_run):
	"""A barcode one substitution from exactly one whitelist barcode is
	written as that barcode, also when that barcode was never seen exactly
	(design "unseen": a single candidate is taken whatever its abundance)."""
	reads = default_run.truth
	for name in design.names["one"]:
		assert _barcode_at(default_run, reads[name]) == design.source[name], reads[name]
	name, cands = design.designs["unseen"]
	assert _barcode_at(default_run, reads[name]) == cands[0]


def test_barcodes_beyond_one_mismatch_are_dropped(design, default_run):
	"""Barcodes two or three substitutions from the whitelist, and random
	ones, cannot be corrected at the default --bc-error-threshold 1: their
	reads are not written."""
	reads = default_run.truth
	for cat in ("two", "three", "far"):
		for name in design.names[cat]:
			assert _barcode_at(default_run, reads[name]) is None, (cat, reads[name])


def test_n_in_barcode(design, default_run):
	"""An N counts as a mismatch: one N in an otherwise whitelisted barcode
	is corrected to it; two Ns, or an N and a substitution, exceed one
	mismatch and the read is dropped."""
	reads = default_run.truth
	for name in design.names["n1"]:
		assert _barcode_at(default_run, reads[name]) == design.source[name], reads[name]
	for cat in ("n2", "n1sub"):
		for name in design.names[cat]:
			assert _barcode_at(default_run, reads[name]) is None, (cat, reads[name])


def test_ambiguous_barcodes_follow_quality_and_abundance(design, default_run):
	"""Barcodes one substitution from several whitelist barcodes (designs in
	build_design): corrected to candidate a when its probability is above 0.9
	(abundance 19:1, a low base quality under a's mismatch, the quality
	outweighing a 1:5 abundance); left uncorrected, so dropped, at a tie, at
	17:2 when both qualities clamp to the same value, and at 17:1:1, where the
	probability is the best score over the sum of all three."""
	reads = default_run.truth
	for label in ("abundance", "quality", "quality_over_abundance"):
		name, cands = design.designs[label]
		assert _barcode_at(default_run, reads[name]) == cands[0], label
	for label in ("tie", "low_clamp", "high_clamp", "three_way"):
		name, cands = design.designs[label]
		assert _barcode_at(default_run, reads[name]) is None, label


def test_barcode_statistics(design, keep_run):
	"""stderr reports "Number of barcodes in whitelist" (reads whose barcode
	is in the whitelist as read) and "Number of corrected barcodes", counted
	over every read before mapping (checked at -t 1, see keep_run)."""
	outcomes = hd.correction_outcomes(design.records, design.whitelist)
	status = collections.Counter(s for s, _ in outcomes.values())
	stats = keep_run.stats()
	assert stats["Number of barcodes in whitelist"] == status[hd.EXACT]
	assert stats["Number of corrected barcodes"] == status[hd.CORRECTED]
	assert status[hd.EXACT] == len(design.names["exact"]) + len(design.names["abundance"])


def test_fragments_are_sorted_by_reference_and_start(index, default_run):
	"""Fragment lines come per reference in FASTA order, by start within one."""
	keys = [(index.genome.rid(f.chrom), f.start) for f in default_run.fragments()]
	assert keys == sorted(keys)


# ------------------------------------------------------------ option runs

def test_output_mappings_not_in_whitelist(design, keep_run):
	"""--output-mappings-not-in-whitelist writes every read: corrected
	barcodes as with the default run, uncorrectable ones as read (an N written
	as A, the base its 2-bit key stores)."""
	want = hd.expected_fragments(design.records, design.whitelist, keep_uncorrectable=True)
	assert _lines(keep_run) == want
	reads = keep_run.truth
	for name in design.names["two"] + design.names["far"]:
		assert _barcode_at(keep_run, reads[name]) == reads[name].barcode
	for name in design.names["n2"]:
		assert _barcode_at(keep_run, reads[name]) == reads[name].barcode.replace("N", "A")


def test_bc_error_threshold_two(design, two_error_run):
	"""--bc-error-threshold 2 corrects barcodes two mismatches from one
	whitelist barcode (two substitutions, two Ns, an N and a substitution)
	and still drops three substitutions and random barcodes; every line
	follows the rule with up to 2 mismatches. A barcode with two Ns has one
	candidate, so it is corrected at the default probability 0.9 whatever
	the whitelist bases under the Ns are."""
	want = hd.expected_fragments(design.records, design.whitelist, max_errors=2)
	assert _lines(two_error_run) == want
	reads = two_error_run.truth
	for cat in ("two", "n2", "n1sub", "one", "n1"):
		for name in design.names[cat]:
			assert _barcode_at(two_error_run, reads[name]) == design.source[name], (cat, reads[name])
	for name in design.names["three"] + design.names["far"]:
		assert _barcode_at(two_error_run, reads[name]) is None


def test_bc_probability_threshold(design, low_probability_run):
	"""--bc-probability-threshold 0.45 corrects the designs whose best
	candidate has probability 0.895 (low_clamp, high_clamp, three_way) to
	candidate a, and the 0.5 tie to the candidate whose mismatch lies later in
	the barcode (candidates are ordered by score, then position: BarcodeWithQual,
	src/utils.h); every line follows the rule at 0.45."""
	want = hd.expected_fragments(design.records, design.whitelist, min_probability=0.45)
	assert _lines(low_probability_run) == want
	reads = low_probability_run.truth
	for label in ("low_clamp", "high_clamp", "three_way", "abundance", "quality"):
		name, cands = design.designs[label]
		assert _barcode_at(low_probability_run, reads[name]) == cands[0], label
	name, cands = design.designs["tie"]
	assert _barcode_at(low_probability_run, reads[name]) == cands[1]


# ------------------------------------------------------------ 5% check

SPARSE_LEN = 12


def _sparse_reads(genome, n_exact, label):
	"""100 reads with 12 bp barcodes: n_exact with whitelisted barcodes, the
	others one substitution from a whitelist barcode (one candidate each)."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_d_barcodes", label))
	rng = hd.stable_rng("test_d_barcodes", label)
	picker = hd.LocusPicker(sim)
	wl = ch.make_whitelist(20, length=SPARSE_LEN, seed=41, min_dist=4)
	recs = []
	for i in range(100):
		w = wl[i % len(wl)]
		bc = w if i < n_exact else hd.mutate(w, [rng.randrange(SPARSE_LEN)])
		recs.append(hd.barcoded_read(sim, picker, bc, hd.random_quality(rng, SPARSE_LEN)))
	return ch.ReadSet(recs), wl


@pytest.fixture(scope="module")
def sparse_whitelist(tmp_path_factory):
	wl = ch.make_whitelist(20, length=SPARSE_LEN, seed=41, min_dist=4)
	return ch.write_lines(tmp_path_factory.mktemp("d_sparse") / "whitelist.txt", wl)


def test_fewer_than_five_percent_whitelisted_exits(index, genome, sparse_whitelist):
	"""With 4 of 100 barcodes in the whitelist chromap stops before mapping
	with status 255 and its "Less than 5% barcodes" message, which names
	--skip-barcode-check."""
	reads, _ = _sparse_reads(genome, 4, "sparse4")
	run = ch.run_chromap(index, reads, ["--barcode-whitelist", str(sparse_whitelist)], check=False)
	assert run.returncode == 255
	assert "Less than 5% barcodes can be found or corrected based on the barcode whitelist." in run.stderr
	assert "--skip-barcode-check" in run.stderr


def test_five_percent_whitelisted_is_enough(index, genome, sparse_whitelist):
	"""5 of 100 barcodes in the whitelist is not fewer than 5%: the run
	maps, and the 95 one-substitution barcodes are corrected (one candidate
	each) and written as 12 bp barcodes."""
	reads, wl = _sparse_reads(genome, 5, "sparse5")
	run = ch.run_chromap(index, reads, ["--barcode-whitelist", str(sparse_whitelist)])
	assert _lines(run) == hd.expected_fragments(reads, wl)
	assert len(run.fragments()) == 100


def test_skip_barcode_check(index, genome, sparse_whitelist):
	"""--skip-barcode-check maps the 4-in-100 input anyway; corrections use
	the abundances of the 4 exact barcodes, and one-candidate barcodes are
	corrected whatever their abundance."""
	reads, wl = _sparse_reads(genome, 4, "sparse4")
	run = ch.run_chromap(index, reads, ["--barcode-whitelist", str(sparse_whitelist), "--skip-barcode-check"])
	assert _lines(run) == hd.expected_fragments(reads, wl)
	assert len(run.fragments()) == 100


# ------------------------------------------------------------ dedup

@pytest.fixture(scope="module")
def dedup_case(index, genome, tmp_path_factory):
	"""Single-end barcoded reads for cell-level deduplication, written with a
	translation table: 8 starts read three times with one whitelist barcode
	(exact + strand, exact - strand, one substitution), 6 starts read with two
	different whitelist barcodes, 4 starts read with an exact barcode and a
	two-substitution one, and 30 single reads (half exact, half one
	substitution). The table is gzipped and tab-separated, one line per 8 bp
	segment of the whitelist barcodes (README: translated barcode first,
	original second; segments of a longer barcode are translated in turn and
	joined with "-")."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_d_barcodes", "dedup_se"))
	rng = hd.stable_rng("test_d_barcodes", "dedup_se")
	picker = hd.LocusPicker(sim)
	wl = ch.make_whitelist(12, length=BC_LEN, seed=43, min_dist=4)
	recs = []
	for g in range(8):
		c, s = picker.take(50)
		w = wl[g % len(wl)]
		recs.append(sim.se_read(c, s, 50, "+", barcode=w))
		recs.append(sim.se_read(c, s, 50, "-", barcode=w))
		recs.append(sim.se_read(c, s, 50, "+", barcode=hd.mutate(w, [rng.randrange(BC_LEN)])))
	for g in range(6):
		c, s = picker.take(50)
		recs.append(sim.se_read(c, s, 50, "+", barcode=wl[g]))
		recs.append(sim.se_read(c, s, 50, "-", barcode=wl[g + 6]))
	for g in range(4):
		c, s = picker.take(50)
		w = wl[g + 3]
		recs.append(sim.se_read(c, s, 50, "-", barcode=w))
		recs.append(sim.se_read(c, s, 50, "-", barcode=hd.mutate(w, rng.sample(range(BC_LEN), 2))))
	for i in range(30):
		w = wl[i % len(wl)]
		recs.append(hd.barcoded_read(sim, picker, w if i % 2 else hd.mutate(w, [rng.randrange(BC_LEN)])))
	reads = sim.shuffled(ch.ReadSet(recs))
	d = tmp_path_factory.mktemp("d_dedup")
	path = ch.write_lines(d / "whitelist.txt", wl)
	segments = sorted(set(w[:8] for w in wl) | set(w[8:] for w in wl))
	ids = dict((x, "seg%02d" % i) for i, x in enumerate(segments))
	table = d / "translate.tsv.gz"
	with gzip.open(str(table), "wt") as f:
		for x in segments:
			f.write("%s\t%s\n" % (ids[x], x))
	run = ch.run_chromap(index, reads, ["--barcode-whitelist", str(path), "--remove-pcr-duplicates",
		"--remove-pcr-duplicates-at-cell-level", "--barcode-translate", str(table)], shared=True)
	return dict(run=run, reads=reads, wl=wl, translate=lambda b: "%s-%s" % (ids[b[:8]], ids[b[8:]]))


def test_cell_level_dedup_of_corrected_barcodes(dedup_case):
	"""--remove-pcr-duplicates --remove-pcr-duplicates-at-cell-level on
	single-end barcoded reads: reads with the same start and the same barcode
	after correction are one line whose count is the group size, whatever
	their strand; the same start with different barcodes stays separate
	lines; a read dropped for its barcode is not counted. Expected lines: the
	reads' corrected (contig, start, end, barcode) counted, the barcode
	translated."""
	want = hd.expected_fragments(dedup_case["reads"], dedup_case["wl"])
	tr = dedup_case["translate"]
	got = collections.Counter((f.chrom, f.start, f.end, f.barcode, f.dups) for f in dedup_case["run"].fragments())
	assert got == collections.Counter((c, s, e, tr(b), n) for (c, s, e, b), n in want.items())
	assert sorted(n for n in want.values() if n > 1) == [3] * 8


def test_barcode_translate_segments(dedup_case):
	"""--barcode-translate with 8 bp segments writes each 16 bp barcode as
	the IDs of its two halves joined with "-", never the sequence itself:
	the barcodes written are exactly the translations of the whitelist
	barcodes the reads were corrected to."""
	want = hd.expected_fragments(dedup_case["reads"], dedup_case["wl"])
	written = set(f.barcode for f in dedup_case["run"].fragments())
	assert written == set(dedup_case["translate"](k[3]) for k in want)
	assert not written & set(dedup_case["wl"])


@pytest.fixture(scope="module")
def pe_design(genome, tmp_path_factory):
	"""Paired-end barcoded fragments: 40 single fragments (20 exact, 14 one
	substitution, 6 two substitutions); 6 fragments read three times (exact +
	strand, exact - strand, one substitution); 4 fragments read with two
	different whitelisted barcodes; 4 starts read with two insert sizes and
	one barcode; 3 fragments read with an exact barcode and with a
	two-substitution one."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_d_barcodes", "pe_design"))
	rng = hd.stable_rng("test_d_barcodes", "pe_design")
	picker = hd.LocusPicker(sim)
	wl = ch.make_whitelist(16, length=BC_LEN, seed=47, min_dist=4)
	recs = []

	def pair(c, s, ins, strand, bc):
		return sim.pair(c, s, ins, strand=strand, barcode=bc, barcode_qual=hd.random_quality(rng, BC_LEN))

	for i in range(40):
		ins = rng.randint(150, 400)
		c, s = picker.take(ins)
		w = wl[i % len(wl)]
		bc = w if i < 20 else hd.mutate(w, rng.sample(range(BC_LEN), 1 if i < 34 else 2))
		recs.append(pair(c, s, ins, "+-"[i % 2], bc))
	for g in range(6):
		ins = rng.randint(150, 400)
		c, s = picker.take(ins)
		w = wl[g]
		recs += [pair(c, s, ins, "+", w), pair(c, s, ins, "-", w),
			pair(c, s, ins, "+", hd.mutate(w, [rng.randrange(BC_LEN)]))]
	for g in range(4):
		c, s = picker.take(300)
		recs += [pair(c, s, 300, "+", wl[g]), pair(c, s, 300, "-", wl[g + 8])]
	for g in range(4):
		c, s = picker.take(337)
		recs += [pair(c, s, 300, "+", wl[g + 4]), pair(c, s, 337, "+", wl[g + 4])]
	for g in range(3):
		c, s = picker.take(250)
		w = wl[g + 10]
		recs += [pair(c, s, 250, "-", w), pair(c, s, 250, "+", hd.mutate(w, rng.sample(range(BC_LEN), 2)))]
	reads = sim.shuffled(ch.ReadSet(recs))
	path = ch.write_lines(tmp_path_factory.mktemp("d_pe_whitelist") / "whitelist.txt", wl)
	return reads, wl, path


@pytest.fixture(scope="module")
def pe_run(index, pe_design):
	reads, wl, path = pe_design
	return ch.run_chromap(index, reads, ["--barcode-whitelist", str(path), "--remove-pcr-duplicates",
		"--remove-pcr-duplicates-at-cell-level", "--cache-size", "2000000"], shared=True)


def test_paired_end_correction_and_cell_level_dedup(pe_design, pe_run):
	"""Paired-end fragments are written as chrom, start, end, corrected
	barcode and count: pairs of one fragment (same start and end) whose
	barcodes agree after correction are one line, whatever read 1's strand;
	the same fragment with two barcodes, or the same start with two insert
	sizes, gives two lines; pairs with uncorrectable barcodes are dropped."""
	reads, wl, _ = pe_design
	want = hd.expected_fragments(reads, wl)
	got = collections.Counter((f.chrom, f.start, f.end, f.barcode, f.dups) for f in pe_run.fragments())
	assert got == collections.Counter(k + (n,) for k, n in want.items())
	assert sorted(n for n in want.values() if n > 1) == [3] * 6


@pytest.fixture(scope="module")
def pe_keep_run(index, pe_design):
	reads, wl, path = pe_design
	"""At -t 1, for the barcode statistics (see keep_run)."""
	return ch.run_chromap(index, reads, ["--barcode-whitelist", str(path), "--remove-pcr-duplicates",
		"--remove-pcr-duplicates-at-cell-level", "--output-mappings-not-in-whitelist", "--cache-size", "2000000"],
		threads=1, shared=True)


def test_paired_end_output_mappings_not_in_whitelist(pe_design, pe_keep_run):
	"""Paired-end --output-mappings-not-in-whitelist also writes the pairs
	whose barcodes cannot be corrected, with the barcode as read; they are
	deduplicated like the others (a fragment read with an exact barcode and
	with a two-substitution one gives two lines)."""
	reads, wl, _ = pe_design
	want = hd.expected_fragments(reads, wl, keep_uncorrectable=True)
	got = collections.Counter((f.chrom, f.start, f.end, f.barcode, f.dups) for f in pe_keep_run.fragments())
	assert got == collections.Counter(k + (n,) for k, n in want.items())
	assert sum(want.values()) == len(reads)


def test_paired_end_barcode_statistics(pe_design, pe_keep_run):
	"""Paired-end runs count each pair's barcode once in "Number of barcodes
	in whitelist" and "Number of corrected barcodes" (checked at -t 1, see
	keep_run)."""
	reads, wl, _ = pe_design
	status = collections.Counter(s for s, _ in hd.correction_outcomes(reads, wl).values())
	stats = pe_keep_run.stats()
	assert stats["Number of barcodes in whitelist"] == status[hd.EXACT]
	assert stats["Number of corrected barcodes"] == status[hd.CORRECTED]


# ------------------------------------------------------------ errors

@pytest.fixture(scope="module")
def small_inputs(genome, tmp_path_factory):
	"""Read and barcode files for the barcode error paths (single-end, which
	stops before allocating the mapping cache)."""
	d = tmp_path_factory.mktemp("d_barcode_errors")
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_d_barcodes", "errors"))
	wl = ch.make_whitelist(6, length=BC_LEN, seed=59)
	reads = sim.se_reads(30, barcodes=wl)
	files = reads.write(d, prefix="ok")
	uneven = list(reads)
	uneven[17] = sim.se_read(uneven[17].contig, uneven[17].start, barcode=wl[0][:15], name=uneven[17].name)
	files_uneven = ch.ReadSet(uneven).write(d, prefix="uneven")
	long_bc = ch.make_whitelist(2, length=33, seed=61)
	long_reads = ch.ReadSet([sim.se_read(r.contig, r.start, barcode=long_bc[i % 2], name=r.name)
		for i, r in enumerate(reads)])
	files_long = long_reads.write(d, prefix="long")
	return dict(dir=d, wl=wl, files=files, uneven=files_uneven, long=files_long, long_bc=long_bc)


def _error_run(index, files, extra):
	return ch.run_chromap(index, {"r1": files["r1"], "bc": files["bc"]}, extra, check=False)


def test_unequal_barcode_lengths_in_input_exit(index, small_inputs):
	"""Barcodes of different lengths among the input's first 1000 stop
	chromap (status 255) with its message."""
	run = _error_run(index, small_inputs["uneven"], [])
	assert run.returncode == 255
	assert "ERROR: barcode lengths are not equal in the sample!" in run.stderr


def test_whitelist_length_differs_from_reads_exit(index, small_inputs):
	"""A whitelist of 15 bp barcodes for 16 bp input barcodes stops chromap."""
	path = ch.write_lines(small_inputs["dir"] / "wl15.txt", [w[:15] for w in small_inputs["wl"]])
	run = _error_run(index, small_inputs["files"], ["--barcode-whitelist", str(path)])
	assert run.returncode == 255
	assert "ERROR: whitelist and input barcode lengths are not equal!" in run.stderr


def test_unequal_whitelist_lengths_exit(index, small_inputs):
	"""A whitelist whose barcodes differ in length (after a first barcode of
	the input's length) stops chromap."""
	wl = list(small_inputs["wl"])
	wl[3] = wl[3][:14]
	path = ch.write_lines(small_inputs["dir"] / "wl_mixed.txt", wl)
	run = _error_run(index, small_inputs["files"], ["--barcode-whitelist", str(path)])
	assert run.returncode == 255
	assert "ERROR: barcode lengths are not equal in the whitelist!" in run.stderr


def test_whitelist_barcodes_of_32_bases_are_accepted(index, genome, tmp_path):
	"""32 bases is the longest barcode a 2-bit key holds (the check is
	"greater than 32"): 32 bp whitelisted barcodes, and barcodes one
	substitution from one of them, are written as 32 bp whitelist
	barcodes."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_d_barcodes", "bc32"))
	rng = hd.stable_rng("test_d_barcodes", "bc32")
	picker = hd.LocusPicker(sim)
	wl = ch.make_whitelist(6, length=32, seed=73, min_dist=4)
	recs = [hd.barcoded_read(sim, picker, w if i < 12 else hd.mutate(w, [rng.randrange(32)]))
		for i, w in enumerate(wl * 3)]
	reads = ch.ReadSet(recs)
	path = ch.write_lines(tmp_path / "wl32.txt", wl)
	run = ch.run_chromap(index, reads, ["--barcode-whitelist", str(path)])
	assert _lines(run) == hd.expected_fragments(reads, wl)
	assert len(run.fragments()) == len(recs)
	assert set(f.barcode for f in run.fragments()) == set(wl)


def test_whitelist_barcodes_over_32_bases_exit(index, small_inputs):
	"""Barcodes are stored as 2-bit keys of at most 32 bases: a 33 bp
	whitelist stops chromap."""
	path = ch.write_lines(small_inputs["dir"] / "wl33.txt", small_inputs["long_bc"])
	run = _error_run(index, small_inputs["long"], ["--barcode-whitelist", str(path)])
	assert run.returncode == 255
	assert "ERROR: barcode length is greater than 32!" in run.stderr


def test_missing_whitelist_file_exits(index, small_inputs):
	"""A whitelist path that does not exist stops chromap with its message."""
	run = _error_run(index, small_inputs["files"], ["--barcode-whitelist",
		str(small_inputs["dir"] / "no_such_whitelist.txt")])
	assert run.returncode == 255
	assert "ERROR: barcode whitelist file does not exist or is truncated!" in run.stderr
