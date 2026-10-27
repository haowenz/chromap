"""Area B, C++ layer: the mapping core's public entry points driven
directly (tests/cpp/b_units.cc, b_probes.cc, b_mmcache.cc).

CT_TESTs check values worked by hand in the driver. The probes print what
ksw_semi_global3, DraftMappingGenerator::GenerateDraftMappings,
CandidateProcessor::GenerateCandidates, Index::GenerateCandidatePositions,
Index::GenerateCandidatePositionsFromRepetitiveReadWithMateInfoOnOneStrand and
CandidateProcessor::ReduceCandidatesForPairedEndRead return, and the expected
values come from the definitions in helpers_b.py (banded edit distance,
affine-gap semi-global alignment, minimizer occurrences of a genome)."""

import collections
import random

import pytest

import chromap_helpers as ch
import helpers_b as hb

pytestmark = pytest.mark.cpp

B_UNITS_TESTS = [
	"candidate_layout_and_order",
	"generating_config_thresholds",
	"draft_mapping_layout",
	"mapping_in_memory_fragments",
	"mapping_parameters_defaults",
	"mapping_metadata_counters_and_reset",
	"mapping_metadata_sort_move_update",
	"mapping_metadata_print_candidates",
	"paired_metadata_best_mappings_by_strand",
	"paired_metadata_acts_on_both_mates",
	"ksw_semi_global_wrappers_agree",
]

B_MMCACHE_TESTS = [
	"cache_stores_after_ten_updates",
	"cache_query_moves_candidates_to_query_positions",
	"cache_query_reverse_complement",
	"cache_needs_a_dominant_fingerprint",
	"cache_replacement_and_saturation",
	"cache_skips_candidates_near_reference_start",
	"cache_update_threshold",
]


def test_b_drivers_list_every_test(cpp):
	"""Every CT_TEST of the area B drivers is run below, and the probes the
	tests use exist."""
	assert cpp("b_units").tests() == B_UNITS_TESTS
	assert cpp("b_mmcache").tests() == B_MMCACHE_TESTS
	assert cpp("b_mmcache").probes() == ["printstats"]
	assert cpp("b_probes").tests() == []
	assert cpp("b_probes").probes() == ["ksw", "draft", "candidates", "positions", "matewindow", "reduce"]


@pytest.mark.parametrize("name", B_UNITS_TESTS)
def test_b_units_cpp(cpp, name):
	cpp("b_units").check(name)


@pytest.mark.parametrize("name", B_MMCACHE_TESTS)
def test_b_mmcache_cpp(cpp, name):
	cpp("b_mmcache").check(name)


def test_mm_cache_print_stats(cpp):
	"""mm_cache::PrintStats (the --debug-cache dump) prints one line per
	slot, and Update with debug on prints its debug lines when it stores an
	entry. The layout of these debugging lines is not asserted."""
	out = cpp("b_mmcache").probe("printstats", 5).splitlines()
	debug = [l for l in out if l.startswith("[DEBUG]")]
	stats = [l for l in out if not l.startswith("[DEBUG]")]
	assert debug
	assert len(stats) == 5


# ------------------------------------------------------------- ksw

def _mutate(rng, seq, n_sub, n_ins, n_del):
	s = list(seq)
	for _ in range(n_sub):
		i = rng.randrange(len(s))
		s[i] = ch.substitute(s[i])
	for _ in range(n_del):
		i = rng.randrange(5, len(s) - 5)
		del s[i:i + rng.randint(1, 3)]
	for _ in range(n_ins):
		i = rng.randrange(5, len(s) - 5)
		s[i:i] = [rng.choice("ACGT") for _ in range(rng.randint(1, 3))]
	return "".join(s)


def _ksw(cpp, window, read, p, w, cigar=True):
	args = [window, read, p.match, p.mismatch, p.o_ins, p.e_ins, p.o_del, p.e_del, w]
	if not cigar:
		args.append("nocigar")
	f = cpp("b_probes").probe("ksw", *args).split()
	return int(f[0]), int(f[1]), int(f[2]), (f[3] if len(f) > 3 else "")


def test_ksw_semi_global3_matches_affine_reference(cpp):
	"""ksw_semi_global3 (the SAM CIGAR of every mapping) aligns the whole read
	inside a reference window of read length + 2e with band w = 2e + 1. Its
	score equals the best affine-gap semi-global score (helpers_b.semi_global:
	match +A, mismatch -B, read gap o + l*e, reference gap o + l*e, gaps open
	only after an aligned base, free window ends, band) on random reads with
	substitutions, insertions and deletions under several scoring schemes;
	the returned CIGAR, start and end describe an alignment with exactly that
	score that consumes the whole read and window[start:end]; and where the
	optimum is unique, it is that alignment. Without CIGAR output the score
	and end are the same. (ksw names the read-gap pair o_del/e_del and the
	window-gap pair o_ins/e_ins, as the probe's argument order shows.)"""
	rng = random.Random(ch.stable_seed("b_ksw"))
	schemes = [hb.DEFAULT_SCORES, hb.ScoreParams(2, 6, 10, 2, 10, 2), hb.ScoreParams(1, 1, 3, 1, 8, 2),
		hb.ScoreParams(1, 9, 2, 1, 2, 1), hb.ScoreParams(3, 2, 6, 3, 4, 1)]
	unique = 0
	for case in range(60):
		e = rng.choice([2, 4, 8])
		L = rng.choice([30, 50])
		w = 2 * e + 1
		core = "".join(rng.choice("ACGT") for _ in range(L + 6))
		read = _mutate(rng, core, rng.randint(0, 3), rng.randint(0, 1), rng.randint(0, 1))[:L]
		lead = rng.randint(0, 2 * e)
		window = ("".join(rng.choice("ACGT") for _ in range(lead)) + core)[:L + 2 * e]
		window += "".join(rng.choice("ACGT") for _ in range(L + 2 * e - len(window)))
		p = schemes[case % len(schemes)]
		ref = hb.semi_global(window, read, p, w)
		score, start, end, cigar = _ksw(cpp, window, read, p, w)
		ctx = (case, window, read, p, ref, (score, start, end, cigar))
		assert score == ref["score"], ctx
		assert hb.cigar_spans(cigar) == (L, end - start), ctx
		assert 0 <= start <= w and end <= len(window), ctx
		assert hb.cigar_score(window, start, read, cigar, p) == score, ctx
		if ref["count"] == 1:
			unique += 1
			assert (start, end, cigar) == ref["alignments"][0], ctx
		score2, _, end2, cigar2 = _ksw(cpp, window, read, p, w, cigar=False)
		assert (score2, end2, cigar2) == (score, end, ""), ctx
	assert unique >= 30


def test_ksw_gap_costs_choose_between_gap_and_mismatches(cpp):
	"""Hand-worked: a 20 bp read that skips one window base after position 17
	(a 1 bp deletion three bases from its end). With -B 4 -O 6 -E 1 the
	deletion costs 7 against 3 mismatches costing 3 * (1 + 4) = 15, so
	17M1D3M (score 20 - 7 = 13) wins; with a reference-gap open of 30 the
	three mismatches win (20M, score 17 - 12 = 5). Read gaps (CIGAR I) are
	priced by the other pair: a read with one extra base gets 1I with an
	insertion open of 6 and mismatches with 30."""
	window = "GG" + "ACGTCAGTTGCATCGATGCA" + "TTTT"
	read = "ACGTCAGTTGCATCGAT" + "CAT"
	assert window[19:22] == "GCA"
	cheap = hb.ScoreParams(1, 4, 6, 1, 6, 1)
	dear = hb.ScoreParams(1, 4, 6, 1, 30, 1)
	assert _ksw(cpp, window[:24], read, cheap, 5) == (13, 2, 23, "17M1D3M")
	assert hb.semi_global(window[:24], read, cheap, 5)["alignments"] == [(2, 23, "17M1D3M")]
	got = _ksw(cpp, window[:24], read, dear, 5)
	assert got[0] == hb.semi_global(window[:24], read, dear, 5)["score"] == 20 - 3 * 5
	assert got[3] == "20M"
	ins_read = "ACGTCAGTTGCATCGATG" + "T" + "C"
	ins_window = "GG" + "ACGTCAGTTGCATCGATGCAAA"
	cheap_i = hb.ScoreParams(1, 4, 6, 1, 6, 1)
	dear_i = hb.ScoreParams(1, 4, 30, 1, 6, 1)
	a = _ksw(cpp, ins_window, ins_read, cheap_i, 5)
	b = _ksw(cpp, ins_window, ins_read, dear_i, 5)
	assert "I" in a[3] and a[0] == hb.semi_global(ins_window, ins_read, cheap_i, 5)["score"]
	assert "I" not in b[3] and b[0] == hb.semi_global(ins_window, ins_read, dear_i, 5)["score"]


# ------------------------------------------------------- draft mappings

def _draft_probe(cpp, genome, reads, cases, e, tmp_path, split=0):
	"""cases: [(read index, number of minimizers, [(strand, contig, position,
	count)])]. Returns per case (S tuple, [mappings]) as printed."""
	fa = tmp_path / ("reads_e%d.fa" % e)
	ch.write_lines(fa, [x for r in reads for x in (">" + r.name, r.seq)])
	lines = []
	for idx, nmin, cands in cases:
		lines.append("R %d %d" % (idx, nmin))
		for strand, contig, pos, count in cands:
			lines.append("C %s %d %d %d" % (strand, genome.rid(contig), pos, count))
		lines.append("E")
	out = cpp("b_probes").probe("draft", genome.fasta_path, fa, e, split, input="\n".join(lines) + "\n")
	results, cur = [], None
	for line in out.splitlines():
		f = line.split()
		if f[0] == "S":
			cur = (tuple(int(x) for x in f[1:]), [], [])
		elif f[0] == "M":
			cur[1].append((f[1], int(f[2]), int(f[3]), int(f[4])))
		elif f[0] == "X":
			cur[2].append((f[1], int(f[2])))
		elif f[0] == "END":
			results.append(cur)
	assert len(results) == len(cases)
	return results


def _summary(errors, e):
	"""(min, second min, #best, #second best) of the mapped error counts,
	e + 1 standing for none."""
	vals = sorted(errors)
	if not vals:
		return (e + 1, e + 1, 0, 0)
	distinct = sorted(set(vals))
	second = distinct[1] if len(distinct) > 1 else e + 1
	return (distinct[0], second, vals.count(distinct[0]), vals.count(second) if len(distinct) > 1 else 0)


@pytest.fixture(scope="module")
def draft_reads(genome):
	"""Reads of 50 and 100 bp from unique loci on both strands with 0-12
	substitutions and indels, each with its candidate positions: the true one,
	the true one shifted by up to 3 (still inside the band), and decoys at
	random unique loci."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("b_units", "draft_reads"))
	rng = sim.rng
	reads = []
	for i in range(24):
		L = 50 if i % 3 else 100
		strand = "+-"[i % 2]
		contig, start = sim.unique_locus(L + 10)
		kind = i % 4
		name = "d%03d" % i
		if kind == 3:
			r = hb.indel_read(genome, contig, start, L, strand, "ID"[i % 2 == 0], L // 2 + rng.randint(-5, 5),
				rng.randint(1, 3), name)
		else:
			nsub = rng.randint(0, 12)
			mm = rng.sample(range(L), nsub)
			r = sim.se_read(contig, start, L, strand, mismatches=mm, name=name)
		reads.append(r)
	return reads


@pytest.mark.parametrize("e", [3, 7, 8, 12, 15])
def test_draft_mappings_match_banded_edit_distance(cpp, genome, draft_reads, tmp_path, e):
	"""GenerateDraftMappings verifies each candidate by the banded edit
	distance of the read against the reference window [p - e, p + L + e)
	(helpers_b.draft_mapping) and keeps those with at most e errors, with the
	end position of the alignment. -e picks the code path (MappingParameters
	::GetNumVPULanes): below 8, eight 16-bit SIMD lanes; 8 to 15, four 32-bit
	lanes; a strand with fewer candidates than lanes is verified one by one.
	Each read gets 3 candidates on one strand (scalar) and 9 or 13 on the
	other (SIMD lanes plus a scalar remainder), all with the same support so
	no candidate is skipped; every path reports the same mappings, and the
	minimum and second minimum error counts with their multiplicities.
	Candidates at the first and last positions whose window fits on the
	contig with a base to spare are included."""
	rng = random.Random(ch.stable_seed("b_units", "draft", e))
	cases, expected = [], []
	for idx, r in enumerate(draft_reads):
		L = len(r.seq)
		contig_seq = genome.contigs[r.contig]
		true_c = r.start if r.strand == "+" else r.start + L - 1
		n = len(contig_seq)

		def bounds(contig, strand):
			# first and last candidate whose window [p - e, p + L + e] fits
			m = len(genome.contigs[contig])
			lo, hi = e, m - L - e - 1
			return (lo, hi) if strand == "+" else (lo + L - 1, hi + L - 1)

		same = [(r.contig, true_c + shift) for shift in (0, -1, 1, 2, -3)]
		same += [(r.contig, b) for b in bounds(r.contig, r.strand)]
		while len(same) < (9 if idx % 2 else 13):
			c = rng.choice(genome.names)
			lo, hi = bounds(c, r.strand)
			same.append((c, rng.randint(lo, hi)))
		other_strand = "-" if r.strand == "+" else "+"
		olo, ohi = bounds(r.contig, other_strand)
		other = [(r.contig, p) for p in (olo, ohi, rng.randint(olo, ohi))]
		cands = [(r.strand, c, p, 2) for c, p in same] + [(other_strand, c, p, 2) for c, p in other]
		cases.append((idx, 30, cands))
		exp = []
		for strand, contig, p, _ in cands:
			res = hb.draft_mapping(genome.contigs[contig], r.seq, strand, p, e)
			assert res is not None
			if res[0] <= e:
				exp.append((strand, res[0], genome.rid(contig), res[1]))
		expected.append(exp)
	results = _draft_probe(cpp, genome, draft_reads, cases, e, tmp_path)
	mapped = 0
	for (summary, got, _), exp, r in zip(results, expected, draft_reads):
		assert sorted(got) == sorted(exp), r
		assert summary == _summary([x[1] for x in exp], e), r
		mapped += bool(exp)
	assert mapped >= 10


def test_draft_mapping_trusts_a_candidate_supported_by_every_minimizer(cpp, genome, tmp_path):
	"""An exact read with exactly one candidate whose support equals its
	number of minimizers is accepted without verification
	(GenerateNonSplitDraftMappingSupportedByAllMinimizers): one mapping with 0
	errors ending at p + L - 1 for a forward candidate (p = start) or at p for
	a reverse one (p = end), summary (0, e + 1, 1, 0). With one supporting
	minimizer fewer, the same candidate is verified and gives the same
	mapping."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("b_units", "direct"))
	reads, cases, expected = [], [], []
	for i in range(6):
		contig, start = sim.unique_locus(50)
		strand = "+-"[i % 2]
		r = sim.se_read(contig, start, 50, strand, name="x%d" % i)
		reads.append(r)
		cand = start if strand == "+" else start + 49
		nmin = 9
		cases.append((i, nmin, [(strand, contig, cand, nmin if i < 4 else nmin - 1)]))
		expected.append(((0, 9, 1, 0), [(strand, 0, genome.rid(contig), start + 49)]))
	results = _draft_probe(cpp, genome, reads, cases, 8, tmp_path)
	for (summary, got, _), (s, m) in zip(results, expected):
		assert (summary, got) == (s, m)


# ------------------------------------------------- candidate generation

@pytest.fixture(scope="module")
def small(index_factory, tmp_path_factory):
	"""The small N-free genome of helpers_b (rep6: six copies on both
	strands, rep2: two copies), its index and its minimizer table."""
	g = hb.make_small_genome().write(tmp_path_factory.mktemp("small_genome"))
	idx = index_factory(genome=g)
	return g, idx, hb.minimizer_table(g)


@pytest.fixture(scope="module")
def small_reads(small):
	"""Reads on the small genome: unique ones (exact, with substitutions,
	with deletions of 2, 3 and 8 bases), ones inside rep6 and rep2, and ones
	straddling a rep6 copy's edge, on both strands."""
	g, _, _ = small
	sim = ch.ReadSimulator(g, seed=ch.stable_seed("b_units", "small_reads"))
	reads = []
	for i in range(8):
		c, s = sim.unique_locus(60, margin=60)
		reads.append(sim.se_read(c, s, 50, "+-"[i % 2], mismatches=[(), (20,), (10, 35)][i % 3]))
	for i in range(6):
		reads.append(sim.segment_read("rep6", 30 + 37 * i, 50, strand="+-"[i % 2], copy=i))
	for i in range(2):
		reads.append(sim.segment_read("rep2", 60 + 50 * i, 50, strand="+-"[i % 2], copy=i))
	for i, (c, s) in enumerate([("s1", 475), ("s1", 775), ("s2", 375), ("s3", 1175)]):
		reads.append(sim.se_read(c, s, 50, "+-"[i % 2]))
	# One minimizer at the copy's edge (2 or 3 occurrences) and the rest
	# inside rep6 (6): with -f 4 the read has one non-repetitive minimizer,
	# so one seed is enough for a candidate.
	reads.append(sim.se_read("s1", 1494, 50, "+", name="onerep_a"))
	reads.append(sim.se_read("s1", 493, 50, "-", name="onerep_b"))
	reads.append(hb.indel_read(g, "s3", 1500, 50, "+", "D", 25, 2, "del"))
	# Deletions of exactly 3 and 8 bases put the seeds on two diagonals
	# exactly e apart for the e = 3 and e = 8 cases: still one candidate.
	reads.append(hb.indel_read(g, "s1", 1200, 50, "+", "D", 24, 3, "del3"))
	reads.append(hb.indel_read(g, "s1", 1000, 50, "-", "D", 24, 8, "del8"))
	return reads


def _write_reads(tmp_path, name, reads):
	fa = tmp_path / name
	ch.write_lines(fa, [x for r in reads for x in (">" + r.name, r.seq)])
	return fa


def _parse_blocks(out):
	blocks = []
	for line in out.splitlines():
		f = line.split()
		if f[0] == "R":
			blocks.append(dict(head=[int(x) for x in f[1:]], P=[], N=[], HP=[], HN=[]))
		else:
			blocks[-1][f[0]].append(tuple(int(x) for x in f[1:]))
	return blocks


@pytest.mark.parametrize("max_freq,rep_freq,heap", [(500, 500, 0), (500, 500, 1), (6, 3, 0), (6, 3, 1), (2, 2, 1),
	(7, 6, 0)])
def test_candidate_positions_from_minimizer_occurrences(cpp, small, small_reads, tmp_path, max_freq, rep_freq, heap):
	"""Index::GenerateCandidatePositions turns every reference occurrence of
	each read minimizer into a candidate position: reference position minus
	read position when the two k-mers are on the same strand (forward
	candidates, the read's start), reference position plus read position
	- k + 1 otherwise (reverse candidates, the read's end), sorted. Minimizers
	occurring at least max_freq times are skipped unless they occur once;
	those occurring at least rep_freq times (and more than once) are counted
	as repetitive, and the repetitive seed length adds k for the first and
	min(k, gap) for each next one within k + w - 1. The heap-merge path gives
	the same hits as the sort path. Occurrences come from the genome's
	minimizers by definition (helpers_b.minimizer_table)."""
	g, idx, table = small
	fa = _write_reads(tmp_path, "reads.fa", small_reads)
	out = cpp("b_probes").probe("positions", idx.path, fa, max_freq, rep_freq, heap)
	blocks = _parse_blocks(out)
	assert len(blocks) == len(small_reads)
	saw_repetitive = False
	for r, b in zip(small_reads, blocks):
		mins = hb.read_minimizers(r.seq)
		pos, neg, nrep, rlen = hb.candidate_hits(table, mins, max_freq, rep_freq)
		assert [h for (h,) in b["HP"]] == pos, r
		assert [h for (h,) in b["HN"]] == neg, r
		assert b["head"][1:] == [nrep, rlen], r
		saw_repetitive = saw_repetitive or nrep > 0
	assert saw_repetitive == (rep_freq <= 6)


@pytest.mark.parametrize("e,s,f0,f1", [(8, 2, 500, 1000), (8, 1, 500, 1000), (8, 5, 500, 1000), (3, 2, 4, 1000),
	(8, 2, 4, 6), (8, 2, 2, 3)])
def test_generate_candidates_groups_hits(cpp, small, small_reads, tmp_path, e, s, f0, f1):
	"""CandidateProcessor::GenerateCandidates groups the sorted hits of each
	strand into candidates (helpers_b.cluster_hits): hits within e of each
	other on one reference form a candidate at their most frequent exact
	position with that frequency as support; a candidate needs min(-s,
	#minimizers - #repetitive minimizers) hits (at least 1). When every
	minimizer is too frequent for the first -f value, a second round uses the
	second value (helpers_b.generate_candidates). Unique exact reads give one
	candidate at the true start (forward) or end (reverse) supported by every
	minimizer; rep6 reads give one per copy."""
	g, idx, table = small
	fa = _write_reads(tmp_path, "reads.fa", small_reads)
	out = cpp("b_probes").probe("candidates", idx.path, fa, e, s, f0, f1)
	blocks = _parse_blocks(out)
	for r, b in zip(small_reads, blocks):
		mins = hb.read_minimizers(r.seq)
		pos, neg, rlen, _ = hb.generate_candidates(table, mins, e, s, f0, f1)
		assert b["head"][1:] == [len(mins), rlen], r
		assert b["P"] == pos, r
		assert b["N"] == neg, r
		unique = r.kind == "genomic" and not g.overlaps_feature(r.contig, r.start, r.end)
		if unique and not r.mismatches and (e, s, f0) == (8, 2, 500):
			rid = g.rid(r.contig)
			want = (rid << 32 | r.start) if r.strand == "+" else (rid << 32 | (r.end - 1))
			got = b["P"] if r.strand == "+" else b["N"]
			assert got == [(want, len(mins))], r
		if r.extra.get("segment") == "rep6" and (e, s, f0, f1) == (8, 2, 500, 1000):
			assert len(b["P"]) + len(b["N"]) == 6, r
		if r.name.startswith("onerep") and f0 == 4:
			# min(-s, #minimizers - #repetitive) = 1 seed: the read keeps a
			# candidate supported by its one non-repetitive minimizer.
			nrep = sum(1 for h, _, _ in mins if len(table.get(h, [])) >= f0)
			assert len(mins) - nrep == 1 and b["P"] + b["N"], (r, b)
		if r.name == "del%d" % e and s == 2:
			# hits on diagonals exactly e apart form one candidate
			assert len(b["P"]) + len(b["N"]) == 1, (r, b)


def _mate_window(cpp, idx, fa, min_seeds, f0, e, queries):
	"""queries: [(read index, strand, range, [(position, count)])] ->
	[(return value, repetitive length, [positions])]."""
	text = "".join("%d %s %d %s\n" % (ri, st, rng_, " ".join("%d:%d" % m for m in mate))
		for ri, st, rng_, mate in queries)
	out = cpp("b_probes").probe("matewindow", idx.path, fa, min_seeds, f0, e, input=text)
	results, cur = [], None
	for line in out.splitlines():
		f = line.split()
		if f[0] == "RET":
			cur = (int(f[1]), int(f[2]), [])
		elif f[0] == "H":
			cur[2].append(int(f[1]))
		else:
			results.append(cur)
	assert len(results) == len(queries)
	return results


def test_mate_window_search_finds_repetitive_read_near_its_mate(cpp, small, small_reads, tmp_path):
	"""Index::GenerateCandidatePositionsFromRepetitiveReadWithMateInfoOnOneStrand
	looks up a repetitive read's minimizers only near the mate's best
	candidates: occurrences whose reference position lies within the search
	range of a best mate candidate (the highest support), on the requested
	strand; a minimizer occurring once is taken wherever it is. The binary
	search may also return the one occurrence just below a window. It returns
	the best mate support and the repetitive seed length (minimizers occurring
	at least F0 times). It gives up, returning minus the best support and no
	positions, when the mate has more than F0 candidates or 300 best ones, or
	at most -s support with 200 best ones."""
	g, idx, table = small
	reads = [r for r in small_reads if r.extra.get("segment") == "rep6"]
	fa = _write_reads(tmp_path, "rep.fa", reads)
	copies = g.segments["rep6"].copies
	queries, expected = [], []
	for ri, r in enumerate(reads):
		mins = hb.read_minimizers(r.seq)
		rep = [p for h, p, s in mins if len(table.get(h, [])) >= 4]
		for ci in (0, 3, 5):
			c = copies[ci]
			mate_pos = (g.rid(c.contig) << 32) | (c.start + 150)
			mate = [(mate_pos - 400, 2), (mate_pos, 9), (mate_pos + 70000, 3)]
			for strand in "+-":
				for rng_ in (300, 1200):
					lo, hi = max(0, mate_pos - rng_), mate_pos + rng_
					required, allowed = set(), set()
					for h, p, s in mins:
						occ = table.get(h, [])
						for o in occ:
							key = (o[0] << 32) | o[1]
							st, cand = hb.candidate_position(o, p, s)
							if st == strand and (len(occ) == 1 or lo <= key <= hi):
								required.add(cand)
						below = [o for o in occ if ((o[0] << 32) | o[1]) < lo]
						for o in below[-1:]:
							st, cand = hb.candidate_position(o, p, s)
							if st == strand:
								allowed.add(cand)
					queries.append((ri, strand, rng_, mate))
					expected.append((r, ci, strand, rng_, required, allowed, hb.repetitive_seed_length(rep)))
	results = _mate_window(cpp, idx, fa, 2, 4, 8, queries)
	found = collections.Counter()
	for (ret, rlen, got), (r, ci, strand, rng_, required, allowed, want_rlen) in zip(results, expected):
		ctx = (r, strand, ci, rng_)
		assert ret == 9, ctx
		assert rlen == want_rlen, ctx
		assert got == sorted(got), ctx
		assert required <= set(got) <= required | allowed, ctx
		found[(r.name, ci)] += len(required)
	# The copy near the mate is found on one of the two strands.
	assert len(found) == 3 * len(reads) and min(found.values()) > 0
	# Giving up: more mate candidates than F0, or many weakly supported ones.
	many = [(1000 + 10 * i, 3) for i in range(5)]
	weak = [(1000 + 10 * i, 2) for i in range(200)]
	assert [(ret, got) for ret, _, got in _mate_window(cpp, idx, fa, 2, 4, 8, [(0, "+", 300, many)])] == [(-3, [])]
	assert [(ret, got) for ret, _, got in _mate_window(cpp, idx, fa, 2, 500, 8, [(0, "+", 300, weak)])] == [(-2, [])]
	assert _mate_window(cpp, idx, fa, 1, 500, 8, [(0, "+", 300, weak)])[0][0] == 2


def test_reduce_candidates_keeps_pairs_within_insert_distance(cpp):
	"""ReduceCandidatesForPairedEndRead keeps a read-1 forward candidate when
	a read-2 reverse candidate lies within the distance of it (and the other
	way round), and likewise read-1 reverse with read-2 forward; with support
	below 6 nothing unpaired is kept. Random sorted candidate lists on two
	references are compared with that rule, list order kept."""
	rng = random.Random(ch.stable_seed("b_units", "reduce"))
	for case in range(40):
		d = rng.choice([50, 300, 1000])
		lists = {}
		for tag in ("p1", "n1", "p2", "n2"):
			vals = sorted({(rng.randint(0, 1) << 32) | rng.randint(0, 6000) for _ in range(rng.randint(0, 8))})
			lists[tag] = [(v, rng.randint(1, 5)) for v in vals]
		text = "\n".join(tag + "".join(" %d:%d" % c for c in lists[tag]) for tag in lists) + "\n"
		out = cpp("b_probes").probe("reduce", d, input=text)
		got = {}
		for line in out.splitlines():
			f = line.split()
			got[f[0]] = [tuple(int(x) for x in item.split(":")) for item in f[1:]]

		def near(c, others):
			return any(abs(c[0] - o[0]) <= d for o in others)

		want = {
			"p1": [c for c in lists["p1"] if near(c, lists["n2"])],
			"n2": [c for c in lists["n2"] if near(c, lists["p1"])],
			"n1": [c for c in lists["n1"] if near(c, lists["p2"])],
			"p2": [c for c in lists["p2"] if near(c, lists["n1"])],
		}
		assert got == want, (case, d, lists)


def _reduce(cpp, d, lists):
	text = "\n".join(tag + "".join(" %d:%d" % c for c in lists.get(tag, [])) for tag in ("p1", "n1", "p2", "n2"))
	out = cpp("b_probes").probe("reduce", d, input=text + "\n")
	got = {}
	for line in out.splitlines():
		f = line.split()
		got[f[0]] = [tuple(int(x) for x in item.split(":")) for item in f[1:]]
	return got


def test_reduce_candidates_keeps_strong_unpaired_candidates(cpp):
	"""Hand-worked, distance 300: an unpaired candidate is still kept when
	its support is at least 6 and at least the support of every paired
	candidate kept before it on its side, it lies on the reference of the
	candidate of the other read it was compared with, and fewer than 5
	unpaired ones were kept on that side."""
	c = lambda pos, n, rid=0: ((rid << 32) | pos, n)
	cases = [
		# read 1 + at 1000 (support 8) and read 2 - at 5000 (3): only the
		# strong one is kept.
		({"p1": [c(1000, 8)], "n2": [c(5000, 3)]}, {"p1": [c(1000, 8)]}),
		# support 6 is enough, 5 is not
		({"p1": [c(1000, 6), c(2000, 5)], "n2": [c(9000, 2)]}, {"p1": [c(1000, 6)]}),
		# the other read's candidate is on another reference: not kept
		({"p1": [c(1000, 9)], "n2": [c(500, 2, 1)]}, {}),
		# at most 5 unpaired candidates per side
		({"p1": [c(1000 + 400 * i, 7) for i in range(7)], "n2": [c(20000, 2)]},
			{"p1": [c(1000 + 400 * i, 7) for i in range(5)]}),
		# read 2's side, and the other direction (read 1 -, read 2 +)
		({"p1": [c(5000, 2)], "n2": [c(1000, 7)]}, {"n2": [c(1000, 7)]}),
		({"n1": [c(1000, 10)], "p2": [c(3000, 1)]}, {"n1": [c(1000, 10)]}),
		# a pair with support 9 raises the bar on read 1's side: 8 is not kept
		({"p1": [c(1000, 9), c(5000, 8)], "n2": [c(1100, 2), c(9000, 2)]},
			{"p1": [c(1000, 9)], "n2": [c(1100, 2)]}),
	]
	for lists, want in cases:
		got = _reduce(cpp, 300, lists)
		assert got == {tag: want.get(tag, []) for tag in ("p1", "n1", "p2", "n2")}, lists
