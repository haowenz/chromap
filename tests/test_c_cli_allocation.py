"""Area C, CLI layer: --allocate-multi-mappings,
--multi-mapping-allocation-distance and --multi-mapping-allocation-seed.

Multi-reads come from the eight-copy repeat rep8 (400 bp, copies on both
strands, tests/README.md) and are mapped with -n 8 -q 0, so all eight
locations are reported with MAPQ 0 (below 4: multi-mappings) and pass the
MAPQ filter. Unique reads (MAPQ 60) straddle the ends of some copies or lie
just outside one. chromap gives each multi-read one of its locations, drawn
with weights equal to the number of unique mappings overlapping the
location widened by the allocation distance, and drops reads whose weights
are all zero (helpers_c.overlap_weight). Each multi-read starts at its own
offset in the repeat, so an output line's reference and start name the read
and the copy."""

import collections

import pytest

import chromap_helpers as ch
import helpers_c as hc

ARGS = "-n 8 -q 0 --allocate-multi-mappings"


@pytest.fixture(scope="module")
def design(genome):
	"""Multi-reads, one per offset into rep8 (labels E: offsets 0-9, M:
	160-189, S: 300-309, E2: 341-350). Unique reads across the left end of
	copies 3 and 5 ([start - 30, start + 20)), across the right end of copy
	6 (+) and copy 2 (-) ([end - 20, end + 30)), and 40 bases past the end
	of copy 7. By the weights, at distance 0 E is drawn among copies 2, 3
	and 5, E2 goes to copy 6, S and M are dropped; at distance 100 E2 and S
	are drawn between copies 6 and 7 and M is still dropped.

	Returns (ReadSet, unique reads, {name: (label, read)}, records) where
	records are the unique mappings and every multi-read location as
	helpers_c records."""
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_c_cli_allocation", "design"))
	copies = genome.segments["rep8"].copies
	uniques = [sim.se_read(copies[i].contig, copies[i].start - 30, 50, "+") for i in (3, 5)]
	uniques += [sim.se_read(copies[i].contig, copies[i].end - 20, 50, "-") for i in (6, 2)]
	uniques.append(sim.se_read(copies[7].contig, copies[7].end + 40, 50, "+"))
	multi = collections.OrderedDict()
	for label, offsets in (("E", range(0, 10)), ("M", range(160, 190)), ("S", range(300, 310)),
			("E2", range(341, 351))):
		for o in offsets:
			r = sim.segment_read("rep8", o, 50, "+")
			multi[r.name] = (label, r)
	records = [hc.record("se", rid=genome.rid(u.contig), id=i, start=u.start, length=50, mapq=60)
		for i, u in enumerate(uniques)]
	for j, (label, r) in enumerate(multi.values()):
		for c, s, e in r.extra["loci"]:
			records.append(hc.record("se", rid=genome.rid(c), id=1000 + j, start=s, length=e - s, mapq=0, uniq=0,
				name=r.name))
	reads = sim.shuffled(ch.ReadSet(uniques + [r for _, r in multi.values()]))
	return reads, uniques, multi, records


def _run(index, design, extra):
	return ch.run_chromap(index, design[0], ARGS + extra, shared=True)


def _allocation(run, design, distance):
	"""Check the output against the weights at `distance` and return
	({read name: (contig, start)}, {name: class})."""
	reads, uniques, multi, records = design
	genome = run.index.genome
	lines = [l.split("\t") for l in run.lines()]
	order = [(genome.rid(f[0]), int(f[1])) for f in lines]
	assert order == sorted(order)
	got_uniques = sorted((f[0], int(f[1]), int(f[2])) for f in lines if f[4] == "60")
	assert got_uniques == sorted((u.contig, u.start, u.end) for u in uniques)
	where = {}
	for label, r in multi.values():
		for c, s, e in r.extra["loci"]:
			assert (c, s) not in where
			where[(c, s)] = r.name
	placed = {}
	for f in lines:
		if f[4] == "60":
			continue
		assert f[4] == "0", f
		name = where[(f[0], int(f[1]))]
		assert name not in placed, "%s placed twice" % name
		placed[name] = (f[0], int(f[1]))
	_, weights = hc.allocation_classes(records, distance)
	names = dict((1000 + j, name) for j, name in enumerate(multi))
	classes = {}
	for read_id, weighted in weights.items():
		name = names[read_id]
		cls = hc.classify(weighted)
		classes[name] = cls
		positive = set((genome.names[loc["rid"]], loc["start"]) for loc, w in weighted if w > 0)
		if cls == "dropped":
			assert name not in placed, name
		else:
			assert placed.get(name) in positive, (name, placed.get(name), positive)
	assert len(placed) == sum(1 for c in classes.values() if c != "dropped")
	return placed, classes


def _by_label(design, classes):
	out = collections.defaultdict(collections.Counter)
	for name, (label, _) in design[2].items():
		out[label][classes[name]] += 1
	return dict((k, dict(v)) for k, v in out.items())


def test_allocation_by_overlap(index, design, genome):
	"""--allocate-multi-mappings (distance 0): each multi-read is written
	once at a location overlapping a unique read, or not at all when none of
	its eight locations does; unique reads are written as they are; lines
	stay sorted. E2 reads all go to copy 6, the only copy where they overlap
	a unique read."""
	run = _run(index, design, "")
	placed, classes = _allocation(run, design, 0)
	assert _by_label(design, classes) == {"E": {"random": 10}, "M": {"dropped": 30}, "S": {"dropped": 10},
		"E2": {"fixed": 10}}
	c6 = genome.segments["rep8"].copies[6]
	for name, (label, r) in design[2].items():
		if label == "E2":
			assert placed[name][0] == c6.contig and c6.start <= placed[name][1] < c6.end


def test_allocation_distance(index, design):
	"""--multi-mapping-allocation-distance 100 widens each location by 100
	bases on both sides before counting overlapping unique reads: S reads,
	50 to 59 bases short of the unique read across copy 6's end and 81 to 90
	short of the one after copy 7, are now placed; E2 reads now reach both;
	M reads, in the middle of the repeat, stay dropped."""
	run = _run(index, design, " --multi-mapping-allocation-distance 100")
	placed, classes = _allocation(run, design, 100)
	assert _by_label(design, classes) == {"E": {"random": 10}, "M": {"dropped": 30}, "S": {"random": 10},
		"E2": {"random": 10}}


def test_allocation_seed(index, design):
	"""--multi-mapping-allocation-seed: a different seed draws differently
	for at least one of the 30 reads with several positive-weight
	locations (the chance that independent draws agree on all of them, the
	product of the reads' sums of squared location probabilities, is checked
	to be below 1e-6); the default seed is 11, so giving 11 reproduces the
	default run."""
	default = _run(index, design, " --multi-mapping-allocation-distance 100")
	eleven = _run(index, design, " --multi-mapping-allocation-distance 100 --multi-mapping-allocation-seed 11")
	other = _run(index, design, " --multi-mapping-allocation-distance 100 --multi-mapping-allocation-seed 12345")
	assert eleven.lines() == default.lines()
	a, classes = _allocation(default, design, 100)
	b, _ = _allocation(other, design, 100)
	_, weights = hc.allocation_classes(design[3], 100)
	agree = 1.0
	for weighted in weights.values():
		if hc.classify(weighted) == "random":
			total = float(sum(w for _, w in weighted))
			agree *= sum((w / total) ** 2 for _, w in weighted)
	assert agree < 1e-6
	assert set(a) == set(b)
	assert [a[n] for n in sorted(a)] != [b[n] for n in sorted(b)]
