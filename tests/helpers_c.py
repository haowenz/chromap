"""Reference functions and read-set designs for the area C tests (output and
post-processing).

Expected values come from the output format definitions (README.md, the SAM
and PAF specifications), chromap's documented rules for sorting,
deduplication and Tn5 shifting (tests/README.md, chromap.1) and the
simulation truth. Nothing here is copied from chromap's output.

Records passed to the C++ probes (tests/cpp/c_pipeline.cc) are dicts; the
`*_line` functions format them the way the probe reads them, and the
`format_*` functions give the lines chromap's writers should produce.
"""

import collections
import math
import random

import chromap_helpers as ch

MAX_DUPS = 255
MIN_UNIQUE_MAPQ = 4  # chromap's min_unique_mapping_mapq_ for multi-mapping allocation

# Fields of each record type, in the probe's input order (rid first, dups
# last), see tests/cpp/c_pipeline.cc.
PROBE_FIELDS = {
	"se": ("rid", "id", "start", "length", "mapq", "dir", "uniq", "dups"),
	"seb": ("rid", "id", "barcode", "start", "length", "mapq", "dir", "uniq", "dups"),
	"pe": ("rid", "id", "start", "length", "mapq", "dir", "uniq", "plen", "nlen", "dups"),
	"peb": ("rid", "id", "barcode", "start", "length", "mapq", "dir", "uniq", "plen", "nlen", "dups"),
	"paf": ("rid", "id", "name", "rlen", "start", "length", "mapq", "dir", "uniq", "dups"),
	"ppaf": ("rid", "id", "name1", "name2", "rlen1", "rlen2", "start", "length", "plen", "nlen", "mapq", "mapq1",
		"mapq2", "dir", "uniq", "dups"),
	"sam": ("rid", "id", "name", "barcode", "pos", "mrid", "mpos", "tlen", "flag", "positive", "mapq", "nm", "cigar",
		"md", "seq", "qual", "dups"),
}

PAIRED_TYPES = ("pe", "peb", "ppaf")
BARCODED_TYPES = ("seb", "peb")


def record(kind, **fields):
	"""A record dict of type `kind` with defaults for what is not given."""
	defaults = {"rid": 0, "mapq": 60, "dir": 1, "uniq": 1, "dups": 1, "barcode": 0, "plen": 50, "nlen": 50,
		"rlen": 50, "rlen1": 50, "rlen2": 50, "mapq1": 60, "mapq2": 60, "mrid": -1, "mpos": 0, "tlen": 0,
		"flag": 0, "positive": 1, "nm": 0, "cigar": "50M", "md": "50", "seq": "-", "qual": "-"}
	rec = dict(defaults)
	rec.update(fields)
	rec["kind"] = kind
	missing = [f for f in PROBE_FIELDS[kind] if f not in rec]
	if missing:
		raise ValueError("record %s misses %s" % (kind, missing))
	return rec


def probe_line(rec, file_index=None):
	vals = [str(rec[f]) for f in PROBE_FIELDS[rec["kind"]]]
	if file_index is not None:
		vals.insert(0, str(file_index))
	return " ".join(vals)


def probe_input(records, file_of=None):
	"""stdin text for a probe: one record per line ("FILE_INDEX record" when
	file_of maps a record index to its temp file)."""
	lines = []
	for i, r in enumerate(records):
		lines.append(probe_line(r, None if file_of is None else file_of(i)))
	return "".join(x + "\n" for x in lines)


def parse_probe_record(kind, tokens):
	"""Inverse of probe_line for the probe's printed records ("rid <fields
	without rid>")."""
	names = PROBE_FIELDS[kind]
	if len(tokens) != len(names):
		raise ValueError("%s record needs %d fields, got %r" % (kind, len(names), tokens))
	rec = {"kind": kind}
	for name, tok in zip(names, tokens):
		rec[name] = tok if name in ("name", "name1", "name2", "cigar", "md", "seq", "qual") else int(tok)
	return rec


# ------------------------------------------------------------ records


def end(rec):
	return rec["start"] + rec["length"]


def tn5_shift(rec):
	"""A copy of `rec` with the Tn5 shift applied (chromap.1 --Tn5-shift:
	forward start positions +4, reverse end positions -5): single-end +
	records start 4 later and keep their length, so their end moves by 4 as
	well (chromap.1 names only the start; this part pins the current
	behaviour), - records lose 5 bases at their end; fragments start 4 later
	and end 5 earlier (read 1's + mate loses 4 bases, the - mate 5). SAM
	records are not shifted (documented)."""
	r = dict(rec)
	kind = r["kind"]
	if kind == "sam":
		return r
	if kind in PAIRED_TYPES:
		r["start"] += 4
		r["length"] -= 9
		r["plen"] -= 4
		r["nlen"] -= 5
	elif r["dir"] == 1:
		r["start"] += 4
	else:
		r["length"] -= 5
	return r


def dedup_key(rec):
	"""Records with equal keys on one reference are PCR duplicates
	(tests/README.md): single-end records share a start whatever their strand
	or length; fragments share start and length; barcoded records also share
	the barcode (cell level); SAM records share position, barcode, read-1 bit
	and the mate's reference and position."""
	kind = rec["kind"]
	if kind in ("se", "paf"):
		return (rec["rid"], rec["start"])
	if kind == "seb":
		return (rec["rid"], rec["barcode"], rec["start"])
	if kind in ("pe", "ppaf"):
		return (rec["rid"], rec["start"], rec["length"])
	if kind == "peb":
		return (rec["rid"], rec["barcode"], rec["start"], rec["length"])
	if kind == "sam":
		return (rec["rid"], rec["pos"], rec["barcode"], rec["flag"] & 64, rec["mrid"], rec["mpos"])
	raise ValueError(kind)


def sort_key(rec):
	"""Output order within a reference: each record class's operator<
	(tests/README.md, *Order*)."""
	r = rec
	kind = r["kind"]
	if kind == "se":
		return (r["start"], r["length"], r["mapq"], r["dir"], r["uniq"], r["id"])
	if kind == "seb":
		return (r["start"], r["length"], r["barcode"], r["mapq"], r["dir"], r["uniq"], r["id"])
	if kind == "pe":
		return (r["start"], r["length"], r["mapq"], r["dir"], r["uniq"], r["id"], r["plen"], r["nlen"])
	if kind == "peb":
		return (r["start"], r["length"], r["barcode"], r["mapq"], r["dir"], r["uniq"], r["id"], r["plen"], r["nlen"])
	if kind == "paf":
		return (r["start"], r["length"], r["mapq"], r["dir"], r["uniq"], r["id"], r["rlen"])
	if kind == "ppaf":
		return (r["start"], r["length"], r["mapq1"], r["mapq2"], r["dir"], r["uniq"], r["id"], r["plen"], r["nlen"])
	if kind == "sam":
		return (r["pos"], r["barcode"], r["mrid"], r["mpos"], r["flag"] & 64, r["mapq"], r["id"])
	raise ValueError(kind)


def group_key_order(members):
	"""Where a duplicate group sorts: the smallest sort key of its members.
	Test designs keep groups whose members differ in sort-relevant fields
	apart from other groups, so the choice does not matter."""
	return min(sort_key(m) for m in members)


# ------------------------------------------------------------ formats


def barcode_string(key, length):
	return ch.decode_seed(key, length)


def format_record(rec, fmt, contigs, barcode_length=0, count=None):
	"""Output lines (a list: paired TagAlign and paired PAF give two) for one
	record in `fmt` ("BED", "TagAlign", "PAF", "SAM"), with duplicate count
	`count` (default the record's own). `contigs` is [(name, length)] in
	reference order.

	Formats (README.md, tests/README.md): bulk BED and single-end TagAlign
	have 7 columns (chrom start end N mapq strand count); barcoded BED is the
	fragment file (chrom start end barcode count); barcoded single-end
	TagAlign has 6 columns; paired-end TagAlign writes read 1's line then
	read 2's, the count only on the second line and only for bulk data; PAF
	has 12 columns (query name, length, start 0, end, strand, target name,
	length, start, end, residue matches, alignment block length, MAPQ); SAM
	has the 11 mandatory fields, NM, MD and CB for barcoded data."""
	kind = rec["kind"]
	n = rec["dups"] if count is None else count
	name, clen = contigs[rec["rid"]]
	if kind == "sam":
		unmapped = rec["flag"] & 4
		rnext = "*" if rec["mrid"] < 0 else ("=" if rec["mrid"] == rec["rid"] else contigs[rec["mrid"]][0])
		pnext = 0 if rec["mrid"] < 0 else rec["mpos"] + 1
		seq = "" if rec["seq"] == "-" else rec["seq"]
		qual = "" if rec["qual"] == "-" else rec["qual"]
		qlen = cigar_query_length(rec["cigar"])
		if not rec["positive"]:
			qual = qual[::-1]
		seq, qual = seq[:qlen], qual[:qlen]
		md = "" if rec["md"] == "-" else rec["md"]
		fields = [rec["name"], str(rec["flag"]), "*" if unmapped else name, str(rec["pos"] + 1), str(rec["mapq"]),
			rec["cigar"], rnext, str(pnext), str(rec["tlen"]), seq, qual, "NM:i:%d" % rec["nm"], "MD:Z:" + md]
		if barcode_length > 0:
			fields.append("CB:Z:" + barcode_string(rec["barcode"], barcode_length))
		return ["\t".join(fields)]
	strand = "+" if rec["dir"] == 1 else "-"
	s, e = rec["start"], end(rec)
	if kind == "paf":
		return ["\t".join(str(x) for x in (rec["name"], rec["rlen"], 0, rec["rlen"], strand, name, clen, s, e,
			rec["rlen"], rec["length"], rec["mapq"]))]
	if kind == "ppaf":
		plus = (s, s + rec["plen"], rec["plen"])
		minus = (e - rec["nlen"], e, rec["nlen"])
		r1 = plus if rec["dir"] == 1 else minus
		r2 = minus if rec["dir"] == 1 else plus
		s1, s2 = ("+", "-") if rec["dir"] == 1 else ("-", "+")
		return ["\t".join(str(x) for x in (rec["name1"], rec["rlen1"], 0, rec["rlen1"], s1, name, clen, r1[0], r1[1],
				rec["rlen1"], r1[2], rec["mapq1"])),
			"\t".join(str(x) for x in (rec["name2"], rec["rlen2"], 0, rec["rlen2"], s2, name, clen, r2[0], r2[1],
				rec["rlen2"], r2[2], rec["mapq2"]))]
	barcoded = kind in BARCODED_TYPES
	if fmt == "BED" or (fmt == "TagAlign" and kind in ("se", "seb")):
		if fmt == "BED" and barcoded:
			return ["%s\t%d\t%d\t%s\t%d" % (name, s, e, barcode_string(rec["barcode"], barcode_length), n)]
		line = "%s\t%d\t%d\tN\t%d\t%s" % (name, s, e, rec["mapq"], strand)
		return [line if barcoded else line + "\t%d" % n]
	if fmt == "TagAlign":
		plus = "%s\t%d\t%d\tN\t%d\t+" % (name, s, s + rec["plen"], rec["mapq"])
		minus = "%s\t%d\t%d\tN\t%d\t-" % (name, e - rec["nlen"], e, rec["mapq"])
		first, second = (plus, minus) if rec["dir"] == 1 else (minus, plus)
		return [first, second if barcoded else second + "\t%d" % n]
	raise ValueError((kind, fmt))


_CIGAR_QUERY = set("MIS=X")
_CIGAR_REF = set("MDN=X")


def parse_cigar(cigar):
	if cigar == "*":
		return []
	out, num = [], ""
	for c in cigar:
		if c.isdigit():
			num += c
		else:
			out.append((int(num), c))
			num = ""
	return out


def cigar_query_length(cigar):
	"""Query bases a CIGAR consumes (SAM specification: M, I, S, =, X)."""
	return sum(n for n, op in parse_cigar(cigar) if op in _CIGAR_QUERY)


def cigar_reference_length(cigar):
	"""Reference bases a CIGAR consumes (M, D, N, =, X)."""
	return sum(n for n, op in parse_cigar(cigar) if op in _CIGAR_REF)


# ------------------------------------------------------------ pipelines


def duplicate_groups(records, apply_tn5_first):
	"""OrderedDict dedup key -> member records (Tn5-shifted when
	`apply_tn5_first`), groups in output order."""
	recs = [tn5_shift(r) if apply_tn5_first else dict(r) for r in records]
	groups = collections.OrderedDict()
	for r in recs:
		groups.setdefault(dedup_key(r), []).append(r)
	order = sorted(groups, key=lambda k: (k[0], group_key_order(groups[k])))
	return collections.OrderedDict((k, groups[k]) for k in order)


def expected_output(records, fmt, contigs, tn5=False, dedup=False, mapq_threshold=0, barcode_length=0,
		tn5_before_dedup=True):
	"""What chromap should write for `records`: a list of (candidate line
	tuples, group members) in output order. Each entry stands for one output
	record; when a duplicate group's members would print differently, any one
	of them may be written (neither mode documents which, tests/README.md), so
	every member's lines are candidates. Groups are dropped when the best
	MAPQ in the group is below the threshold (test designs give a group's
	members equal lengths when they differ in MAPQ, so both modes keep a
	member with that MAPQ). Without dedup every record is its own group with
	count 1."""
	if dedup:
		groups = duplicate_groups(records, tn5 and tn5_before_dedup)
		if tn5 and not tn5_before_dedup:
			groups = collections.OrderedDict((k, [tn5_shift(m) for m in v]) for k, v in groups.items())
	else:
		recs = [tn5_shift(r) if tn5 else dict(r) for r in records]
		recs.sort(key=lambda r: (r["rid"], sort_key(r)))
		groups = collections.OrderedDict((i, [r]) for i, r in enumerate(recs))
	out = []
	for members in groups.values():
		if max(m["mapq"] for m in members) < mapq_threshold:
			continue
		count = min(MAX_DUPS, len(members)) if dedup else members[0]["dups"]
		best = max(m["mapq"] for m in members)
		cands = set(tuple(format_record(m, fmt, contigs, barcode_length, count)) for m in members if m["mapq"] == best)
		out.append((cands, members))
	return out


def check_output(lines, expected, header=()):
	"""Assert `lines` (output without header lines) match `expected` from
	expected_output: one candidate per entry, in order."""
	lines = [l for l in lines if l not in header]
	width = [len(next(iter(c))) for c, _ in expected]
	assert len(lines) == sum(width), "%d output lines, expected %d" % (len(lines), sum(width))
	i = 0
	for (cands, members), w in zip(expected, width):
		got = tuple(lines[i:i + w])
		assert got in cands, "output %r is none of %r (members %r)" % (got, sorted(cands), members)
		i += w


# ------------------------------------------------------------ allocation


def overlap_weight(rec, uniques, distance):
	"""Number of unique mappings on the record's reference that overlap
	[start - distance, end + distance) (start clipped at 0), half-open
	intervals: the weight chromap gives one location of a multi-read
	(--multi-mapping-allocation-distance)."""
	lo = max(0, rec["start"] - distance)
	hi = end(rec) + distance
	return sum(1 for u in uniques.get(rec["rid"], ()) if lo < end(u) and u["start"] < hi)


def allocation_classes(records, distance, min_unique_mapq=MIN_UNIQUE_MAPQ):
	"""Split single-end records into unique ones (MAPQ >= min_unique_mapq,
	kept as they are) and multi-reads (grouped by read id) with each
	location's weight. Returns (uniques, {read id: [(location, weight)]})."""
	uniques = collections.defaultdict(list)
	multi = collections.OrderedDict()
	for r in records:
		if r["mapq"] >= min_unique_mapq:
			uniques[r["rid"]].append(r)
	for r in records:
		if r["mapq"] < min_unique_mapq:
			multi.setdefault(r["id"], []).append(r)
	weights = collections.OrderedDict(
		(rid, [(loc, overlap_weight(loc, uniques, distance)) for loc in locs]) for rid, locs in multi.items())
	return uniques, weights


def classify(weighted):
	"""'dropped' (every weight 0), 'fixed' (exactly one location with a
	positive weight) or 'random'."""
	positive = [loc for loc, w in weighted if w > 0]
	if not positive:
		return "dropped"
	return "fixed" if len(positive) == 1 else "random"


# ------------------------------------------------------------ summary


def inverse_logit(x):
	return 1.0 / (1.0 + math.exp(-x))


def estimated_frip(coeffs, fric, dup, unmapped, lowmapq):
	"""README.md, *Estimating FRiP*: inverse logit of a linear model in fric,
	duplicate, unmapped and lowmapq; 0 when fric is 0 (no cache hits)."""
	if fric == 0:
		return 0.0
	return inverse_logit(coeffs[0] + coeffs[1] * fric + coeffs[2] * dup + coeffs[3] * unmapped + coeffs[4] * lowmapq)


SUMMARY_HEADER = ["barcode", "total", "duplicate", "unmapped", "lowmapq", "cachehit", "fric", "estfrip"]


def summary_row_set(summary):
	"""Rows of a parsed summary as a set of tuples (row order is not
	specified)."""
	return set((bc,) + tuple(row[k] for k in summary.header[1:]) for bc, row in summary.rows.items())


# ------------------------------------------------------------ read designs


def interval_key(read):
	return (read.contig, read.start)


class Group(object):
	"""Reads designed to form one duplicate group: `reads`, a label `kind`
	and `key`, what every member shares and the output line shows: (contig,
	start) for single-end reads, (contig, start, end) for fragments; None
	where the location is not known in advance (multimappers)."""

	def __init__(self, kind, reads, key=None):
		self.kind = kind
		self.reads = list(reads)
		if key is None and kind != "multimapper":
			r = self.reads[0]
			key = (r.contig, r.start, r.end) if r.paired else (r.contig, r.start)
		self.key = key

	@property
	def size(self):
		return len(self.reads)

	@property
	def names(self):
		return set(r.name for r in self.reads)

	def __repr__(self):
		return "Group(%s, %d reads, %r)" % (self.kind, self.size, self.reads[0])


def single_end_dedup_design(genome, seed, identical_sizes=(2, 3, 4, 5, 6, 7, 9, 12), singletons=40):
	"""Single-end reads for the deduplication tests.

	Groups (each at its own unique locus, 50 bp unless stated):
	- `singletons` unique reads (groups of one), both strands;
	- identical duplicates: for each size in identical_sizes, that many copies
	  of one read (strand alternating between groups);
	- mixed strand: four groups of one + read and two - reads sharing their
	  leftmost base (duplicates whatever the strand);
	- mixed length: three groups of + reads of 50, 60 and 70 bp sharing a
	  start;
	- multimappers: three groups of 2, 3 and 4 identical reads from the
	  eight-copy repeat rep8 (MAPQ 0, so dropped at the default -q 30;
	  identical reads land on the same copy, tests/README.md);
	- unmapped: five random reads, and three genomic reads of 25 bp, shorter
	  than --min-read-length 30.

	Returns (ReadSet, groups, unmapped reads)."""
	sim = ch.ReadSimulator(genome, seed=seed)
	groups = []
	for i in range(singletons):
		contig, start = sim.unique_locus(50)
		groups.append(Group("singleton", [sim.se_read(contig, start, 50, "+-"[i % 2])]))
	for i, n in enumerate(identical_sizes):
		contig, start = sim.unique_locus(50)
		r = sim.se_read(contig, start, 50, "+-"[i % 2])
		groups.append(Group("identical", [r] + [sim.se_read(contig, start, 50, r.strand) for _ in range(n - 1)]))
	for _ in range(4):
		contig, start = sim.unique_locus(50)
		groups.append(Group("mixed_strand", [sim.se_read(contig, start, 50, s) for s in "+--"]))
	for _ in range(3):
		contig, start = sim.unique_locus(70)
		groups.append(Group("mixed_length", [sim.se_read(contig, start, n, "+") for n in (50, 60, 70)]))
	for n, offset in zip((2, 3, 4), (60, 160, 260)):
		r = sim.segment_read("rep8", offset, 50, "+")
		groups.append(Group("multimapper", [r] + [sim.segment_read("rep8", offset, 50, "+") for _ in range(n - 1)]))
	unmapped = [sim.random_read(50) for _ in range(5)]
	for _ in range(3):
		contig, start = sim.unique_locus(25)
		unmapped.append(sim.se_read(contig, start, 25, "+"))
	reads = [r for g in groups for r in g.reads] + unmapped
	return sim.shuffled(ch.ReadSet(reads)), groups, unmapped


def single_end_group_lines(group, count, shift=False):
	"""Candidate BED lines (7 columns) for a single-end duplicate group with
	MAPQ 60: any member's interval and strand with the group's count,
	Tn5-shifted when `shift`."""
	out = set()
	for r in group.reads:
		s, e = r.start, r.end
		if shift:
			s, e = (s + 4, e + 4) if r.strand == "+" else (s, e - 5)
		out.add("%s\t%d\t%d\tN\t60\t%s\t%d" % (r.contig, s, e, r.strand, count))
	return out


def match_groups(records, groups, record_key, candidates):
	"""Assert every output record (a line, or a tuple of lines) belongs to
	exactly one designed group, found by record_key(record) == group.key,
	and is one of candidates(group); returns {group index: record}."""
	by_key = {}
	for i, g in enumerate(groups):
		assert g.key not in by_key, "two designed groups share %r" % (g.key,)
		by_key[g.key] = i
	seen = {}
	for rec in records:
		k = record_key(rec)
		assert k in by_key, "output %r belongs to no designed group" % (rec,)
		i = by_key[k]
		assert i not in seen, "group %r written twice: %r and %r" % (groups[i], seen[i], rec)
		cands = candidates(groups[i])
		assert rec in cands, "output %r is not one of %r" % (rec, sorted(cands))
		seen[i] = rec
	return seen


def pair_tagalign_lines(pair, count, barcoded=False):
	"""The two TagAlign lines of a pair with MAPQ 60: read 1's interval and
	strand first, then read 2's with the count (bulk data only)."""
	first = "%s\t%d\t%d\tN\t60\t%s" % (pair.contig, pair.r1.start, pair.r1.end, pair.r1.strand)
	second = "%s\t%d\t%d\tN\t60\t%s" % (pair.contig, pair.r2.start, pair.r2.end, pair.r2.strand)
	return (first, second if barcoded else second + "\t%d" % count)


def contig_rank(genome, contig):
	return genome.names.index(contig)


def random_pair(sim, name=None, barcode=None):
	"""A read pair of random sequence (should not map)."""
	name = name or sim.next_name("u")
	return ch.SimPair(name, sim.random_read(50, name=name), sim.random_read(50, name=name), None, None, None, None,
		kind="random", barcode=barcode)


def repeat_pair(sim, offset, insert=200, strand="+", barcode=None):
	"""A read pair inside one copy of the eight-copy repeat rep8 (both mates
	at eight equally good places: MAPQ below 30)."""
	seg = sim.g.segments["rep8"]
	c = seg.copies[0]
	return sim.pair(c.contig, c.start + offset, insert, strand=strand, barcode=barcode, kind="repeat")


def paired_dedup_design(genome, seed, identical_sizes=(2, 3, 4, 5), mixed_groups=3, uniques=50, repeats=2,
		unmapped=3):
	"""Read pairs for the paired-end deduplication tests: `uniques` pairs at
	unique loci (groups of one); identical duplicate groups of the given
	sizes (read 1's strand alternating between groups); `mixed_groups`
	groups of two pairs from one fragment with read 1 on + in one and on - in
	the other; `repeats` pairs inside rep8 (MAPQ below 30) and `unmapped`
	random pairs. Returns (ReadSet, groups, repeat pairs, unmapped pairs)."""
	sim = ch.ReadSimulator(genome, seed=seed)
	groups = [Group("singleton", [p]) for p in sim.pairs(uniques)]
	for i, n in enumerate(identical_sizes):
		ins = 180 + 20 * i
		contig, start = sim.unique_locus(ins)
		groups.append(Group("identical", [sim.pair(contig, start, ins, strand="+-"[i % 2]) for _ in range(n)]))
	for i in range(mixed_groups):
		ins = 230 + 10 * i
		contig, start = sim.unique_locus(ins)
		groups.append(Group("mixed_strand", [sim.pair(contig, start, ins, strand=s) for s in "+-"]))
	reps = [repeat_pair(sim, 20 + 130 * i, strand="+-"[i % 2]) for i in range(repeats)]
	unm = [random_pair(sim) for _ in range(unmapped)]
	reads = [p for g in groups for p in g.reads] + reps + unm
	return sim.shuffled(ch.ReadSet(reads)), groups, reps, unm


def single_end_cap_design(genome, seed, uniques=60):
	"""Single-end reads for the count-cap test: `uniques` reads at unique
	loci on both strands, identical groups of 2, 3 and 4 reads, a group of
	300 identical + reads and one of 256 identical - reads (counts capped at
	255). Returns (ReadSet, groups)."""
	sim = ch.ReadSimulator(genome, seed=seed)
	groups = [Group("singleton", [r]) for r in sim.se_reads(uniques)]
	for n, strand in ((2, "+"), (3, "-"), (4, "+"), (300, "+"), (256, "-")):
		contig, start = sim.unique_locus(50)
		groups.append(Group("identical", [sim.se_read(contig, start, 50, strand) for _ in range(n)]))
	reads = [r for g in groups for r in g.reads]
	return sim.shuffled(ch.ReadSet(reads)), groups


def non_whitelisted_barcode(whitelist, seed, min_dist=3):
	"""A barcode at least `min_dist` mismatches from every whitelist entry,
	so a one-mismatch correction cannot move it into the whitelist."""
	rng = random.Random(seed)
	while True:
		b = "".join(rng.choice("ACGT") for _ in range(len(whitelist[0])))
		if all(ch.hamming(b, w) >= min_dist for w in whitelist):
			return b


def barcoded_pair_design(genome, seed, whitelist, tie_breaks=False):
	"""Barcoded read pairs (whitelisted barcodes unless stated):
	- four fragments, each read three times with whitelist[0], once with
	  whitelist[1] and twice with whitelist[2] (read 1's strand alternating
	  between fragments): six pairs per fragment, three cell-level groups;
	- ten pairs at unique loci with barcodes from whitelist[3:8];
	- one pair whose barcode is three or more mismatches from every
	  whitelisted barcode (dropped, counted as non-whitelist);
	- one pair inside rep8 with whitelist[8] (MAPQ below 30) and one random
	  pair with whitelist[9] (unmapped);
	- with tie_breaks, two more fragment groups for bulk-level duplicate
	  removal (appended to the fragment groups, kinds "most_pairs" and
	  "abundance_tie"): a fragment read twice with whitelist[10] and once
	  with the more abundant whitelist[0], and a fragment read once with a
	  barcode used nowhere else and once with whichever of whitelist[0] and
	  whitelist[2] has the larger 2-bit key (more abundant, and later in
	  chromap's barcode order).
	Returns (ReadSet, fragment groups (all pairs of a fragment), singleton
	groups, the non-whitelisted pair, the repeat pair, the unmapped pair)."""
	sim = ch.ReadSimulator(genome, seed=seed)
	fragments = []
	for i in range(4):
		ins = 200 + 15 * i
		contig, start = sim.unique_locus(ins)
		members = []
		for bc, n in ((whitelist[0], 3), (whitelist[1], 1), (whitelist[2], 2)):
			members += [sim.pair(contig, start, ins, strand="+-"[i % 2], barcode=bc) for _ in range(n)]
		fragments.append(Group("fragment", members))
	singles = [Group("singleton", [p]) for p in sim.pairs(10, barcodes=whitelist[3:8])]
	contig, start = sim.unique_locus(250)
	bad = sim.pair(contig, start, 250, barcode=non_whitelisted_barcode(whitelist, seed))
	rep = repeat_pair(sim, 90, barcode=whitelist[8])
	unm = random_pair(sim, barcode=whitelist[9])
	if tie_breaks:
		contig, start = sim.unique_locus(260)
		members = [sim.pair(contig, start, 260, barcode=whitelist[10]) for _ in range(2)]
		members.append(sim.pair(contig, start, 260, barcode=whitelist[0]))
		fragments.append(Group("most_pairs", members))
		heavy = max((whitelist[0], whitelist[2]), key=ch.encode_seed)
		light = min(whitelist[11:], key=ch.encode_seed)
		if ch.encode_seed(light) > ch.encode_seed(heavy):
			raise ValueError("no unused whitelist barcode sorts before %s" % heavy)
		contig, start = sim.unique_locus(270)
		fragments.append(Group("abundance_tie", [sim.pair(contig, start, 270, barcode=b) for b in (light, heavy)]))
	reads = [p for g in fragments + singles for p in g.reads] + [bad, rep, unm]
	return sim.shuffled(ch.ReadSet(reads)), fragments, singles, bad, rep, unm
