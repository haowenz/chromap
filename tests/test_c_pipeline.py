"""Area C, C++ layer: records, sorting, deduplication, Tn5 shift, writers,
temp files, the low-memory merge, multi-mapping allocation and the summary,
driven through tests/cpp/c_records.cc, c_pipeline.cc and c_minhash.cc.

The CT_TESTs check hand-worked values. The probes feed hand-made records to
the classes chromap calls (MappingProcessor, MappingWriter,
TempMappingFileHandle) and the expected output is computed here with the
reference functions in helpers_c.py (format definitions, chromap's
documented sort and duplicate rules)."""

import collections
import random

import pytest

import chromap_helpers as ch
import helpers_c as hc

pytestmark = pytest.mark.cpp

RECORDS_CPP_TESTS = [
	"single_end_record_order_and_duplicates",
	"single_end_barcoded_record",
	"paired_end_record_order_and_duplicates",
	"paired_end_barcoded_record",
	"paf_records",
	"mapping_interface_dispatch",
	"sam_cigar_strings_and_lengths",
	"sam_constructor_orients_and_truncates",
	"sam_record_order_and_duplicates",
	"sam_flag_bits_and_cigar_table",
	"paf_write_load_round_trip",
	"paired_paf_write_load_round_trip",
	"sam_write_load_round_trip",
]
PIPELINE_CPP_TESTS = ["move_buffers_to_container", "mapping_statistics_counts"]
MINHASH_CPP_TESTS = ["minhash_keeps_k_smallest_distinct", "minhash_below_k_is_zero", "minhash_estimate_arithmetic"]
PIPELINE_PROBES = ["inmemory", "lowmem", "temp", "allocate", "summary"]

CONTIGS = [("chrA", 1000), ("chrB", 2000), ("chrC", 500), ("chrD", 3000)]
BC_LEN = 16


@pytest.fixture(scope="module")
def fasta(tmp_path_factory):
	"""A four-contig FASTA (chrA 1000, chrB 2000, chrC 500, chrD 3000) for
	the writers' reference names and lengths."""
	rng = random.Random(ch.stable_seed("test_c_pipeline", "fasta"))
	path = tmp_path_factory.mktemp("c_pipeline") / "ref.fa"
	lines = []
	for name, n in CONTIGS:
		lines += [">%s test" % name, "".join(rng.choice("ACGT") for _ in range(n))]
	return ch.write_lines(path, lines)


def test_drivers_list_every_test(cpp):
	"""Every CT_TEST of the three area C drivers is run below, and the probes
	used here exist."""
	assert cpp("c_records").tests() == RECORDS_CPP_TESTS
	assert cpp("c_pipeline").tests() == PIPELINE_CPP_TESTS
	assert cpp("c_pipeline").probes() == PIPELINE_PROBES
	assert cpp("c_minhash").tests() == MINHASH_CPP_TESTS


@pytest.mark.parametrize("name", RECORDS_CPP_TESTS)
def test_records_cpp(cpp, name):
	cpp("c_records").check(name)


@pytest.mark.parametrize("name", PIPELINE_CPP_TESTS)
def test_pipeline_cpp(cpp, name):
	cpp("c_pipeline").check(name)


@pytest.mark.parametrize("name", MINHASH_CPP_TESTS)
def test_minhash_cpp(cpp, name):
	cpp("c_minhash").check(name)


# ----------------------------------------------------------- record designs


def _bc(seq):
	return ch.encode_seed(seq)


BARCODES = [_bc(s) for s in ("ACGTACGTACGTACGT", "AAAAAAAAAAAAAAAC", "TTTTGGGGCCCCAAAA", "CATCATCATCATCATC")]


def _seq(rng, n):
	return "".join(rng.choice("ACGT") for _ in range(n))


def _qual(rng, n):
	return "".join(chr(rng.randint(35, 74)) for _ in range(n))


def format_records(kind, rng):
	"""Unique records (no two share a duplicate key) on chrA, chrB and chrD
	(chrC has none): ties in start broken by length, MAPQ, strand and
	barcode, so the sort order is visible in the output; two records with
	MAPQ below 30. Except for SAM (whose order ignores the strand), a + record
	(read id 100) and a - record (read id 101) on chrB tie on every field
	before the strand, so the - record must come first against the read ids;
	both have MAPQ exactly 30, which passes -q 30."""
	recs = []
	rid_starts = [(0, [10, 10, 10, 10, 200, 950]), (1, [5, 5, 700, 1990]), (3, [0, 2999])]
	n = 0
	for rid, starts in rid_starts:
		for j, s in enumerate(starts):
			n += 1
			length = 50 + (j % 3) * 5
			mapq = (60, 31, 10, 60)[j % 4]
			d = (j // 2) % 2
			bc = BARCODES[j % 4]
			if kind == "se":
				recs.append(hc.record(kind, rid=rid, id=n, start=s, length=length, mapq=mapq, dir=d))
			elif kind == "seb":
				recs.append(hc.record(kind, rid=rid, id=n, barcode=bc, start=s, length=length, mapq=mapq, dir=d))
			elif kind in ("pe", "peb"):
				frag = 200 + 7 * j
				recs.append(hc.record(kind, rid=rid, id=n, barcode=bc, start=s, length=frag, mapq=mapq, dir=d,
					plen=50 - j, nlen=45 + j))
			elif kind == "paf":
				recs.append(hc.record(kind, rid=rid, id=n, name="q%d" % n, rlen=length + 1, start=s, length=length,
					mapq=mapq, dir=d))
			elif kind == "ppaf":
				frag = 250 + j
				recs.append(hc.record(kind, rid=rid, id=n, name1="p%d/1" % n, name2="p%d/2" % n, rlen1=50, rlen2=76,
					start=s, length=frag, plen=48, nlen=76, mapq=mapq, mapq1=mapq, mapq2=max(0, mapq - 1), dir=d))
			elif kind == "sam":
				cigar = ("50M", "20M1I29M", "30M2D20M", "5S45M")[j % 4]
				qlen = 50
				seq, qual = _seq(rng, qlen), _qual(rng, qlen)
				mrid = (-1, rid, (rid + 1) % 4)[j % 3]
				flag = 0 if mrid < 0 else (99, 147, 163, 83)[j % 4]
				recs.append(hc.record(kind, rid=rid, id=n, name="s%d" % n, barcode=bc, pos=s, mrid=mrid,
					mpos=0 if mrid < 0 else s + 150 + j, tlen=0 if mrid < 0 else (1 - 2 * (j % 2)) * (200 + j),
					flag=flag, positive=1 - d, mapq=mapq, nm=j % 3, cigar=cigar, md="%d" % (40 + j), seq=seq,
					qual=qual))
	if kind != "sam":
		for d, read_id in ((1, 100), (0, 101)):
			fields = dict(rid=1, id=read_id, start=1500, length=60, mapq=30, dir=d, barcode=BARCODES[1])
			if kind in ("pe", "peb"):
				fields.update(length=230)
			elif kind == "paf":
				fields.update(name="q%d" % read_id, rlen=60)
			elif kind == "ppaf":
				fields.update(length=230, name1="p%d/1" % read_id, name2="p%d/2" % read_id, mapq1=30, mapq2=30)
			recs.append(hc.record(kind, **fields))
	return recs


FORMAT_CASES = [("se", "BED"), ("se", "TagAlign"), ("seb", "BED"), ("seb", "TagAlign"), ("pe", "BED"),
	("pe", "TagAlign"), ("peb", "BED"), ("peb", "TagAlign"), ("paf", "PAF"), ("ppaf", "PAF"), ("sam", "SAM")]


def _run_inmemory(cpp, kind, fmt, fasta, out, records, ops, barcode_length=0):
	args = [kind, fmt, fasta, out] + list(ops)
	if barcode_length:
		args.append("bc:%d" % barcode_length)
	cpp("c_pipeline").probe("inmemory", *args, input=hc.probe_input(records))
	return out.read_text().splitlines()


@pytest.mark.parametrize("kind,fmt", FORMAT_CASES)
def test_writer_formats_and_sort_order(cpp, fasta, tmp_path, kind, fmt):
	"""Without deduplication every record is written once, sorted per
	reference (FASTA order) by its class's order, in the line format of its
	output type (helpers_c.format_record: 7-column BED and single-end
	TagAlign, the 5-column fragment file, 6-column barcoded TagAlign, two
	TagAlign lines per fragment with the count on the second for bulk data,
	12-column PAF, SAM with RNEXT '*', '=' or the mate's reference name and
	CB for barcoded data); records with MAPQ below -q 30 are not written.
	The SAM header lists every reference with its length."""
	rng = random.Random(ch.stable_seed("formats", kind, fmt))
	records = format_records(kind, rng)
	rng.shuffle(records)
	barcoded = kind in ("seb", "peb", "sam")
	# ParallelSortOutputMappings on 3 threads, with 0 threads (chromap's
	# low-memory flush: the caller's threads) and SortOutputMappings.
	ops = [{"se": "psort:3", "pe": "psort:3", "sam": "psort:3", "seb": "psort:0", "peb": "psort:0",
		"paf": "psort:0"}.get(kind, "sort"), "mapq:30"]
	lines = _run_inmemory(cpp, kind, fmt, fasta, tmp_path / "out.txt", records, ops,
		BC_LEN if barcoded else 0)
	header = []
	if fmt == "SAM":
		header = ["@SQ\tSN:%s\tLN:%d" % c for c in CONTIGS]
		assert lines[:4] == header
		lines = lines[4:]
	want = hc.expected_output(records, fmt, CONTIGS, mapq_threshold=30, barcode_length=BC_LEN if barcoded else 0)
	assert len(want) == len(records) - 2
	hc.check_output(lines, want)


def test_sam_writer_unmapped_and_no_barcode(cpp, fasta, tmp_path):
	"""SAM records of bulk data carry no CB tag; a record with the unmapped
	flag (4) has RNAME '*' (SAM specification); a record without a CIGAR
	writes '*'."""
	recs = [
		hc.record("sam", rid=1, id=1, name="u1", pos=10, flag=4, cigar="*", md="-", seq="-", qual="-", mapq=0),
		hc.record("sam", rid=1, id=2, name="m2", pos=20, flag=16, positive=0, seq="ACGTA" * 10,
			qual="".join(chr(40 + i) for i in range(50))),
	]
	lines = _run_inmemory(cpp, "sam", "SAM", fasta, tmp_path / "out.sam", recs, ["sort", "mapq:0"])[4:]
	assert lines == [
		"u1\t4\t*\t11\t0\t*\t*\t0\t0\t\t\tNM:i:0\tMD:Z:",
		"m2\t16\tchrB\t21\t60\t50M\t*\t0\t0\t%s\t%s\tNM:i:0\tMD:Z:50" % ("ACGTA" * 10,
			"".join(chr(40 + i) for i in range(50))[::-1]),
	]


# ----------------------------------------------------------- deduplication


def dedup_records(kind, rng):
	"""Duplicate groups for `kind` on chrA, chrB and chrD:
	- identical groups of 1 to 6 records;
	- a group of 300 (count capped at 255);
	- groups whose members differ in strand (and, for single-end records, in
	  length) with one MAPQ;
	- a group mixing MAPQ 10 and 60 with equal lengths (kept: its best MAPQ
	  passes -q 30, the count covers every member);
	- a group whose MAPQs are all below 30 (dropped);
	- barcoded types: groups at one position with different barcodes (cell
	  level: separate groups).
	Every group sits at its own start (and length for fragments), so the
	group order does not depend on which member is written."""
	recs = []
	n = [0]

	def add(rid, start, size, length=60, mapqs=(60,), dirs=(1,), lengths=None, barcode=None):
		for i in range(size):
			n[0] += 1
			L = lengths[i % len(lengths)] if lengths else length
			mapq = mapqs[i % len(mapqs)]
			d = dirs[i % len(dirs)]
			bc = BARCODES[0] if barcode is None else barcode
			if kind in ("se", "seb"):
				recs.append(hc.record(kind, rid=rid, id=n[0], barcode=bc, start=start, length=L, mapq=mapq, dir=d))
			elif kind in ("pe", "peb"):
				recs.append(hc.record(kind, rid=rid, id=n[0], barcode=bc, start=start, length=L + 150, mapq=mapq,
					dir=d, plen=50, nlen=50))
			elif kind == "paf":
				recs.append(hc.record(kind, rid=rid, id=n[0], name="q%d" % n[0], rlen=L, start=start, length=L,
					mapq=mapq, dir=d))
			elif kind == "ppaf":
				recs.append(hc.record(kind, rid=rid, id=n[0], name1="p%d/1" % n[0], name2="p%d/2" % n[0],
					start=start, length=L + 150, mapq=mapq, mapq1=mapq, mapq2=mapq, dir=d))
			elif kind == "sam":
				recs.append(hc.record(kind, rid=rid, id=n[0], name="s%d" % n[0], barcode=bc, pos=start, mrid=rid,
					mpos=start + 200, tlen=250, flag=99 if d else 83, positive=d, mapq=mapq, seq=_seq(rng, 50),
					qual=_qual(rng, 50)))

	start = 10
	for size in (1, 2, 3, 4, 5, 6):
		add(0, start, size, dirs=(size % 2,))
		start += 40
	add(1, 100, 300)
	add(1, 300, 4, dirs=(1, 0))
	if kind in ("se", "seb", "paf"):
		add(1, 400, 3, lengths=(50, 60, 70))
	add(3, 50, 4, mapqs=(10, 60, 10, 10))
	add(3, 150, 3, mapqs=(10, 20, 0))
	add(3, 2500, 2)
	if kind in ("seb", "peb", "sam"):
		for b in BARCODES:
			add(3, 1000, 2, barcode=b)
	return recs


DEDUP_TYPES = ["se", "seb", "pe", "peb", "paf", "ppaf", "sam"]
FORMAT_OF = {"se": "BED", "seb": "BED", "pe": "BED", "peb": "BED", "paf": "PAF", "ppaf": "PAF", "sam": "SAM"}


@pytest.mark.parametrize("kind", DEDUP_TYPES)
def test_in_memory_deduplication(cpp, fasta, tmp_path, kind):
	"""RemovePCRDuplicate (as chromap runs it in memory, sorting on 2
	threads) leaves one record per duplicate group (helpers_c.dedup_key) with
	the group size as its count, capped at 255; the MAPQ threshold applies
	after deduplication (a group with a MAPQ 60 member passes with all its
	members counted; a group below 30 throughout is dropped); groups come out
	in sort order."""
	rng = random.Random(ch.stable_seed("dedup", kind))
	records = dedup_records(kind, rng)
	rng.shuffle(records)
	barcoded = kind in ("seb", "peb", "sam")
	fmt = FORMAT_OF[kind]
	lines = _run_inmemory(cpp, kind, fmt, fasta, tmp_path / "out.txt", records, ["dedup:2", "mapq:30"],
		BC_LEN if barcoded else 0)
	if fmt == "SAM":
		lines = lines[4:]
	want = hc.expected_output(records, fmt, CONTIGS, dedup=True, mapq_threshold=30,
		barcode_length=BC_LEN if barcoded else 0)
	hc.check_output(lines, want)
	if fmt == "BED":
		counts = sorted(int(l.split("\t")[-1]) for l in lines)
		assert counts[-1] == 255
		assert counts.count(4) >= 1


@pytest.mark.parametrize("kind,fmt", [("se", "BED"), ("seb", "BED"), ("pe", "BED"), ("peb", "BED"), ("paf", "PAF"),
	("ppaf", "PAF"), ("pe", "TagAlign"), ("peb", "TagAlign")])
def test_in_memory_tn5_shift_then_dedup(cpp, fasta, tmp_path, kind, fmt):
	"""ApplyTn5ShiftOnMappings then RemovePCRDuplicate (chromap's in-memory
	order): every written interval is shifted (+4 at the start of + records
	and fragments, -5 at the end of - records and fragments; in paired
	TagAlign and PAF lines the + read loses 4 bases at its start and the -
	read 5 at its end) and duplicate groups keep their counts. The records
	place no + and - record where shifting would merge or split groups, so
	the result does not depend on whether the shift comes before or after
	deduplication."""
	rng = random.Random(ch.stable_seed("tn5", kind))
	records = dedup_records(kind, rng)
	if kind in ("se", "seb", "paf"):
		# A single-end mixed-strand group splits when shifted first; leave it
		# out (see the known-bug note c-single-end-tn5-dedup-order).
		records = [r for r in records if not (r["rid"] == 1 and r["start"] == 300)]
	bc_len = BC_LEN if kind in ("seb", "peb") else 0
	before = hc.expected_output(records, fmt, CONTIGS, tn5=True, dedup=True, mapq_threshold=30,
		barcode_length=bc_len)
	after = hc.expected_output(records, fmt, CONTIGS, tn5=True, dedup=True, mapq_threshold=30,
		tn5_before_dedup=False, barcode_length=bc_len)
	assert [c for c, _ in before] == [c for c, _ in after]
	lines = _run_inmemory(cpp, kind, fmt, fasta, tmp_path / "out.txt", records, ["tn5", "dedup:1", "mapq:30"],
		bc_len)
	hc.check_output(lines, before)


# ----------------------------------------------------------- low memory


def lowmem_records(rng):
	"""Single-end PAF records spread over four temp files (file 2 empty):
	duplicate groups split across files, a 300-record group (count 255), a
	group mixing MAPQ 10 and 60, a group below MAPQ 30, a reference with no
	records (chrC), the first group of chrB at the start of the last group of
	chrA (not duplicates: different references) and the largest record of
	the last reference alone in its own file (file 3)."""
	recs, files = [], []
	n = [0]

	def add(rid, start, size, mapqs=(60,), d=1, length=50, file_of=None):
		for i in range(size):
			n[0] += 1
			recs.append(hc.record("paf", rid=rid, id=n[0], name="q%d" % n[0], rlen=length, start=start,
				length=length, mapq=mapqs[i % len(mapqs)], dir=d))
			files.append(file_of(i) if file_of else rng.choice((0, 1)))

	for k, s in enumerate(range(5, 900, 37)):
		add(0, s, 1 + k % 4, d=k % 2)
	add(1, 893, 2)
	add(1, 1000, 300)
	add(1, 1300, 4, mapqs=(10, 60, 10, 10))
	add(1, 1500, 3, mapqs=(0, 20, 10))
	add(3, 7, 2, d=0)
	add(3, 2900, 1, d=1, length=60, file_of=lambda i: 3)
	return recs, files


@pytest.mark.parametrize("dedup,tn5,mapq", [(1, 1, 30), (0, 0, 0)])
def test_low_memory_merge_of_several_temp_files(cpp, fasta, tmp_path, dedup, tn5, mapq):
	"""ProcessAndOutputMappingsInLowMemory merges sorted temp files (one per
	OutputTempMappings call, here four, one of them empty) into one sorted
	output. With deduplication, duplicates in different files form one group
	with the whole group's count (capped at 255), the MAPQ threshold applies
	to the group (best MAPQ), and the Tn5 shift applies to what is written;
	without, every record is written once in sort order (ties broken by read
	id). The temp files are removed afterwards. PAF records are used because
	their merge buffers total 1 GiB (BED records get 10 GiB)."""
	rng = random.Random(ch.stable_seed("lowmem", dedup, tn5, mapq))
	records, files = lowmem_records(rng)
	out = tmp_path / "out.paf"
	text = "".join("%d %s\n" % (f, hc.probe_line(r)) for r, f in zip(records, files))
	cpp("c_pipeline").probe("lowmem", "paf", "PAF", fasta, out, mapq, tn5, dedup, 4, input=text)
	lines = out.read_text().splitlines()
	want = hc.expected_output(records, "PAF", CONTIGS, tn5=bool(tn5), dedup=bool(dedup), mapq_threshold=mapq,
		tn5_before_dedup=False)
	hc.check_output(lines, want)
	assert sorted(p.name for p in tmp_path.iterdir()) == ["out.paf"]
	if not dedup:
		assert len(lines) == len(records)
	else:
		assert len(lines) == len(want) and len(want) < len(records)


# ----------------------------------------------------------- temp files


def temp_records(kind, rng, layout):
	"""Records per reference for `layout`, a list of counts per reference
	(e.g. [0, 3, 0, 5]); fields vary so the round trip is visible."""
	recs = []
	n = 0
	for rid, count in enumerate(layout):
		for i in range(count):
			n += 1
			d = (n // 2) % 2
			mapq = (n * 7) % 61
			if kind == "sam":
				seq, qual = _seq(rng, 50), _qual(rng, 50)
				recs.append(hc.record("sam", rid=rid, id=n, name="read_%d/%d" % (n, 1 + n % 2),
					barcode=rng.getrandbits(32), pos=n * 13, mrid=(-1, rid)[n % 2], mpos=n * 13 + 200,
					tlen=-n, flag=(0, 99, 147, 256)[n % 4], positive=1, mapq=mapq % 61, nm=n % 5,
					cigar=("50M", "10M2I38M", "25M3D25M")[n % 3], md="%d^AC%d" % (n % 20, 25), seq=seq,
					qual=qual, dups=1 + n % 3))
			else:
				recs.append(hc.record(kind, rid=rid, id=n * 3, barcode=rng.getrandbits(32), start=n * 101 + i,
					length=50 + n % 30, mapq=mapq, dir=d, uniq=n % 2, plen=40 + n % 9, nlen=30 + n % 11,
					name="q%d" % n, name1="a%d/1" % n, name2="a%d/2" % n, rlen=50 + n % 7, rlen1=50, rlen2=60 + n % 5,
					mapq1=mapq, mapq2=(mapq + 5) % 61, dups=1 + n % 200))
	return recs


@pytest.mark.parametrize("kind", ["se", "seb", "pe", "peb", "paf", "ppaf", "sam"])
@pytest.mark.parametrize("block", [1, 2, 3, 1000])
def test_temp_file_round_trip_in_blocks(cpp, tmp_path, kind, block):
	"""A temp file written by OutputTempMappings and read back with
	TempMappingFileHandle in blocks of 1, 2, 3 or 1000 records yields every
	record once, per reference in order (rid), in the order written, with
	every field intact: blocks smaller than one reference's records,
	references with no records (first, middle and last) and the last
	reference are all handled."""
	rng = random.Random(ch.stable_seed("temp", kind, block))
	for layout in ([0, 3, 0, 5, 1, 0], [4], [1, 1, 1, 1], [0, 0, 7]):
		records = temp_records(kind, rng, layout)
		out = tmp_path / ("out_%s.txt" % "_".join(map(str, layout)))
		text = cpp("c_pipeline").probe("temp", kind, block, len(layout), out, input=hc.probe_input(records))
		got = []
		lines = text.splitlines()
		assert lines[0].startswith("path ") and lines[-1] == "end"
		for line in lines[1:-1]:
			tokens = line.split()
			got.append(hc.parse_probe_record(kind, tokens))
		fields = hc.PROBE_FIELDS[kind]
		want = [dict((f, r[f]) for f in fields) for r in sorted(records, key=lambda r: r["rid"])]
		got = [dict((f, r[f]) for f in fields) for r in got]
		assert got == want, "layout %r" % (layout,)


# ----------------------------------------------------------- allocation


def allocation_records(rng):
	"""Unique reads (MAPQ 60) in [300, 10000) of references 0-2 and
	multi-reads (MAPQ 0-3, eight locations at most) designed to fall in every
	class at every distance tested:
	- edge cases in the empty region [15000, 20000) of reference 1 around a
	  unique read at [16000, 16050): a location abutting it on the right
	  [16050, 16100) and on the left [15950, 16000) (weight 0 at distance 0,
	  1 from distance 1 on), and one overlapping it by one base [16049, 16099);
	  each read's other locations are far away;
	- a location near position 0 ([10, 60)) and a unique read at [0, 5) on
	  reference 2, reached only from distance 6 on (start - distance clips at
	  0);
	- a read whose locations are all far from unique reads (always dropped)
	  and one overlapping unique reads at [15000, 15050) of references 0 and
	  2 (always drawn between them);
	- the abutting and one-base cases again on reference 3, which has only
	  three unique reads ([1000, 1050), [2000, 2050), [3000, 3050)), so its
	  interval tree is searched leaf by leaf rather than through inner nodes;
	- 40 random multi-reads with two to five locations, some overlapping
	  unique reads, some in the empty region."""
	recs = []
	n = [0]

	def uniq(rid, start, length):
		n[0] += 1
		recs.append(hc.record("se", rid=rid, id=n[0], start=start, length=length, mapq=60, uniq=1))

	def multi(locs):
		n[0] += 1
		for rid, start in locs:
			recs.append(hc.record("se", rid=rid, id=n[0], start=start, length=50, mapq=rng.randint(0, 3), uniq=0))

	for _ in range(300):
		uniq(rng.randrange(3), rng.randrange(300, 9900), rng.randint(30, 80))
	uniq(1, 16000, 50)
	uniq(2, 0, 5)
	uniq(0, 15000, 50)
	uniq(2, 15000, 50)
	for start in (1000, 2000, 3000):
		uniq(3, start, 50)
	far = iter(range(17000, 400000, 600))
	multi([(1, 16050), (0, next(far))])
	multi([(1, 15950), (2, next(far))])
	multi([(1, 16049), (0, next(far)), (2, next(far))])
	multi([(2, 10), (1, next(far))])
	multi([(0, next(far)), (2, next(far))])
	multi([(0, 15010), (2, 15020), (1, next(far))])
	multi([(3, 1050), (0, next(far))])
	multi([(3, 1950), (0, next(far))])
	multi([(3, 2999), (0, next(far))])
	for _ in range(40):
		k = rng.randint(2, 5)
		locs = []
		for _ in range(k):
			if rng.random() < 0.6:
				locs.append((rng.randrange(3), rng.randrange(0, 9900)))
			else:
				locs.append((rng.randrange(3), next(far) if rng.random() < 0.3 else rng.randrange(10000, 14000)))
		multi(locs)
	return recs


def _allocate(cpp, records, seed, distance):
	text = cpp("c_pipeline").probe("allocate", 4, seed, distance, input=hc.probe_input(records))
	out = []
	for line in text.splitlines():
		tokens = line.split()
		out.append(hc.parse_probe_record("se", tokens))
	return out


@pytest.fixture(scope="module")
def allocation_case():
	rng = random.Random(ch.stable_seed("test_c_pipeline", "allocation"))
	return allocation_records(rng)


@pytest.mark.parametrize("distance", [0, 1, 25, 200])
def test_multi_mapping_allocation_classes(cpp, allocation_case, distance):
	"""AllocateMultiMappings keeps every unique mapping (MAPQ >= 4) and gives
	each multi-read (MAPQ < 4) one of its locations, drawn with weights equal
	to the number of unique mappings overlapping the location widened by the
	allocation distance (half-open intervals): a read whose weights are all 0
	is dropped, a read with one positive weight goes there, and otherwise it
	goes to one of its positive-weight locations. The abutting and clipped
	edge cases switch class exactly at distance 1 and 6, on a large and a
	small interval tree."""
	records = allocation_case
	uniques, weights = hc.allocation_classes(records, distance)
	got = _allocate(cpp, records, 11, distance)
	key = lambda r: (r["rid"], r["start"], r["length"], r["mapq"], r["id"])
	got_uniques = sorted(key(r) for r in got if r["mapq"] >= hc.MIN_UNIQUE_MAPQ)
	assert got_uniques == sorted(key(r) for v in uniques.values() for r in v)
	placed = collections.defaultdict(list)
	for r in got:
		if r["mapq"] < hc.MIN_UNIQUE_MAPQ:
			placed[r["id"]].append(r)
	classes = collections.Counter()
	for read_id, weighted in weights.items():
		cls = hc.classify(weighted)
		classes[cls] += 1
		if cls == "dropped":
			assert read_id not in placed, read_id
			continue
		assert len(placed[read_id]) == 1, (read_id, placed[read_id])
		loc = placed[read_id][0]
		allowed = [key(l) for l, w in weighted if w > 0]
		assert key(loc) in allowed, (read_id, loc, weighted)
	assert classes["dropped"] and classes["fixed"] and classes["random"], classes
	edge = list(weights)[:9]
	abut = ["dropped", "dropped", "fixed"] if distance == 0 else ["fixed", "fixed", "fixed"]
	near_zero = "fixed" if distance >= 6 else "dropped"
	assert [hc.classify(weights[i]) for i in edge] == abut + [near_zero, "dropped", "random"] + abut
	for i, start in zip(edge, (16050, 15950, 16049, 10, None, None, 1050, 1950, 2999)):
		if start is not None and hc.classify(weights[i]) == "fixed":
			assert placed[i][0]["start"] == start
	assert placed[edge[5]][0]["start"] in (15010, 15020)


def test_multi_mapping_allocation_seed(cpp, allocation_case):
	"""--multi-mapping-allocation-seed: the same seed gives the same
	allocation; another seed draws differently for at least one of the
	multi-reads with several positive-weight locations (the chance that two
	independent draws agree on all of them, the product over those reads of
	the sum of squared location probabilities, is checked to be below 1e-6),
	and never changes reads with one or no positive-weight location."""
	records = allocation_case
	_, weights = hc.allocation_classes(records, 200)
	key = lambda r: (r["id"], r["rid"], r["start"])
	a = sorted(key(r) for r in _allocate(cpp, records, 11, 200))
	b = sorted(key(r) for r in _allocate(cpp, records, 11, 200))
	c = sorted(key(r) for r in _allocate(cpp, records, 987654, 200))
	assert a == b
	random_ids = set(i for i, w in weights.items() if hc.classify(w) == "random")
	agree = 1.0
	for i in random_ids:
		total = float(sum(x for _, x in weights[i]))
		agree *= sum((x / total) ** 2 for _, x in weights[i])
	assert agree < 1e-6
	assert [x for x in a if x[0] not in random_ids] == [x for x in c if x[0] not in random_ids]
	assert [x for x in a if x[0] in random_ids] != [x for x in c if x[0] in random_ids]


@pytest.mark.parametrize("seed", [11, 12345, 987654])
def test_multi_mapping_allocation_follows_weights(cpp, seed):
	"""Each multi-read goes to one location drawn with probability
	proportional to its weight: 400 multi-reads, each with one location
	under 9 unique mappings and one under 1 (distance 0), should go to the
	first about 90% of the time (binomial: 360, standard deviation 6). The
	count must lie in [330, 385] at each seed; a draw that ignored the
	weights' sizes (uniform among positive weights) would give about 200,
	one that always took the heaviest 400."""
	recs = []
	n = 0
	for i in range(9):
		n += 1
		recs.append(hc.record("se", rid=0, id=n, start=5000 + i, length=50, mapq=60))
	n += 1
	recs.append(hc.record("se", rid=1, id=n, start=8000, length=50, mapq=60))
	for _ in range(400):
		n += 1
		recs.append(hc.record("se", rid=0, id=n, start=5010, length=50, mapq=0, uniq=0))
		recs.append(hc.record("se", rid=1, id=n, start=8010, length=50, mapq=0, uniq=0))
	_, weights = hc.allocation_classes(recs, 0)
	assert len(weights) == 400
	assert all(sorted(w for _, w in v) == [1, 9] for v in weights.values())
	got = _allocate(cpp, recs, seed, 0)
	multi = [r for r in got if r["mapq"] < hc.MIN_UNIQUE_MAPQ]
	assert len(multi) == 400 and len(set(r["id"] for r in multi)) == 400
	heavy = sum(1 for r in multi if r["rid"] == 0)
	assert 330 <= heavy <= 385, heavy


# ----------------------------------------------------------- summary


def _summary(cpp, tmp_path, commands, fmt="BED", barcode_length=BC_LEN, whitelist=0, not_in_whitelist=0):
	out, csv = tmp_path / "out.txt", tmp_path / "summary.csv"
	cpp("c_pipeline").probe("summary", fmt, barcode_length, whitelist, not_in_whitelist, out, csv,
		input="".join(c + "\n" for c in commands))
	return ch.parse_summary(csv)


DEFAULT_COEFFS = (-1.0996, 4.2391, 3.0164e-05, -2.1087e-04, -5.5825e-05)


def test_summary_columns_and_estimated_frip(cpp, tmp_path):
	"""One row per barcode (its sequence, BARCODE_LENGTH bases), columns
	total, duplicate, unmapped = total - mapped, lowmapq, cachehit, fric =
	cachehit / mapped, estfrip = inverse logit of c0 + c1 fric + c2 duplicate
	+ c3 unmapped + c4 lowmapq (README.md), 0 when fric is 0, numcacheslots;
	fractions printed with five decimals."""
	bc1, bc2, bc3 = "ACGTACGTACGTACGT", "TTTTTTTTTTTTTTTA", "GATTACAGATTACAGA"
	counts = {
		bc1: dict(total=1000, dup=120, mapped=900, lowmapq=40, cachehit=450, cardinality=77),
		bc2: dict(total=10, dup=0, mapped=4, lowmapq=1, cachehit=0, cardinality=0),
		bc3: dict(total=5000, dup=2500, mapped=4999, lowmapq=0, cachehit=1, cardinality=3),
	}
	commands = []
	for bc, c in counts.items():
		for field, value in c.items():
			# split each count over two updates to check accumulation
			commands.append("count %d %s %d" % (ch.encode_seed(bc), field, value - value // 3))
			commands.append("count %d %s %d" % (ch.encode_seed(bc), field, value // 3))
	coeffs = (0.5, -2.0, 0.001, -0.0004, 0.01)
	commands.append("output %s 1" % " ".join(repr(x) for x in coeffs))
	s = _summary(cpp, tmp_path, commands)
	assert s.header == hc.SUMMARY_HEADER + ["numcacheslots"]
	assert set(s.rows) == set(counts)
	for bc, c in counts.items():
		row = s[bc]
		fric = c["cachehit"] / float(c["mapped"])
		unmapped = c["total"] - c["mapped"]
		assert (row["total"], row["duplicate"], row["unmapped"], row["lowmapq"], row["cachehit"],
			row["numcacheslots"]) == (c["total"], c["dup"], unmapped, c["lowmapq"], c["cachehit"], c["cardinality"])
		assert abs(row["fric"] - fric) <= 5e-6
		assert abs(row["estfrip"] - hc.estimated_frip(coeffs, fric, c["dup"], unmapped, c["lowmapq"])) <= 5e-6
	assert s[bc2]["estfrip"] == 0.0
	line = [l for l in s.lines if l.startswith(bc1)][0]
	assert line.split(",")[6:8] == ["%.5f" % (450 / 900.0), "%.5f" % hc.estimated_frip(coeffs, 0.5, 120, 100, 40)]


def test_summary_default_coefficients_and_bulk_row(cpp, tmp_path):
	"""Bulk data has one row with an empty barcode (barcode length 0); with
	chromap's default --frip-est-params coefficients the estimate follows the
	same formula; without the cache-slot column the header has 8 columns."""
	commands = ["count 0 total 2000", "count 0 mapped 1500", "count 0 dup 300", "count 0 lowmapq 25",
		"count 0 cachehit 600", "output %s 0" % " ".join(repr(x) for x in DEFAULT_COEFFS)]
	s = _summary(cpp, tmp_path, commands, barcode_length=0)
	assert s.header == hc.SUMMARY_HEADER
	assert list(s.rows) == [""]
	row = s[""]
	assert (row["total"], row["duplicate"], row["unmapped"], row["lowmapq"], row["cachehit"]) == (2000, 300, 500, 25,
		600)
	assert abs(row["estfrip"] - hc.estimated_frip(DEFAULT_COEFFS, 0.4, 300, 500, 25)) <= 5e-6


@pytest.mark.parametrize("whitelist,not_in_whitelist,has_row", [(1, 0, True), (1, 1, False), (0, 0, False)])
def test_summary_non_whitelist_row(cpp, tmp_path, whitelist, not_in_whitelist, has_row):
	"""With a whitelist (and without --output-mappings-not-in-whitelist) the
	summary ends with a 'non-whitelist' row holding the reads whose barcodes
	were not corrected into the whitelist; otherwise there is no such row."""
	commands = ["count %d total 10" % ch.encode_seed("ACGTACGTACGTACGT"), "count %d mapped 8" %
		ch.encode_seed("ACGTACGTACGTACGT"), "nonwl total 7", "nonwl total 2", "output 0 0 0 0 0 1"]
	s = _summary(cpp, tmp_path, commands, whitelist=whitelist, not_in_whitelist=not_in_whitelist)
	assert ("non-whitelist" in s) == has_row
	if has_row:
		assert s.lines[-1].startswith("non-whitelist,")
		row = s["non-whitelist"]
		assert (row["total"], row["duplicate"], row["unmapped"], row["lowmapq"]) == (9, 0, 9, 0)
	assert s["ACGTACGTACGTACGT"]["unmapped"] == 2


@pytest.mark.parametrize("fmt,halved", [("SAM", True), ("BED", False)])
def test_summary_paired_end_sam_halving(cpp, tmp_path, fmt, halved):
	"""Paired-end SAM writes two records per pair, so chromap halves the
	duplicate, mapped and lowmapq counts (integer division) before writing
	the summary (AdjustSummaryPairedEndOverCount); total and cachehit are
	per pair already. Other formats are left alone."""
	key = ch.encode_seed("CCCCAAAAGGGGTTTT")
	commands = ["count %d total 50" % key, "count %d dup 14" % key, "count %d mapped 81" % key,
		"count %d lowmapq 7" % key, "count %d cachehit 3" % key, "adjust", "output 0 0 0 0 0 1"]
	row = _summary(cpp, tmp_path, commands, fmt=fmt)["CCCCAAAAGGGGTTTT"]
	if halved:
		assert (row["total"], row["duplicate"], row["unmapped"], row["lowmapq"], row["cachehit"]) == (50, 7, 10, 3, 3)
		assert abs(row["fric"] - 3 / 40.0) <= 5e-6
	else:
		# mapped 81 > total 50 here, so only the counts that stay defined
		assert (row["total"], row["duplicate"], row["lowmapq"], row["cachehit"]) == (50, 14, 7, 3)
		assert abs(row["fric"] - 3 / 81.0) <= 5e-6
