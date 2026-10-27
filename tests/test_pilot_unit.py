"""Pilot of the C++ unit layer: the pattern for tests that reach chromap's
C++ classes directly through a driver in tests/cpp/.

Two kinds of check: CT_TEST cases inside the driver assert values worked by
hand (one pytest item each, listed below so none is skipped silently), and
CT_PROBE functions print what a chromap function returns so the expected
values can be computed here, independently, in Python."""

import random

import pytest

import chromap_helpers as ch

pytestmark = pytest.mark.cpp

PILOT_CPP_TESTS = [
	"seed_encoding_hand_worked",
	"hit_encoding_round_trip",
	"tn5_shift_single_end_record",
	"tn5_shift_paired_end_record",
]


def test_pilot_driver_lists_every_test(cpp):
	"""Every CT_TEST in cpp/pilot.cc is run by test_pilot_cpp, and the probes
	used below exist."""
	driver = cpp("pilot")
	assert driver.tests() == PILOT_CPP_TESTS
	assert driver.probes() == ["minimizers"]


@pytest.mark.parametrize("name", PILOT_CPP_TESTS)
def test_pilot_cpp(cpp, name):
	cpp("pilot").check(name)


def _minimizers_from_probe(cpp, sequences, k, w, tmp_path):
	fasta = tmp_path / ("seqs_k%d_w%d.fa" % (k, w))
	ch.write_lines(fasta, [x for i, s in enumerate(sequences) for x in (">s%d" % i, s)])
	out = cpp("pilot").probe("minimizers", fasta, k, w)
	got = {}
	for line in out.splitlines():
		i, h, hit = (int(x) for x in line.split("\t"))
		got.setdefault(i, []).append((h, hit))
	return got


@pytest.mark.parametrize("k,w", [(17, 7), (15, 10), (21, 11), (11, 3), (5, 1)])
def test_minimizer_generator_matches_definition(cpp, tmp_path, k, w):
	"""MinimizerGenerator returns exactly the k-mers that are the smallest
	hash in some window of w consecutive k-mers (chromap_helpers
	.reference_minimizers), with the hit encoding sequence index, end position
	and strand, in position order. Sequences are N-free random sequence of
	lengths around and well above w + k - 1, including one shorter than a
	full window (the single smallest k-mer is returned) and one exactly a
	window long."""
	rng = random.Random(ch.stable_seed("minimizers", k, w))
	lengths = [k, w + k - 2, w + k - 1, w + k, 2 * (w + k), 150, 1000]
	seqs = ["".join(rng.choice("ACGT") for _ in range(n)) for n in lengths]
	got = _minimizers_from_probe(cpp, seqs, k, w, tmp_path)
	for i, s in enumerate(seqs):
		expected = ch.reference_minimizers(s, k, w, seq_index=i)
		assert got.get(i, []) == expected, "sequence %d (length %d)" % (i, len(s))
