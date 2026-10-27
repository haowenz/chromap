"""Fixtures for the chromap unit tests; tests/README.md documents the API."""

import os
import subprocess

import pytest

import chromap_build
import chromap_helpers as ch


def pytest_report_header(config):
	if os.environ.get("CHROMAP_BIN"):
		where = "CHROMAP_BIN=%s" % os.environ["CHROMAP_BIN"]
	else:
		where = "built from src/ into %s" % chromap_build.cache_root()
	return ["chromap binary: %s" % where, "CHROMAP_TEST_THREADS: %s" % ch.default_threads()]


@pytest.fixture(scope="session", autouse=True)
def _run_root(tmp_path_factory):
	"""Runs go under pytest's base temp directory (local disk by default)."""
	ch.set_run_root(tmp_path_factory.mktemp("runs"))


@pytest.fixture(scope="session")
def chromap_bin():
	"""Path of the chromap binary ($CHROMAP_BIN, or built from src/ on first
	use). A build failure is an error, not a skip."""
	return ch.find_chromap()


@pytest.fixture(scope="session")
def chromap_threads():
	"""The default -t of this session ($CHROMAP_TEST_THREADS or 1)."""
	return ch.default_threads()


@pytest.fixture(scope="session")
def genome(tmp_path_factory):
	"""The default SyntheticGenome with genome.fa written."""
	return ch.make_synthetic_genome().write(tmp_path_factory.mktemp("genome"))


@pytest.fixture(scope="session")
def index_factory(genome, chromap_bin, tmp_path_factory):
	"""make(k=None, w=None, genome=None, extra_args=()) -> ChromapIndex, built
	once per distinct argument set for the session. `genome` is a written
	SyntheticGenome (default: the session genome)."""
	cache = {}
	default_genome = genome

	def make(k=None, w=None, genome=None, extra_args=()):
		g = genome or default_genome
		key = (id(g), k, w, tuple(extra_args))
		if key not in cache:
			cache[key] = ch.build_index(g, tmp_path_factory.mktemp("index"), k=k, w=w, extra_args=extra_args)
		return cache[key]

	return make


@pytest.fixture(scope="session")
def index(index_factory):
	"""Index of the default genome with chromap's default k=17, w=7."""
	return index_factory()


@pytest.fixture
def sim(genome, request):
	"""ReadSimulator seeded from the test's module and name, so a test draws
	the same reads however the suite is invoked or ordered."""
	module = os.path.splitext(os.path.basename(request.node.module.__file__))[0]
	return ch.ReadSimulator(genome, seed=ch.stable_seed(module, request.node.name))


@pytest.fixture(scope="session")
def whitelist(tmp_path_factory):
	"""(barcodes, path): 32 cell barcodes of 16 bp, pairwise >= 3 mismatches
	apart, written one per line."""
	barcodes = ch.make_whitelist(32, seed=11)
	path = ch.write_lines(tmp_path_factory.mktemp("whitelist") / "whitelist.txt", barcodes)
	return barcodes, path


class CppDriver(object):
	"""A built C++ test driver (tests/cpp/<name>.cc).

	tests(), probes(): names it defines; check(test): run one CT_TEST and fail
	with its messages unless it passes; probe(name, *args, input=None) -> the
	probe's stdout (fails on a non-zero exit); run(argv, input_text=None) ->
	the CompletedProcess, for tests that expect a failure."""

	def __init__(self, name, path):
		self.name = name
		self.path = path

	def run(self, args, input_text=None, timeout=120):
		return subprocess.run([str(self.path)] + [str(a) for a in args],
			input=None if input_text is None else input_text.encode(), stdout=subprocess.PIPE,
			stderr=subprocess.PIPE, env=ch.subprocess_env(ch.default_threads()), timeout=timeout)

	def tests(self):
		return self.run(["--list"]).stdout.decode().split()

	def probes(self):
		return self.run(["--list-probes"]).stdout.decode().split()

	def check(self, test):
		proc = self.run(["--run", test])
		if proc.returncode != 0:
			pytest.fail("C++ test %s::%s failed (exit %d):\n%s%s" % (self.name, test, proc.returncode,
				proc.stderr.decode(errors="replace")[-4000:], proc.stdout.decode(errors="replace")[-2000:]),
				pytrace=False)

	def probe(self, name, *args, **kwargs):
		input_text = kwargs.pop("input", None)
		if kwargs:
			raise TypeError("unexpected arguments %s" % sorted(kwargs))
		proc = self.run(["--probe", name] + list(args), input_text=input_text)
		if proc.returncode != 0:
			pytest.fail("C++ probe %s::%s failed (exit %d):\n%s" % (self.name, name, proc.returncode,
				proc.stderr.decode(errors="replace")[-4000:]), pytrace=False)
		return proc.stdout.decode()


@pytest.fixture(scope="session")
def cpp():
	"""cpp(name) -> CppDriver for tests/cpp/<name>.cc, compiled with the
	Makefile's compiler and flags against src/ objects. Every driver is built
	in parallel on first use and cached across sessions; a build failure
	fails the tests of that driver only."""
	cache = {}
	try:
		chromap_build.prebuild_drivers()
	except chromap_build.BuildError:
		pass

	def get(name):
		if name not in cache:
			cache[name] = CppDriver(name, chromap_build.cpp_driver(name))
		return cache[name]

	return get
