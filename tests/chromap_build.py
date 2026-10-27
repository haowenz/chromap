"""Building chromap and the C++ unit-test drivers for the test suite.

The compiler and flags are the Makefile's own: `make` is asked to print CXX,
CXXFLAGS and LDFLAGS, with any of CHROMAP_TEST_CXX, CHROMAP_TEST_CXXFLAGS and
CHROMAP_TEST_LDFLAGS passed on its command line as overrides, exactly as one
would pass `make CXX=clang++`. Builds go to a cache directory keyed by a hash
of src/, the Makefile and the toolchain, so a source change always rebuilds
and an unchanged tree never does. Nothing is written into the repository.

Environment:
	CHROMAP_BIN            an already built chromap; the CLI tests use it as is
	CHROMAP_TEST_SRC_DIR   build chromap and the C++ drivers from this src/
	                       directory instead of the checkout's (e.g. a mutant
	                       copy); flags still come from the checkout's Makefile
	CHROMAP_TEST_BUILD_DIR cache directory (default <tempdir>/chromap-tests-<user>)
	CHROMAP_TEST_BUILD_JOBS parallel compile jobs (default 4)
	CHROMAP_TEST_CXX, CHROMAP_TEST_CXXFLAGS, CHROMAP_TEST_LDFLAGS  make overrides
"""

import getpass
import hashlib
import os
import shlex
import shutil
import subprocess
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


def _src_dir():
	"""The src/ directory to build: $CHROMAP_TEST_SRC_DIR when set (it must
	hold chromap's sources), else this checkout's src/."""
	env = os.environ.get("CHROMAP_TEST_SRC_DIR")
	if not env:
		return REPO_ROOT / "src"
	p = Path(env).expanduser().resolve()
	if not (p / "chromap.h").is_file():
		raise RuntimeError("CHROMAP_TEST_SRC_DIR=%s does not hold chromap's sources (no chromap.h)" % env)
	return p


SRC_DIR = _src_dir()
TESTS_DIR = Path(__file__).resolve().parent
CPP_DIR = TESTS_DIR / "cpp"
MAKEFILE = REPO_ROOT / "Makefile"

# Sources linked into the C++ unit-test drivers: every .cc in src/ except the
# one holding main(). chromap_driver.cc also instantiates every mapping
# template and takes ~15 s at -O3, so leaving it out keeps unit builds fast.
DRIVER_SOURCE = "chromap_driver.cc"

_LOCK = threading.Lock()
_CACHE = {}


class BuildError(AssertionError):
	"""A compile or link step failed; the message carries the command and its
	output."""


def build_jobs():
	n = int(os.environ.get("CHROMAP_TEST_BUILD_JOBS", "4"))
	return max(1, n)


def cache_root():
	root = os.environ.get("CHROMAP_TEST_BUILD_DIR")
	if root:
		return Path(root).expanduser().resolve()
	try:
		user = getpass.getuser()
	except Exception:
		user = "user"
	return Path(tempfile.gettempdir()) / ("chromap-tests-" + user)


def make_overrides():
	"""VAR=value arguments for make from the CHROMAP_TEST_* variables."""
	out = []
	for var in ("CXX", "CXXFLAGS", "LDFLAGS"):
		val = os.environ.get("CHROMAP_TEST_" + var)
		if val is not None:
			out.append("%s=%s" % (var, val))
	return out


def _run(cmd, cwd=None, what="command"):
	proc = subprocess.run(cmd, cwd=None if cwd is None else str(cwd), stdout=subprocess.PIPE,
		stderr=subprocess.STDOUT)
	if proc.returncode != 0:
		raise BuildError("%s failed (exit %d)\n$ %s\n%s" % (what, proc.returncode,
			" ".join(shlex.quote(str(c)) for c in cmd), proc.stdout.decode(errors="replace")[-6000:]))
	return proc.stdout.decode(errors="replace")


def toolchain():
	"""dict(cxx, cxxflags, ldflags, version): the values the repository's
	Makefile uses, after CHROMAP_TEST_* overrides. Cached per session."""
	with _LOCK:
		if "toolchain" in _CACHE:
			return _CACHE["toolchain"]
	make = shutil.which("make")
	if make is None:
		raise BuildError("`make` is not on PATH; set CHROMAP_BIN to a built chromap")
	fd, helper = tempfile.mkstemp(suffix=".mk")
	try:
		with os.fdopen(fd, "w") as f:
			f.write("_chromap_test_print:\n"
				"\t@echo 'CXX=$(CXX)'\n\t@echo 'CXXFLAGS=$(CXXFLAGS)'\n\t@echo 'LDFLAGS=$(LDFLAGS)'\n")
		out = _run([make, "-s", "--no-print-directory", "-f", str(MAKEFILE), "-f", helper,
			"_chromap_test_print"] + make_overrides(), cwd=REPO_ROOT, what="make (print flags)")
	finally:
		os.unlink(helper)
	vals = {}
	for line in out.splitlines():
		if "=" in line:
			k, v = line.split("=", 1)
			vals[k.strip()] = v.strip()
	cxx = shlex.split(vals.get("CXX", "g++"))
	version = ""
	try:
		version = subprocess.run(cxx + ["--version"], stdout=subprocess.PIPE,
			stderr=subprocess.STDOUT).stdout.decode(errors="replace").splitlines()[0]
	except (OSError, IndexError):
		pass
	tc = dict(cxx=cxx, cxxflags=shlex.split(vals.get("CXXFLAGS", "")),
		ldflags=shlex.split(vals.get("LDFLAGS", "")), version=version)
	with _LOCK:
		_CACHE["toolchain"] = tc
	return tc


def _digest(paths, extra=()):
	h = hashlib.sha1()
	for p in paths:
		h.update(str(Path(p).name).encode())
		h.update(b"\0")
		h.update(Path(p).read_bytes())
		h.update(b"\0")
	for e in extra:
		h.update(str(e).encode())
		h.update(b"\0")
	return h.hexdigest()[:16]


def source_files():
	return sorted(p for p in SRC_DIR.iterdir() if p.is_file() and p.suffix in (".cc", ".h", ".hpp"))


def _tool_key(tc):
	return [" ".join(tc["cxx"]), " ".join(tc["cxxflags"]), " ".join(tc["ldflags"]), tc["version"]]


def _in_tmp(final, prefix, build):
	"""Run build(tmp) in a fresh directory next to `final` and move it into
	place; the directory is removed if the build fails."""
	final.parent.mkdir(parents=True, exist_ok=True)
	tmp = Path(tempfile.mkdtemp(prefix=prefix, dir=str(final.parent)))
	try:
		build(tmp)
	except BaseException:
		shutil.rmtree(str(tmp), ignore_errors=True)
		raise
	_publish(tmp, final)


def _publish(tmp, final):
	"""Move a finished build directory into place; another session may have
	won the race, in which case its result is kept."""
	try:
		os.rename(str(tmp), str(final))
	except OSError:
		shutil.rmtree(str(tmp), ignore_errors=True)
		if not final.exists():
			raise


def find_chromap():
	"""Path of the chromap binary under test: $CHROMAP_BIN when set, else a
	build of this repository's sources in the cache directory (built on first
	use with the Makefile's flags)."""
	env = os.environ.get("CHROMAP_BIN")
	if env:
		p = Path(env).expanduser()
		if not (p.is_file() and os.access(str(p), os.X_OK)):
			raise BuildError("CHROMAP_BIN=%s is not an executable file" % env)
		return p.resolve()
	return _cached("chromap", _build_chromap)


def _build_chromap():
	tc = toolchain()
	key = _digest(source_files() + [MAKEFILE], _tool_key(tc))
	final = cache_root() / ("bin-" + key)
	exe = final / "chromap"
	if not exe.exists():
		_in_tmp(final, "bin-tmp-", lambda tmp: _run([shutil.which("make"), "-s", "--no-print-directory",
			"-f", str(MAKEFILE), "-C", str(tmp), "src_dir=%s" % SRC_DIR, "objs_dir=%s" % (tmp / "objs"),
			"exec=%s" % (tmp / "chromap"), "-j%d" % build_jobs()] + make_overrides(), what="building chromap"))
	return exe


def unit_objects():
	"""Object files of src/*.cc (all but chromap_driver.cc), compiled with the
	Makefile's CXX and CXXFLAGS into the cache. Returns a list of paths."""
	return _cached("objects", _build_unit_objects)


def _build_unit_objects():
	tc = toolchain()
	sources = [p for p in source_files() if p.suffix == ".cc" and p.name != DRIVER_SOURCE]
	key = _digest(source_files(), _tool_key(tc) + ["units"])
	final = cache_root() / ("units-" + key)
	objs = [final / (p.stem + ".o") for p in sources]
	if not all(o.exists() for o in objs):
		def build(tmp):
			def compile_one(src):
				_run(tc["cxx"] + tc["cxxflags"] + ["-c", str(src), "-o", str(tmp / (src.stem + ".o"))],
					what="compiling %s" % src.name)

			with ThreadPoolExecutor(max_workers=build_jobs()) as pool:
				for fut in [pool.submit(compile_one, s) for s in sources]:
					fut.result()

		_in_tmp(final, "units-tmp-", build)
	return objs


def _cached(key, build):
	"""build() once per session; a failure is remembered and re-raised, so a
	broken build is reported once per test without recompiling."""
	with _LOCK:
		if key in _CACHE:
			val = _CACHE[key]
			if isinstance(val, BuildError):
				raise BuildError(str(val))
			return val
	try:
		val = build()
	except BuildError as e:
		with _LOCK:
			_CACHE[key] = e
		raise
	with _LOCK:
		_CACHE[key] = val
	return val


def cpp_driver(name):
	"""Build tests/cpp/<name>.cc against the unit objects and return the
	executable's path. The link line is the Makefile's: CXX CXXFLAGS objs -o
	exe LDFLAGS, with -I src -I tests/cpp added."""
	return _cached(("driver", name), lambda: _build_driver(name))


def _build_driver(name):
	src = CPP_DIR / (name + ".cc")
	if not src.exists():
		raise BuildError("no C++ test driver %s" % src)
	tc = toolchain()
	objs = unit_objects()
	headers = sorted(p for p in CPP_DIR.iterdir() if p.suffix in (".h", ".hpp"))
	key = _digest([src] + headers + source_files(), _tool_key(tc))
	final = cache_root() / ("driver-%s-%s" % (name, key))
	exe = final / name
	if not exe.exists():
		_in_tmp(final, "driver-tmp-", lambda tmp: _run(tc["cxx"] + tc["cxxflags"]
			+ ["-I", str(SRC_DIR), "-I", str(CPP_DIR), str(src)] + [str(o) for o in objs]
			+ ["-o", str(tmp / name)] + tc["ldflags"], what="building C++ driver %s" % name))
	return exe


def driver_names():
	"""Names of every driver in tests/cpp (the .cc files)."""
	return sorted(p.stem for p in CPP_DIR.glob("*.cc"))


def prebuild_drivers():
	"""Build every driver in tests/cpp in parallel. A driver that fails to
	build is remembered and fails only the tests that use it."""
	names = driver_names()
	unit_objects()

	def build(name):
		try:
			cpp_driver(name)
		except BuildError:
			pass

	with ThreadPoolExecutor(max_workers=build_jobs()) as pool:
		list(pool.map(build, names))
