"""Area D: the command-line front end (ChromapDriver::ParseArgsAndRun):
help, version, option precedence and every error path that does not need
barcodes (those are in test_d_barcodes.py, the pairs and --chr-order ones
in test_d_hic.py, --read-format in test_d_read_format.py).

Error paths exit through ExitWithMessage: the message on stderr and status
255. An option chromap does not register, or a value that does not parse,
makes the option parser throw and the process abort, so those tests assert
a non-zero status and the message only. Each error-path command line is a
complete, valid mapping command except for the one fault under test, so a
missing check would show as a mapping run that succeeds."""

import re

import pytest

import chromap_build
import chromap_helpers as ch
import helpers_d as hd

# Options chromap.1 documents (its OPTIONS section and synopsis), which the
# help must list.
DOCUMENTED_OPTIONS = [
	"-i", "-k", "-w", "--min-frag-length",
	"--preset", "--split-alignment", "-e", "-s", "-f", "-l", "-q", "--min-read-length", "--trim-adapters",
	"--remove-pcr-duplicates", "--remove-pcr-duplicates-at-bulk-level", "--remove-pcr-duplicates-at-cell-level",
	"--Tn5-shift", "--low-mem", "--bc-error-threshold", "--bc-probability-threshold", "-t",
	"-r", "-x", "-1", "-2", "-b", "--barcode-whitelist", "--read-format",
	"-o", "--output-mappings-not-in-whitelist", "--chr-order", "--BED", "--TagAlign", "--SAM", "--pairs",
	"--pairs-natural-chr-order", "--barcode-translate", "--summary", "-v",
]


@pytest.fixture(scope="module")
def inputs(genome, tmp_path_factory):
	"""Small single-end and paired-end read files and paths for the error
	paths."""
	d = tmp_path_factory.mktemp("d_cli")
	sim = ch.ReadSimulator(genome, seed=ch.stable_seed("test_d_cli", "inputs"))
	se = sim.se_reads(5).write(d, prefix="se")
	pe = sim.pairs(5).write(d, prefix="pe")
	return dict(dir=d, se=se["r1"], r1=pe["r1"], r2=pe["r2"])


def _mapping_args(index, inputs, out="out.bed"):
	return ["-x", index.path, "-r", index.fasta, "-1", inputs["se"], "-o", inputs["dir"] / out]


def _help_options(text):
	"""Option names the help lists ("-k, --kmer INT", "    --low-mem")."""
	names = set()
	for line in text.splitlines():
		m = re.match(r"^\s+(?:(-\w),\s+)?(--[\w-]+)", line)
		if m:
			names.update(x for x in m.groups() if x)
	return names


@pytest.fixture(scope="module")
def help_run():
	return ch.run_cli(["-h"])


def test_help_lists_documented_options(help_run):
	"""-h prints the usage to stderr (stdout stays empty) and exits 0; it
	lists every option chromap.1 documents, and -h itself."""
	assert help_run.returncode == 0
	assert help_run.stdout == ""
	assert "Usage:" in help_run.stderr
	listed = _help_options(help_run.stderr)
	missing = [o for o in DOCUMENTED_OPTIONS + ["-h", "--help"] if o not in listed]
	assert missing == []


@pytest.mark.parametrize("args", [["--help"], [], ["-r", "ref.fa", "-o", "out", "-t", "2"]])
def test_help_without_mode(help_run, args):
	"""--help, no arguments, or options without -i or -1 all print the same
	usage on stderr and exit 0."""
	run = ch.run_cli(args)
	assert run.returncode == 0
	assert run.stderr == help_run.stderr
	assert run.stdout == ""


def _source_version():
	text = (chromap_build.SRC_DIR / "chromap.h").read_text()
	return re.search(r'#define CHROMAP_VERSION "([^"]+)"', text).group(1)


@pytest.mark.parametrize("flag", ["-v", "--version"])
def test_version(flag):
	"""-v and --version print the version defined in src/chromap.h
	(CHROMAP_VERSION) and nothing else, and exit 0."""
	run = ch.run_cli([flag])
	assert run.returncode == 0
	assert (run.stdout + run.stderr).strip() == _source_version()


@pytest.mark.parametrize("missing", ["-r", "-o", "-x"])
def test_mapping_without_required_option_exits(index, inputs, missing):
	"""Mapping (-1) needs -r, -o and -x: leaving one out stops chromap with
	status 255 and the matching message; no output is written."""
	args = ["-1", inputs["se"]]
	if missing != "-r":
		args += ["-r", index.fasta]
	if missing != "-o":
		args += ["-o", inputs["dir"] / "req.bed"]
	if missing != "-x":
		args += ["-x", index.path]
	run = ch.run_cli(args)
	message = {"-r": "No reference specified!", "-o": "No output file specified!",
		"-x": "No index file specified!"}[missing]
	assert run.returncode == 255
	assert message in run.stderr
	assert not (inputs["dir"] / "req.bed").exists()


@pytest.mark.parametrize("missing", ["-r", "-o"])
def test_index_without_required_option_exits(index, inputs, missing):
	"""Index building (-i) needs -r and -o."""
	args = ["-i"]
	if missing != "-r":
		args += ["-r", index.fasta]
	if missing != "-o":
		args += ["-o", inputs["dir"] / "idx_req"]
	run = ch.run_cli(args)
	assert run.returncode == 255
	assert {"-r": "No reference specified!", "-o": "No output file specified!"}[missing] in run.stderr
	assert not (inputs["dir"] / "idx_req").exists()


@pytest.mark.parametrize("extra,message", [
	(["--preset", "rnaseq"], "Unrecognized preset parameters rnaseq"),
	(["--cache-size", "1999999"], "cache size is not in appropriate range"),
	(["--cache-size", "15000001"], "cache size is not in appropriate range"),
	(["--cache-update-param", "1.5"], "cache update param is not approriate"),
	(["--cache-update-param", "-0.5"], "cache update param is not approriate"),
	(["--k-for-minhash", "0"], "Invalid paramter for size of MinHash sketch (--k-for-minhash)"),
	(["--k-for-minhash", "2000"], "Invalid paramter for size of MinHash sketch (--k-for-minhash)"),
	(["--barcode-whitelist", "whitelist.txt"],
		"No barcode file specified but the barcode whitelist file is given!"),
])
def test_bad_option_value_exits(index, inputs, extra, message):
	"""Values chromap checks itself: an unknown preset, --cache-size outside
	[2000000, 15000000], --cache-update-param outside [0, 1], --k-for-minhash
	outside [1, 1999], and --barcode-whitelist without -b each stop chromap
	with status 255 and the message, before any output is written."""
	run = ch.run_cli(_mapping_args(index, inputs, "bad.bed") + extra)
	assert run.returncode == 255
	assert message in run.stderr
	assert run.stdout == ""
	assert not (inputs["dir"] / "bad.bed").exists()


@pytest.mark.parametrize("value,message", [
	("1;2", "Invalid number of parameters, expecting 5 parameters but found 2 parameters"),
	("1;2;3;4;5;6", "Invalid number of parameters, expecting 5 parameters but found 6 parameters"),
	("a;b;c;d;e", "Exception occurred while processing chromap score parameters"),
])
def test_malformed_frip_est_params_exit(index, inputs, value, message):
	"""Paired-end runs parse --frip-est-params as five numbers separated by
	semicolons; another count, or a token that is not a number, stops chromap
	with status 255 and the message."""
	run = ch.run_cli(["-x", index.path, "-r", index.fasta, "-1", inputs["r1"], "-2", inputs["r2"],
		"-o", inputs["dir"] / "frip.bed", "--frip-est-params", value])
	assert run.returncode == 255
	assert message in run.stderr


@pytest.mark.parametrize("extra", [
	["--cell-by-bin"],
	["--bin-size", "5000"],
	["--depth-cutoff", "3"],
	["--peak-min-length", "30"],
	["--peak-merge-max-length", "30"],
	["-p", "matrix"],
	["--matrix-output-prefix", "matrix"],
	["--no-such-option"],
])
def test_unregistered_option_rejected(index, inputs, extra):
	"""The peak and matrix options are defined but never registered
	(AddPeakOptions is not called, -p is commented out): chromap rejects
	them, as any unknown option, with a non-zero status and a message naming
	the option, before mapping."""
	run = hd.run_cli_expect_abort(_mapping_args(index, inputs, "unreg.bed") + extra)
	assert run.returncode != 0
	name = extra[0].lstrip("-")
	assert re.search(r"Option .%s. does not exist" % re.escape(name), run.stderr), run.stderr
	assert not (inputs["dir"] / "unreg.bed").exists()


@pytest.mark.parametrize("extra", [
	["-q", "x"],
	["-e", "1.5"],
	["-l", "12x"],
	["-t", "four"],
	["--bc-error-threshold", "two"],
	["--bc-probability-threshold", "high"],
	["--min-read-length", "thirty"],
])
def test_malformed_value_rejected(index, inputs, extra):
	"""A value that does not parse as the option's type is rejected with a
	non-zero status and a message quoting it, before mapping."""
	run = hd.run_cli_expect_abort(_mapping_args(index, inputs, "malformed.bed") + extra)
	assert run.returncode != 0
	assert re.search(r"Argument .%s. failed to parse" % re.escape(extra[1]), run.stderr), run.stderr
	assert not (inputs["dir"] / "malformed.bed").exists()


# ------------------------------------------------------------ precedence

FORMAT_LINES = {
	"BED": "Output mappings in BED/BEDPE format.",
	"TagAlign": "Output mappings in TagAlign/PairedTagAlign format.",
	"PAF": "Output mappings in PAF format.",
	"SAM": "Output mappings in SAM format.",
	"pairs": "Output mappings in pairs format.",
}


def _settings(index, inputs, extra):
	"""chromap's settings report for a paired-end command whose reference
	file does not exist: the driver reports the settings it parsed, then
	mapping stops at loading the reference (status 255)."""
	run = ch.run_cli(["-x", index.path, "-r", inputs["dir"] / "no_such_ref.fa", "-1", inputs["r1"],
		"-2", inputs["r2"], "-o", inputs["dir"] / "prec.out"] + extra)
	assert run.returncode == 255
	assert "Cannot find sequence file" in run.stderr
	return run.stderr


@pytest.mark.parametrize("flags,fmt", [
	(["--TagAlign", "--BED"], "TagAlign"),
	(["--BED", "--TagAlign"], "TagAlign"),
	(["--PAF", "--TagAlign", "--BED"], "PAF"),
	(["--BED", "--TagAlign", "--PAF"], "PAF"),
	(["--pairs", "--PAF"], "pairs"),
	(["--PAF", "--pairs"], "pairs"),
	(["--SAM", "--pairs", "--BED"], "SAM"),
	(["--pairs", "--SAM"], "SAM"),
	([], "BED"),
	(["--preset", "hic"], "pairs"),
	(["--preset", "hic", "--SAM"], "SAM"),
	(["--SAM", "--preset", "hic"], "SAM"),
	(["--preset", "atac", "--TagAlign"], "TagAlign"),
	(["--preset", "chip"], "BED"),
])
def test_output_format_precedence(index, inputs, flags, fmt):
	"""Output format flags rank SAM > pairs > PAF > TagAlign > BED whatever
	their order, BED is the default, and a preset's format yields to an
	explicit flag given before or after it (presets are applied first,
	chromap.1 --preset)."""
	err = _settings(index, inputs, flags)
	assert FORMAT_LINES[fmt] in err
	assert [k for k, v in FORMAT_LINES.items() if v in err] == [fmt]


@pytest.mark.parametrize("flags,error,mapq", [
	(["--preset", "hic"], 4, 1),
	(["--preset", "hic", "-e", "6", "-q", "20"], 6, 20),
	(["-e", "6", "-q", "20", "--preset", "hic"], 6, 20),
	(["-q", "20"], 8, 20),
])
def test_preset_applied_before_options(index, inputs, flags, error, mapq):
	"""--preset hic sets -e 4 and -q 1 (chromap.1); -e and -q given with it
	override the preset whether they come before or after it."""
	err = _settings(index, inputs, flags)
	assert "Parameters: error threshold: %d," % error in err
	assert "MAPQ-threshold: %d," % mapq in err
