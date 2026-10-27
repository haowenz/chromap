// Minimal self-contained test framework for chromap's C++ unit-test drivers.
//
// A driver is one .cc file under tests/cpp/ that includes this header once
// and defines tests and probes:
//
//   CT_TEST(name) { CT_CHECK(cond); CT_CHECK_EQ(a, b); }
//       A self-checking test. Run with `driver --run name`; exit 0 on pass,
//       1 with file:line messages on stderr on failure.
//   CT_PROBE(name) { ...; return 0; }
//       Exposes a C++ function to pytest: `driver --probe name args...` with
//       `args` (std::vector<std::string>) and std::cin/std::cout for data.
//       Python computes the expected values and compares.
//
// `driver --list` prints every test name, one per line; `driver --list-probes`
// every probe. Without arguments every test runs. Only the C++11 standard
// library is used, so a driver builds with the repository's own toolchain.

#ifndef CHROMAP_TEST_H_
#define CHROMAP_TEST_H_

// Several src/ headers use standard names without including their headers
// (bed_mapping.h uses std::tie, mapping_processor.h std::mt19937 and
// std::discrete_distribution) and rely on what was included before them.
// Include this header first in every driver.
#include <algorithm>
#include <cassert>
#include <cinttypes>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <functional>
#include <iostream>
#include <limits>
#include <queue>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace chromap_test {

typedef void (*TestFunction)();
typedef int (*ProbeFunction)(const std::vector<std::string> &args);

struct Registry {
  std::vector<std::pair<std::string, TestFunction> > tests;
  std::vector<std::pair<std::string, ProbeFunction> > probes;
  int failures = 0;
};

inline Registry &GetRegistry() {
  static Registry registry;
  return registry;
}

struct TestRegistrar {
  TestRegistrar(const char *name, TestFunction f) {
    GetRegistry().tests.push_back(std::make_pair(std::string(name), f));
  }
};

struct ProbeRegistrar {
  ProbeRegistrar(const char *name, ProbeFunction f) {
    GetRegistry().probes.push_back(std::make_pair(std::string(name), f));
  }
};

inline void ReportFailure(const char *file, int line, const std::string &what) {
  ++GetRegistry().failures;
  std::cerr << file << ":" << line << ": FAILED: " << what << std::endl;
}

template <typename T>
std::string Show(const T &value) {
  std::ostringstream os;
  os << value;
  return os.str();
}

inline std::string Show(const unsigned char &value) {
  return Show(static_cast<unsigned>(value));
}

inline std::string Show(const signed char &value) {
  return Show(static_cast<int>(value));
}

inline std::string Show(const std::string &value) { return "\"" + value + "\""; }

inline int RunTest(const std::string &name, TestFunction f) {
  int before = GetRegistry().failures;
  try {
    f();
  } catch (const std::exception &e) {
    ReportFailure("<exception>", 0, std::string("uncaught exception: ") + e.what());
  } catch (...) {
    ReportFailure("<exception>", 0, "uncaught non-standard exception");
  }
  bool ok = GetRegistry().failures == before;
  std::cerr << (ok ? "PASS " : "FAIL ") << name << std::endl;
  return ok ? 0 : 1;
}

inline int Main(int argc, char **argv) {
  Registry &r = GetRegistry();
  std::vector<std::string> args(argv + 1, argv + argc);
  if (args.size() == 1 && args[0] == "--list") {
    for (size_t i = 0; i < r.tests.size(); ++i) std::cout << r.tests[i].first << "\n";
    return 0;
  }
  if (args.size() == 1 && args[0] == "--list-probes") {
    for (size_t i = 0; i < r.probes.size(); ++i) std::cout << r.probes[i].first << "\n";
    return 0;
  }
  if (args.size() == 2 && args[0] == "--run") {
    for (size_t i = 0; i < r.tests.size(); ++i) {
      if (r.tests[i].first == args[1]) return RunTest(args[1], r.tests[i].second);
    }
    std::cerr << "no test named " << args[1] << std::endl;
    return 2;
  }
  if (args.size() >= 2 && args[0] == "--probe") {
    for (size_t i = 0; i < r.probes.size(); ++i) {
      if (r.probes[i].first == args[1]) {
        std::vector<std::string> rest(args.begin() + 2, args.end());
        return r.probes[i].second(rest);
      }
    }
    std::cerr << "no probe named " << args[1] << std::endl;
    return 2;
  }
  if (!args.empty()) {
    std::cerr << "usage: driver [--list | --list-probes | --run NAME | --probe NAME ARGS...]" << std::endl;
    return 2;
  }
  int failed = 0;
  for (size_t i = 0; i < r.tests.size(); ++i) failed += RunTest(r.tests[i].first, r.tests[i].second);
  return failed == 0 ? 0 : 1;
}

}  // namespace chromap_test

#define CT_TEST(name)                                                      \
  static void name();                                                      \
  static chromap_test::TestRegistrar name##_test_registrar(#name, name);   \
  static void name()

#define CT_PROBE(name)                                                     \
  static int name(const std::vector<std::string> &args);                   \
  static chromap_test::ProbeRegistrar name##_probe_registrar(#name, name); \
  static int name(const std::vector<std::string> &args)

#define CT_CHECK(cond)                                                     \
  do {                                                                     \
    if (!(cond)) chromap_test::ReportFailure(__FILE__, __LINE__, #cond);   \
  } while (0)

#define CT_CHECK_EQ(a, b)                                                  \
  do {                                                                     \
    const auto &ct_a_ = (a);                                               \
    const auto &ct_b_ = (b);                                               \
    if (!(ct_a_ == ct_b_))                                                 \
      chromap_test::ReportFailure(__FILE__, __LINE__,                      \
                                  std::string(#a " == " #b ": ") +         \
                                      chromap_test::Show(ct_a_) + " vs " + \
                                      chromap_test::Show(ct_b_));          \
  } while (0)

#define CT_MAIN \
  int main(int argc, char **argv) { return chromap_test::Main(argc, argv); }

#endif  // CHROMAP_TEST_H_
