// Area C: K_MinHash (src/chromap.h), the sketch behind the summary's
// numcacheslots column. It keeps the k smallest distinct values added and
// estimates the number of distinct values as k * range / (largest kept) - 1,
// or 0 while fewer than k values are kept. Values below are worked by hand.

#include "chromap_test.h"

#include "chromap.h"

using namespace chromap;

// k = 3, range 100: 50, 10, 30, 10 (repeat), 20 keep {10, 20, 30};
// 3 * 100 / 30 - 1 = 9.
CT_TEST(minhash_keeps_k_smallest_distinct) {
  K_MinHash sketch(3, 100);
  sketch.add(50);
  sketch.add(10);
  sketch.add(30);
  sketch.add(10);
  CT_CHECK_EQ(sketch.compute_cardinality(), (size_t)(3 * 100 / 50 - 1));  // {10, 30, 50}
  sketch.add(20);
  CT_CHECK_EQ(sketch.compute_cardinality(), (size_t)9);
  // A value larger than every kept one changes nothing; an evicted value
  // added again is evicted again.
  sketch.add(99);
  sketch.add(50);
  CT_CHECK_EQ(sketch.compute_cardinality(), (size_t)9);
}

// Fewer than k distinct values: 0, however often they repeat.
CT_TEST(minhash_below_k_is_zero) {
  K_MinHash sketch(3, 1000);
  CT_CHECK_EQ(sketch.compute_cardinality(), (size_t)0);
  sketch.add(5);
  sketch.add(5);
  sketch.add(7);
  sketch.add(7);
  CT_CHECK_EQ(sketch.compute_cardinality(), (size_t)0);
  sketch.add(900);
  // {5, 7, 900}: 3 * 1000 / 900 - 1 = 3 - 1 = 2 (integer division).
  CT_CHECK_EQ(sketch.compute_cardinality(), (size_t)2);
}

// Integer arithmetic: k = 2, range 1000, {40, 80}: 2000 / 80 - 1 = 24; with
// k = 250 and the values 1..250 the estimate is 250 * R / 250 - 1 = R - 1.
CT_TEST(minhash_estimate_arithmetic) {
  K_MinHash a(2, 1000);
  a.add(80);
  a.add(40);
  CT_CHECK_EQ(a.compute_cardinality(), (size_t)24);
  K_MinHash b(250, 4000003);
  for (size_t v = 250; v >= 1; --v) b.add(v);
  CT_CHECK_EQ(b.compute_cardinality(), (size_t)4000002);
  // Adding 251..1000 keeps the 250 smallest: the estimate does not move.
  for (size_t v = 251; v <= 1000; ++v) b.add(v);
  CT_CHECK_EQ(b.compute_cardinality(), (size_t)4000002);
}

CT_MAIN
