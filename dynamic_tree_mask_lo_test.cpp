#include "gtest/gtest.h"

#include <bitset>

#include "dynamic_tree_mask_lo.hpp"
#include "two_state_markov_process.hpp"

namespace dtl {
using bitmap = boost::dynamic_bitset<$u32>;
} // namespace dtl

//===----------------------------------------------------------------------===//
dtl::bitmap
gen(u64 n, $f64 f, $f64 d) {

  // init bitset
  f64 f_min = d >= 1.0 ? n : d/(1-d);
  f64 f_actual = std::max(f, f_min);
  two_state_markov_process mp(f_actual, d);
  dtl::bitmap bs(n);
  for ($u64 i = 0; i < n; i++) {
    bs[i] = mp.next();
  }

  $f64 d_actual = (bs.count() * 1.0) / n;
  if (std::abs(d - d_actual) > 1
      || std::abs(f - f_actual) > 0.25) {
    throw std::invalid_argument("Failed to construct a random bitmap with the given parameters.");
  }
  return bs;
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
TEST(dynamic_tree_mask_lo, all_bits_set) {
  dtl::bitmap b(64);
  b.set();
  dtl::dynamic_tree_mask_lo tm(b);
  std::cout << tm << std::endl;
  auto it = tm.it();
  ASSERT_TRUE(!it.end());
  ASSERT_EQ(it.pos(), 0);
  ASSERT_EQ(it.length(), 64);
  it.next();
  ASSERT_TRUE(it.end());
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
TEST(dynamic_tree_mask_lo, no_bits_set) {
  dtl::bitmap b(64);
  dtl::dynamic_tree_mask_lo tm(b);
  std::cout << tm << std::endl;
  auto it = tm.it();
  ASSERT_TRUE(it.end());
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
TEST(dynamic_tree_mask_lo, single_1fill) {
  dtl::bitmap b(64, 0b11110000);
  dtl::dynamic_tree_mask_lo tm(b);
  std::cout << tm << std::endl;
  auto it = tm.it();
  ASSERT_TRUE(!it.end());
  ASSERT_EQ(it.pos(), 4);
  ASSERT_EQ(it.length(), 4);
  it.next();
  ASSERT_TRUE(it.end());
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
TEST(dynamic_tree_mask_lo, multiple_1fills_at_different_levels) {
  dtl::bitmap b(64, 0b0011000011110000);
  dtl::dynamic_tree_mask_lo tm(b);
  std::cout << tm << std::endl;
  auto it = tm.it();
  ASSERT_TRUE(!it.end());
  ASSERT_EQ(it.pos(), 4);
  ASSERT_EQ(it.length(), 4);
  it.next();
  ASSERT_EQ(it.pos(), 12);
  ASSERT_EQ(it.length(), 2);
  ASSERT_TRUE(!it.end());
  it.next();
  ASSERT_TRUE(it.end());
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
void
skip_test(u64 n, u64 bitmap, u64 skip_to_pos, u64 expected_pos, u64 expected_len) {
  assert(n <= 64);
  dtl::bitmap b(n, bitmap);
  dtl::dynamic_tree_mask_lo tm(b);
  std::cout << tm << std::endl;
  auto it = tm.it();
  it.skip_to(skip_to_pos);
  if (expected_len == 0) {
    ASSERT_TRUE(it.end());
  }
  else {
    ASSERT_TRUE(!it.end());
  }
  ASSERT_EQ(it.pos(), expected_pos);
  if (expected_pos != n) {
    ASSERT_EQ(it.length(), expected_len);
  }
}

TEST(dynamic_tree_mask_lo, skip_to_1fill) {
  skip_test(8, 0b00001001, 3, 3, 1);
  skip_test(8, 0b00000101, 2, 2, 1);
}

TEST(dynamic_tree_mask_lo, skip_to_1fill_passing_the_root_node) {
  skip_test(8, 0b01000001, 6, 6, 1);
  skip_test(8, 0b10000001, 7, 7, 1);
}

TEST(dynamic_tree_mask_lo, skip_before_1fill) {
  skip_test(8, 0b01000001, 5, 6, 1);
  skip_test(8, 0b01000001, 4, 6, 1);
  skip_test(8, 0b01000001, 3, 6, 1);
  skip_test(8, 0b01000001, 2, 6, 1);
  skip_test(8, 0b01000001, 1, 6, 1);
}

TEST(dynamic_tree_mask_lo, skip_beyond_last_1fill) {
  skip_test(8, 0b01000001, 7, 8, 0);
}

TEST(dynamic_tree_mask_lo, skip_into_a_1fill) {
  skip_test(8, 0b11111101, 6, 6, 2);
}

TEST(dynamic_tree_mask_lo, skip_multiple_times) {
  dtl::bitmap b(8, 0b11001101);
  dtl::dynamic_tree_mask_lo tm(b);
  std::cout << tm << std::endl;
  auto it = tm.it();
  it.skip_to(2);
  ASSERT_EQ(it.pos(), 2);
  ASSERT_EQ(it.length(), 2);
  it.skip_to(6);
  ASSERT_EQ(it.pos(), 6);
  ASSERT_EQ(it.length(), 2);
}
//===----------------------------------------------------------------------===//

