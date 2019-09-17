#include "gtest/gtest.h"
#include <dtl/dtl.hpp>
#include <dtl/bitmap.hpp>
#include <dtl/bitmap/util/convert.hpp>
#include <dtl/bitmap/util/random.hpp>
#include "api_types.hpp"
//===----------------------------------------------------------------------===//
// Typed API tests for loss-less compressed bitmaps.
// Iterate runs.
//===----------------------------------------------------------------------===//
constexpr std::size_t LEN = 8;
//===----------------------------------------------------------------------===//
// Fixture for the parameterized test case.
template<typename T>
class api_run_iterator_skip_test : public ::testing::Test {};
TYPED_TEST_CASE(api_run_iterator_skip_test, types_under_test);
// TODO test skip iterator and scan iterator explicitly
//===----------------------------------------------------------------------===//
template<typename T>
void
skip_test(u64 n, u64 bitmap, u64 skip_to_pos,
          u64 expected_pos, u64 expected_len) {
  assert(n <= 64);
  std::stringstream info;
  info << "skip test: n=" << n
       << ", bitmap=" << boost::dynamic_bitset<$u32>(n, bitmap)
       << " (=" << bitmap << ")"
       << ", skip_to=" << skip_to_pos
       << ", expected_pos=" << expected_pos
       << ", expected_len=" << expected_len
                            << std::endl;
  std::cout << info.str();
  dtl::bitmap b(n, bitmap);
  T t(b);
  auto it = t.it();
  it.skip_to(skip_to_pos);
  if (expected_len == 0) {
    ASSERT_TRUE(it.end()) << info.str();
  }
  else {
    ASSERT_TRUE(!it.end()) << info.str();
  }
  ASSERT_EQ(it.pos(), expected_pos) << info.str();
  if (expected_pos != n) {
//    ASSERT_EQ(it.length(), expected_len);
  }
}

TYPED_TEST(api_run_iterator_skip_test, skip_to_1fill) {
  using T = TypeParam;
  skip_test<T>(8, 0b00000101, 2, 2, 1);
  skip_test<T>(8, 0b00001001, 3, 3, 1);
  skip_test<T>(8, 0b00010001, 4, 4, 1);
  skip_test<T>(8, 0b00100001, 5, 5, 1);
  skip_test<T>(8, 0b01000001, 6, 6, 1);
  skip_test<T>(8, 0b10000001, 7, 7, 1);
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_run_iterator_skip_test, skip_to_1fill_passing_the_root_node) {
  using T = TypeParam;
  skip_test<T>(8, 0b01000001, 6, 6, 1);
  skip_test<T>(8, 0b10000001, 7, 7, 1);
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_run_iterator_skip_test, skip_before_1fill) {
  using T = TypeParam;
  skip_test<T>(8, 0b01000001, 5, 6, 1);
  skip_test<T>(8, 0b01000001, 4, 6, 1);
  skip_test<T>(8, 0b01000001, 3, 6, 1);
  skip_test<T>(8, 0b01000001, 2, 6, 1);
  skip_test<T>(8, 0b01000001, 1, 6, 1);
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_run_iterator_skip_test, skip_beyond_last_1fill) {
  using T = TypeParam;
  skip_test<T>(8, 0b01000001, 7, 8, 0);
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_run_iterator_skip_test, skip_into_a_1fill) {
  using T = TypeParam;
  skip_test<T>(8, 0b11111101, 6, 6, 2);
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_run_iterator_skip_test, skip_multiple_times) {
  using T = TypeParam;
  dtl::bitmap b(8, 0b11001101);
  T tm(b);
  auto it = tm.it();
  it.skip_to(2);
  ASSERT_EQ(it.pos(), 2);
//  ASSERT_EQ(it.length(), 2);
  it.skip_to(6);
  ASSERT_EQ(it.pos(), 6);
//  ASSERT_EQ(it.length(), 2);
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_run_iterator_skip_test, skip_into_a_1fill_length_one) {
  using T = TypeParam;
  skip_test<T>(16, 0b1100000110001111, 6, 7, 1);
}
//===----------------------------------------------------------------------===//
template<typename T>
void
skip_next_test(u64 n, u64 bitmap, u64 skip_to_pos,
               u64 expected_pos, u64 expected_len) {
  assert(n <= 64);
  std::stringstream info;
  info << "skip/next test: n=" << n
       << ", bitmap=" << boost::dynamic_bitset<$u32>(n, bitmap)
       << " (=" << bitmap << ")"
       << ", skip_to=" << skip_to_pos
       << ", expected_pos=" << expected_pos
       << ", expected_len=" << expected_len
       << std::endl;
  std::cout << info.str();
  dtl::bitmap b(n, bitmap);
  T t(b);
  std::cout << t << std::endl;
  auto it = t.it();
  it.skip_to(skip_to_pos);
  it.next();
  if (expected_len == 0) {
    ASSERT_TRUE(it.end()) << info.str();
  }
  else {
    ASSERT_TRUE(!it.end()) << info.str();
  }
  ASSERT_EQ(it.pos(), expected_pos) << info.str();
  if (expected_pos != n) {
//    ASSERT_EQ(it.length(), expected_len);
  }
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_run_iterator_skip_test, skip_next) {
  using T = TypeParam;
  skip_next_test<T>(8, 0b00010101, 2, 4, 1);
  skip_next_test<T>(8, 0b00100101, 2, 5, 1);
  skip_next_test<T>(8, 0b01000101, 2, 6, 1);
  skip_next_test<T>(8, 0b10000101, 2, 7, 1);
  skip_next_test<T>(8, 0b11110101, 2, 4, 4);
  skip_next_test<T>(8, 0b11100101, 2, 5, 3);
  skip_next_test<T>(8, 0b11000101, 2, 6, 2);
  skip_next_test<T>(8, 0b10000101, 2, 7, 1);
}
//===----------------------------------------------------------------------===//
