#include "api_types.hpp"
#include "gtest/gtest.h"

#include <dtl/bitmap.hpp>
#include <dtl/bitmap/util/convert.hpp>
#include <dtl/bitmap/util/random.hpp>
#include <dtl/dtl.hpp>
//===----------------------------------------------------------------------===//
// Typed API tests for loss-less compressed bitmaps.
// Iterate runs.
//===----------------------------------------------------------------------===//
constexpr std::size_t LEN = 8;
//===----------------------------------------------------------------------===//
// Fixture for the parameterized test case.
template<typename T>
class api_run_iterator_test : public ::testing::Test {};
TYPED_TEST_CASE(api_run_iterator_test, types_under_test);
// TODO test skip iterator and scan iterator explicitly
//===----------------------------------------------------------------------===//
TYPED_TEST(api_run_iterator_test, all_bits_set) {
  using T = TypeParam;
  dtl::bitmap b(64);
  b.set();
  T t(b);
  auto actual = dtl::to_bitmap_using_iterator(t);
  ASSERT_EQ(b, actual);
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_run_iterator_test, no_bits_set) {
  using T = TypeParam;
  dtl::bitmap b(64);
  T t(b);
  auto actual = dtl::to_bitmap_using_iterator(t);
  ASSERT_EQ(b, actual);
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_run_iterator_test, single_1fill) {
  using T = TypeParam;
  dtl::bitmap b(64, 0b11110000);
  T t(b);
  auto actual = dtl::to_bitmap_using_iterator(t);
  ASSERT_EQ(b, actual);
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_run_iterator_test, multiple_1fills_at_different_levels) {
  using T = TypeParam;
  dtl::bitmap b(64, 0b0011000011110000);
  T t(b);
  auto actual = dtl::to_bitmap_using_iterator(t);
  ASSERT_EQ(b, actual);
}
//===----------------------------------------------------------------------===//
