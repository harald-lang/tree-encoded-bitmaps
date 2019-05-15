#include "gtest/gtest.h"

#include <dtl/dtl.hpp>
#include <dtl/bitmap.hpp>
#include <dtl/bitmap/dynamic_bitmap.hpp>
#include <dtl/bitmap/dynamic_roaring_bitmap.hpp>
#include <dtl/bitmap/dynamic_wah.hpp>
#include <dtl/bitmap/partitioned_position_list.hpp>
#include <dtl/bitmap/partitioned_range_list.hpp>
#include <dtl/bitmap/position_list.hpp>
#include <dtl/bitmap/range_list.hpp>
#include <dtl/bitmap/teb.hpp>
#include <dtl/bitmap/teb_scan.hpp>
#include <dtl/bitmap/util/convert.hpp>
#include <dtl/bitmap/util/random.hpp>

//===----------------------------------------------------------------------===//
// Typed API tests for loss-less compressed bitmaps.
// Iterate ranges.
//===----------------------------------------------------------------------===//

constexpr std::size_t LEN = 8;

// Types under test.
using teb_o0 = dtl::teb<0>;
using teb_o1 = dtl::teb<1>;
using teb_o2 = dtl::teb<2>;
using teb_scan_o0 = dtl::teb_scan<0>;
using teb_scan_o1 = dtl::teb_scan<1>;
using teb_scan_o2 = dtl::teb_scan<2>;
// Competitors
using plain_bitmap_32 = dtl::dynamic_bitmap<$u32>;
using position_list_32 = dtl::position_list<$u32>;
using partitioned_position_list_32_u8 = dtl::partitioned_position_list<$u32, $u8>;
using partitioned_position_list_32_u16 = dtl::partitioned_position_list<$u32, $u16>;
using range_list_32 = dtl::range_list<$u32>;
using partitioned_range_list_32_u8 = dtl::partitioned_range_list<$u32, $u8>;
using partitioned_range_list_32_u16 = dtl::partitioned_range_list<$u32, $u16>;
using roaring_bitmap = dtl::dynamic_roaring_bitmap;
using wah = dtl::dynamic_wah32;

// Fixture for the parameterized test case.
template<typename T>
class api_range_iterator_skip_test : public ::testing::Test {};

// Specify the types for which we want to run the API tests.
using types_under_test = ::testing::Types<
    teb_o0,
    teb_o1,
    teb_o2,
//    teb_scan_o0,
//    teb_scan_o1,
//    teb_scan_o2,
    plain_bitmap_32,
    position_list_32,
    partitioned_position_list_32_u8,
    partitioned_position_list_32_u16,
    range_list_32,
    partitioned_range_list_32_u8,
    partitioned_range_list_32_u16,
    roaring_bitmap,
    wah
>;
TYPED_TEST_CASE(api_range_iterator_skip_test, types_under_test);


//===----------------------------------------------------------------------===//
TYPED_TEST(api_range_iterator_skip_test, all_bits_set) {
  using T = TypeParam;
  dtl::bitmap b(64);
  b.set();
  T t(b);
//  std::cout << t << std::endl;
  auto actual = dtl::to_bitmap_using_iterator(t);
  ASSERT_EQ(b, actual);
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
TYPED_TEST(api_range_iterator_skip_test, no_bits_set) {
  using T = TypeParam;
  dtl::bitmap b(64);
  T t(b);
  auto actual = dtl::to_bitmap_using_iterator(t);
  ASSERT_EQ(b, actual);
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
TYPED_TEST(api_range_iterator_skip_test, single_1fill) {
  using T = TypeParam;
  dtl::bitmap b(64, 0b11110000);
  T t(b);
  auto actual = dtl::to_bitmap_using_iterator(t);
  ASSERT_EQ(b, actual);
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
TYPED_TEST(api_range_iterator_skip_test, multiple_1fills_at_different_levels) {
  using T = TypeParam;
  dtl::bitmap b(64, 0b0011000011110000);
  T t(b);
  auto actual = dtl::to_bitmap_using_iterator(t);
  ASSERT_EQ(b, actual);
}
//===----------------------------------------------------------------------===//


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

TYPED_TEST(api_range_iterator_skip_test, skip_to_1fill) {
  using T = TypeParam;
  skip_test<T>(8, 0b00000101, 2, 2, 1);
  skip_test<T>(8, 0b00001001, 3, 3, 1);
  skip_test<T>(8, 0b00010001, 4, 4, 1);
  skip_test<T>(8, 0b00100001, 5, 5, 1);
  skip_test<T>(8, 0b01000001, 6, 6, 1);
  skip_test<T>(8, 0b10000001, 7, 7, 1);
}

TYPED_TEST(api_range_iterator_skip_test, skip_to_1fill_passing_the_root_node) {
  using T = TypeParam;
  skip_test<T>(8, 0b01000001, 6, 6, 1);
  skip_test<T>(8, 0b10000001, 7, 7, 1);
}

TYPED_TEST(api_range_iterator_skip_test, skip_before_1fill) {
  using T = TypeParam;
  skip_test<T>(8, 0b01000001, 5, 6, 1);
  skip_test<T>(8, 0b01000001, 4, 6, 1);
  skip_test<T>(8, 0b01000001, 3, 6, 1);
  skip_test<T>(8, 0b01000001, 2, 6, 1);
  skip_test<T>(8, 0b01000001, 1, 6, 1);
}

TYPED_TEST(api_range_iterator_skip_test, skip_beyond_last_1fill) {
  using T = TypeParam;
  skip_test<T>(8, 0b01000001, 7, 8, 0);
}

TYPED_TEST(api_range_iterator_skip_test, skip_into_a_1fill) {
  using T = TypeParam;
  skip_test<T>(8, 0b11111101, 6, 6, 2);
}

TYPED_TEST(api_range_iterator_skip_test, skip_multiple_times) {
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

TYPED_TEST(api_range_iterator_skip_test, skip_into_a_1fill_length_one) {
  using T = TypeParam;
  skip_test<T>(16, 0b1100000110001111, 6, 7, 1);
}
//===----------------------------------------------------------------------===//



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

TYPED_TEST(api_range_iterator_skip_test, skip_next) {
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
