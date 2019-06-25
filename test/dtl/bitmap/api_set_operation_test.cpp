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
#include <navin/experiments/util/gen.hpp>

#include "util/bitwise_operations.hpp"

//===----------------------------------------------------------------------===//
// API tests for set operations.
//===----------------------------------------------------------------------===//

constexpr std::size_t LEN = 8;

constexpr std::size_t RANDOM_REPEAT = 50;
constexpr std::size_t RANDOM_LENGTH = 1ull << 10;


// Types under test.
using teb_o0 = dtl::teb<0>;
using teb_o1 = dtl::teb<1>;
using teb_o2 = dtl::teb<2>;
using teb_o3 = dtl::teb<3>;
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
class api_set_operation_test : public ::testing::Test {};

// Specify the types for which we want to run the API tests.
using types_under_test = ::testing::Types<
    teb_o0,
//    teb_o1,
    teb_o2,
    teb_o3,

//    teb_scan_o0,
//    teb_scan_o1,
//    teb_scan_o2,

//    plain_bitmap_32,
//    position_list_32,
//    partitioned_position_list_32_u8,
//    partitioned_position_list_32_u16,
//    range_list_32,
//    partitioned_range_list_32_u8,
//    partitioned_range_list_32_u16,
    roaring_bitmap
//    ,
//    wah
>;
TYPED_TEST_CASE(api_set_operation_test, types_under_test);


//===----------------------------------------------------------------------===//
// Intersection / bitwise and.
TYPED_TEST(api_set_operation_test, bitwise_and) {
  using T = TypeParam;

  for (std::size_t a = 1; a < (1u << LEN); a++) {
    dtl::bitmap bm_a(LEN, a);
    T tm_a(bm_a);

    for (std::size_t b = 0; b < (1u << LEN); b++) {
      dtl::bitmap bm_b(LEN, b);
      dtl::bitmap bm_expected = bm_a & bm_b;
      std::cout << "a=" << bm_a << " (" << a <<  ")"
                << ", b=" << bm_b << " (" << b <<  ")"
                << ", r=" << bm_expected << std::endl;

      T tm_b(bm_b);
      const auto bm_actual = bitwise_and(tm_a, tm_b);

      ASSERT_EQ(bm_actual, bm_expected)
          << "Test (a & b) failed for a=" << bm_a << " (" << a << ")"
          << " and b=" << bm_b << "(" << b << ")" << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_set_operation_test, bitwise_and_iter) {
  using T = TypeParam;

  for (std::size_t a = 1; a < (1u << LEN); a++) {
    dtl::bitmap bm_a(LEN, a);
    T tm_a(bm_a);

    for (std::size_t b = 0; b < (1u << LEN); b++) {
      dtl::bitmap bm_b(LEN, b);
      dtl::bitmap bm_expected = bm_a & bm_b;
      std::cout << "a=" << bm_a << " (" << a <<  ")"
                << ", b=" << bm_b << " (" << b <<  ")"
                << ", r=" << bm_expected << std::endl;

      T tm_b(bm_b);
      const auto bm_actual = bitwise_and_iter(tm_a, tm_b);

      ASSERT_EQ(bm_actual, bm_expected)
          << "Test (a & b) failed for a=" << bm_a << " (" << a << ")"
          << " and b=" << bm_b << "(" << b << ")" << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_set_operation_test, bitwise_and_random) {
  using T = TypeParam;
  static u64 n = RANDOM_LENGTH;
  for (std::size_t r = 0; r < RANDOM_REPEAT; ++r) {
    auto bm_a = gen_random_bitmap_uniform(n, std::min(0.01 * r, 1.0));
    auto bm_b = gen_random_bitmap_uniform(n, std::min(0.025 * r, 1.0));

    auto bm_expected = bm_a & bm_b;
    T tm_expected(bm_expected);

    T tm_a(bm_a);
    T tm_b(bm_b);
    const auto bm_actual = bitwise_and(tm_a, tm_b);

    ASSERT_EQ(bm_actual, bm_expected);
  }
}
//===----------------------------------------------------------------------===//
// Intersection / bitwise or.
TYPED_TEST(api_set_operation_test, bitwise_or) {
  using T = TypeParam;

  for (std::size_t a = 1; a < (1u << LEN); a++) {
    dtl::bitmap bm_a(LEN, a);
    T tm_a(bm_a);

    for (std::size_t b = 0; b < (1u << LEN); b++) {
      dtl::bitmap bm_b(LEN, b);
      dtl::bitmap bm_expected = bm_a | bm_b;
      std::cout << "a=" << bm_a << " (" << a <<  ")"
                << ", b=" << bm_b << " (" << b <<  ")"
                << ", r=" << bm_expected << std::endl;

      T tm_b(bm_b);
      const auto bm_actual = bitwise_or(tm_a, tm_b);

      ASSERT_EQ(bm_actual, bm_expected)
          << "Test (a | b) failed for a=" << bm_a << " (" << a << ")"
          << " and b=" << bm_b << "(" << b << ")" << std::endl;

    }
  }
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_set_operation_test, bitwise_or_iter) {
  using T = TypeParam;

  for (std::size_t a = 1; a < (1u << LEN); a++) {
    dtl::bitmap bm_a(LEN, a);
    T tm_a(bm_a);

    for (std::size_t b = 0; b < (1u << LEN); b++) {
      dtl::bitmap bm_b(LEN, b);
      dtl::bitmap bm_expected = bm_a | bm_b;
      std::cout << "a=" << bm_a << " (" << a <<  ")"
                << ", b=" << bm_b << " (" << b <<  ")"
                << ", r=" << bm_expected << std::endl;

      T tm_b(bm_b);
      const auto bm_actual = bitwise_or_iter(tm_a, tm_b);

      ASSERT_EQ(bm_actual, bm_expected)
          << "Test (a | b) failed for a=" << bm_a << " (" << a << ")"
          << " and b=" << bm_b << "(" << b << ")" << std::endl;

    }
  }
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_set_operation_test, bitwise_or_random) {
  using T = TypeParam;
  static u64 n = RANDOM_LENGTH;
  for (std::size_t r = 0; r < RANDOM_REPEAT; ++r) {
    auto bm_a = gen_random_bitmap_uniform(n, std::min(0.01 * r, 1.0));
    auto bm_b = gen_random_bitmap_uniform(n, std::min(0.025 * r, 1.0));

    auto bm_expected = bm_a | bm_b;

    T tm_a(bm_a);
    T tm_b(bm_b);
    const auto bm_actual = bitwise_or(tm_a, tm_b);

    ASSERT_EQ(bm_actual, bm_expected);
  }
}
//===----------------------------------------------------------------------===//
// Difference / bitwise xor.
TYPED_TEST(api_set_operation_test, bitwise_xor) {
  using T = TypeParam;

  for (std::size_t a = 1; a < (1u << LEN); a++) {
    dtl::bitmap bm_a(LEN, a);
    T tm_a(bm_a);

    for (std::size_t b = 0; b < (1u << LEN); b++) {
      dtl::bitmap bm_b(LEN, b);
      dtl::bitmap bm_expected = bm_a ^ bm_b;
      std::cout << "a=" << bm_a << " (" << a <<  ")"
                << ", b=" << bm_b << " (" << b <<  ")"
                << ", r=" << bm_expected << std::endl;

      T tm_b(bm_b);
      const auto bm_actual = bitwise_xor(tm_a, tm_b);

      ASSERT_EQ(bm_actual, bm_expected)
          << "Test (a ^ b) failed for a=" << bm_a << " (" << a << ")"
          << " and b=" << bm_b << "(" << b << ")" << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_set_operation_test, bitwise_xor_iter) {
  using T = TypeParam;

  for (std::size_t a = 1; a < (1u << LEN); a++) {
    dtl::bitmap bm_a(LEN, a);
    T tm_a(bm_a);

    for (std::size_t b = 0; b < (1u << LEN); b++) {
      dtl::bitmap bm_b(LEN, b);
      dtl::bitmap bm_expected = bm_a ^ bm_b;
      std::cout << "a=" << bm_a << " (" << a <<  ")"
                << ", b=" << bm_b << " (" << b <<  ")"
                << ", r=" << bm_expected << std::endl;

      T tm_b(bm_b);
      const auto bm_actual = bitwise_xor_iter(tm_a, tm_b);

      ASSERT_EQ(bm_actual, bm_expected)
          << "Test (a ^ b) failed for a=" << bm_a << " (" << a << ")"
          << " and b=" << bm_b << "(" << b << ")" << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_set_operation_test, bitwise_xor_random) {
  using T = TypeParam;
  static u64 n = RANDOM_LENGTH;

  for (std::size_t r = 0; r < RANDOM_REPEAT; ++r) {
    auto bm_a = gen_random_bitmap_uniform(n, std::min(0.01 * r, 1.0));
    auto bm_b = gen_random_bitmap_uniform(n, std::min(0.025 * r, 1.0));

    auto bm_expected = bm_a ^ bm_b;
    for (std::size_t i = 0; i < bm_a.size(); ++i) {
      assert(bm_expected[i] == bm_a[i] ^ bm_b[i]);
    }

    T tm_a(bm_a);
    T tm_b(bm_b);

    auto bm_actual = bitwise_xor(tm_a, tm_b);
    ASSERT_EQ(bm_actual, bm_expected)
        << "Test (a ^ b) failed for"
        << "\na=" << bm_a
        << "\nb=" << bm_b
        << "\nd=" << (bm_actual ^ bm_expected)
        << "\nf=" << (bm_actual ^ bm_expected).find_first()
        << std::endl;
  }
}
//===----------------------------------------------------------------------===//
// Special case of bitwise xor which occurs with range encoding (RE). In
// range-encoded indexes, the following holds: Each bit that is set in the i-th
// bitmap, then these bits are also set in the (i+1)-th bitmap.
TYPED_TEST(api_set_operation_test, bitwise_xor_re) {
  using T = TypeParam;

  for (std::size_t a = 1; a < (1u << LEN); a++) {
    dtl::bitmap bm_a(LEN, a);
    T tm_a(bm_a);

    for (std::size_t b = 0; b < (1u << LEN); b++) {
      dtl::bitmap bm_b(LEN, b);
      // consider only RE cases
      if ((bm_a & bm_b) != bm_a) continue;
      dtl::bitmap bm_expected = bm_a ^ bm_b;
      std::cout << "a=" << bm_a << " (" << a <<  ")"
                << ", b=" << bm_b << " (" << b <<  ")"
                << ", r=" << bm_expected << std::endl;

      T tm_b(bm_b);
      const auto bm_actual = bitwise_xor_re(tm_a, tm_b);

      ASSERT_EQ(bm_actual, bm_expected)
          << "Test (a ^ b) failed for a=" << bm_a << " (" << a << ")"
          << " and b=" << bm_b << "(" << b << ")" << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_set_operation_test, bitwise_xor_re_iter) {
  using T = TypeParam;

  for (std::size_t a = 1; a < (1u << LEN); a++) {
    dtl::bitmap bm_a(LEN, a);
    T tm_a(bm_a);

    for (std::size_t b = 0; b < (1u << LEN); b++) {
      dtl::bitmap bm_b(LEN, b);
      // consider only RE cases
      if ((bm_a & bm_b) != bm_a) continue;
      dtl::bitmap bm_expected = bm_a ^ bm_b;
      std::cout << "a=" << bm_a << " (" << a <<  ")"
                << ", b=" << bm_b << " (" << b <<  ")"
                << ", r=" << bm_expected << std::endl;

      T tm_b(bm_b);
      const auto bm_actual = bitwise_xor_re_iter(tm_a, tm_b);

      ASSERT_EQ(bm_actual, bm_expected)
          << "Test (a ^ b) failed for a=" << bm_a << " (" << a << ")"
          << " and b=" << bm_b << "(" << b << ")" << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_set_operation_test, bitwise_xor_re_random) {
  using T = TypeParam;
  static u64 n = RANDOM_LENGTH;

  for (std::size_t r = 0; r < RANDOM_REPEAT; ++r) {
    auto bm_a = gen_random_bitmap_uniform(n, std::min(0.01 * r, 1.0));
    auto bm_b = gen_random_bitmap_uniform(n, std::min(0.025 * r, 1.0));
    bm_b |= bm_a; // RE case

    auto bm_expected = bm_a ^ bm_b;
    for (std::size_t i = 0; i < bm_a.size(); ++i) {
      assert(bm_expected[i] == bm_a[i] ^ bm_b[i]);
    }

    T tm_a(bm_a);
    T tm_b(bm_b);

    auto bm_actual = bitwise_xor_re(tm_a, tm_b);
    ASSERT_EQ(bm_actual, bm_expected)
        << "Test (a ^ b) failed for"
        << "\na=" << bm_a
        << "\nb=" << bm_b
        << "\nd=" << (bm_actual ^ bm_expected)
        << "\nf=" << (bm_actual ^ bm_expected).find_first()
        << std::endl;
  }
}
//===----------------------------------------------------------------------===//
