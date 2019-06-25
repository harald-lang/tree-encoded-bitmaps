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
// API tests for range indices and bitmap indices.
// Random access / point lookup
//===----------------------------------------------------------------------===//

constexpr std::size_t LEN = 8;

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
class api_random_access_test : public ::testing::Test {};

// Specify the types for which we want to run the API tests.
using types_under_test = ::testing::Types<
//    teb_o0,
//    teb_o1,
//    teb_o2,
    teb_o3,
//    teb_scan_o0,
//    teb_scan_o1,
//    teb_scan_o2,
    plain_bitmap_32,
//    position_list_32,
//    partitioned_position_list_32_u8,
//    partitioned_position_list_32_u16,
//    range_list_32,
//    partitioned_range_list_32_u8,
//    partitioned_range_list_32_u16,
    roaring_bitmap,
    wah
>;
TYPED_TEST_CASE(api_random_access_test, types_under_test);

//===----------------------------------------------------------------------===//
TYPED_TEST(api_random_access_test, empty_bitmap) {
  using T = TypeParam;

  for (auto n_log2 = 10; n_log2 <= 20; ++n_log2) {
    const auto n = 1ull << n_log2;
    dtl::bitmap bs(n);
    T t(bs);
    for (std::size_t i = 0; i < n; ++i) {
      ASSERT_FALSE(t.test(i)) << "Point lookup failed at index i=" << i
          << ". - Expected 'false' but got 'true'."
          << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_random_access_test, sparse_uniform_bitmap) {
  using T = TypeParam;

  for (auto n_log2 = 10; n_log2 <= 20; ++n_log2) {
    const auto n = 1ull << n_log2;
    dtl::bitmap bs = dtl::gen_random_bitmap_uniform(n, 0.0001);
    T t(bs);
    for (std::size_t i = 0; i < n; ++i) {
      u1 a = bs[i];
      u1 b = t.test(i);
      if (a != b) {
        std::cout << t.info() << std::endl;
        u1 b = t.test(i);
        std::cout << b << std::endl;
      }
      ASSERT_EQ(a, b) << "Point lookup failed at index i=" << i
          << ". - Expected '" << (a ? "true" : "false")
          << "' but got '" << (b ? "true" : "false") << "'.\n"
          << "Bitmap info:\n"
          << t.info()
          << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_random_access_test, dense_uniform_bitmap) {
  using T = TypeParam;

  for (auto n_log2 = 10; n_log2 <= 20; ++n_log2) {
    const auto n = 1ull << n_log2;
    dtl::bitmap bs = dtl::gen_random_bitmap_uniform(n, 0.25);
    T t(bs);
    for (std::size_t i = 0; i < n; ++i) {
      u1 a = bs[i];
      u1 b = t.test(i);
      ASSERT_EQ(a, b) << "Point lookup failed at index i=" << i
          << ". - Expected '" << (a ? "true" : "false")
          << "' but got '" << (b ? "true" : "false") << "'."
          << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//
