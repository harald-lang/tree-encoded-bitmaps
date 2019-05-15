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
#include <dtl/bitmap/util/random.hpp>

//===----------------------------------------------------------------------===//
// API tests for range indices and bitmap indices.
// Encoding/Decoding using the iterator interface.
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
class api_encode_decode_test : public ::testing::Test {};

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
TYPED_TEST_CASE(api_encode_decode_test, types_under_test);

//===----------------------------------------------------------------------===//
// Encode, decode, and compare the results.
TYPED_TEST(api_encode_decode_test, info) {
  using T = TypeParam;

  dtl::bitmap bs(LEN, 1u << 8);
  T t(bs);
  std::cout << t.info() << std::endl;
}
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
// Encode, decode, and compare the results.
TYPED_TEST(api_encode_decode_test, encode_decode_0_to_255) {
  using T = TypeParam;

  for (auto i = 0; i < (1u << LEN); ++i) {
    std::cout << std::bitset<LEN>(i) << std::endl;
    dtl::bitmap bs(LEN, i);
    T t(bs);
    dtl::bitmap dec = t.to_bitset();
    ASSERT_EQ(bs, dec) << "Decoding failed for i=" << i
      << ". - '" << bs << "' -> '" << t
      << "' -> '" << dec << "'"
      << std::endl;
  }
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
// Encode, decode, and compare the results (varying bitmap sizes).
TYPED_TEST(api_encode_decode_test, encode_decode_varying_bitmap_sizes) {
  using T = TypeParam;

  for (auto len = 128; len <= 1ull << 4; ++len) {
    dtl::bitmap bs(len);
    // all bits zero
    bs.clear();
    {
      T t(bs);
      dtl::bitmap dec = t.to_bitset();
      ASSERT_EQ(bs, dec)
        << "Decoding failed: "
        << "'" << bs << "' -> '" << t
        << "' -> '" << dec << "'"
        << std::endl;
    }
    // all bits one
    bs.flip();
    {
      T t(bs);
      dtl::bitmap dec = t.to_bitset();
      ASSERT_EQ(bs, dec)
        << "Decoding failed: "
        << "'" << bs << "' -> '" << t
        << "' -> '" << dec << "'"
        << std::endl;
    }
    // random bitmap
    {
      for (std::size_t rep = 0; rep < 10; ++rep) {
        bs = dtl::gen_random_bitmap_markov(len, 4.0, 0.2);
        T t(bs);
        dtl::bitmap dec = t.to_bitset();
        ASSERT_EQ(bs, dec)
          << "Decoding failed: "
          << "'" << bs << "' -> '" << t
          << "' -> '" << dec << "'"
          << std::endl;
      }
    }

  }
}
//===----------------------------------------------------------------------===//
