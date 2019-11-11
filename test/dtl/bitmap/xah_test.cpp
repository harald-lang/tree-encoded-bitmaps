#include "gtest/gtest.h"

#include <dtl/bitmap.hpp>
#include <dtl/bitmap/traits.hpp>
#include <dtl/bitmap/util/convert.hpp>
#include <dtl/bitmap/util/plain_bitmap_iter.hpp>
#include <dtl/bitmap/util/random.hpp>
#include <dtl/bitmap/uah.hpp>
#include <dtl/bitmap/uah_skip.hpp>
#include <dtl/bitmap/xah.hpp>
#include <dtl/bitmap/xah_skip.hpp>
//===----------------------------------------------------------------------===//
// Tests for XAH.
//===----------------------------------------------------------------------===//
using types_under_test = ::testing::Types<
    dtl::uah8,
    dtl::uah16,
    dtl::uah32,
    dtl::uah64,
    dtl::uah_skip<u8, 1>,
    dtl::uah_skip<u16, 1>,
    dtl::uah_skip<u32, 1>,
    dtl::uah_skip<u64, 1>,
    dtl::xah8,
    dtl::xah16,
    dtl::xah32,
    dtl::xah64,
    dtl::xah_skip<u8, 1>,
    dtl::xah_skip<u16, 1>,
    dtl::xah_skip<u32, 1>,
    dtl::xah_skip<u64, 1>
    >;
//===----------------------------------------------------------------------===//
// Fixture for the parameterized test case.
template<typename T>
class xah_test : public ::testing::Test {};
TYPED_TEST_CASE(xah_test, types_under_test);
//===----------------------------------------------------------------------===//
/// Test encoding followed by decoding.
TYPED_TEST(xah_test, encode_decode_test) {
  using T = TypeParam;
  const auto LEN = 8;
  for (auto i = 0; i < (1u << LEN); ++i) {
    dtl::plain_bitmap<$u64> b(LEN);

    // Populate the bitmap.
    std::array<$u32, LEN> pos;
    auto pos_cnt = dtl::bitmap_fun<$u32>::to_positions(i, &pos[0]);
    ASSERT_EQ(dtl::bits::pop_count(i), pos_cnt);

    for (std::size_t i = 0; i < pos_cnt; ++i) {
      b.set(pos[i], true);
    }

    dtl::bitmap t(LEN, i);
    // Compress
    T enc(t);
    // Decompress and validate it.
    auto dec = enc.to_plain_bitmap();

    ASSERT_EQ(b, dec)
        << "Decoding failed for i=" << i
        << ". - '" << b
        << "' -> '" << dec << "'"
        << std::endl;
  }
}
//===----------------------------------------------------------------------===//
/// Test encoding followed by decoding using random bitmaps with lengths that
/// are NOT a power of two.
TYPED_TEST(xah_test, encode_decode_non_pow2_length_test) {
  using T = TypeParam;
  for (auto len = 10; len < 1000; len *= 1.2) {
    for (auto d = 0.0; d <= 1.0; d += 0.01) {
      auto b = dtl::gen_random_bitmap_uniform(len, d);
      std::cout << b << std::endl;
      // Compress
      T enc(b);

      // Decompress and validate it.
      auto dec = dtl::to_bitmap_using_iterator(enc);

      ASSERT_EQ(b, dec)
          << "Decoding failed for n=" << len << " and d=" << d
          << ". - '" << b
          << "' -> '" << dec << "'"
          << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//
/// Test compression with long bitmaps that are either entirely 0 or 1.
TYPED_TEST(xah_test, encode_decode_long_bitmaps_consisting_of_a_single_run_test) {
  using T = TypeParam;
  auto len = 1ull << 12;
  dtl::bitmap b(len);
  {
    // Compress
    T enc(b);
    // Decompress and validate it.
    auto dec = dtl::to_bitmap_using_iterator(enc);
    ASSERT_EQ(b, dec) << "\n" << enc;
  }

  b.flip();
  {
    // Compress
    T enc(b);
    // Decompress and validate it.
    auto dec = dtl::to_bitmap_using_iterator(enc);
    ASSERT_EQ(b, dec) << "\n" << enc;
  }
}
//===----------------------------------------------------------------------===//
/// Test random access (bit test).
TYPED_TEST(xah_test, random_access_test) {
  using T = TypeParam;
  const auto LEN = 20;
  for (auto i = 0; i < (1u << LEN); ++i) {
    dtl::bitmap b(LEN, i);
    T enc(b);
    for (std::size_t k = 0; k < LEN; ++k) {
      ASSERT_EQ(b.test(k), enc.test(k))
          << "Bit test failed for i=" << i
          << " at pos k=" << k << " - '" << b
          << "' -> \n" << enc
          << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//
/// Test 1-run iterator.
TYPED_TEST(xah_test, run_iterator_test) {
  using T = TypeParam;
  const auto LEN = 20;
  for (auto i = 0; i < (1u << LEN); ++i) {
    dtl::bitmap b(LEN, i);
    T enc(b);
    auto dec = dtl::to_bitmap_using_iterator(enc);
    ASSERT_EQ(b, dec)
        << "Decoding failed for i=" << i
        << ". - '" << b
        << "' -> '" << dec << "'\n"
        << enc
        << std::endl;
  }
}
//===----------------------------------------------------------------------===//
