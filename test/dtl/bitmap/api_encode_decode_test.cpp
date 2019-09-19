#include "api_types.hpp"
#include "gtest/gtest.h"

#include <dtl/bitmap.hpp>
#include <dtl/bitmap/util/convert.hpp>
#include <dtl/bitmap/util/random.hpp>
#include <dtl/dtl.hpp>
//===----------------------------------------------------------------------===//
// API tests for range indices and bitmap indices.
// Encoding/Decoding.
//===----------------------------------------------------------------------===//
constexpr std::size_t LEN = 8;
//===----------------------------------------------------------------------===//
// Fixture for the parameterized test case.
template<typename T>
class api_encode_decode_test : public ::testing::Test {};
TYPED_TEST_CASE(api_encode_decode_test, types_under_test);
//===----------------------------------------------------------------------===//
// Encode, decode, and compare the results.
TYPED_TEST(api_encode_decode_test, encode_decode_0_to_255) {
  using T = TypeParam;

  for (auto i = 0; i < (1u << LEN); ++i) {
    dtl::bitmap bs(LEN, i);
    std::cout << bs << " (" << i << ")" << std::endl;
    T t(bs);
    dtl::bitmap dec = dtl::to_bitmap_using_iterator(t);
    ASSERT_EQ(bs, dec)
        << "Decoding failed for i=" << i
        << ". - '" << bs << "' -> '" << t
        << "' -> '" << dec << "'"
        << std::endl;
  }
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_encode_decode_test, encode_decode_0_to_65536) {
  using T = TypeParam;
  constexpr std::size_t LEN = 16;

  for (auto i = 0; i < (1u << LEN); ++i) {
    dtl::bitmap bs(LEN, i);
    std::cout << bs << std::endl;
    T t(bs);
    dtl::bitmap dec = dtl::to_bitmap_using_iterator(t);
    ASSERT_EQ(bs, dec)
        << "Decoding failed for i=" << i
        << ". - '" << bs << "' -> '" << t
        << "' -> '" << dec << "'"
        << std::endl;
  }
}
//===----------------------------------------------------------------------===//
// Encode, decode, and compare the results (varying bitmap sizes).
TYPED_TEST(api_encode_decode_test,
    encode_decode_varying_bitmap_sizes) {
  using T = TypeParam;

  for (auto len = 128; len <= 1024; len *= 2) {
    dtl::bitmap bs(len);
    // all bits zero
    bs.reset();
    {
      T t(bs);
      dtl::bitmap dec = dtl::to_bitmap_using_iterator(t);
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
      dtl::bitmap dec = dtl::to_bitmap_using_iterator(t);
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
        dtl::bitmap dec = dtl::to_bitmap_using_iterator(t);
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
