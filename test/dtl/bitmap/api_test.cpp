#include "gtest/gtest.h"

#include <bitset>

#include <dtl/bitmap/static/roaring_bitmap.hpp>
#include <dtl/bitmap/static/tree_mask_lo.hpp>
#include <dtl/bitmap/static/tree_mask_po.hpp>
#include <dtl/bitmap/static/wah.hpp>

//===----------------------------------------------------------------------===//
// Typed API tests for loss-less compressed bitmaps.
//===----------------------------------------------------------------------===//

// Instantiate the type templates to a fixed length.
constexpr std::size_t LEN = 8;
using roaring_bitmap = dtl::roaring_bitmap<LEN>;
using tree_mask_lo = dtl::tree_mask_lo<LEN>;
using tree_mask_po = dtl::tree_mask_po<LEN>;
using wah32 = dtl::wah32<LEN>;
using wah64 = dtl::wah64<LEN>;

// Fixture for the parameterized test case.
template<typename T>
class api_test : public ::testing::Test {};

// Specify the types for which we want to run the API tests.
using types_under_test = ::testing::Types<
//    roaring_bitmap,
    tree_mask_lo //,
//    tree_mask_po,
//    wah32
//    wah64 // FIXME: Position iterator does not work with 64-bit impl.
>;
TYPED_TEST_CASE(api_test, types_under_test);


//===----------------------------------------------------------------------===//
// Encode, decode, compare the results.
TYPED_TEST(api_test, encode_decode) {
  using T = TypeParam;

  for (auto i = 0; i < (1u << LEN); ++i) {
    std::bitset<LEN> bs(i);
    T t(bs);
    std::bitset<LEN> dec = t.to_bitset();
    ASSERT_EQ(bs, dec) << "Decoding failed for i=" << i << ". - '" << bs << "' -> '" << t << "' -> '" << dec << "'"
                       << std::endl;
  }
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
// Bitwise AND.
TYPED_TEST(api_test, bitwise_and) {
  using T = TypeParam;

  for (std::size_t a = 0; a < (1u << LEN); a++) {
    std::bitset<LEN> bm_a(a);
    T tm_a(bm_a);

    for (std::size_t b = 0; b < (1u << LEN); b++) {
      std::bitset<LEN> bm_b(b);
      std::bitset<LEN> bm_expected = bm_a & bm_b;

      T tm_b(bm_b);
      T tm_c = tm_a & tm_b;

      std::bitset<LEN> bm_actual = tm_c.to_bitset();

      ASSERT_EQ(bm_actual, bm_expected) << "Test (a & b) failed for a=" << a << " and b=" << b << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
// Bitwise AND (with range encoding).
// Note: For range encoding the following must hold: (a[i] == true) => (b[i] == true)
TYPED_TEST(api_test, bitwise_and_range_encoding) {
  using T = TypeParam;

  for (std::size_t a = 0; a < (1u << LEN); a++) {
    std::bitset<LEN> bm_a(a);
    T tm_a(bm_a);

    for (std::size_t b = 0; b < (1u << LEN); b++) {
      // make sure, that all bit that are set in a are also set in b (as guaranteed in range encoding)
      if ((a & b) != a) continue;

      std::bitset<LEN> bm_b(b);
      std::bitset<LEN> bm_expected = bm_a & bm_b;

      T tm_b(bm_b);
      T tm_c = tm_a.and_re(tm_b);

      std::bitset<LEN> bm_actual = tm_c.to_bitset();

      ASSERT_EQ(bm_actual, bm_expected) << "Test (a and_re b) failed for a=" << a << " and b=" << b << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
// Bitwise XOR.
TYPED_TEST(api_test, bitwise_xor) {
  using T = TypeParam;

  for (std::size_t a = 0; a < (1u << LEN); a++) {
    std::bitset<LEN> bm_a(a);
    T tm_a(bm_a);

    for (std::size_t b = 0; b < (1u << LEN); b++) {
      std::bitset<LEN> bm_b(b);
      std::bitset<LEN> bm_expected = bm_a ^ bm_b;

      T tm_b(bm_b);
      T tm_c = tm_a ^ tm_b;

      std::bitset<LEN> bm_actual = tm_c.to_bitset();

      ASSERT_EQ(bm_actual, bm_expected) << "Test (a ^ b) failed for a=" << a << " and b=" << b << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
// Bitwise XOR (with range encoding).
// Note: For range encoding the following must hold: (a[i] == true) => (b[i] == true)
TYPED_TEST(api_test, bitwise_xor_range_encoding) {
  using T = TypeParam;

  for (std::size_t a = 0; a < (1u << LEN); a++) {
    std::bitset<LEN> bm_a(a);
    T tm_a(bm_a);

    for (std::size_t b = 0; b < (1u << LEN); b++) {
      // make sure, that all bit that are set in a are also set in b (as guaranteed in range encoding)
      if ((a & b) != a) continue;

      std::bitset<LEN> bm_b(b);
      std::bitset<LEN> bm_expected = bm_a ^ bm_b;

      T tm_b(bm_b);
      T tm_c = tm_a.xor_re(tm_b);

      std::bitset<LEN> bm_actual = tm_c.to_bitset();

      ASSERT_EQ(bm_actual, bm_expected) << "Test (a xor_re b) failed for a=" << a << " and b=" << b << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
// Fused XOR-AND.
// Note: Fused XOR-AND is used with multi-dimensional predicates,
//       e.g. (10 < a < 100) & (20 < b < 30), in combination with range encoding.
//       XOR is used to evaluate the predicate on a column (i.e., a or b) and
//       the AND operation is used to combine the results of multiple
//       predicates.
//       The idea is to fuse both operations together and be (possibly)
//       faster than evaluating them one-by-one.
TYPED_TEST(api_test, bitwise_fused_xor_and_range_encoding) {
  using T = TypeParam;

  for (std::size_t a = 0; a < (1u << LEN); a++) {
    std::bitset<LEN> bm_a(a);
    T tm_a(bm_a);

    for (std::size_t b = 0; b < (1u << LEN); b++) {
      // make sure, that all bit that are set in a are also set in b (as guaranteed in range encoding)
      if ((a & b) != a) continue;

      std::bitset<LEN> bm_b(b);
      T tm_b(bm_b);

      std::bitset<LEN> bm_a_xor_b = bm_a ^ bm_b;
      T tm_a_xor_b = tm_a.xor_re(tm_b);

      for (std::size_t c = 0; c < (1u << LEN); c++) {
        std::bitset<LEN> bm_c(c);
        T tm_c(bm_c);

        T tm_res = tm_a_xor_b & tm_c;

        std::bitset<LEN> bm_expected = bm_a_xor_b & bm_c;
        std::bitset<LEN> bm_actual = tm_res.to_bitset();

        ASSERT_EQ(bm_actual, bm_expected) << "Test ((a ^_re b) & c) failed for a=" << a << ", b=" << b << ", c=" << c << std::endl;
      }

    }
  }
}
//===----------------------------------------------------------------------===//