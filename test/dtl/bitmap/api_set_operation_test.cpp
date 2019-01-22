#include "gtest/gtest.h"

#include <bitset>

#include <dtl/dtl.hpp>
#include <dtl/bitmap.hpp>
#include <dtl/bitmap/position_list.hpp>
#include <dtl/bitmap/partitioned_position_list.hpp>
#include <dtl/bitmap/util/random.hpp>

//===----------------------------------------------------------------------===//
// API tests for range indices and bitmap indices.
// Set operations.
//===----------------------------------------------------------------------===//

constexpr std::size_t LEN = 8;
using teb = dtl::teb<$u32>;
using plain_bitmap_32 = dtl::dynamic_bitmap<$u32>;
using plain_bitmap_64 = dtl::dynamic_bitmap<$u64>;
using position_list_32 = dtl::position_list<$u32>;
using partitioned_position_list_32_u8 = dtl::partitioned_position_list<$u32>;
using roaring_bitmap = dtl::dynamic_roaring_bitmap;
using wah = dtl::dynamic_wah32;

// Fixture for the parameterized test case.
template<typename T>
class api_set_operation_test : public ::testing::Test {};

// Specify the types for which we want to run the API tests.
using types_under_test = ::testing::Types<
//    plain_bitmap_32,
//    position_list_32,
    partitioned_position_list_32_u8,
//    roaring_bitmap,
//    teb,
    wah
>;
TYPED_TEST_CASE(api_set_operation_test, types_under_test);


//===----------------------------------------------------------------------===//
// Bitwise AND.
TYPED_TEST(api_set_operation_test, bitwise_and) {
  using T = TypeParam;

  for (std::size_t a = 1; a < (1u << LEN); a++) {
    dtl::bitmap bm_a(LEN, a);
    T tm_a(bm_a);

    for (std::size_t b = 0; b < (1u << LEN); b++) {
      dtl::bitmap bm_b(LEN, b);
      dtl::bitmap bm_expected = bm_a & bm_b;

      T tm_b(bm_b);
      T tm_c = tm_a & tm_b;

      dtl::bitmap bm_actual = tm_c.to_bitset();

      ASSERT_EQ(bm_actual, bm_expected) << "Test (a & b) failed for a=" << a
                                        << " and b=" << b << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
// Bitwise AND (with range encoding).
// Note: For range encoding the following must hold:
//       (a[i] == true) => (b[i] == true)
TYPED_TEST(api_set_operation_test, bitwise_and_range_encoding) {
  using T = TypeParam;

  for (std::size_t a = 0; a < (1u << LEN); a++) {
    dtl::bitmap bm_a(LEN, a);
    T tm_a(bm_a);

    for (std::size_t b = 0; b < (1u << LEN); b++) {
      // make sure, that all bit that are set in a are also set in b (as
      // guaranteed in range encoding)
      if ((a & b) != a) continue;

      dtl::bitmap bm_b(LEN, b);
      dtl::bitmap bm_expected = bm_a & bm_b;

      T tm_b(bm_b);
      T tm_c = tm_a.and_re(tm_b);

      dtl::bitmap bm_actual = tm_c.to_bitset();

      ASSERT_EQ(bm_actual, bm_expected) << "Test (a and_re b) failed for a="
                                        << a << " and b=" << b << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
// Bitwise XOR.
TYPED_TEST(api_set_operation_test, bitwise_xor) {
  using T = TypeParam;

  for (std::size_t a = 0; a < (1u << LEN); a++) {
    dtl::bitmap bm_a(LEN, a);
    T tm_a(bm_a);

    for (std::size_t b = 0; b < (1u << LEN); b++) {
      dtl::bitmap bm_b(LEN, b);
      dtl::bitmap bm_expected = bm_a ^ bm_b;

      T tm_b(bm_b);
      T tm_c = tm_a ^ tm_b;

      dtl::bitmap bm_actual = tm_c.to_bitset();

      ASSERT_EQ(bm_actual, bm_expected) << "Test (a ^ b) failed for a=" << a
                                        << " and b=" << b << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//


////===----------------------------------------------------------------------===//
//// Bitwise XOR (with range encoding).
//// Note: For range encoding the following must hold:
////       (a[i] == true) => (b[i] == true)
//TYPED_TEST(api_encode_decode_test, bitwise_xor_range_encoding) {
//  using T = TypeParam;
//
//  for (std::size_t a = 0; a < (1u << LEN); a++) {
//    dtl::bitmap bm_a(LEN, a);
//    T tm_a(bm_a);
//
//    for (std::size_t b = 0; b < (1u << LEN); b++) {
//      // make sure, that all bit that are set in a are also set in b (as
//      // guaranteed in range encoding)
//      if ((a & b) != a) continue;
//
//      dtl::bitmap bm_b(LEN, b);
//      dtl::bitmap bm_expected = bm_a ^ bm_b;
//
//      T tm_b(bm_b);
//      T tm_c = tm_a.xor_re(tm_b);
//
//      dtl::bitmap bm_actual = tm_c.to_bitset();
//
//      ASSERT_EQ(bm_actual, bm_expected) << "Test (a xor_re b) failed for a="
//                                        << a << " and b=" << b << std::endl;
//    }
//  }
//}
////===----------------------------------------------------------------------===//


////===----------------------------------------------------------------------===//
//// Fused XOR-AND.
//// Note: Fused XOR-AND is used with multi-dimensional predicates,
////       e.g. (10 < a < 100) & (20 < b < 30), in combination with range encoding.
////       XOR is used to evaluate the predicate on a column (i.e., a or b) and
////       the AND operation is used to combine the results of multiple
////       predicates.
////       The idea is to fuse both operations together and be (possibly)
////       faster than evaluating them one-by-one.
//TYPED_TEST(api_encode_decode_test, bitwise_fused_xor_and_range_encoding) {
//  using T = TypeParam;
//
//  for (std::size_t a = 0; a < (1u << LEN); a++) {
//    dtl::bitmap bm_a(LEN, a);
//    T tm_a(bm_a);
//
//    for (std::size_t b = 0; b < (1u << LEN); b++) {
//      // make sure, that all bit that are set in a are also set in b (as guaranteed in range encoding)
//      if ((a & b) != a) continue;
//
//      dtl::bitmap bm_b(LEN, b);
//      T tm_b(bm_b);
//
//      dtl::bitmap bm_a_xor_b = bm_a ^ bm_b;
//      T tm_a_xor_b = tm_a.xor_re(tm_b);
//
//      for (std::size_t c = 0; c < (1u << LEN); c++) {
//        dtl::bitmap bm_c(LEN, c);
//        T tm_c(bm_c);
//
//        T tm_res = tm_a_xor_b & tm_c;
//
//        dtl::bitmap bm_expected = bm_a_xor_b & bm_c;
//        dtl::bitmap bm_actual = tm_res.to_bitset();
//
//        ASSERT_EQ(bm_actual, bm_expected)
//          << "Test ((a ^_re b) & c) failed for a=" << a
//          << ", b=" << b << ", c=" << c
//          << std::endl;
//      }
//
//    }
//  }
//}
////===----------------------------------------------------------------------===//
