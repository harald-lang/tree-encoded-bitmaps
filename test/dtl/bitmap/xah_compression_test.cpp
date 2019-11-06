#include "gtest/gtest.h"

#include <dtl/bitmap.hpp>
#include <dtl/bitmap/dynamic_wah.hpp>
#include <dtl/bitmap/util/convert.hpp>
#include <dtl/bitmap/util/plain_bitmap_iter.hpp>
#include <dtl/bitmap/util/random.hpp>
#include <dtl/bitmap/xah.hpp>
//===----------------------------------------------------------------------===//
// AB tests for XAH regarding compression ratios.
//===----------------------------------------------------------------------===//
template<
    typename A,
    typename B>
struct ab {
  using a_type = A;
  using b_type = B;
};
//===----------------------------------------------------------------------===//
using types_under_test = ::testing::Types<
    ab<dtl::xah32, dtl::dynamic_wah32>,
    ab<dtl::xah64, dtl::dynamic_wah64>
    >;
//===----------------------------------------------------------------------===//
// Fixture for the parameterized test case.
template<typename T>
class xah_compression_test : public ::testing::Test {};
TYPED_TEST_CASE(xah_compression_test, types_under_test);
//===----------------------------------------------------------------------===//
/// Test the compression ratio of XAH (must not be larger than WAH).
TYPED_TEST(xah_compression_test, ensure_xah_size_is_le_to_wah_size) {
  using A = typename TypeParam::a_type;
  using B = typename TypeParam::b_type;
  for (auto n = 1ull << 10; n <= 1ull << 16; n <<= 1) {
    for (auto d = 0.01; d <= 0.75; d += 0.01) {
      auto bm = dtl::gen_random_bitmap_markov(n, 4.0, d);
      // Compress
      A a(bm);
      i64 a_size = a.size_in_byte();
      B b(bm);
      i64 b_size = b.size_in_byte();
      std::cout
          << std::setw(10) << n << ", "
          << std::setw(5) << d << ", "
          << std::setw(10) << a_size << ", "
          << std::setw(10) << b_size << ", "
          << std::setw(10) << std::abs(a_size - b_size)
          << std::setw(10) << std::setprecision(3)
          << ((std::abs(a_size - b_size)/ (n/8.0)) * 100) << "%"
          << std::endl;
      ASSERT_TRUE(a_size <= b_size);
    }
  }
}
//===----------------------------------------------------------------------===//
