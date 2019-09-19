#include "api_types.hpp"
#include "gtest/gtest.h"

#include <dtl/bitmap.hpp>
#include <dtl/bitmap/util/convert.hpp>
#include <dtl/bitmap/util/random.hpp>
#include <dtl/dtl.hpp>
//===----------------------------------------------------------------------===//
// API tests for range indices and bitmap indices.
// Random access / point lookup
//===----------------------------------------------------------------------===//
constexpr std::size_t LEN = 8;
//===----------------------------------------------------------------------===//
// Fixture for the parameterized test case.
template<typename T>
class api_random_access_test : public ::testing::Test {};
TYPED_TEST_CASE(api_random_access_test, types_under_test);
//===----------------------------------------------------------------------===//
TYPED_TEST(api_random_access_test, empty_bitmap) {
  using T = TypeParam;

  for (auto n_log2 = 10; n_log2 <= 10; ++n_log2) {
    std::cout << "n_log2=" << n_log2 << std::endl;
    const auto n = 1ull << n_log2;
    dtl::bitmap bs(n);
    T t(bs);
    for (std::size_t i = 0; i < n; ++i) {
      ASSERT_FALSE(t.test(i))
          << "Point lookup failed at index i=" << i
          << ". - Expected 'false' but got 'true'."
          << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//
TYPED_TEST(api_random_access_test, sparse_uniform_bitmap) {
  using T = TypeParam;

  for (auto n_log2 = 10; n_log2 <= 10; ++n_log2) {
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
      ASSERT_EQ(a, b)
          << "Point lookup failed at index i=" << i
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

  for (auto n_log2 = 10; n_log2 <= 10; ++n_log2) {
    const auto n = 1ull << n_log2;
    dtl::bitmap bs = dtl::gen_random_bitmap_uniform(n, 0.25);
    T t(bs);
    for (std::size_t i = 0; i < n; ++i) {
      u1 a = bs[i];
      u1 b = t.test(i);
      ASSERT_EQ(a, b)
          << "Point lookup failed at index i=" << i
          << ". - Expected '" << (a ? "true" : "false")
          << "' but got '" << (b ? "true" : "false") << "'."
          << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//
