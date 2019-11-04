#include "gtest/gtest.h"

#include <dtl/bitmap.hpp>
#include <dtl/bitmap/util/convert.hpp>
#include <dtl/bitmap/util/plain_bitmap.hpp>
#include <dtl/bitmap/util/plain_bitmap_iter.hpp>
//===----------------------------------------------------------------------===//
// Tests for the plain bitmap iterator.
//===----------------------------------------------------------------------===//
// The iterator is supposed to be working with boost::dynamic_bitset and
// with dtl::plain_bitmap.
using types_under_test = ::testing::Types<
    boost::dynamic_bitset<$u32>,
    boost::dynamic_bitset<$u64>,
    dtl::plain_bitmap<$u32>,
    dtl::plain_bitmap<$u64>>;
//===----------------------------------------------------------------------===//
// Fixture for the parameterized test case.
template<typename T>
class plain_bitmap_iter_test : public ::testing::Test {};
TYPED_TEST_CASE(plain_bitmap_iter_test, types_under_test);
//===----------------------------------------------------------------------===//
/// Test the 1-run iterator by making a copy of the input bitmap.
TYPED_TEST(plain_bitmap_iter_test, make_copy) {
  using T = TypeParam;
  constexpr auto LEN = 8;
  for (auto i = 0; i < (1u << LEN); ++i) {
    // Start with an empty bitmap.
    T b(LEN);

    // Populate the bitmap.
    std::array<$u32, LEN> pos;
    auto pos_cnt = dtl::bitmap_fun<$u32>::to_positions(i, &pos[0]);
    ASSERT_EQ(dtl::bits::pop_count(i), pos_cnt);

    for (std::size_t i = 0; i < pos_cnt; ++i) {
      b.set(pos[i], true);
    }

    // Decompress and validate it.
    dtl::plain_bitmap_iter<T> it(b);
    dtl::bitmap dec = dtl::to_bitmap_from_iterator(it, b.size());

    for (std::size_t i = 0; i < LEN; ++i) {
      ASSERT_EQ(b.test(i), dec.test(i));
    }
  }
}
//===----------------------------------------------------------------------===//
