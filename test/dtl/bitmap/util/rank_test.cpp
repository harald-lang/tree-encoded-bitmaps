#include "gtest/gtest.h"

#include <dtl/bitmap/util/plain_bitmap.hpp>
#include <dtl/bitmap/util/random.hpp>
#include <dtl/bitmap/util/rank1.hpp>
#include <dtl/bitmap/util/rank1_logic_linear.hpp>
#include <dtl/bitmap/util/rank1_logic_surf.hpp>
//===----------------------------------------------------------------------===//
// Tests for the rank1 implementations.
//===----------------------------------------------------------------------===//
static constexpr u1 inclusive = true;
static constexpr u1 exclusive = false;
using word_type = $u64;
using rank_types_under_test = ::testing::Types<
    dtl::rank1_logic_surf<word_type, inclusive, 4096>,
    dtl::rank1_logic_surf<word_type, exclusive, 4096>,
    dtl::rank1_logic_surf<word_type, inclusive, 2048>,
    dtl::rank1_logic_surf<word_type, exclusive, 2048>,
    dtl::rank1_logic_surf<word_type, inclusive, 1024>,
    dtl::rank1_logic_surf<word_type, exclusive, 1024>,
    dtl::rank1_logic_surf<word_type, inclusive, 512>,
    dtl::rank1_logic_surf<word_type, exclusive, 512>,
    dtl::rank1_logic_surf<word_type, inclusive, 256>,
    dtl::rank1_logic_surf<word_type, exclusive, 256>,
    dtl::rank1_logic_surf<word_type, inclusive, 128>,
    dtl::rank1_logic_surf<word_type, exclusive, 128>,
    dtl::rank1_logic_surf<word_type, inclusive, 64>,
    dtl::rank1_logic_surf<word_type, exclusive, 64>,
    dtl::rank1_logic_linear<word_type, inclusive>,
    dtl::rank1_logic_linear<word_type, exclusive>,
    >;
//===----------------------------------------------------------------------===//
// Fixture for the parameterized test case.
template<typename T>
class rank_test : public ::testing::Test {};
TYPED_TEST_CASE(rank_test, rank_types_under_test);
//===----------------------------------------------------------------------===//
std::size_t
linear_rank(std::size_t idx, dtl::plain_bitmap<word_type>& b, u1 inclusive) {
  std::size_t cnt = 0;
  for (std::size_t i = 0; i < idx; ++i) {
    cnt += b.test(i);
  }
  if (inclusive) {
    cnt += b.test(idx);
  }
  return cnt;
}
//===----------------------------------------------------------------------===//
TYPED_TEST(rank_test, empty_bitmap) {
  using R = dtl::rank1<TypeParam>;

  dtl::plain_bitmap<word_type> b(4096);
  R rank;
  rank.init(b.data_begin(), b.data_end());

  for (std::size_t i = 0; i < 4096; ++i) {
    const auto r = rank(i, b.data());
    ASSERT_EQ(r, 0)
        << "test failed for i=" << i << "."
        << std::endl;
  }
}
//===----------------------------------------------------------------------===//
TYPED_TEST(rank_test, fully_populated_bitmap) {
  using R = dtl::rank1<TypeParam>;

  dtl::plain_bitmap<word_type> b(4096);
  b.set(0ul, 4096ul);
  R rank;
  rank.init(b.data_begin(), b.data_end());

  for (std::size_t i = 0; i < 4096; ++i) {
    const auto r = rank(i, b.data());
    ASSERT_EQ(r, i + rank.is_inclusive())
        << "test failed for i=" << i << "."
        << std::endl;
  }
}
//===----------------------------------------------------------------------===//
TYPED_TEST(rank_test, random_bitmaps) {
  using R = dtl::rank1<TypeParam>;

  for (std::size_t rep = 0; rep < 5; ++rep) {
    dtl::plain_bitmap<word_type> b(4096);
    two_state_markov_process mp(2.0, 0.25);
    for (std::size_t i = 0; i < 4096; ++i) {
      b.set(i, mp.next());
    }

    R rank;
    rank.init(b.data_begin(), b.data_end());

    for (std::size_t i = 0; i < 4096; ++i) {
      const auto r = rank(i, b.data());
      const auto l = linear_rank(i, b, rank.is_inclusive());
      ASSERT_EQ(r, l)
          << "bitmap: " << b << "\n"
          << "test failed for i=" << i << "."
          << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//
