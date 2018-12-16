#include "gtest/gtest.h"

#include <dtl/dtl.hpp>

#include "dynamic_bitmap.hpp"
#include "dynamic_partitioned_tree_mask.hpp"
#include "dynamic_tree_mask_lo.hpp"
#include "two_state_markov_process.hpp"
#include "roaring_bitmap.hpp"
#include "tree_mask_lo.hpp"
#include "tree_mask_po.hpp"
#include "wah.hpp"
#include "test_utils.hpp"

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
class iterate_test : public ::testing::Test {};

// Specify the types for which we want to run the API tests.
using types_under_test = ::testing::Types<
    dtl::dynamic_tree_mask_lo,
    dtl::dynamic_partitioned_tree_mask
//    dtl::dynamic_bitmap<$u32>
//    roaring_bitmap,
//    tree_mask_lo,
//    tree_mask_po,
//    wah32
//    wah64 // FIXME: Position iterator does not work with 64-bit impl.
>;
TYPED_TEST_CASE(iterate_test, types_under_test);


//===----------------------------------------------------------------------===//
template<typename T>
dtl::bitmap
decode_using_iterator(T bitmap) {
  dtl::bitmap ret_val(bitmap.size());
  auto it = bitmap.it();
  while (!it.end()) {
    for (std::size_t i = it.pos(); i < it.pos() + it.length(); ++i) {
      ret_val.set(i);
    }
    it.next();
  }
  return ret_val;
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
template<typename T>
dtl::bitmap
bitwise_and_using_iterator(T bitmap_a, T bitmap_b) {
  dtl::bitmap ret_val(bitmap_a.size());
  auto it_a = bitmap_a.it();
  auto it_b = bitmap_b.it();
  while (!(it_a.end() && it_b.end())) {
    const auto a_begin = it_a.pos();
    const auto a_end   = it_a.pos() + it_a.length();
    const auto b_begin = it_b.pos();
    const auto b_end   = it_b.pos() + it_b.length();

    const auto begin_max = std::max(a_begin, b_begin);
    const auto end_min   = std::min(a_end,   b_end);

    for (std::size_t i = begin_max; i < end_min; ++i) {
      ret_val[i] = true;
    }

    if (a_end == b_end) {
      it_a.next();
      it_b.next();
    }
    else if (a_end <= b_end) {
      it_a.next();
    }
    else {
      it_b.next();
    }
  }
  return ret_val;
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
template<typename T>
dtl::bitmap
bitwise_and_using_skip_iterator(T bitmap_a, T bitmap_b) {
  dtl::bitmap ret_val(bitmap_a.size());
  auto it_a = bitmap_a.it();
  auto it_b = bitmap_b.it();
  while (!(it_a.end() && it_b.end())) {
    const auto a_begin = it_a.pos();
    const auto a_end   = it_a.pos() + it_a.length();
    const auto b_begin = it_b.pos();
    const auto b_end   = it_b.pos() + it_b.length();

    const auto begin_max = std::max(a_begin, b_begin);
    const auto end_min   = std::min(a_end,   b_end);

    for (std::size_t i = begin_max; i < end_min; ++i) {
      ret_val[i] = true;
    }

    u1 overlap = begin_max < end_min;

    if (overlap) {
      if (a_end == b_end) {
        it_a.next();
        it_b.next();
      }
      else if (a_end <= b_end) {
        it_a.next();
      }
      else {
        it_b.next();
      }
    }
    else {
      // no overlap
      if (a_end < b_end) {
        it_a.skip_to(b_begin);
      }
      else {
        it_b.skip_to(a_begin);
      }
    }
  }
  return ret_val;
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
TYPED_TEST(iterate_test, decode) {
  using T = TypeParam;

  for ($u32 i = 0; i < 256; ++i) {
    dtl::bitmap bm(8, i);
    T enc_bm(bm);
    auto dec_bm = decode_using_iterator(enc_bm);
    ASSERT_EQ(bm, dec_bm);
  }
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
TYPED_TEST(iterate_test, bitwise_and) {
  using T = TypeParam;

  for ($u32 i = 0; i < 256; ++i) {
    for ($u32 j = 0; j < 256; ++j) {
      dtl::bitmap bm_a(8, i);
      dtl::bitmap bm_b(8, j);
      T enc_bm_a(bm_a);
      T enc_bm_b(bm_b);
      auto result = bitwise_and_using_iterator(enc_bm_a, enc_bm_b);
      ASSERT_EQ(bm_a & bm_b, result);
    }
  }
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
TYPED_TEST(iterate_test, bitwise_and_skip) {
  using T = TypeParam;

  for ($u32 i = 0; i < 256; ++i) {
    for ($u32 j = 0; j < 256; ++j) {
//      std::cout
//          << "a=" << std::bitset<8>(i)
//          << ", b=" << std::bitset<8>(j)
//          << std::endl;
      dtl::bitmap bm_a(8, i);
      dtl::bitmap bm_b(8, j);
      T enc_bm_a(bm_a);
      T enc_bm_b(bm_b);
      auto result = bitwise_and_using_skip_iterator(enc_bm_a, enc_bm_b);
      ASSERT_EQ(bm_a & bm_b, result);
    }
  }
}
//===----------------------------------------------------------------------===//

