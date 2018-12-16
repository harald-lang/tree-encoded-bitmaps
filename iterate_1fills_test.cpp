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
class iterate_1fills_test : public ::testing::Test {};

// Specify the types for which we want to run the API tests.
using types_under_test = ::testing::Types<
    dtl::dynamic_tree_mask_lo,
    dtl::dynamic_bitmap<$u32>
//    dtl::dynamic_partitioned_tree_mask
//    roaring_bitmap,
//    tree_mask_lo,
//    tree_mask_po,
//    wah32
//    wah64 // FIXME: Position iterator does not work with 64-bit impl.
>;
TYPED_TEST_CASE(iterate_1fills_test, types_under_test);



//===----------------------------------------------------------------------===//
TYPED_TEST(iterate_1fills_test, all_bits_set) {
  using T = TypeParam;
  dtl::bitmap b(64);
  b.set();
  T t(b);
  std::cout << t << std::endl;
  auto it = t.it();
  ASSERT_TRUE(!it.end());
  ASSERT_EQ(it.pos(), 0);
  ASSERT_EQ(it.length(), 64);
  it.next();
  ASSERT_TRUE(it.end());
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
TYPED_TEST(iterate_1fills_test, no_bits_set) {
  using T = TypeParam;
  dtl::bitmap b(64);
  T t(b);
  std::cout << t << std::endl;
  auto it = t.it();
  ASSERT_TRUE(it.end());
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
TYPED_TEST(iterate_1fills_test, single_1fill) {
  using T = TypeParam;
  dtl::bitmap b(64, 0b11110000);
  T t(b);
  std::cout << t << std::endl;
  auto it = t.it();
  ASSERT_TRUE(!it.end());
  ASSERT_EQ(it.pos(), 4);
  ASSERT_EQ(it.length(), 4);
  it.next();
  ASSERT_TRUE(it.end());
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
TYPED_TEST(iterate_1fills_test, multiple_1fills_at_different_levels) {
  using T = TypeParam;
  dtl::bitmap b(64, 0b0011000011110000);
  T tm(b);
  std::cout << tm << std::endl;
  auto it = tm.it();
  ASSERT_TRUE(!it.end());
  ASSERT_EQ(it.pos(), 4);
  ASSERT_EQ(it.length(), 4);
  it.next();
  ASSERT_EQ(it.pos(), 12);
  ASSERT_EQ(it.length(), 2);
  ASSERT_TRUE(!it.end());
  it.next();
  ASSERT_TRUE(it.end());
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
template<typename T>
void
skip_test(u64 n, u64 bitmap, u64 skip_to_pos, u64 expected_pos, u64 expected_len) {
  assert(n <= 64);
  dtl::bitmap b(n, bitmap);
  T tm(b);
  std::cout << tm << std::endl;
  auto it = tm.it();
  it.skip_to(skip_to_pos);
  if (expected_len == 0) {
    ASSERT_TRUE(it.end());
  }
  else {
    ASSERT_TRUE(!it.end());
  }
  ASSERT_EQ(it.pos(), expected_pos);
  if (expected_pos != n) {
//    ASSERT_EQ(it.length(), expected_len);
  }
}

TYPED_TEST(iterate_1fills_test, skip_to_1fill) {
  using T = TypeParam;
  skip_test<T>(8, 0b00000101, 2, 2, 1);
  skip_test<T>(8, 0b00001001, 3, 3, 1);
  skip_test<T>(8, 0b00010001, 4, 4, 1);
  skip_test<T>(8, 0b00100001, 5, 5, 1);
  skip_test<T>(8, 0b01000001, 6, 6, 1);
  skip_test<T>(8, 0b10000001, 7, 7, 1);
}

TYPED_TEST(iterate_1fills_test, skip_to_1fill_passing_the_root_node) {
  using T = TypeParam;
  skip_test<T>(8, 0b01000001, 6, 6, 1);
  skip_test<T>(8, 0b10000001, 7, 7, 1);
}

TYPED_TEST(iterate_1fills_test, skip_before_1fill) {
  using T = TypeParam;
  skip_test<T>(8, 0b01000001, 5, 6, 1);
  skip_test<T>(8, 0b01000001, 4, 6, 1);
  skip_test<T>(8, 0b01000001, 3, 6, 1);
  skip_test<T>(8, 0b01000001, 2, 6, 1);
  skip_test<T>(8, 0b01000001, 1, 6, 1);
}

TYPED_TEST(iterate_1fills_test, skip_beyond_last_1fill) {
  using T = TypeParam;
  skip_test<T>(8, 0b01000001, 7, 8, 0);
}

TYPED_TEST(iterate_1fills_test, skip_into_a_1fill) {
  using T = TypeParam;
  skip_test<T>(8, 0b11111101, 6, 6, 2);
}

TYPED_TEST(iterate_1fills_test, skip_multiple_times) {
  using T = TypeParam;
  dtl::bitmap b(8, 0b11001101);
  T tm(b);
  std::cout << tm << std::endl;
  auto it = tm.it();
  it.skip_to(2);
  ASSERT_EQ(it.pos(), 2);
//  ASSERT_EQ(it.length(), 2);
  it.skip_to(6);
  ASSERT_EQ(it.pos(), 6);
//  ASSERT_EQ(it.length(), 2);
}

TYPED_TEST(iterate_1fills_test, skip_into_a_1fill_length_one) {
  using T = TypeParam;
  skip_test<T>(16, 0b1100000110001111, 6, 7, 1);
}

TEST(iterate_1fills_full_8bit_test, skip_to_test_all_combinations) {
  using T = dtl::dynamic_tree_mask_lo;
  using range = dtl::dynamic_tree_mask_lo::range;
  #define N 8

  for($u64 i = 0; i < 256; i++){

    assert(N <= 64);
    dtl::bitmap b(N, i);

    std::bitset<N> b_set(static_cast<u8>(i));

    std::cout << std::endl << "bitmap: " << b_set << std::endl;

    std::vector<dtl::dynamic_tree_mask_lo::range> r;
    T tm(b);
    std::cout << tm << std::endl;

    auto it_1 = tm.it();
    while(!it_1.end()){
      range tmp(it_1.pos(), it_1.length());
      r.push_back(tmp);
      it_1.next();
    }

    std::cout << "All possible ranges: ";

    for(range e : r){
      e.print();
      std::cout << " ";
    }
    std::cout << std::endl;


    for($u8 p = 0; p < N; p++){
      auto it_2 = tm.it();
      std::cout << "pos=" << it_2.pos() << std::endl;

      if(!it_2.end()){
        it_2.skip_to(p);

        std::cout << "Pos: " << it_2.pos() << " Len: " << it_2.length() << std::endl;

        if (it_2.length() == 0 || it_2.pos() == N) {
          ASSERT_TRUE(it_2.end());
        }       else {
          ASSERT_TRUE(!it_2.end());
        }

        if(!it_2.end()){
          bool found_range = false;

          for(range e: r){
            if(e.first == it_2.pos() && (e.last-e.first) == it_2.length()){
              found_range = true;
              break;
            }
          }

          ASSERT_TRUE(found_range);
        }
      }
    }
  }
}

//===----------------------------------------------------------------------===//

