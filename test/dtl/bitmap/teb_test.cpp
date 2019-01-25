#include "gtest/gtest.h"

#include <dtl/dtl.hpp>
#include <dtl/bitmap.hpp>
#include <dtl/bitmap/util/random.hpp>

//===----------------------------------------------------------------------===//
// Test the compression levels of TEBs.
TEST(teb,
     ensure_that_higher_compression_do_not_increase_the_size) {
  for (auto n = 1ull << 10; n <= 1ull << 10; n <<= 1) {
    for (auto f = 1; f <= 16; f *= 2) {
      for (auto d : {0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5} ) {
        for (std::size_t rep = 0; rep < 1; ++rep) {
          dtl::bitmap bs = dtl::gen_random_bitmap(n, f, d);
          std::cout << "plain: " << bs << std::endl;
          dtl::teb<0> teb_o0(bs);
          std::cout << "o0:\n" << teb_o0 << std::endl;
          dtl::teb<1> teb_o1(bs);
          std::cout << "o1:\n" << teb_o1 << std::endl;
          dtl::teb<2> teb_o2(bs);
          std::cout << "o2:\n" << teb_o2 << std::endl;

          ASSERT_LE(teb_o1.size_in_byte(), teb_o0.size_in_byte())
              << "\n"
              << "plain:  " << bs
              << "\n"
              << "TEB-o0: " << teb_o0
              << "\n"
              << "TEB-o1: " << teb_o1
              << std::endl;
          ASSERT_LE(teb_o2.size_in_byte(), teb_o1.size_in_byte())
              << "\n"
              << "plain:  " << bs
              << "\n"
              << "TEB-o1: " << teb_o1
              << "\n"
              << "TEB-o2: " << teb_o2
              << std::endl;
        }
      }
    }
  }
}
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
TEST(teb,
     ensure_the_correct_number_of_perfect_levels) {
  using teb = dtl::teb<2>;
  //                      implicit_inner_node_cnt, expected # of perfect levels
  ASSERT_EQ(teb::determine_perfect_tree_levels( 0), 1);

  ASSERT_EQ(teb::determine_perfect_tree_levels( 1), 2);
  ASSERT_EQ(teb::determine_perfect_tree_levels( 2), 2);

  ASSERT_EQ(teb::determine_perfect_tree_levels( 3), 3);
  ASSERT_EQ(teb::determine_perfect_tree_levels( 4), 3);
  ASSERT_EQ(teb::determine_perfect_tree_levels( 5), 3);
  ASSERT_EQ(teb::determine_perfect_tree_levels( 6), 3);

  ASSERT_EQ(teb::determine_perfect_tree_levels( 7), 4);
  ASSERT_EQ(teb::determine_perfect_tree_levels( 8), 4);
  ASSERT_EQ(teb::determine_perfect_tree_levels( 9), 4);
  ASSERT_EQ(teb::determine_perfect_tree_levels(10), 4);
  ASSERT_EQ(teb::determine_perfect_tree_levels(11), 4);
  ASSERT_EQ(teb::determine_perfect_tree_levels(12), 4);
  ASSERT_EQ(teb::determine_perfect_tree_levels(13), 4);
  ASSERT_EQ(teb::determine_perfect_tree_levels(14), 4);

  ASSERT_EQ(teb::determine_perfect_tree_levels(15), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(16), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(17), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(18), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(19), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(20), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(21), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(22), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(23), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(24), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(25), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(26), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(27), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(28), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(29), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(30), 5);

  ASSERT_EQ(teb::determine_perfect_tree_levels(31), 6);
  ASSERT_EQ(teb::determine_perfect_tree_levels(32), 6);
}