#include "gtest/gtest.h"

#include <dtl/dtl.hpp>
#include <dtl/bitmap.hpp>
#include <dtl/bitmap/util/bitmap_tree.hpp>

//===----------------------------------------------------------------------===//
// Set the optimization level to 1, which means that we enable implicit tree
// nodes, but restricted pruning is not used.
static constexpr auto opt_level = 1;
TEST(bitmap_tree,
    test) {
//  dtl::bitmap bm(4, 0b0000);
//  dtl::bitmap bm(4, 0b0100);
//  dtl::bitmap bm(4, 0b1110);
//  dtl::bitmap bm(4, 0b0111);
//  dtl::bitmap bm(4, 0b0001);
//  dtl::bitmap_tree<opt_level> bt(bm);
//  std::cout << "\n" << bt << std::endl;
//  std::cout << "\n" << bt.info() << std::endl;

  for (std::size_t i = 1; i < 16; ++i) {
    std::cout << "testing " << std::bitset<4>(i) << std::endl;
    dtl::bitmap bm(4, i);
    dtl::bitmap_tree<opt_level> bt(bm);
  }
}
//===----------------------------------------------------------------------===//
