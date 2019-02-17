#include "gtest/gtest.h"

#include <chrono>

#include <dtl/dtl.hpp>
#include <dtl/bitmap/util/bitmap_fun.hpp>

//===----------------------------------------------------------------------===//
TEST(bitmap_fun,
     fetch_bits) {
  using bmf = dtl::bitmap_fun<$u64>;
  std::vector<$u64> bitmap;
  bitmap.push_back(0x0101010101010101);
  bitmap.push_back(0x1010101010101010);
  bitmap.push_back(0xFFFFFFFFFFFFFFFF);
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 0,    0+64),
      0x0101010101010101);
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 64,   64+64),
      0x1010101010101010);
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 128,  128+64),
      0xFFFFFFFFFFFFFFFF);
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 64+1, 64+1+64),
      (0x1010101010101010 >> 1) | (1ul << 63));
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 1,    1+64),
      (0x0101010101010101 >> 1));
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 8,    8+64),
      (0x0101010101010101 >> 8) | (0x1010101010101010 << (64 - 8)));
}
//===----------------------------------------------------------------------===//
