#include "gtest/gtest.h"

#include <dtl/bitmap/teb_scan_util.hpp>
#include <dtl/dtl.hpp>
//===----------------------------------------------------------------------===//
TEST(teb_scan_util,
     fetch_bits) {
  std::vector<$u64> bitmap;
  bitmap.push_back(0x0101010101010101);
  bitmap.push_back(0x1010101010101010);
  bitmap.push_back(0xFFFFFFFFFFFFFFFF);
  dtl::data_view<u64> bitmap_view {
      bitmap.data(), bitmap.data() + 3,
  };
  ASSERT_EQ(fetch_bits(bitmap_view, -3,   -1),
            0b0000000000011);
  ASSERT_EQ(fetch_bits(bitmap_view, -1,   2),
            0b0000000000011);
  ASSERT_EQ(fetch_bits(bitmap_view, 130,  133),
            0b111);
  ASSERT_EQ(fetch_bits(bitmap_view, 192,  200),
            0);
}
//===----------------------------------------------------------------------===//
TEST(teb_scan_util,
     fetch_bits_from_empty_bitmap) {
  std::vector<$u64> bitmap;
  dtl::data_view<u64> bitmap_view { bitmap.data(), bitmap.data() };
  ASSERT_EQ(fetch_bits(bitmap_view, -3,   3),
            0b000111);
  ASSERT_EQ(fetch_bits(bitmap_view, -63,  1),
            (~0ul) >> 1);
}
//===----------------------------------------------------------------------===//
