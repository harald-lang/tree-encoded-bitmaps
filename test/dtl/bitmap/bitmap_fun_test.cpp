#include "gtest/gtest.h"

#include <chrono>

#include <dtl/dtl.hpp>
#include <dtl/bitmap/util/bitmap_fun.hpp>
#include <boost/dynamic_bitset.hpp>
#include <dtl/iterator.hpp>
#include <dtl/bitmap/util/bitmap_view.hpp>

//===----------------------------------------------------------------------===//
TEST(bitmap_fun,
     fetch_words) {
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
TEST(bitmap_fun,
     fetch_bits) {
  using bmf = dtl::bitmap_fun<$u64>;
  std::vector<$u64> bitmap;
  bitmap.push_back(0x0101010101010101);
  bitmap.push_back(0x1010101010101010);
  bitmap.push_back(0xFFFFFFFFFFFFFFFF);
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 0, 0+1),
            0x0000000000000001);
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 0, 0+9),
            0x0000000000000101);
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 1, 1+8),
            0b00000010000000);
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 8, 64),
            0x0001010101010101);
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 63,  63+8),
            0x0000000000000020);
}
//===----------------------------------------------------------------------===//
TEST(bitmap_fun,
    find_first_find_last) {
  boost::dynamic_bitset<$u64> bs_fun(128);
  dtl::data_view<u64> bs_data {
      bs_fun.m_bits.data(),
      bs_fun.m_bits.data() + bs_fun.m_bits.size()
  };
  dtl::bitmap_view<$u64> bs_view(bs_data);
  ASSERT_EQ(bs_view.find_first(), 128);
  ASSERT_EQ(bs_view.find_last(), 128);
  bs_fun[0] = true;
  ASSERT_EQ(bs_view.find_first(), 0);
  ASSERT_EQ(bs_view.find_last(), 0);
  bs_fun[0] = false;
  bs_fun[127] = true;
  ASSERT_EQ(bs_view.find_first(), 127);
  ASSERT_EQ(bs_view.find_last(), 127);
  bs_fun[127] = false;
  for (std::size_t i = 0; i < 128; ++i) {
    bs_fun[i] = true;
    ASSERT_EQ(bs_view.find_first(), i);
    ASSERT_EQ(bs_view.find_last(), i);
    bs_fun[i] = false;
  }
  bs_fun[2] = true;
  bs_fun[3] = true;
  ASSERT_EQ(bs_view.find_first(), 2);
  ASSERT_EQ(bs_view.find_last(), 3);
}
//===----------------------------------------------------------------------===//
