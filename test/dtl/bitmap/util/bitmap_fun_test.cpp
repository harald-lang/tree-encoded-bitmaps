#include "gtest/gtest.h"

#include <dtl/bitmap/util/bitmap_fun.hpp>
#include <dtl/bitmap/util/bitmap_view.hpp>
#include <dtl/dtl.hpp>
#include <dtl/iterator.hpp>

#include <boost/dynamic_bitset.hpp>

#include <dtl/bitmap/util/plain_bitmap.hpp>
#include <chrono>
//===----------------------------------------------------------------------===//
TEST(bitmap_fun,
    fetch_words) {
  using bmf = dtl::bitmap_fun<$u64>;
  std::vector<$u64> bitmap;
  bitmap.push_back(0x0101010101010101);
  bitmap.push_back(0x1010101010101010);
  bitmap.push_back(0xFFFFFFFFFFFFFFFF);
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 0, 0 + 64),
      0x0101010101010101);
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 64, 64 + 64),
      0x1010101010101010);
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 128, 128 + 64),
      0xFFFFFFFFFFFFFFFF);
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 64 + 1, 64 + 1 + 64),
      (0x1010101010101010 >> 1) | (1ul << 63));
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 1, 1 + 64),
      (0x0101010101010101 >> 1));
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 8, 8 + 64),
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
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 0, 0 + 1),
      0x0000000000000001);
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 0, 0 + 9),
      0x0000000000000101);
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 1, 1 + 8),
      0b00000010000000);
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 8, 64),
      0x0001010101010101);
  ASSERT_EQ(bmf::fetch_bits(bitmap.data(), 63, 63 + 8),
      0x0000000000000020);
}
//===----------------------------------------------------------------------===//
TEST(bitmap_fun,
    find_first_find_last) {
  boost::dynamic_bitset<$u64> bs_fun(128);
  dtl::data_view<$u64> bs_data {
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
TEST(bitmap_fun,
    find_first_in_range__word_aligned) {
  // We use the plain bitmap here for simplicity. Internally, the call is
  // forwarded to bitmap_fun.
  const auto n = 128;
  dtl::plain_bitmap<$u32> bs(n);
  for (std::size_t i = 0; i < n; ++i) {
    bs.set(i);
    ASSERT_EQ(bs.find_first(0, n), i);
    bs.clear(i);
  }
}
//===----------------------------------------------------------------------===//
TEST(bitmap_fun,
    find_first_in_range__with_offset_at_the_begin) {
  // We use the plain bitmap here for simplicity. Internally, the call is
  // forwarded to bitmap_fun.
  const auto n = 128;
  dtl::plain_bitmap<$u32> bs(n);
  for (std::size_t offset = 1; offset < 64; ++offset) {
    for (std::size_t i = offset; i < n - offset; ++i) {
      // Set all bit outside of the range. Those should be ignored.
      bs.set(0, offset);

      bs.set(i);
      ASSERT_EQ(bs.find_first(offset, n), i)
        << " offset was " << offset;

      bs.clear(0, offset);
      bs.clear(i);
    }
  }
}
//===----------------------------------------------------------------------===//
TEST(bitmap_fun,
    find_first_in_range__with_offset_at_the_end) {
  // We use the plain bitmap here for simplicity. Internally, the call is
  // forwarded to bitmap_fun.
  const std::size_t n = 128;
  dtl::plain_bitmap<$u32> bs(n);
  for (std::size_t offset = 1; offset < 64; ++offset) {
    for (std::size_t i = 0; i < n - offset; ++i) {
      // Set all bit outside of the range. Those should be ignored.
      bs.set(n - offset, n);

      bs.set(i);
      ASSERT_EQ(bs.find_first(0, n - offset), i)
        << " offset was " << offset;

      bs.clear(n - offset, n);
      bs.clear(i);
    }
  }
}
//===----------------------------------------------------------------------===//
TEST(bitmap_fun,
    find_first_in_range__with_offsets) {
  // We use the plain bitmap here for simplicity. Internally, the call is
  // forwarded to bitmap_fun.
  const std::size_t n = 256;
  dtl::plain_bitmap<$u32> bs(n);
  for (std::size_t offset_b = 1; offset_b < 64; ++offset_b) {
    for (std::size_t offset_e = 1; offset_e < 64; ++offset_e) {
      for (std::size_t i = offset_b; i < n - offset_b - offset_e; ++i) {
        // Set all bit outside of the range. Those should be ignored.
        bs.set(0, offset_b);
        bs.set(n - offset_e, n);

        bs.set(i);
        ASSERT_EQ(bs.find_first(offset_b, n - offset_e), i)
          << " offsets were " << offset_b << "/" << offset_e;

        bs.clear(0, offset_b);
        bs.clear(n - offset_e, n);
        bs.clear(i);
      }
    }
  }
}
//===----------------------------------------------------------------------===//
TEST(bitmap_fun,
    find_first_in_range__full) {
  // We use the plain bitmap here for simplicity. Internally, the call is
  // forwarded to bitmap_fun.
  const std::size_t n = 128;
  dtl::plain_bitmap<$u32> bs(n);
  for (std::size_t range_offset = 0; range_offset < n; ++range_offset) {
    for (std::size_t range_len = 0; range_len < n; ++range_len) {
      if (range_offset + range_len >= n) continue;
      for (std::size_t i = range_offset; i < range_offset + range_len; ++i) {
        // Set all bit outside of the range. Those should be ignored.
        bs.set(0, range_offset);
        bs.set(range_offset + range_len, n);

        bs.set(i);
        ASSERT_EQ(bs.find_first(range_offset, range_offset + range_len), i)
          << " range offset was " << range_offset << ", length "
          << range_len;

        bs.clear(0, range_offset);
        bs.clear(range_offset + range_len, n);
        bs.clear(i);
      }
    }
  }
}
//===----------------------------------------------------------------------===//
TEST(bitmap_fun,
    to_positions__single_32_bit_word) {
  using word_type = $u32;
  constexpr auto word_bitwidth = sizeof(word_type) * 8;
  std::vector<$u32> positions(word_bitwidth, 0);
  std::size_t result_cnt = 0;
  u32 off = 42;

  auto test = [&](word_type w) {
    result_cnt =
        dtl::bitmap_fun<word_type>::to_positions(w, positions.data(), off);
    ASSERT_EQ(result_cnt, dtl::bits::pop_count(w))
        << "w=" << std::bitset<word_bitwidth>(w);

    auto* reader = positions.data();
    for (std::size_t i = 0; i < word_bitwidth; ++i) {
      if (dtl::bits::bit_test(w, i)) {
        ASSERT_EQ(*reader, i + off)
            << "w=" << std::bitset<word_bitwidth>(w);
        reader++;
      }
    }
  };

  // Testing the lower 16 bits.
  for (word_type w = 0; w < 1u << 16 - 1; ++w) { test(w); }
  // Testing the higher 16 bits.
  for (word_type w = 0; w < 1u << 16 - 1; ++w) { test(w << 16); }
  test(~word_type(0));
}
//===----------------------------------------------------------------------===//
TEST(bitmap_fun,
    to_positions__single_64_bit_word) {
  using word_type = $u64;
  constexpr auto word_bitwidth = sizeof(word_type) * 8;
  std::vector<$u32> positions(word_bitwidth, 0);
  std::size_t result_cnt = 0;
  u32 off = 42;

  auto test = [&](word_type w) {
    result_cnt =
        dtl::bitmap_fun<word_type>::to_positions(w, positions.data(), off);
    ASSERT_EQ(result_cnt, dtl::bits::pop_count(w))
        << "w=" << std::bitset<word_bitwidth>(w);

    auto* reader = positions.data();
    for (std::size_t i = 0; i < word_bitwidth; ++i) {
      if (dtl::bits::bit_test(w, i)) {
        ASSERT_EQ(*reader, i + off)
            << "w=" << std::bitset<word_bitwidth>(w);
        reader++;
      }
    }
  };

  for (word_type w = 0; w < 1u << 16 - 1; ++w) { test(w); }
  for (word_type w = 0; w < 1u << 16 - 1; ++w) { test(w << 16); }
  for (word_type w = 0; w < 1u << 16 - 1; ++w) { test(w << 32); }
  for (word_type w = 0; w < 1u << 16 - 1; ++w) { test(w << 48); }
  test(~word_type(0));
}
//===----------------------------------------------------------------------===//
TEST(bitmap_fun,
    to_positions__word_aligned_bitmap) {
  using word_type = $u64;
  constexpr auto word_bitwidth = sizeof(word_type) * 8;

  dtl::plain_bitmap<word_type> bm(word_bitwidth * 10);
  std::vector<$u32> positions(bm.size(), 0);

  std::size_t result_cnt = 0;
  u32 off = 42;

  for (std::size_t distance = 1; distance < bm.size(); distance <<= 1) {
    bm.clear(0, bm.size());
    std::size_t popcnt = 0;
    for (std::size_t i = 0; i < bm.size(); i += distance) {
      bm.set(i);
      ++popcnt;
    }

    dtl::bitmap_fun<word_type>::to_positions(
        bm.data(),
        bm.data() + ((bm.size() + word_bitwidth - 1) / word_bitwidth),
        positions.data(),
        off);

    auto* reader = positions.data();
    for (std::size_t i = 0; i < bm.size(); ++i) {
      if (bm[i] == true) {
        ASSERT_EQ(*reader, i + off);
        reader++;
      }
    }
  }
}
//===----------------------------------------------------------------------===//
TEST(bitmap_fun,
    to_positions__range) {
  using word_type = $u64;
  constexpr auto word_bitwidth = sizeof(word_type) * 8;

  dtl::plain_bitmap<word_type> bm(word_bitwidth * 10);
  std::vector<$u32> positions(bm.size(), 0);

  std::size_t result_cnt = 0;

  bm.set(0, bm.size());

  // Vary the range size.
  for (std::size_t range_length = 1; range_length <= bm.size(); ++range_length) {
    // Slide the range over the bitmap.
    for (std::size_t b = 0; b < bm.size() - range_length; ++b) {
      const std::size_t e = b + range_length;

      // Set all bits in the bitmap.
      {
        bm.set(0, bm.size());

        dtl::bitmap_fun<word_type>::to_positions(
            bm.data(), b, e, positions.data());

        auto* reader = positions.data();
        for (std::size_t i = b; i < e; ++i) {
          if (bm[i] == true) {
            ASSERT_EQ(*reader, i)
                << "range_length=" << range_length
                << ", b=" << b << ", e=" << e;
            reader++;
          }
        }
      }

      // Set only the bits within the given range.
      {
        bm.clear(0, bm.size());
        bm.set(b, e);

        auto cnt = dtl::bitmap_fun<word_type>::to_positions(
            bm.data(), b, e, positions.data());

        ASSERT_EQ(cnt, range_length)
            << "range_length=" << range_length
            << ", b=" << b << ", e=" << e;

        auto* reader = positions.data();
        for (std::size_t i = b; i < e; ++i) {
          if (bm[i] == true) {
            ASSERT_EQ(*reader, i)
                << "range_length=" << range_length
                << ", b=" << b << ", e=" << e;
            reader++;
          }
        }
      }

      // Set the bits outside the given range.
      {
        bm.clear(0, bm.size());
        bm.set(0, b);
        bm.set(e, bm.size());

        auto cnt = dtl::bitmap_fun<word_type>::to_positions(
            bm.data(), b, e, positions.data());

        ASSERT_EQ(cnt, 0)
            << "range_length=" << range_length
            << ", b=" << b << ", e=" << e;
      }
    }
  }
}
//===----------------------------------------------------------------------===//
TEST(bitmap_fun,
    find_next_zero) {
  using word_type = $u64;
  constexpr auto word_bitwidth = sizeof(word_type) * 8;

  const std::size_t n = word_bitwidth * 3;
  dtl::plain_bitmap<word_type> bm(n);

  for (std::size_t b = 0; b < n; ++b) {
    for (std::size_t e = b; e < n; ++e) {
      bm.clear(0, bm.size());
      bm.set(b, e);
//      std::cout << bm << std::endl;
      if (b == e) {
        ASSERT_EQ(bm.find_next_zero(b, n), std::min(b + 1, n))
                      << " with b = " << b << " and e = " << e << std::endl;
      }
      else {
        ASSERT_TRUE(bm.test(b));
        ASSERT_TRUE(bm.test(e - 1));
        ASSERT_EQ(bm.find_next_zero(b, n), std::min(e, n))
                      << " with b = " << b << " and e = " << e << std::endl;
      }
    }
  }

  bm.set(1, bm.size());
  ASSERT_EQ(bm.size(), bm.find_next_zero(0, bm.size()));
  bm.clear(1);
  ASSERT_EQ(1, bm.find_next_zero(0, bm.size()));
  bm.set(0);
  ASSERT_EQ(1, bm.find_next_zero(0, bm.size()));
  bm.set(1);
  bm.clear(2);
  ASSERT_EQ(2, bm.find_next_zero(0, bm.size()));
}
//===----------------------------------------------------------------------===//
