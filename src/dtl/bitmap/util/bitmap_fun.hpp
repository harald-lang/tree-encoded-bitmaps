#pragma once
//===----------------------------------------------------------------------===//
#include <cassert>

#include <dtl/dtl.hpp>
#include <dtl/bits.hpp>

#include "config.hpp"
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// Static functions to work with plain bitmaps.
template<typename _word_type>
struct bitmap_fun {

  using word_type = _word_type;
  static constexpr std::size_t word_bitlength = sizeof(word_type) * 8;

  /// Test the bit at position i.
  static u1 __teb_inline__
  test(const word_type* bitmap, std::size_t i) noexcept {
    const auto block_idx = i / word_bitlength;
    const auto word = bitmap[block_idx];
    const auto bit_idx = i % word_bitlength;
    return (word & (word_type(1) << i)) != 0;
  }

  /// Set the bit at position i.
  static void __teb_inline__
  set(word_type* bitmap, std::size_t i, u1 val) noexcept {
    if (val) {
      set(bitmap, i);
    }
    else {
      clear(bitmap, i);
    }
  }

  /// Set the bit at position i.
  static void __teb_inline__
  set(word_type* bitmap, std::size_t i) noexcept {
    const auto block_idx = i / word_bitlength;
    const auto bit_idx = i % word_bitlength;
    bitmap[block_idx] |= word_type(1) << bit_idx;
  }

  /// Set the bits in [b,e).
  static void __teb_inline__
  set(word_type* bitmap, std::size_t b, std::size_t e) noexcept {
    assert(b < e);
    // The algorithm below has been adapted from the paper "Consistently faster
    // and smaller compressed bitmaps with Roaring" by Lemire et al.
    const auto x = b / word_bitlength;
    const auto y = (e - 1) / word_bitlength;
    const word_type Z = ~word_type(0);

    const word_type X = Z << (b % word_bitlength);
    const word_type Y = Z >> ((64 - (e % word_bitlength)) % word_bitlength);
    if (x == y) {
      bitmap[x] |= (X & Y);
    }
    else {
      bitmap[x] |= X;
      for (std::size_t k = x + 1; k < y; ++k) {
        bitmap[k] |= Z;
      }
      bitmap[y] |= Y;
    }
  }

  /// Clear the bit at position i.
  static void __teb_inline__
  clear(word_type* bitmap, std::size_t i) noexcept {
    const auto block_idx = i / word_bitlength;
    const auto bit_idx = i % word_bitlength;
    bitmap[block_idx] &= ~(word_type(1) << bit_idx);
  }

  /// Clear the bits in [b,e).
  static void __teb_inline__
  clear(word_type* bitmap, std::size_t b, std::size_t e) noexcept {
    assert(b < e);
    // The algorithm below has been adapted from the paper "Consistently faster
    // and smaller compressed bitmaps with Roaring" by Lemire et al.
    const auto x = b / word_bitlength;
    const auto y = (e - 1) / word_bitlength;
    const word_type Z = ~word_type(0);

    const word_type X = Z << (b % word_bitlength);
    const word_type Y = Z >> ((64 - (e % word_bitlength)) % word_bitlength);
    if (x == y) {
      bitmap[x] &= ~(X & Y);
    }
    else {
      bitmap[x] &= ~X;
      for (std::size_t k = x + 1; k < y; ++k) {
        bitmap[k] &= ~Z;
      }
      bitmap[y] &= ~Y;
    }
  }

  /// Fetch up to size(word_type)*8 consecutive bits.
  static word_type __teb_inline__
  fetch_bits(const word_type* bitmap,
      u64 bit_idx_begin,
      u64 bit_idx_end /* non-inclusive */) {
    assert(bit_idx_end > bit_idx_begin);
    static constexpr u64 word_bitlength = sizeof(word_type) * 8;
    const auto word_idx_begin = bit_idx_begin / word_bitlength;
    const auto word_idx_end = (bit_idx_end - 1) / word_bitlength;
    assert(word_idx_end - word_idx_begin <= 1);
    u64 cnt = bit_idx_end - bit_idx_begin;
    if (word_idx_begin == word_idx_end) {
      const auto word_idx = word_idx_begin;
      word_type bitmap_word = bitmap[word_idx];
      bitmap_word >>= (bit_idx_begin % word_bitlength);
      bitmap_word &= (~word_type(0)) >> (word_bitlength - cnt);
      return bitmap_word;
    }
    else {
      word_type bitmap_word_0 = bitmap[word_idx_begin];
      word_type bitmap_word_1 = bitmap[word_idx_end];
      bitmap_word_0 >>= (bit_idx_begin % word_bitlength);
      bitmap_word_1 &= (~word_type(0)) >> (word_bitlength -
          ((bit_idx_end % word_bitlength)));
      return bitmap_word_0 | (bitmap_word_1 << (word_bitlength -
          (bit_idx_begin % word_bitlength)));
    }
  }

  /// Find the first set bit.  Returns the index of the first set bit. If no
  /// bits are set, the length of the bitmap is returned.
  static std::size_t
  find_first(const word_type* bitmap_begin, const word_type* bitmap_end) {
    const std::size_t word_cnt = bitmap_end - bitmap_begin;
    std::size_t word_idx = 0;
    while (word_idx < word_cnt
      && dtl::bits::pop_count(bitmap_begin[word_idx]) == 0) {
      ++word_idx;
    }
    if (word_idx < word_cnt) {
      return dtl::bits::tz_count(bitmap_begin[word_idx])
          + (word_idx * sizeof(word_type) * 8);
    }
    return word_idx * sizeof(word_type) * 8;
  }

  /// Find the last set bit.  Returns the index of the last set bit. If no bits
  /// are set, the length of the bitmap is returned.
  static std::size_t
  find_last(const word_type* bitmap_begin, const word_type* bitmap_end) {
    const std::size_t word_cnt = bitmap_end - bitmap_begin;
    std::size_t word_idx = word_cnt;
    while (word_idx > 0
      && dtl::bits::pop_count(bitmap_begin[word_idx - 1]) == 0) {
      --word_idx;
    }
    if (word_idx > 0) {
      return (word_idx * sizeof(word_type) * 8)
          - (dtl::bits::lz_count(bitmap_begin[word_idx - 1]) + 1);
    }
    return word_cnt * sizeof(word_type) * 8;
  }

  // Construction not allowed.
  bitmap_fun() = delete;

};
//===----------------------------------------------------------------------===//
} // namespace dtl
