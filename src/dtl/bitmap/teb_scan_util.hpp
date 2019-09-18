#pragma once
//===----------------------------------------------------------------------===//
#include "util/bitmap_fun.hpp"

#include <dtl/dtl.hpp>
#include <dtl/iterator.hpp>

#include <cassert>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
static u64 __forceinline__
fetch_bits(const dtl::data_view<u64>& bitmap,
           i64 bit_idx_begin,
           i64 bit_idx_end /* non-inclusive */) {
  using word_type = $u64;
  static constexpr u64 word_bitlength = sizeof(word_type) * 8;

  assert(bit_idx_end > bit_idx_begin);
  const auto cnt = bit_idx_end - bit_idx_begin;
  assert(cnt > 0);
  assert(cnt <= word_bitlength);

  i64 bitmap_bitlength = bitmap.size() * word_bitlength;
  if (bitmap_bitlength == 0) {
    return (~word_type(0)) >> (word_bitlength - std::abs(bit_idx_begin));
  }

  if (bit_idx_end <= 0) {
    return (~word_type(0)) >> (word_bitlength - cnt);
  }
  else if (bit_idx_begin < 0) {
    const auto offset = std::abs(bit_idx_begin);
    word_type ret_val = (~word_type(0)) >> (word_bitlength - offset);
    ret_val |= dtl::bitmap_fun<$u64>::fetch_bits(bitmap.begin(), 0,
                                                 cnt - offset) << offset;
    return ret_val;
  }
  else {
    if (bit_idx_begin >= bitmap_bitlength) {
      return word_type(0);
    }
    return dtl::bitmap_fun<$u64>::fetch_bits(bitmap.begin(),
        bit_idx_begin, std::min(bit_idx_end, bitmap_bitlength));
  }
}
//===----------------------------------------------------------------------===//
static u64 __forceinline__
fetch_n_bits(const dtl::data_view<u64>& bitmap,
             i64 bit_idx_begin,
             i64 cnt) {
  return fetch_bits(bitmap, bit_idx_begin, bit_idx_begin + cnt);
}
//===----------------------------------------------------------------------===//
static u64 __forceinline__
fetch_bits(const dtl::data_view<u64>& bitmap,
    i64 bit_idx_begin,
    i64 bit_idx_end, /* non-inclusive */
    u1 lo_bit, u1 hi_bit) {
  using word_type = $u64;
  static constexpr u64 word_bitlength = sizeof(word_type) * 8;
  const auto lo_word = lo_bit ? ~word_type(0) : word_type(0);
  const auto hi_word = hi_bit ? ~word_type(0) : word_type(0);

  assert(bit_idx_end > bit_idx_begin);
  const auto cnt = bit_idx_end - bit_idx_begin;
  assert(cnt > 0);
  assert(cnt <= word_bitlength);

  i64 bitmap_bitlength = bitmap.size() * word_bitlength;
  if (bitmap_bitlength == 0) {
    return (lo_word) >> (word_bitlength - std::abs(bit_idx_begin));
  }

  if (bit_idx_end <= 0) {
    return (lo_word) >> (word_bitlength - cnt);
  }
  else if (bit_idx_begin < 0) {
    const auto offset = std::abs(bit_idx_begin);
    word_type ret_val = (lo_word) >> (word_bitlength - offset);
    ret_val |= dtl::bitmap_fun<$u64>::fetch_bits(bitmap.begin(), 0,
        cnt - offset) << offset;
    return ret_val;
  }
  else {
    if (bit_idx_begin >= bitmap_bitlength) {
      return hi_word;
    }
    return dtl::bitmap_fun<$u64>::fetch_bits(bitmap.begin(),
        bit_idx_begin, std::min(bit_idx_end, bitmap_bitlength));
  }
}
//===----------------------------------------------------------------------===//
static u64 __forceinline__
fetch_n_bits(const dtl::data_view<u64>& bitmap,
    i64 bit_idx_begin,
    i64 cnt,
    u1 lo_bit, u1 hi_bit) {
  return fetch_bits(bitmap, bit_idx_begin, bit_idx_begin + cnt, lo_bit, hi_bit);
}
//===----------------------------------------------------------------------===//
} // namespace dtl