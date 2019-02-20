#pragma once

#include <dtl/dtl.hpp>
#include <dtl/iterator.hpp>
#include <dtl/bitmap/util/bitmap_fun.hpp>

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
} // namespace dtl