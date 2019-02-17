#pragma once

#include <dtl/dtl.hpp>
#include <dtl/bits.hpp>

namespace dtl {

template<typename _word_type>
struct bitmap_fun {

  using word_type = _word_type;
  static constexpr std::size_t word_bitlength = sizeof(word_type) * 8;

  static u1 __forceinline__
  test(const word_type* b, std::size_t i) noexcept {
    const auto block_idx = i / word_bitlength;
    const auto bit_idx = i % word_bitlength;
    return dtl::bits::bit_test(b[block_idx], bit_idx);
  }

  static void __forceinline__
  set(const word_type* b, std::size_t i) noexcept {
    const auto block_idx = i / word_bitlength;
    const auto bit_idx = i % word_bitlength;
    b[block_idx] ^= word_type(1) << bit_idx;
  }

  /// Fetch consecutive bits.
  static u64 __forceinline__
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

};

} // namespace dtl
