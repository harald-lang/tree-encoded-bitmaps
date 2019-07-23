#pragma once

#if !defined(__teb_inline__)
#if defined(NDEBUG)
// Release build.
#define __teb_inline__ inline __attribute__((always_inline))
//#define __teb_inline__ __attribute__((noinline))
#else
#define __teb_inline__
#endif
#endif


#include <dtl/dtl.hpp>
#include <dtl/bits.hpp>

namespace dtl {

template<typename _word_type>
struct bitmap_fun {

  using word_type = _word_type;
  static constexpr std::size_t word_bitlength = sizeof(word_type) * 8;

  static u1 __teb_inline__
  test(const word_type* b, std::size_t i) noexcept {
    const auto block_idx = i / word_bitlength;
    const auto word = b[block_idx];
    const auto bit_idx = i % word_bitlength;
    return (word & (word_type(1) << i)) != 0;
  }

  static void __teb_inline__
  set(const word_type* b, std::size_t i) noexcept {
    const auto block_idx = i / word_bitlength;
    const auto bit_idx = i % word_bitlength;
    b[block_idx] ^= word_type(1) << bit_idx;
  }

  /// Fetch consecutive bits.
  static u64 __teb_inline__
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

  /// Find the first set bit.
  /// Returns the index of the first set bit. If no bits are set, the length
  /// of the bitmap is returned.
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

  /// Find the last set bit.
  /// Returns the index of the last set bit. If no bits are set, the length
  /// of the bitmap is returned.
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

};

} // namespace dtl
