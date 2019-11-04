#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/dtl.hpp>
#include <dtl/math.hpp>

#include <string>
#include <type_traits>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// Rank implementation that uses a precomputed lookup table (aka a dictionary).
/// Adapted from https://github.com/efficient/SuRF/blob/master/include/popcount.h
template<
    /// The word type used to store bitmaps.
    typename _word_type = $u64,
    /// Inclusive [b,e] vs exclusive [b,e)
    u1 _inclusive = false,
    /// The granularity of the rank lookup table.
    u64 _block_bitlength = 512>
struct rank1_logic_surf {
  using word_type = typename std::remove_cv<_word_type>::type;
  using size_type = $u32;

  static constexpr u64 block_bitlength = _block_bitlength;
  static constexpr u64 word_bitlength = sizeof(word_type) * 8;
  static constexpr u64 words_per_block = block_bitlength / word_bitlength;
  static constexpr u64 is_inclusive = _inclusive ? 1 : 0;

  static_assert(dtl::is_power_of_two(block_bitlength),
      "The block size must be a power of two.");
  static_assert((block_bitlength % word_bitlength) == 0,
      "The block size must be a multiple of word size.");

  // Pure static.
  rank1_logic_surf() = delete;
  ~rank1_logic_surf() = delete;

  /// Computes the popcount within the range [x, x + nbits).
  static size_type __forceinline__
  popcount_linear(const uint64_t* bits, uint64_t x, uint64_t nbits) noexcept {
    if (nbits == 0) { return 0; }
    assert(bits != nullptr);

    static constexpr u64 popcount_size = 64ull;
    static constexpr u64 popcount_mask = (popcount_size - 1);

    uint64_t last_word_idx = (nbits - 1) / popcount_size;
    uint64_t p = 0;

    for (uint64_t i = 0; i < last_word_idx; i++) {
      p += dtl::bits::pop_count(bits[x + i]);
    }

    // 'nbits' may or may not fall on a multiple of 64 boundary,
    // so we may need to zero out the right side of the last word
    // (accomplished by shifting it right, since we're just popcounting)
    uint64_t lastshifted =
        bits[x + last_word_idx] << (63 - ((nbits - 1) & popcount_mask));
    p += dtl::bits::pop_count(lastshifted);
    return p;
  }

public:
  /// Initializes the rank LuT in place. The function estimate_size_in_bytes()
  /// allows to predetermine the required memory.
  static void
  init_inplace(
      const word_type* const bitmap_begin,
      const word_type* const bitmap_end,
      size_type* lut) noexcept {
    u64 bitmap_word_cnt = bitmap_end - bitmap_begin;
    u64 bitmap_bitlength = bitmap_word_cnt * word_bitlength;
    u64 block_cnt = (bitmap_bitlength + block_bitlength - 1) / block_bitlength;
    u64 lut_entry_cnt = block_cnt + 1;

    size_type bit_cntr = 0;
    for ($u64 i = 0; i < block_cnt; ++i) {
      lut[i] = bit_cntr;
      const auto word_cnt_in_current_block =
          (i + 1) * words_per_block <= bitmap_word_cnt
          ? words_per_block
          : bitmap_word_cnt % words_per_block;
      const auto nbits = word_cnt_in_current_block * word_bitlength;
      bit_cntr += popcount_linear(
          bitmap_begin, i * words_per_block, nbits);
    }
    lut[lut_entry_cnt - 1] = bit_cntr;
  }

  /// Returns the number of LuT entries for a bitmap of the given size.
  static constexpr u64 __forceinline__
  lut_entry_cnt(u64 bitmap_size) noexcept {
    u64 bitmap_bitlength = bitmap_size;
    u64 block_cnt = (bitmap_bitlength + block_bitlength - 1) / block_bitlength;
    u64 lut_entry_cnt = block_cnt + 1;
    return lut_entry_cnt;
  }

  /// Returns the size of the rank LuT in bytes for a bitmap of the given size.
  static constexpr u64 __forceinline__
  estimate_size_in_bytes(u64 bitmap_size) noexcept {
    return lut_entry_cnt(bitmap_size) * sizeof(size_type);
  }

  /// Computes the rank1 of the bit at position 'idx'.
  static size_type __forceinline__
  get(const size_type* lut, u64 idx, const word_type* bitmap_ptr) noexcept {
    const auto block_id = idx / block_bitlength;
    const auto offset = idx & (block_bitlength - 1);
    return (lut[block_id]
        + popcount_linear(bitmap_ptr, block_id * words_per_block,
            offset + is_inclusive));
  }

  /// Returns the name of the implementation.
  static std::string
  name() noexcept {
    return std::string("SuRF");
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
