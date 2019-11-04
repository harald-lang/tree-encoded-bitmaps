#pragma once
//===----------------------------------------------------------------------===//
#include "rank1_logic_surf.hpp"

#include <dtl/dtl.hpp>
#include <dtl/math.hpp>

#include <boost/dynamic_bitset.hpp>

#include <string>
#include <type_traits>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// Very fast rank1 implementation. The LuT has 64-bit granularity which
/// results in 50% additional space consumption.
template<
    /// The word type used to store bitmaps.
    typename _word_type = $u64,
    /// Inclusive [b,e] vs exclusive [b,e)
    u1 _inclusive = false>
struct rank1_logic_word_blocked {
  using word_type = typename std::remove_cv<_word_type>::type;
  using size_type = $u32;
  static constexpr u64 is_inclusive = _inclusive ? 1 : 0;

  static constexpr u64 word_bitlength = sizeof(word_type) * 8;
  static constexpr u64 block_bitlength = 64;
  static_assert(word_bitlength == block_bitlength,
      "A word and a block must be equal in size."); // TODO remove limitation

  using rank_helper =
      dtl::rank1_logic_surf<_word_type, _inclusive, block_bitlength>;

  // Pure static.
  rank1_logic_word_blocked() = delete;
  ~rank1_logic_word_blocked() = delete;

  /// Initializes the rank LuT in place. The function estimate_size_in_bytes()
  /// allows to predetermine the required memory.
  static void
  init_inplace(
      const word_type* const bitmap_begin,
      const word_type* const bitmap_end,
      size_type* lut) noexcept {
    rank_helper::init_inplace(bitmap_begin, bitmap_end, lut);
  }

  /// Returns the number of LuT entries for a bitmap of the given size.
  static constexpr u64
  lut_entry_cnt(u64 bitmap_size) noexcept {
    return rank_helper::lut_entry_cnt(bitmap_size);
  }

  /// Returns the size of the rank LuT in bytes for a bitmap of the given size.
  static constexpr u64
  estimate_size_in_bytes(u64 bitmap_size) noexcept {
    return lut_entry_cnt(bitmap_size) * sizeof(size_type);
  }

  /// Computes the rank1 of the bit at position 'idx'.
  static size_type __forceinline__
  get(const size_type* lut, u64 idx, const word_type* bitmap_ptr) noexcept {
    const word_type mask63 = ~word_type(0) >> 1;

    const auto bit_idx = idx + is_inclusive;
    const auto block_idx = bit_idx / block_bitlength;
    const auto lut_entry = lut[block_idx];
    const auto bitmap_word = bitmap_ptr[block_idx];

    const auto offset = bit_idx & (block_bitlength - 1);
    // The following works, because we never have to popcount a full word.
    // In cases where 'offset' is 0, we only need the value from the LuT.
    const auto bitmap_word_masked = (bitmap_word << (63 - offset)) & mask63;
    const auto popcnt = dtl::bits::pop_count(bitmap_word_masked);
    return lut_entry + popcnt;
  }

  /// Returns the name of the implementation.
  static std::string
  name() noexcept {
    return std::string("WordBlocked");
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
