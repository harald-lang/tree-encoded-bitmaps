#pragma once
//===----------------------------------------------------------------------===//
#include "rank1_logic_surf.hpp" // popcount_linear()

#include <dtl/dtl.hpp>
#include <dtl/math.hpp>

#include <boost/dynamic_bitset.hpp>

#include <string>
#include <type_traits>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// Rank1 implementation which does not require additional space but has linear
/// time complexity.
template<
    /// The word type used to store bitmaps.
    typename _word_type = $u64,
    /// Inclusive [b,e] vs exclusive [b,e)
    u1 _inclusive = false>
struct rank1_logic_linear {
  using word_type = typename std::remove_cv<_word_type>::type;
  using size_type = $u32;
  static constexpr u64 is_inclusive = _inclusive ? 1 : 0;

  // Pure static.
  rank1_logic_linear() = delete;
  ~rank1_logic_linear() = delete;

  // Does nothing. Exists for compatibility reasons.
  static void
  init_inplace(
      const word_type* const bitmap_begin,
      const word_type* const bitmap_end,
      size_type* lut) {
  }

  // Always returns 0. Exists for compatibility reasons.
  static constexpr u64
  lut_entry_cnt(u64 bitmap_size) {
    return 0;
  }

  // Always returns 0. Exists for compatibility reasons.
  static constexpr u64
  estimate_size_in_bytes(u64 bitmap_size) {
    return lut_entry_cnt(bitmap_size) * sizeof(size_type);
  }

  /// Computes the rank1 of the bit at position 'idx'.
  static size_type __forceinline__
  get(const size_type* lut, u64 idx, const word_type* bitmap_ptr) noexcept {
    return rank1_logic_surf<_word_type, _inclusive>
        ::popcount_linear(bitmap_ptr, 0, idx + is_inclusive);
  }

  /// Returns the name of the implementation.
  static std::string
  name() noexcept {
    return std::string("\"Linear\"");
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
