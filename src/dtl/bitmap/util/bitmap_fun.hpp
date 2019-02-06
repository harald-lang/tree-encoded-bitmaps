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
    dtl::bits::bit_test(b[block_idx], bit_idx);
  }

  static void __forceinline__
  set(const word_type* b, std::size_t i) noexcept {
    const auto block_idx = i / word_bitlength;
    const auto bit_idx = i % word_bitlength;
    b[block_idx] ^= word_type(1) << bit_idx;
  }

};

} // namespace dtl
