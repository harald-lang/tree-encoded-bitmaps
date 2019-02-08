#pragma once

#include <vector>

#include <dtl/dtl.hpp>
#include <dtl/math.hpp>
#include <boost/dynamic_bitset.hpp>

namespace dtl {

//===----------------------------------------------------------------------===//
// Taken and adapted from
//  https://github.com/efficient/SuRF/blob/master/include/popcount.h
//===----------------------------------------------------------------------===//

template<typename _word_type = $u64, u1 _inclusive = false>
struct rank1_super_fast {

  using word_type = typename std::remove_cv<_word_type>::type;

  // Lookup table for the pre-computed rank values on block level.
  using size_type = $u32;
  std::vector<size_type> lut;

  static constexpr u64 block_bitlength = 64;
  static constexpr u64 word_bitlength = sizeof(word_type) * 8;
  static constexpr u64 words_per_block = block_bitlength / word_bitlength;

  static constexpr u64 is_inclusive = _inclusive ? 1 : 0;
  ~rank1_super_fast() = default;

#define popcountsize 64ULL
#define popcountmask (popcountsize - 1)

  uint64_t __forceinline__
  popcount_last_word(const uint64_t bits, uint64_t nbits) const {
    if (nbits == 0) { return 0; }
    assert(nbits <= 64);
    const uint64_t bits_shifted = bits << (63 - ((nbits - 1)));
    return dtl::bits::pop_count(bits_shifted);
  }

  void
  init(const boost::dynamic_bitset<word_type>& bitmap) {
    u64 bitmap_bitlength = bitmap.size() * word_bitlength;
    u64 block_cnt = (bitmap_bitlength + block_bitlength - 1) / block_bitlength;
    u64 lut_entry_cnt = block_cnt + 1;
    lut.resize(lut_entry_cnt, 0);

    size_type bit_cntr = 0;
    for ($u64 i = 0; i < block_cnt; ++i) {
      lut[i] = bit_cntr;
      bit_cntr += popcount_last_word(bitmap.m_bits.data()[i], block_bitlength);
    }
    lut[lut_entry_cnt - 1] = bit_cntr;
  }

  size_type __forceinline__
  operator()(u64 idx, const word_type* bitmap_ptr) const {
    const auto block_id = idx / block_bitlength;
    const auto offset = idx & (block_bitlength - 1);
    return lut[block_id]
        + popcount_last_word(bitmap_ptr[block_id], offset + is_inclusive);
  }

  u64
  size_in_bytes() const {
    return lut.size() * sizeof(size_type); // lut size
  }

  static u64
  estimate_size_in_bytes(u64 bitmap_size) {
    u64 bitmap_bitlength = bitmap_size;
    u64 block_cnt = (bitmap_bitlength + block_bitlength - 1) / block_bitlength;
    u64 lut_entry_cnt = block_cnt + 1;
    return lut_entry_cnt * sizeof(size_type);
  }

  /// Returns the important properties in JSON.
  std::string
  info() const noexcept {
    return "{\"name\":" + std::string("\"super_fast\"")
        + ",\"size\":" + std::to_string(size_in_bytes())
        + ",\"block_size\":" + std::to_string(block_bitlength / 8)
        + "}";
  }

};
//===----------------------------------------------------------------------===//
} // namespace dtl
