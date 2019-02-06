#pragma once

#include <vector>

#include <dtl/dtl.hpp>
#include <dtl/math.hpp>

namespace dtl {

//===----------------------------------------------------------------------===//
struct rank1 {

  // Lookup table for the pre-computed rank values on block level.
  using size_type = $u32;
  std::vector<size_type> lut;
  static constexpr u64 block_bitlength = 512;

  // A copy of the bitmap // TODO remove the copy
  using word_t = $u64;
  std::vector<word_t> bitmap_;
  static constexpr u64 word_bitlength = sizeof(word_t) * 8;
  static constexpr u64 words_per_block = block_bitlength / word_bitlength;

  ~rank1() = default;

  void
  init(std::vector<$u1>& bitmap) {
    copy(bitmap);
    init();
  }

  void
  init(boost::dynamic_bitset<$u32>& bitmap) {
    copy(bitmap);
    init();
  }

  /// Counts the 1-bits within the range [bit_idx_from, bit_idx_to). The range
  /// must start at a word boundary.
  static size_type __forceinline__
  popcount(const word_t* bitmap, u64 bit_idx_from, u64 bit_idx_to) {
    assert(bit_idx_from % word_bitlength == 0);
    size_type bit_cntr = 0;
    const word_t* words_begin = &bitmap[bit_idx_from / word_bitlength];
    const word_t* words_end = words_begin +
        ((bit_idx_to - bit_idx_from) / word_bitlength);
    for (word_t const* wi = words_begin; wi != words_end; wi++) {
      bit_cntr += dtl::bits::pop_count(*wi);
    }
    u64 rem = bit_idx_to % word_bitlength;
    if (rem != 0) {
      word_t w = *words_end << (word_bitlength - rem);
      bit_cntr += dtl::bits::pop_count(w);
    }
    return bit_cntr;
  }

  /// Counts the 1-bits within a block.
  static size_type __forceinline__
  popcount_block(const word_t* block) {
    size_type bit_cntr = 0;
    for (std::size_t i = 0; i < words_per_block; ++i) {
      bit_cntr += dtl::bits::pop_count(block[i]);
    }
    return bit_cntr;
  }

  void
  init() {
    u64 bitmap_bitlength = bitmap_.size() * word_bitlength;
    u64 block_cnt = (bitmap_bitlength + block_bitlength - 1) / block_bitlength;
    u64 lut_entry_cnt = block_cnt + 1;
    lut.resize(lut_entry_cnt, 0);

    size_type bit_cntr = 0;
    for ($u64 i = 0; i < block_cnt; ++i) {
      lut[i] = bit_cntr;
      u64 bit_idx_from = i * block_bitlength;
      u64 bit_idx_to = std::min(bit_idx_from + block_bitlength, bitmap_bitlength);
      bit_cntr += popcount(bitmap_.data(), bit_idx_from, bit_idx_to);
    }
    lut[lut_entry_cnt - 1] = bit_cntr;
  }

  size_type __forceinline__
  operator()(u64 idx) const {
    u64 lut_idx = idx >> dtl::ct::log_2<block_bitlength>::value;
    u64 bit_idx_from = lut_idx << dtl::ct::log_2<block_bitlength>::value;
    u64 bit_idx_to = idx;
    return lut[lut_idx] + popcount(bitmap_.data(), bit_idx_from, bit_idx_to);
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

  // TODO remove redundancy (for testing purposes only)
  void
  copy(std::vector<$u1>& bitmap) {
    bitmap_.clear();
    bitmap_.resize((bitmap.size() + word_bitlength - 1) / word_bitlength, 0);
    for (std::size_t i = 0; i < bitmap.size(); i++) {
      if (bitmap[i]) {
        u64 word_idx = i / 64;
        u64 bit_idx = i % 64;
        assert(word_idx < bitmap_.size());
        bitmap_[word_idx] |= u64(1) << bit_idx;
      }
    }
  }

  // TODO remove redundancy (for testing purposes only)
  void
  copy(boost::dynamic_bitset<$u32>& bitmap) {
    bitmap_.clear();
    bitmap_.resize((bitmap.size() + word_bitlength - 1) / word_bitlength, 0);
    for (std::size_t i = 0; i < bitmap.size(); ++i) {
      if (bitmap[i]) {
        u64 word_idx = i / 64;
        u64 bit_idx = i % 64;
        assert(word_idx < bitmap_.size());
        bitmap_[word_idx] |= u64(1) << bit_idx;
      }
    }
  }

  /// Returns the important properties in JSON.
  std::string
  info() const noexcept {
    return "{\"size\":" + std::to_string(size_in_bytes())
        + ",\"block_size\":" + std::to_string(block_bitlength / 8)
        + "}";
  }

};
//===----------------------------------------------------------------------===//


} // namespace dtl

