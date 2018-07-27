#pragma once

#include <vector>

#include <dtl/dtl.hpp>
#include <dtl/math.hpp>

namespace dtl {

//===----------------------------------------------------------------------===//
struct rank1 {

  using size_type = $u32;
  std::vector<size_type> lut;
//  static constexpr u64 lut_entry_span = 4096; // bits
  static constexpr u64 lut_entry_span = 512; // bits

  std::vector<$u64> bitmap_;

  // TODO remove redundancy (for testing purposes only)
  void
  copy(std::vector<$u1>& bitmap) {
    std::vector<$u64> clone;
    clone.resize((bitmap.size() + 63) / 64);
    for (std::size_t i = 0; i < bitmap.size(); i++) {
      if (bitmap[i]) {
        u64 word_idx = i / 64;
        u64 bit_idx = i % 64;
        clone[word_idx] |= u64(1) << bit_idx;
      }
    }
    std::swap(bitmap_, clone);
  }

  // TODO remove redundancy (for testing purposes only)
  void
  copy(boost::dynamic_bitset<$u32>& bitmap) {
    std::vector<$u64> clone;
    clone.resize((bitmap.size() + 63) / 64);
    for (std::size_t i = 0; i < bitmap.size(); i++) {
      if (bitmap[i]) {
        u64 word_idx = i / 64;
        u64 bit_idx = i % 64;
        clone[word_idx] |= u64(1) << bit_idx;
      }
    }
    std::swap(bitmap_, clone);
  }

  // TODO remove the first LuT entry (to save 4 byte)
  void
  init(std::vector<$u1>& bitmap) {
    copy(bitmap);
    u64 bitmap_size = bitmap.size();
    u64 lut_entry_cnt = (bitmap_size + (lut_entry_span - 1)) / lut_entry_span;
    lut.resize(lut_entry_cnt, 0);
    size_type bit_cntr = 0;
    $u64 lut_write_pos = 0;
    for ($u64 i = 0; i < bitmap_size; i++) {
      if (i % lut_entry_span == 0) {
        assert(lut_write_pos < lut_entry_cnt);
        lut[lut_write_pos] = bit_cntr;
        lut_write_pos++;
      }
      bit_cntr += bitmap[i];
    }
  }

  void
  init(boost::dynamic_bitset<$u32>& bitmap) {
    copy(bitmap);
    u64 bitmap_size = bitmap.size();
    u64 lut_entry_cnt = (bitmap_size + (lut_entry_span - 1)) / lut_entry_span;
    lut.resize(lut_entry_cnt, 0);
    size_type bit_cntr = 0;
    $u64 lut_write_pos = 0;
    for ($u64 i = 0; i < bitmap_size; i++) {
      if (i % lut_entry_span == 0) {
        assert(lut_write_pos < lut_entry_cnt);
        lut[lut_write_pos] = bit_cntr;
        lut_write_pos++;
      }
      bit_cntr += bitmap[i];
    }
  }

  __forceinline__
  size_type
  operator()(u64 idx) const {
    u64 lut_idx = idx >> dtl::ct::log_2<lut_entry_span>::value;
    u64 bit_idx_from = lut_idx << dtl::ct::log_2<lut_entry_span>::value;
    u64 bit_idx_to = idx;
    size_type bit_cntr = lut[lut_idx];

    u64* words_begin = &bitmap_[bit_idx_from / 64];
    u64* words_end = words_begin + ((bit_idx_to - bit_idx_from) / 64);
    for ($u64 const* word_iter = words_begin; word_iter != words_end; word_iter++) {
      bit_cntr += dtl::bits::pop_count(*word_iter);
    }
    u64 rem = bit_idx_to % 64;
    if (rem != 0) {
      u64 w = *words_end << (64 - rem);
      bit_cntr += dtl::bits::pop_count(w);
    }
    return bit_cntr;
  }

  u64 size_in_bytes() const {
    return lut.size() * sizeof(size_type) // lut size
           + 8;   // pointer to bitmap
  }

};
//===----------------------------------------------------------------------===//


} // namespace dtl

