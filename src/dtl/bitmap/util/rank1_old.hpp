#pragma once

#include <vector>

#include <dtl/dtl.hpp>
#include <dtl/math.hpp>

namespace dtl {

//===----------------------------------------------------------------------===//
struct rank1 {

  using size_type = $u32;
  std::vector<size_type> lut;
//  static constexpr u64 lut_entry_span = 4096 * 16; // bits
//  static constexpr u64 lut_entry_span = 4096; // bits
  static constexpr u64 lut_entry_span = 512; // bits
//  static constexpr u64 lut_entry_span = 64; // bits

  std::vector<$u64> bitmap_;

  ~rank1() = default;

  // TODO remove redundancy (for testing purposes only)
  void
  copy(std::vector<$u1>& bitmap) {
    bitmap_.clear();
    bitmap_.resize((bitmap.size() + 63) / 64);
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
    bitmap_.resize((bitmap.size() + 63) / 64);
    for (std::size_t i = 0; i < bitmap.size(); i++) {
      if (bitmap[i]) {
        u64 word_idx = i / 64;
        u64 bit_idx = i % 64;
        assert(word_idx < bitmap_.size());
        bitmap_[word_idx] |= u64(1) << bit_idx;
      }
    }
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
    $u64 lut_entry_cnt = (bitmap_size + (lut_entry_span - 1)) / lut_entry_span;
    if (lut_entry_cnt == 0) lut_entry_cnt = 1;
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

//  __forceinline__
//  size_type
//  operator()(u64 idx) const {
//    u64 lut_idx = idx >> dtl::ct::log_2<lut_entry_span>::value;
//    u64 bit_idx_from = lut_idx << dtl::ct::log_2<lut_entry_span>::value;
//    u64 bit_idx_to = idx;
//    size_type bit_cntr = 0;
//
//    u64 word_idx_begin = bit_idx_from / 64;
//    u64 word_idx_end = word_idx_begin + ((bit_idx_to - bit_idx_from) / 64);
//    u64 to_cnt = word_idx_end - word_idx_begin;
//    switch (to_cnt) {
//      case 7: bit_cntr += dtl::bits::pop_count(bitmap_[word_idx_begin + 7]);
//      case 6: bit_cntr += dtl::bits::pop_count(bitmap_[word_idx_begin + 6]);
//      case 5: bit_cntr += dtl::bits::pop_count(bitmap_[word_idx_begin + 5]);
//      case 4: bit_cntr += dtl::bits::pop_count(bitmap_[word_idx_begin + 4]);
//      case 3: bit_cntr += dtl::bits::pop_count(bitmap_[word_idx_begin + 3]);
//      case 2: bit_cntr += dtl::bits::pop_count(bitmap_[word_idx_begin + 2]);
//      case 1: bit_cntr += dtl::bits::pop_count(bitmap_[word_idx_begin + 1]);
//    }
//    u64 rem = bit_idx_to % 64;
//    if (rem != 0) {
//      u64 w = bitmap_[word_idx_end] << (64 - rem);
//      bit_cntr += dtl::bits::pop_count(w);
//    }
//    return bit_cntr + lut[lut_idx];
//  }

//  __forceinline__
//  size_type
//  operator()(u64 idx) const {
//    u64 lut_idx = idx >> dtl::ct::log_2<lut_entry_span>::value;
//    u64 bit_idx_from = lut_idx << dtl::ct::log_2<lut_entry_span>::value;
//    u64 bit_idx_to = idx;
//    size_type bit_cntr = 0;
//
//    u64 word_idx_begin = bit_idx_from / 64;
//    u64 word_idx_end = word_idx_begin + ((bit_idx_to - bit_idx_from) / 64);
//    u64 to_cnt = word_idx_end - word_idx_begin;
//    size_type bit_cntr_0 = 0;
//    size_type bit_cntr_1 = 0;
//    size_type bit_cntr_2 = 0;
//    size_type bit_cntr_3 = 0;
//    size_type bit_cntr_4 = 0;
//    size_type bit_cntr_5 = 0;
//    size_type bit_cntr_6 = 0;
//    switch (to_cnt) {
//      case 7: bit_cntr_0 += dtl::bits::pop_count(bitmap_[word_idx_begin + 7]);
//      case 6: bit_cntr_1 += dtl::bits::pop_count(bitmap_[word_idx_begin + 6]);
//      case 5: bit_cntr_2 += dtl::bits::pop_count(bitmap_[word_idx_begin + 5]);
//      case 4: bit_cntr_3 += dtl::bits::pop_count(bitmap_[word_idx_begin + 4]);
//      case 3: bit_cntr_4 += dtl::bits::pop_count(bitmap_[word_idx_begin + 3]);
//      case 2: bit_cntr_5 += dtl::bits::pop_count(bitmap_[word_idx_begin + 2]);
//      case 1: bit_cntr_6 += dtl::bits::pop_count(bitmap_[word_idx_begin + 1]);
//    }
//    bit_cntr_0 += bit_cntr_1;
//    bit_cntr_2 += bit_cntr_3;
//    bit_cntr_4 += bit_cntr_5;
//    bit_cntr_0 += bit_cntr_2;
//    bit_cntr_4 += bit_cntr_6;
//    bit_cntr += bit_cntr_0 + bit_cntr_4;
//    u64 rem = bit_idx_to % 64;
//    if (rem != 0) {
//      u64 w = bitmap_[word_idx_end] << (64 - rem);
//      bit_cntr += dtl::bits::pop_count(w);
//    }
//    return bit_cntr + lut[lut_idx];
//  }

  __forceinline__
  size_type //__attribute__ ((noinline))
  operator()(u64 idx) const {
    u64 lut_idx = idx >> dtl::ct::log_2<lut_entry_span>::value;
    u64 bit_idx_from = lut_idx << dtl::ct::log_2<lut_entry_span>::value;
    u64 bit_idx_to = idx;
    size_type bit_cntr = 0;

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
    return bit_cntr + lut[lut_idx];
  }

  u64
  size_in_bytes() const {
    return lut.size() * sizeof(size_type) // lut size
           + 8;   // pointer to bitmap
  }

};
//===----------------------------------------------------------------------===//


} // namespace dtl

