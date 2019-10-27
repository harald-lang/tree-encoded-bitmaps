#pragma once
//===----------------------------------------------------------------------===//
#include "bitmap_fun.hpp"

#include <dtl/dtl.hpp>

#include <cassert>
#include <cstddef>
#include <memory>
#include <type_traits>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// A buffered bitmap reader.
template<typename _word_type>
class bitmap_seq_reader {
  using word_type = typename std::remove_cv<_word_type>::type;
  static_assert(std::is_integral<word_type>::value,
      "The word type must be an integral type.");
  static_assert(!std::is_same<word_type, bool>::value,
      "The word type must not be a boolean.");

  static constexpr auto word_bitlength = sizeof(_word_type) * 8;
  using fn = dtl::bitmap_fun<word_type>;

  const word_type* bitmap_;
  std::size_t read_idx_;

  word_type buffer_ = 0;
  std::size_t buffer_read_idx_ = word_bitlength;

public:
  inline bitmap_seq_reader(const word_type* bitmap, std::size_t start_idx)
      : bitmap_(bitmap), read_idx_(start_idx) {
  }

  /// Returns the next bit.
  inline u1
  next() {
    if (buffer_read_idx_ == word_bitlength) {
      fill_buffer();
    }
    assert(buffer_read_idx_ < word_bitlength);
    u1 ret_val = dtl::bits::bit_test(buffer_, buffer_read_idx_);
    ++buffer_read_idx_;
    return ret_val;
  }

private:
  inline void
  fill_buffer() {
    const auto word_idx = read_idx_ / word_bitlength;
    buffer_ = bitmap_[word_idx];
    const auto offset = read_idx_ % word_bitlength;
    buffer_read_idx_ = offset;
    read_idx_ += word_bitlength - offset;
  }
};
//===----------------------------------------------------------------------===//
/// A buffered bitmap reader. When the reader exceeds the end of the bitmap
/// the next() function returns 0.
template<typename _word_type>
class bitmap_limit_seq_reader {
  using word_type = typename std::remove_cv<_word_type>::type;
  static_assert(std::is_integral<word_type>::value,
      "The word type must be an integral type.");
  static_assert(!std::is_same<word_type, bool>::value,
      "The word type must not be a boolean.");

  static constexpr auto word_bitlength = sizeof(_word_type) * 8;
  using fn = dtl::bitmap_fun<word_type>;

  const word_type* bitmap_;
  std::size_t read_idx_;
  const std::size_t read_idx_end_;

  word_type buffer_ = 0;
  std::size_t buffer_read_idx_ = word_bitlength;

public:
  inline bitmap_limit_seq_reader(const word_type* bitmap,
      std::size_t start_idx, std::size_t end_idx)
      : bitmap_(bitmap), read_idx_(start_idx), read_idx_end_(end_idx) {
    assert(start_idx <= end_idx);
  }

  /// Returns the next bit.
  inline u1
  next() {
    if (buffer_read_idx_ == word_bitlength) {
      fill_buffer();
    }
    assert(buffer_read_idx_ < word_bitlength);
    u1 ret_val = dtl::bits::bit_test(buffer_, buffer_read_idx_);
    ++buffer_read_idx_;
    return ret_val;
  }

private:
  inline void
  fill_buffer() {
    if (read_idx_ >= read_idx_end_) {
      buffer_ = 0;
      buffer_read_idx_ = 0;
      read_idx_ += word_bitlength;
      return;
    }
    const auto word_idx = read_idx_ / word_bitlength;
    buffer_ = bitmap_[word_idx];
    const auto offset = read_idx_ % word_bitlength;
    const auto last_word_idx = (read_idx_end_ - 1) / word_bitlength;
    buffer_read_idx_ = offset;
    read_idx_ += word_bitlength - offset;

    if (word_idx == last_word_idx) {
      const auto offset_end = read_idx_end_ % word_bitlength;
      if (offset_end > 0) {
        const auto mask = ~word_type(0) >> (word_bitlength - offset_end);
        buffer_ &= mask;
      }
    }
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
