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
/// A buffered bitmap writer.
template<typename _word_type>
class bitmap_writer {
  using word_type = typename std::remove_cv<_word_type>::type;
  static_assert(std::is_integral<word_type>::value,
      "The word type must be an integral type.");
  static_assert(!std::is_same<word_type, bool>::value,
      "The word type must not be a boolean.");

  static constexpr auto word_bitlength = sizeof(_word_type) * 8;
  using fn = dtl::bitmap_fun<word_type>;

  word_type* bitmap_;
  std::size_t write_idx_;

  word_type buffer_ = 0;
  std::size_t buffer_write_idx_ = 0;

public:
  inline bitmap_writer(word_type* bitmap, std::size_t start_idx) noexcept
      : bitmap_(bitmap), write_idx_(start_idx) {
  }

  /// Write the 'bit_cnt' low order bits from the given word to the bitmap.
  inline void
  write(word_type word, std::size_t bit_cnt) noexcept {
    if (bit_cnt == 0) return;
    assert(bit_cnt <= word_bitlength);
    //    assert(write_idx_ + bit_cnt <= bitmap_.size());
    const std::size_t buf_free = word_bitlength - buffer_write_idx_;
    if (bit_cnt < buf_free) {
      buffer_ |= (word << buffer_write_idx_);
      buffer_write_idx_ += bit_cnt;
    }
    else {
      flush();
      buffer_ = word;
      buffer_write_idx_ = bit_cnt;
    }
  }

  /// Flush the buffer.
  inline void
  flush() noexcept {
    if (buffer_write_idx_ == 0) return;
    fn::store_bits(bitmap_, write_idx_, write_idx_ + buffer_write_idx_, buffer_);
#ifndef NDEBUG
    // Validation code.
    for (std::size_t i = 0; i < buffer_write_idx_; ++i) {
      u1 bit_to_write = dtl::bits::bit_test(buffer_, i);
      assert(fn::test(bitmap_, write_idx_ + i) == bit_to_write);
    }
#endif
    write_idx_ += buffer_write_idx_;
    buffer_ = 0;
    buffer_write_idx_ = 0;
  }
};
//===----------------------------------------------------------------------===//
/// A buffered bitmap writer. When the writer exceeds the end of the bitmap
/// further calls to 'write' have no effect.
template<typename _word_type>
class bitmap_limit_writer {
  using word_type = typename std::remove_cv<_word_type>::type;
  static_assert(std::is_integral<word_type>::value,
      "The word type must be an integral type.");
  static_assert(!std::is_same<word_type, bool>::value,
      "The word type must not be a boolean.");

  static constexpr auto word_bitlength = sizeof(_word_type) * 8;
  using fn = dtl::bitmap_fun<word_type>;

  word_type* bitmap_;
  std::size_t write_idx_;
  std::size_t limit_;

  word_type buffer_ = 0;
  std::size_t buffer_write_idx_ = 0;

  inline u1
  exceeded() noexcept {
    return write_idx_ == limit_;
  }

public:
  inline bitmap_limit_writer(word_type* bitmap,
      std::size_t start_idx, std::size_t limit_idx) noexcept
      : bitmap_(bitmap),
        write_idx_(start_idx),
        limit_(limit_idx) {
    assert(start_idx <= limit_idx);
  }

  /// Write the 'bit_cnt' low order bits from the given word to the bitmap.
  inline void
  write(word_type word, std::size_t bit_cnt) noexcept {
    if (bit_cnt == 0) return;
    if (exceeded()) return;
    assert(bit_cnt <= word_bitlength);
    //    assert(write_idx_ + bit_cnt <= bitmap_.size());
    const std::size_t buf_free = word_bitlength - buffer_write_idx_;
    if (bit_cnt < buf_free) {
      buffer_ |= (word << buffer_write_idx_);
      buffer_write_idx_ += bit_cnt;
    }
    else {
      flush();
      buffer_ = word;
      buffer_write_idx_ = bit_cnt;
    }
  }

  /// Flush the buffer.
  inline void
  flush() noexcept {
    std::size_t end_idx = std::min(limit_, write_idx_ + buffer_write_idx_);
    if (write_idx_ != end_idx) {
      fn::store_bits(bitmap_, write_idx_, end_idx, buffer_);
#ifndef NDEBUG
      // Validation code.
      for (std::size_t i = 0; i < buffer_write_idx_; ++i) {
        if (write_idx_ + i == end_idx) break;
        u1 bit_to_write = dtl::bits::bit_test(buffer_, i);
        assert(fn::test(bitmap_, write_idx_ + i) == bit_to_write);
      }
#endif
    }
    if (end_idx != limit_) {
      write_idx_ += buffer_write_idx_;
      buffer_ = 0;
      buffer_write_idx_ = 0;
    }
    else {
      write_idx_ = limit_;
      buffer_ = 0;
      buffer_write_idx_ = 0;
    }
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
