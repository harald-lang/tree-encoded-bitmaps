#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/dtl.hpp>

#include <cassert>
#include <cstddef>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// 1-run iterator, with skip support.
template<typename _bitmap_type>
class plain_bitmap_iter {
  const _bitmap_type& bitmap_;

  /// Points to the beginning of a 1-run.
  $u64 pos_;
  /// The length of the current 1-run.
  $u64 length_;

public:
  /// C'tor
  explicit __forceinline__
  plain_bitmap_iter(const _bitmap_type& bitmap)
      : bitmap_(bitmap), pos_(bitmap_.find_first()), length_(0) {
    // determine the length of the current 1fill
    if (pos_ < bitmap_.size()) {
      length_ = 1;
      while (pos_ + length_ < bitmap_.size()
          && bitmap_.test(pos_ + length_)) {
        length_++;
      }
    }
    else {
      pos_ = bitmap_.size();
      length_ = 0;
    }
  }

  /// Forwards the iterator to the next 1-run.
  void __forceinline__
  next() {
    pos_ += length_;
    length_ = 0;
    pos_ = bitmap_.find_next(pos_);
    // determine the length of the current 1fill
    if (pos_ < bitmap_.size()) {
      length_ = 1;
      while (pos_ + length_ < bitmap_.size()
          && bitmap_.test(pos_ + length_)) {
        length_++;
      }
    }
    else {
      pos_ = bitmap_.size();
      length_ = 0;
    }
  }

  /// Forwards the iterator to the given position. If the given position is a
  /// 0 bit, the iterator is forwarded to the next 1-run.
  void __forceinline__
  skip_to(const std::size_t to_pos) {
    if (to_pos >= bitmap_.size()) {
      pos_ = bitmap_.size();
      length_ = 0;
      return;
    }
    pos_ = to_pos;
    length_ = 0;
    if (!bitmap_.test(pos_)) {
      pos_ = bitmap_.find_next(pos_);
    }
    if (pos_ < bitmap_.size()) {
      length_ = 1;
      while (pos_ + length_ < bitmap_.size() && bitmap_.test(pos_ + length_)) {
        length_++;
      }
    }
    else {
      pos_ = bitmap_.size();
      length_ = 0;
    }
  }

  /// Returns true when the iterator reached the end of the bitmap.
  u1
  end() const noexcept {
    return length_ == 0;
  }

  /// Returns the start index of the current 1-run.
  u64
  pos() const noexcept {
    return pos_;
  }

  /// Returns the length of the current 1-run.
  u64
  length() const noexcept {
    return length_;
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
