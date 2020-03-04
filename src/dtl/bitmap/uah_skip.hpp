#pragma once
//===----------------------------------------------------------------------===//
#include "uah.hpp"

#include <dtl/bitmap/util/plain_bitmap.hpp>
#include <dtl/bitmap/util/plain_bitmap_iter.hpp>
#include <dtl/dtl.hpp>

#include <boost/dynamic_bitset.hpp>

#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <vector>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// Un-Aligned Hybrid: An RLE compressed representation of a bitmap of length N.
/// Unlike WAH or BBC, the encoding is not word or byte aligned.
/// This implementation maintains a small index that allows for faster skips.
template<typename _word_type = u32, std::size_t _skip_distance = 1024>
class uah_skip : public uah<_word_type> {
  using super = uah<_word_type>;
  using word_type = typename super::word_type;

  using uah<_word_type>::uah;
  using uah<_word_type>::is_fill_word;
  using uah<_word_type>::extract_fill_value;
  using uah<_word_type>::extract_fill_length;
  using uah<_word_type>::word_bitlength;
  using uah<_word_type>::payload_bit_cnt;

  std::vector<$u32> offsets_;

  /// Initialize the skip offsets.
  void
  init_skip_offsets() {
    const auto word_cnt = this->data_.size();
    offsets_.reserve((word_cnt + (_skip_distance - 1)) / _skip_distance);

    std::size_t i = 0;
    if (word_cnt > 1) {
      for (std::size_t word_idx = 0; word_idx < word_cnt - 1; ++word_idx) {
        if (word_idx % _skip_distance == 0) {
          offsets_.push_back(i);
        }
        const auto w = this->data_[word_idx];
        i += super::is_fill_word(w)
            ? super::extract_fill_length(w)
            : super::payload_bit_cnt;
      }
    }

    if (offsets_.empty()) {
      offsets_.push_back(0);
    }
  }

public:
  uah_skip() = default;

  explicit uah_skip(const boost::dynamic_bitset<$u32>& in) : super(in), offsets_() {
    init_skip_offsets();
  }

  ~uah_skip() = default;
  uah_skip(const uah_skip& other) = default;
  uah_skip(uah_skip&& other) noexcept = default;
  uah_skip& operator=(const uah_skip& other) = default;
  uah_skip& operator=(uah_skip&& other) noexcept = default;

  /// Return the size in bytes.
  std::size_t __forceinline__
  size_in_bytes() const {
    return super::size_in_bytes()
        + offsets_.size() * sizeof(u32);
  }

  static std::string
  name() {
    return "uah_skip" + std::to_string(super::word_bitlength);
  }

  /// Returns the value of the bit at the position pos.
  u1 __forceinline__
  test(const std::size_t pos) const {
    std::size_t word_idx = 0;
    std::size_t i = 0;

    auto search = std::upper_bound(offsets_.begin(), offsets_.end(), $u32(pos));
    const std::size_t offset_idx = std::distance(offsets_.begin(), search) - 1;

    word_idx = offset_idx * _skip_distance;
    i = offsets_[offset_idx];

    // Find the corresponding word.
    for (; word_idx < this->data_.size(); ++word_idx) {
      auto& w = this->data_[word_idx];
      if (this->is_literal_word(w)) {
        if (pos >= i + super::payload_bit_cnt) {
          i += super::payload_bit_cnt;
          continue;
        }
        else {
          return dtl::bits::bit_test(w, pos - i + 1); // TODO optimize
        }
      }
      else {
        auto fill_len = this->extract_fill_length(w);
        if (pos >= i + fill_len) {
          i += fill_len;
          continue;
        }
        else {
          return this->extract_fill_value(w);
        }
      }
    }
    return false;
  }

  /// Try to reduce the memory consumption. This function is supposed to be
  /// called after the bitmap has been modified.
  __forceinline__ void
  shrink() {
    super::shrink();
    offsets_.shrink_to_fit();
  }

  //===--------------------------------------------------------------------===//
  /// 1-run iterator
  class iter {
    const uah_skip& outer_;

    std::size_t word_idx_;
    std::size_t in_word_idx_;

    //===------------------------------------------------------------------===//
    // Iterator state
    //===------------------------------------------------------------------===//
    /// Points to the beginning of a 1-run.
    $u64 pos_;
    /// The length of the current 1-run.
    $u64 length_;
    //===------------------------------------------------------------------===//

  public:
    explicit __forceinline__
    iter(const uah_skip& outer)
        : outer_(outer),
          word_idx_(0),
          in_word_idx_(0),
          pos_(0), length_(0) {
      const auto word_cnt = outer_.data_.size();
      word_type w;
      while (word_idx_ < word_cnt) {
        w = outer_.data_[word_idx_];
        if (is_fill_word(w)) {
          if (extract_fill_value(w) == false) {
            pos_ += extract_fill_length(w);
          }
          else {
            length_ = extract_fill_length(w);
            in_word_idx_ = 0;
            break;
          }
        }
        else {
          const word_type payload = w >> 1;
          if (payload == 0) {
            pos_ += payload_bit_cnt;
          }
          else {
            const std::size_t b = dtl::bits::tz_count(payload);
            std::size_t e = b + 1;
            for (; e < payload_bit_cnt; ++e) {
              u1 is_set = dtl::bits::bit_test(payload, e);
              if (!is_set) break;
            }
            pos_ += b;
            length_ = e - b;
            in_word_idx_ = e;
            break;
          }
        }
        ++word_idx_;
      }
      if (word_idx_ == word_cnt) {
        pos_ = outer_.encoded_bitmap_length_;
        length_ = 0;
      }
      else {
        if (is_fill_word(w)) {
          ++word_idx_;
          in_word_idx_ = 0;
        }
      }
    }

    /// Forward the iterator to the next 1-run.
    void __forceinline__
    next() {
      pos_ += length_;
      length_ = 0;
      const auto word_cnt = outer_.data_.size();
      word_type w;
      while (word_idx_ < word_cnt) {
        w = outer_.data_[word_idx_];
        if (is_fill_word(w)) {
          if (extract_fill_value(w) == false) {
            pos_ += extract_fill_length(w);
          }
          else {
            length_ = extract_fill_length(w);
            break;
          }
        }
        else {
          if (in_word_idx_ < payload_bit_cnt) { // TODO decode the entire literal word at once.
            const word_type payload = w >> (1 + in_word_idx_);
            if (payload == 0) {
              pos_ += payload_bit_cnt - in_word_idx_;
            }
            else {
              const std::size_t b = dtl::bits::tz_count(payload);
              std::size_t e = b + 1;
              for (; e < (payload_bit_cnt - in_word_idx_); ++e) {
                u1 is_set = dtl::bits::bit_test(payload, e);
                if (!is_set) break;
              }
              pos_ += b;
              length_ = e - b;
              in_word_idx_ += e;
              break;
            }
          }
        }
        ++word_idx_;
        in_word_idx_ = 0;
      }
      if (word_idx_ == word_cnt) {
        pos_ = outer_.encoded_bitmap_length_;
        length_ = 0;
      }
      else {
        if (is_fill_word(w)) {
          ++word_idx_;
          in_word_idx_ = 0;
        }
      }
    }

    /// Forward the iterator to the desired position.
    void __forceinline__
    skip_to(const std::size_t to_pos) {
      assert(pos_ <= to_pos);
      if (to_pos >= outer_.encoded_bitmap_length_) {
        pos_ = outer_.encoded_bitmap_length_;
        length_ = 0;
        return;
      }

      // Skip close to the desired position.
      auto search = std::upper_bound(
          outer_.offsets_.begin(), outer_.offsets_.end(), $u32(to_pos));
      const std::size_t offset_idx =
          std::distance(outer_.offsets_.begin(), search) - 1;

      word_idx_ = offset_idx * _skip_distance;
      in_word_idx_ = 0;
      pos_ = outer_.offsets_[offset_idx];
      length_ = 0;

      // Call next until the desired position has been reached.
      while (!end() && pos() + length() <= to_pos) {
        next();
      }

      // Adjust the current position and run length.
      if (!end() && pos() < to_pos) {
        length_ -= to_pos - pos_;
        pos_ = to_pos;
      }
    }

    u1 __forceinline__
    end() const noexcept {
      return pos_ >= outer_.encoded_bitmap_length_;
    }

    u64 __forceinline__
    pos() const noexcept {
      return pos_;
    }

    u64 __forceinline__
    length() const noexcept {
      return length_;
    }
  };
  //===--------------------------------------------------------------------===//

  using skip_iter_type = iter;
  using scan_iter_type = iter;

  /// Returns a 1-run iterator.
  skip_iter_type __forceinline__
  it() const {
    return skip_iter_type(*this);
  }

  /// Returns a 1-run iterator.
  scan_iter_type __forceinline__
  scan_it() const {
    return scan_iter_type(*this);
  }

  /// Returns the name of the instance including the most important parameters
  /// in JSON.
  std::string
  info() const {
    return "{\"name\":\"" + name() + "\""
        + ",\"n\":" + std::to_string(this->size())
        + ",\"size\":" + std::to_string(size_in_bytes())
        + ",\"word_size\":" + std::to_string(sizeof(_word_type))
        + ",\"skip_distance\":" + std::to_string(_skip_distance)
        + ",\"skip_offsets_size\":" + std::to_string(offsets_.size() * sizeof($u32))
        + "}";
  }

  // For debugging purposes.
  void
  print(std::ostream& os) const {
    super::print(os);
    os << "     idx | offset " << std::endl;
    for (std::size_t i = 0; i < offsets_.size(); ++i) {
      os << std::setw(8) << i << " | ";
      os << std::setw(8) << offsets_[i];
      os << std::endl;
    }
  }
};
//===----------------------------------------------------------------------===//
/// UAH compressed representation of a bitmap of length N using 8-bit words.
using uah8_skip = uah_skip<u8>;
/// UAH compressed representation of a bitmap of length N using 16-bit words.
using uah16_skip = uah_skip<u16>;
/// UAH compressed representation of a bitmap of length N using 32-bit words.
using uah32_skip = uah_skip<u32>;
/// UAH compressed representation of a bitmap of length N using 64-bit words.
using uah64_skip = uah_skip<u64>;
//===----------------------------------------------------------------------===//
} // namespace dtl
