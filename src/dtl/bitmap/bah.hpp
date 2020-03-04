#pragma once
//===----------------------------------------------------------------------===//
#include "bah_patterns.hpp"
#include "bah_types.hpp"

#include <dtl/bitmap/util/plain_bitmap.hpp>
#include <dtl/dtl.hpp>

#include <boost/dynamic_bitset.hpp>

#include <cstddef>
#include <iomanip>
#include <type_traits>
#include <vector>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// Implementation of the Byte Aligned Hybrid bitmap compression technique as
/// described in the paper 'BAH: A Bitmap Index Compression Algorithm for Fast
/// Data Retrieval' by Li et al.
class bah {
  static constexpr std::size_t word_bitlength = sizeof(bah_word_t) * 8;

  /// The encoded bitmap.
  std::vector<bah_byte_t> main_;
  std::vector<bah_word_t> data_;
  std::vector<bah_byte_t> index_;
  std::vector<bah_word_t> counter_;
  /// The length of the encoded bitmap.
  std::size_t encoded_bitmap_length_ = 0;

  /// The different word kinds.
  enum class word_kind {
    ZERO,
    ENCODABLE,
    LITERAL
  };

  /// Determines the word kind for a single word of the plain bitmap.
  static word_kind
  get_word_kind(bah_word_t word) {
    if (word == 0) {
      return word_kind::ZERO;
    }
    if (get_oep_code(word) != dtl::bah_code_not_found) {
      return word_kind::ENCODABLE;
    }
    if (get_tep_code(word) != dtl::bah_code_not_found) {
      return word_kind::ENCODABLE;
    }
    return word_kind::LITERAL;
  }

public:
  bah() = default;

  explicit bah(const boost::dynamic_bitset<$u32>& in) {
    auto* bitmap_begin = in.m_bits.data();
    auto* bitmap_end = in.m_bits.data() + in.m_bits.size();
    const auto word_cnt = bitmap_end - bitmap_begin; // TODO handle the case when the bitmap size is not a multiple of the word size

    // The plain bitmap is segmented into words of the same kind.
    // Keep track of the first and last word of the current segment.
    std::size_t seg_begin = 0;
    std::size_t seg_end = 0;

    // Encode the bitmap segment by segment.
    while (seg_begin < word_cnt) {
      auto& w = bitmap_begin[seg_begin];
      auto seg_kind = get_word_kind(w);
      seg_end = seg_begin + 1;
      // Find the end of the current 'segment'.
      while (seg_end < word_cnt
          && get_word_kind(bitmap_begin[seg_end]) == seg_kind) {
        auto seg_len = seg_end - seg_begin;
        if ((seg_kind == word_kind::LITERAL && seg_len == 63)
            || (seg_kind == word_kind::ENCODABLE)) {
          break;
        }
        ++seg_end;
      }
      // Encode the segment.
      if (seg_kind == word_kind::ZERO) {
        auto len = seg_end - seg_begin;
        if (len > 252) {
          main_.push_back(bah_byte_t(0));
          counter_.push_back(len);
        }
        else {
          while (len > 0) {
            auto l = std::min(len, 63ul);
            main_.push_back(static_cast<bah_byte_t>(l));
            len -= l;
          }
        }
      }
      else if (seg_kind == word_kind::LITERAL) {
        auto len = seg_end - seg_begin;
        main_.push_back(static_cast<bah_byte_t>(len) | (bah_byte_t(0b01) << 6));
        for (std::size_t i = seg_begin; i < seg_end; ++i) {
          data_.push_back(bitmap_begin[i]);
        }
      }
      else { /* ENCODABLE */
        auto code = get_oep_code(w);
        if (code != dtl::bah_code_not_found) {
          main_.push_back(
              static_cast<bah_byte_t>(code) | (bah_byte_t(0b10) << 6));
        }
        else {
          code = get_tep_code(w);
          bah_byte_t code_lo = code & 0b111111;
          bah_byte_t code_hi = code >> 6;
          main_.push_back(code_lo | (bah_byte_t(0b11) << 6));
          index_.push_back(code_hi);
        }
      }

      // Next.
      seg_begin = seg_end;
    }
    encoded_bitmap_length_ = in.size();

    // Try to reduce the memory consumption.
    shrink();
  }

  ~bah() = default;
  bah(const bah& other) = default;
  bah(bah&& other) noexcept = default;
  bah& operator=(const bah& other) = default;
  bah& operator=(bah&& other) noexcept = default;

  /// Return the size in bytes.
  std::size_t __forceinline__
  size_in_bytes() const {
    return main_.size() * sizeof(bah_byte_t)
        + data_.size() * sizeof(bah_word_t)
        + index_.size() * sizeof(bah_byte_t)
        + counter_.size() * sizeof(bah_word_t)
        + 3 * sizeof(void*) /* pointers/offsets to the three other arrays */
        + sizeof(encoded_bitmap_length_); /* bit-length of the original bitmap */
  }

  /// Returns the size of the bitmap.
  std::size_t __forceinline__
  size() const {
    return encoded_bitmap_length_;
  }

  /// Conversion to an plain bitmap. // TODO remove
  dtl::plain_bitmap<$u64>
  to_plain_bitmap() {
    dtl::plain_bitmap<$u64> ret(encoded_bitmap_length_, false);
    // The write position.
    std::size_t i = 0;

    // The read position(s).
    std::size_t data_idx = 0;
    std::size_t index_idx = 0;
    std::size_t counter_idx = 0;
    for (std::size_t main_idx = 0; main_idx < main_.size(); ++main_idx) {
      const auto& w = main_[main_idx];
      const auto type = w >> 6;
      const auto n = w & 0b111111;
      switch (type) {
        case 0b00: { /* ZERO */
          std::size_t count = 0;
          if (n == 0) {
            count = counter_[counter_idx++];
          }
          else {
            count = n;
          }
          auto l = count * word_bitlength;
          ret.clear(i, std::min(i + l, encoded_bitmap_length_));
          i += l;
          break;
        }
        case 0b01: { /* LITERAL */
          for (std::size_t j = 0; j < n; ++j) {
            auto lit = data_[data_idx++];
            ret.store_bits(i,
                std::min(i + word_bitlength, encoded_bitmap_length_), lit);
            i += word_bitlength;
          }
          break;
        }
        case 0b10: { /* ENCODABLE (one byte) */
          auto pattern = get_oep(n);
          ret.store_bits(i,
              std::min(i + word_bitlength, encoded_bitmap_length_), pattern);
          i += word_bitlength;
          break;
        }
        case 0b11: { /* ENCODABLE (two bytes) */
          auto code = bah_word_t(n) | (bah_word_t(index_[index_idx++]) << 6);
          auto pattern = get_tep(code);
          ret.store_bits(i,
              std::min(i + word_bitlength, encoded_bitmap_length_), pattern);
          i += word_bitlength;
          break;
        }
      }
    }
    if (i < encoded_bitmap_length_) {
      ret.clear(i, encoded_bitmap_length_);
    }
    return std::move(ret);
  }

  static std::string
  name() {
    return "bah";
  }

  /// Returns the value of the bit at the position pos.
  u1 __forceinline__
  test(const std::size_t pos) const {
    // The current range spanned by the current encoded word.
    std::size_t b = 0;
    std::size_t e = 0;

    // The read position(s).
    std::size_t data_idx = 0;
    std::size_t index_idx = 0;
    std::size_t counter_idx = 0;
    for (std::size_t main_idx = 0; main_idx < main_.size(); ++main_idx) {
      const auto& w = main_[main_idx];
      const auto type = w >> 6;
      const auto n = w & 0b111111;
      switch (type) {
        case 0b00: { /* ZERO */
          std::size_t count = 0;
          if (n == 0) {
            count = counter_[counter_idx++];
          }
          else {
            count = n;
          }
          auto l = count * word_bitlength;
          e = std::min(b + l, encoded_bitmap_length_);
          if (b <= pos && pos < e) {
            return false;
          }
          break;
        }
        case 0b01: { /* LITERAL */
          e = std::min(b + n * word_bitlength, encoded_bitmap_length_);
          if (b <= pos && pos < e) {
            auto offset = pos - b;
            auto lit_word = data_[data_idx + (offset / word_bitlength)];
            return dtl::bits::bit_test(lit_word, offset % word_bitlength);
          }
          data_idx += n;
          break;
        }
        case 0b10: { /* ENCODABLE (one byte) */
          e = std::min(b + word_bitlength, encoded_bitmap_length_);
          if (b <= pos && pos < e) {
            auto pattern = get_oep(n);
            return dtl::bits::bit_test(pattern, pos % word_bitlength);
          }
          break;
        }
        case 0b11: { /* ENCODABLE (two bytes) */
          e = std::min(b + word_bitlength, encoded_bitmap_length_);
          auto code = bah_word_t(n) | (bah_word_t(index_[index_idx++]) << 6);
          if (b <= pos && pos < e) {
            auto pattern = get_tep(code);
            return dtl::bits::bit_test(pattern, pos % word_bitlength);
          }
          break;
        }
      }
      // Next.
      b = e;
    }
    return false;
  }

  /// Try to reduce the memory consumption. This function is supposed to be
  /// called after the bitmap has been modified.
  __forceinline__ void
  shrink() {
    main_.shrink_to_fit();
    data_.shrink_to_fit();
    index_.shrink_to_fit();
    counter_.shrink_to_fit();
  }

  //===--------------------------------------------------------------------===//
  /// 1-run iterator
  class iter {
    const bah& outer_;

    /// The read positions.
    std::size_t main_idx_;
    std::size_t data_idx_;
    std::size_t index_idx_;
    std::size_t counter_idx_;

    //    word_kind segment_kind_;
    bah_word_t literal_;
    std::size_t literal_words_remaining_;
    std::size_t literal_in_word_idx_;

    /// Points to the beginning of a 1-run.
    $u64 pos_;
    /// The length of the current 1-run.
    $u64 length_;

  public:
    explicit __forceinline__
    iter(const bah& outer)
        : outer_(outer),
          main_idx_(0),
          data_idx_(0),
          index_idx_(0),
          counter_idx_(0),
          literal_(0),
          literal_words_remaining_(0),
          literal_in_word_idx_(0),
          pos_(0), length_(0) {
      // Find the first 1-run.
      const auto main_cnt = outer_.main_.size();
      while (main_idx_ < main_cnt) {
        pos_ += length_;
        length_ = 0;

        const auto& w = outer_.main_[main_idx_++];
        const auto type = w >> 6;
        const auto n = w & 0b111111;

        switch (type) {
          case 0b00: { /* ZERO */
            std::size_t count = 0;
            if (n == 0) {
              count = outer_.counter_[counter_idx_++];
            }
            else {
              count = n;
            }
            auto l = count * word_bitlength;
            length_ = l;
            break; // Skip over zeros.
          }
          case 0b01: { /* LITERAL */
            literal_ = outer_.data_[data_idx_++];
            assert(literal_ != 0);
            literal_words_remaining_ = n - 1;
            assert(literal_words_remaining_ < n);
            const std::size_t b = dtl::bits::tz_count(literal_);
            std::size_t e = b + 1;
            for (; e < word_bitlength; ++e) {
              u1 is_set = dtl::bits::bit_test(literal_, e);
              if (!is_set) break;
            }
            literal_in_word_idx_ = e;
            pos_ += b;
            length_ = e - b;
            return;
          }
          case 0b10: { /* ENCODABLE (one byte) */
            literal_ = get_oep(n);
            assert(literal_ != 0);
            literal_words_remaining_ = 0;
            const std::size_t b = dtl::bits::tz_count(literal_);
            std::size_t e = b + 1;
            for (; e < word_bitlength; ++e) {
              u1 is_set = dtl::bits::bit_test(literal_, e);
              if (!is_set) break;
            }
            literal_in_word_idx_ = e;
            pos_ += b;
            length_ = e - b;
            return;
          }
          case 0b11: { /* ENCODABLE (two bytes) */
            auto code =
                bah_word_t(n) | (bah_word_t(outer_.index_[index_idx_++]) << 6);
            literal_ = get_tep(code);
            assert(literal_ != 0);
            literal_words_remaining_ = 0;
            const std::size_t b = dtl::bits::tz_count(literal_);
            std::size_t e = b + 1;
            for (; e < word_bitlength; ++e) {
              u1 is_set = dtl::bits::bit_test(literal_, e);
              if (!is_set) break;
            }
            literal_in_word_idx_ = e;
            pos_ += b;
            length_ = e - b;
            return;
          }
        }
      }
      pos_ = outer_.encoded_bitmap_length_;
      length_ = 0;
    }

    /// Forward the iterator to the next 1-run.
    void __forceinline__
    next() {
      const auto main_cnt = outer_.main_.size();
      while (main_idx_ < main_cnt || literal_ != 0) {
        pos_ += length_;
        length_ = 0;

        // Check, whether we are currently in a literal segment.
        if (literal_ != 0) {
          if (literal_in_word_idx_ < word_bitlength) {
            // Find the next set bit.
            std::size_t b = literal_in_word_idx_;
            for (; b < word_bitlength; b++) { // TODO optimize
              u1 is_set = dtl::bits::bit_test(literal_, b);
              if (is_set) break;
            }
            if (b < word_bitlength) {
              // Determine the length of the current 1-run and return.
              std::size_t e = b + 1;
              for (; e < word_bitlength; ++e) { // TODO optimize
                u1 is_set = dtl::bits::bit_test(literal_, e);
                if (!is_set) break;
              }
              pos_ += b - literal_in_word_idx_;
              length_ = e - b;
              literal_in_word_idx_ = e;
              return;
            }
            pos_ += b - literal_in_word_idx_;
          }
          // Fetch the next literal word, if any.
          if (literal_words_remaining_ > 0) {
            literal_ = outer_.data_[data_idx_++];
            --literal_words_remaining_;
            literal_in_word_idx_ = 0;
            continue;
          }
          else {
            // Reached the end of the literal segment.
            literal_ = 0;
            literal_in_word_idx_ = 0;
          }
        }

        assert(literal_ == 0);

        if (main_idx_ >= main_cnt) break;

        const auto& w = outer_.main_[main_idx_++];
        const auto type = w >> 6;
        const auto n = w & 0b111111;
        switch (type) {
          case 0b00: { /* ZERO */
            std::size_t count = 0;
            if (n == 0) {
              count = outer_.counter_[counter_idx_++];
            }
            else {
              count = n;
            }
            auto l = count * word_bitlength;
            length_ = l;
            break; // Skip over zeros.
          }
          case 0b01: { /* LITERAL */
            literal_ = outer_.data_[data_idx_++];
            assert(literal_ != 0);
            literal_words_remaining_ = n - 1;
            assert(literal_words_remaining_ < n);
            const std::size_t b = dtl::bits::tz_count(literal_);
            std::size_t e = b + 1;
            for (; e < word_bitlength; ++e) {
              u1 is_set = dtl::bits::bit_test(literal_, e);
              if (!is_set) break;
            }
            literal_in_word_idx_ = e;
            pos_ += b;
            length_ = e - b;
            return;
          }
          case 0b10: { /* ENCODABLE (one byte) */
            literal_ = get_oep(n);
            assert(literal_ != 0);
            literal_words_remaining_ = 0;
            const std::size_t b = dtl::bits::tz_count(literal_);
            std::size_t e = b + 1;
            for (; e < word_bitlength; ++e) { // TODO optimize
              u1 is_set = dtl::bits::bit_test(literal_, e);
              if (!is_set) break;
            }
            literal_in_word_idx_ = e;
            pos_ += b;
            length_ = e - b;
            return;
          }
          case 0b11: { /* ENCODABLE (two bytes) */
            auto code =
                bah_word_t(n) | (bah_word_t(outer_.index_[index_idx_++]) << 6);
            literal_ = get_tep(code);
            assert(literal_ != 0);
            literal_words_remaining_ = 0;
            const std::size_t b = dtl::bits::tz_count(literal_);
            std::size_t e = b + 1;
            for (; e < word_bitlength; ++e) { // TODO optimize
              u1 is_set = dtl::bits::bit_test(literal_, e);
              if (!is_set) break;
            }
            literal_in_word_idx_ = e;
            pos_ += b;
            length_ = e - b;
            return;
          }
        }
      }
      pos_ = outer_.encoded_bitmap_length_;
      length_ = 0;
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
      return length_ == 0;
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
        + ",\"n\":" + std::to_string(encoded_bitmap_length_)
        + ",\"size\":" + std::to_string(size_in_bytes())
        + ",\"word_size\":" + std::to_string(sizeof(bah_word_t))
        + "}";
  }

  // For debugging purposes.
  void
  print(std::ostream& os) const {
    os << "main idx | word type | content" << std::endl;
    os << "---------|-----------|---------------------------------" << std::endl;
    // The read position(s).
    std::size_t data_idx = 0;
    std::size_t index_idx = 0;
    std::size_t counter_idx = 0;
    for (std::size_t main_idx = 0; main_idx < main_.size(); ++main_idx) {
      os << std::setw(8) << main_idx << " | ";
      auto& w = main_[main_idx];
      const auto type = w >> 6;
      const auto n = w & 0b111111;
      switch (type) {
        case 0b00: { /* ZERO */
          std::size_t count = 0;
          if (n == 0) {
            count = counter_[counter_idx++];
          }
          else {
            count = n;
          }
          os << count << " x ZERO";
          break;
        }
        case 0b01: { /* LITERAL */
          os << n << " x LITERAL";
          for (std::size_t i = 0; i < n; ++i) {
            os << std::endl;
            os << "         |           | "
               << std::bitset<word_bitlength>(data_[data_idx + i]);
          }
          data_idx += n;
          break;
        }
        case 0b10: { /* ENCODABLE (one byte) */
          os << "1"
             << " x ENCODABLE (one byte)";
          os << std::endl;
          os << "         |           | "
              << std::bitset<word_bitlength>(get_oep(n));
          break;
        }
        case 0b11: { /* ENCODABLE (two bytes) */
          os << "1"
             << " x ENCODABLE (two bytes)";
          auto code = bah_word_t(n) | (bah_word_t(index_[index_idx++]) << 6);
          os << std::endl;
          os << "         |           | "
              << std::bitset<word_bitlength>(get_tep(code));
          break;
        }
      }
      os << std::endl;
    }
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
