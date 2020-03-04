#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/bitmap/util/plain_bitmap.hpp>
#include <dtl/bitmap/util/plain_bitmap_iter.hpp>
#include <dtl/dtl.hpp>

#include <boost/dynamic_bitset.hpp>

#include <cstddef>
#include <type_traits>
#include <vector>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// Byte-aligned Bitmap Code (BBC)
///
/// The implementation is based on the paper "Notes on Design and Implementation
/// of Compressed Bit Vectors" by Wu et al.
///
/// Note: The implementation is NOT yet optimized for performance.
class bbc {

protected:
  using byte_t = $u8;

  static constexpr byte_t all_zero = 0b00000000;
  static constexpr byte_t all_one = 0b11111111;

  /// The maximum number of fill bytes per header byte.
  static constexpr std::size_t max_fill_byte_cnt = 1u << (7 * 4);
  /// The maximum number of tail bytes per header byte.
  static constexpr std::size_t max_tail_byte_cnt = 15;


  /// The encoded bitmap.
  std::vector<byte_t> data_;
  /// The length of the encoded bitmap.
  std::size_t encoded_bitmap_length_ = 0;

  //===--------------------------------------------------------------------===//
  // Helper functions to work with the code words.
  //===--------------------------------------------------------------------===//
  /// The different kind of input bytes to distinguish.
  enum byte_kind {
    ALL_ZERO,
    ALL_ONE,
    ALMOST_ALL_ZERO,
    ALMOST_ALL_ONE,
    MIXED
  };

  /// Determine the kind the byte at index i.
  static inline constexpr byte_kind
  determine_byte_kind(const byte_t* byte_array, const std::size_t byte_array_len,
      const std::size_t i) noexcept {
    assert(i < byte_array_len);
    const byte_t b = byte_array[i];
    if (b == all_zero)
      return ALL_ZERO;
    if (b == all_one)
      return ALL_ONE;
    if (dtl::bits::pop_count(b) == 1) {
      // Could be a dirty 1-bit.
      if (i == 0) return MIXED;
      if (byte_array[i - 1] != all_zero) return MIXED;
      if (((i + 1) < byte_array_len)
        && !(byte_array[i - 1] == all_zero || byte_array[i - 1] == all_one)) return MIXED;
      // Is a dirty 1-bit.
      return ALMOST_ALL_ZERO;
    }
    if (dtl::bits::pop_count(b) == 7) {
      // Could be a dirty 0-bit.
      if (i == 0) return MIXED;
      if (byte_array[i - 1] != all_one) return MIXED;
      if (((i + 1) < byte_array_len)
          && !(byte_array[i - 1] == all_zero || byte_array[i - 1] == all_one)) return MIXED;
      // Is a dirty 0-bit.
      return ALMOST_ALL_ONE;
    }
    return MIXED;
  }
  //===--------------------------------------------------------------------===//

  /// Compress the given input bitmap.
  void
  compress(const byte_t* bitmap_begin, const byte_t* bitmap_end) {
    const byte_t* byte_array = bitmap_begin;
    const std::size_t byte_cnt = bitmap_end - bitmap_begin;

    // The current range of input bytes that are inspected.
    std::size_t b = 0;
    std::size_t e = 0;

    while (e < byte_cnt) {
      // Inspect the input bitmap byte by byte.

      // Counter for the fill bytes.
      std::size_t fill_byte_cnt = 0;
      // The fill value.
      $u1 fill_val = 0;

      // Determine the fill value and length, if there is one.
      if (e < byte_cnt) {
        const auto kind = determine_byte_kind(byte_array, byte_cnt, e);
        if (kind == ALL_ZERO || kind == ALL_ONE) {
          ++e;
          fill_val = kind == ALL_ZERO ? 0 : 1;
          while ((e < byte_cnt)
              && ((e - b) < max_fill_byte_cnt)
              && (determine_byte_kind(byte_array, byte_cnt, e) == kind)) {
            ++e;
          }
          fill_byte_cnt = e - b;
          assert(fill_byte_cnt <= max_fill_byte_cnt);
        }
      }

      // Indicates whether a dirty bit was found in the tail.
      $u1 dirty_tail = false;
      // The position of the dirty bit.
      std::size_t dirty_bit_pos = 0;

      // Check, whether the run has a tail byte with a dirty bit.
      if (e < byte_cnt) {
        const auto kind = determine_byte_kind(byte_array, byte_cnt, e);
        if (kind == ALMOST_ALL_ZERO || kind == ALMOST_ALL_ONE) {
          // Found a literal with a dirty bit.
          dirty_tail = true;
          dirty_bit_pos = dtl::bits::tz_count(
              (kind == ALMOST_ALL_ZERO)
              ? u32(byte_array[e])
              : ~u32(byte_array[e]));
          ++e;
        }
      }

      // Counter for the literal tail bytes.
      std::size_t literal_byte_cnt = 0;
      // The index of the first literal byte in the input bitmap.
      std::size_t literal_bytes_begin = e;
      // The index + 1 of the last literal byte in the input bitmap.
      std::size_t literal_bytes_end = e;

      // Determine the number of tail bytes, if there are any.
      if (!dirty_tail && e < byte_cnt) {
        while ((e < byte_cnt)
            && ((e - b) < max_tail_byte_cnt)
            && (determine_byte_kind(byte_array, byte_cnt, e) == MIXED)) {
          ++e;
        }
        literal_bytes_end = e;
        literal_byte_cnt = literal_bytes_end - literal_bytes_begin;
        assert(literal_byte_cnt <= max_tail_byte_cnt);
      }

      // Write the header.
      byte_t hdr = 0;
      if (fill_byte_cnt < 4) {
        if (!dirty_tail) {
          // BBC case 1
          hdr = (1 << 7)
              | (byte_t(fill_val) << 6)
              | (fill_byte_cnt << 4)
              | literal_byte_cnt;
        }
        else {
          // BBC case 2
          hdr = (1 << 6)
              | (byte_t(fill_val) << 5)
              | (fill_byte_cnt << 3)
              | dirty_bit_pos;
        }
        data_.push_back(hdr);
      }
      else {
        if (!dirty_tail) {
          // BBC case 3
          hdr = (1 << 5)
              | (byte_t(fill_val) << 4)
              | literal_byte_cnt;
        }
        else {
          // BBC case 4
          hdr = (1 << 4)
              | (byte_t(fill_val) << 3)
              | dirty_bit_pos;
        }
        data_.push_back(hdr);
        // Write the variable-byte (VB) counter.
        const auto var_byte_cnt =
            (((sizeof(fill_byte_cnt) * 8) - dtl::bits::lz_count(fill_byte_cnt))
            + 6) / 7;
        for (std::size_t k = 0; k < var_byte_cnt; ++k) {
          const auto o = var_byte_cnt - k - 1;
          byte_t c = (fill_byte_cnt >> (o * 7)) & 0b01111111;
          if (o != 0) c |= 0b10000000;
          data_.push_back(c);
        }
      }

      // Write the literals.
      if (!dirty_tail) {
        for (std::size_t i = literal_bytes_begin; i < literal_bytes_end; ++i) {
          data_.push_back(byte_array[i]);
        }
      }

      b = e;
    }
  }


public:
  bbc() = default;

  explicit bbc(const boost::dynamic_bitset<$u32>& in) {
    // Reinterpret the input bitmap as a byte array.
    const byte_t* bitmap_begin =
        reinterpret_cast<const byte_t*>(in.m_bits.data());
    const byte_t* bitmap_end = bitmap_begin + ((in.size() + 7) / 8);
    // Compress the bitmap.
    compress(bitmap_begin, bitmap_end);
    encoded_bitmap_length_ = in.size();
    // Try to reduce the memory consumption.
    shrink();
  }

  ~bbc() = default;
  bbc(const bbc& other) = default;
  bbc(bbc&& other) noexcept = default;
  bbc& operator=(const bbc& other) = default;
  bbc& operator=(bbc&& other) noexcept = default;

  /// Return the size in bytes.
  std::size_t __forceinline__
  size_in_bytes() const {
    return data_.size() * sizeof(byte_t) /* size of the compressed bitmap */
        + sizeof(encoded_bitmap_length_); /* bit-length of the original bitmap */
  }

  /// Returns the size of the bitmap.
  std::size_t __forceinline__
  size() const {
    return encoded_bitmap_length_;
  }

  /// Conversion to an plain bitmap. // TODO remove
  dtl::plain_bitmap<$u64>
  to_plain_bitmap() const {
    dtl::plain_bitmap<$u64> bitmap(encoded_bitmap_length_);
    byte_t* byte_array = reinterpret_cast<byte_t*>(bitmap.data());
    std::size_t o = 0;
    std::size_t i = 0;
    if (data_.empty()) return bitmap;

    while (i < data_.size()) {
      const auto hdr = data_[i];
      const auto bbc_case = dtl::bits::lz_count(u32(hdr)) - 23;
      switch (bbc_case) {
        case 1: {
          u1 fill_val = (hdr & 0b01000000) != 0;
          const std::size_t fill_byte_cnt = (hdr & 0b00110000) >> 4;
          const std::size_t tail_byte_cnt = hdr & 0b00001111;

          // Write the fill.
          const byte_t fill_byte = fill_val ? all_one : all_zero;
          for (std::size_t k = 0; k < fill_byte_cnt; ++k) {
            byte_array[o + k] = fill_byte;
          }
          o += fill_byte_cnt;

          // Write the literals.
          for (std::size_t k = 0; k < tail_byte_cnt; ++k) {
            byte_array[o + k] = data_[i + 1 + k];
          }
          o += tail_byte_cnt;

          i += 1 + tail_byte_cnt;
          break;
        }

        case 2: {
          u1 fill_val = (hdr & 0b00100000) != 0;
          const std::size_t fill_byte_cnt = (hdr & 0b00011000) >> 3;
          const std::size_t dirty_bit_pos = hdr & 0b00000111;

          // Write the fill.
          const byte_t fill_byte = fill_val ? all_one : all_zero;
          for (std::size_t k = 0; k < fill_byte_cnt; ++k) {
            byte_array[o + k] = fill_byte;
          }
          o += fill_byte_cnt;

          // Write the tail.
          const byte_t dirty_tail_byte = fill_byte ^ (1 << dirty_bit_pos);
          byte_array[o] = dirty_tail_byte;
          o += 1;

          i += 1;
          break;
        }

        case 3: {
          u1 fill_val = (hdr & 0b00010000) != 0;
          // Read the variable-byte (VB) counter.
          std::size_t fill_byte_cnt = 0;
          std::size_t fill_byte_cntr_len = 1;
          for (; fill_byte_cntr_len <= 4; ++fill_byte_cntr_len) {
            const auto x = data_[i + fill_byte_cntr_len];
            fill_byte_cnt <<= 7;
            fill_byte_cnt |= (x & 0b01111111);
            if ((x & 0b10000000) == 0) break;
          }
          const std::size_t tail_byte_cnt = hdr & 0b00001111;

          // Write the fill.
          const byte_t fill_byte = fill_val ? all_one : all_zero;
          for (std::size_t k = 0; k < fill_byte_cnt; ++k) {
            byte_array[o + k] = fill_byte;
          }
          o += fill_byte_cnt;

          // Write the literals.
          for (std::size_t k = 0; k < tail_byte_cnt; ++k) {
            byte_array[o + k] = data_[i + 1 + k + fill_byte_cntr_len];
          }
          o += tail_byte_cnt;

          i += 1 + fill_byte_cntr_len + tail_byte_cnt;
          break;
        }

        case 4: {
          u1 fill_val = (hdr & 0b00001000) != 0;
          // Read the variable-byte (VB) counter.
          std::size_t fill_byte_cnt = 0;
          std::size_t fill_byte_cntr_len = 1;
          for (; fill_byte_cntr_len <= 4; ++fill_byte_cntr_len) {
            const auto x = data_[i + fill_byte_cntr_len];
            fill_byte_cnt <<= 7;
            fill_byte_cnt |= (x & 0b01111111);
            if ((x & 0b10000000) == 0) break;
          }
          const std::size_t dirty_bit_pos = hdr & 0b00000111;

          // Write the fill.
          const byte_t fill_byte = fill_val ? all_one : all_zero;
          for (std::size_t k = 0; k < fill_byte_cnt; ++k) {
            byte_array[o + k] = fill_byte;
          }
          o += fill_byte_cnt;

          // Write the tail.
          const byte_t dirty_tail_byte = fill_byte ^ (1 << dirty_bit_pos);
          byte_array[o] = dirty_tail_byte;
          o += 1;

          i += 1 + fill_byte_cntr_len;
          break;
        }

      }
    }
    return bitmap;
  }

  static std::string
  name() {
    return "bbc";
  }

  /// Returns the value of the bit at the position pos.
  u1 __forceinline__
  test(const std::size_t pos) const {
    const auto dec = to_plain_bitmap(); // FIXME HACK
    return dec.test(pos);
  }

  /// Try to reduce the memory consumption. This function is supposed to be
  /// called after the bitmap has been modified.
  __forceinline__ void
  shrink() {
    data_.shrink_to_fit();
  }

  //===--------------------------------------------------------------------===//
  /// 1-run iterator
  class iter {
    const bbc& outer_;

    //===------------------------------------------------------------------===//
    // Iterator state
    //===------------------------------------------------------------------===//
    using plain_bitmap_t = dtl::plain_bitmap<$u64>;
    std::unique_ptr<plain_bitmap_t> plain_bitmap_; // FIXME HACK
    using nested_iter_t = dtl::plain_bitmap_iter<dtl::plain_bitmap<$u64>>;
    std::unique_ptr<nested_iter_t> plain_bitmap_iter_;
    //===------------------------------------------------------------------===//

  public:
    explicit __forceinline__
    iter(const bbc& outer)
        : outer_(outer) {
      plain_bitmap_ = std::make_unique<plain_bitmap_t>(outer.to_plain_bitmap());
      plain_bitmap_iter_ = std::make_unique<nested_iter_t>(*plain_bitmap_);
    }

    /// Forward the iterator to the next 1-run.
    void __forceinline__
    next() {
      plain_bitmap_iter_->next();
    }

    /// Forward the iterator to the desired position.
    void __forceinline__
    skip_to(const std::size_t to_pos) {
      plain_bitmap_iter_->skip_to(to_pos);
    }

    u1 __forceinline__
    end() const noexcept {
      return plain_bitmap_iter_->end();
    }

    u64 __forceinline__
    pos() const noexcept {
      return plain_bitmap_iter_->pos();
    }

    u64 __forceinline__
    length() const noexcept {
      return plain_bitmap_iter_->length();
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
        + ",\"word_size\":" + std::to_string(sizeof(byte_t))
        + "}";
  }

  // For debugging purposes.
  void
  print(std::ostream& os) const {
    os << "\nword idx | content" << std::endl;
    for (std::size_t i = 0; i < data_.size(); ++i) {
      os << std::setw(8) << i << " | ";
      auto& w = data_[i];
      os << "'" << std::bitset<8>(w);
      os << std::endl;
    }
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
