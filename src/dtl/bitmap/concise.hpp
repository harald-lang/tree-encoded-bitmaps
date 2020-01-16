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
/// CONCISE -- a variant of WAH that improves compression by dealing with dirty
/// bits. These dirty bits are piggybacked by the fill words.
///
/// The implementation is based on the paper "CONCISE: Compressed 'n' Composable
/// Integer Set" of Colantonio et al.
///
/// Note: The implementation is NOT yet optimized for performance.
class concise {

protected:
  using word_type = $u32;
  static constexpr std::size_t word_bitlength = sizeof(word_type) * 8;
  /// The number of bits required to store the position of a dirty bit.
  static constexpr std::size_t dirty_position_bit_cnt = 5; // log2(32)
  /// The maximum repetition value that can be stored in a fill word.
  static constexpr std::size_t max_fill_repetitions =
      (1ull << (word_bitlength - 2 - dirty_position_bit_cnt)) - 1;
  /// The number of bits that can be stored in a literal word.
  static constexpr std::size_t payload_bit_cnt = word_bitlength - 1;
  static constexpr word_type all_ones = word_type(~word_type(0));

  /// The encoded bitmap.
  std::vector<word_type> data_;
  /// The length of the encoded bitmap.
  std::size_t encoded_bitmap_length_ = 0;
  std::size_t max_ = 0;

  std::size_t remaining_bit_cnt_in_last_literal_word_ = 0;

  //===--------------------------------------------------------------------===//
  // Helper functions to work with the code words.
  //===--------------------------------------------------------------------===//
  /// Returns true if the given word is a literal word.
  static inline constexpr u1
  is_literal_word(const word_type w) noexcept {
    return (w & word_type(1)) == 0;
  }

  /// Returns true if the given word is a fill word.
  static inline constexpr u1
  is_fill_word(const word_type w) noexcept {
    return (w & word_type(1)) == 1;
  }

  /// Extract the payload from a literal word.
  static inline constexpr word_type
  extract_payload(const word_type w) noexcept {
    assert(is_literal_word(w));
    return w >> 1;
  }

  /// Extract the fill value from a fill word.
  static inline constexpr u1
  extract_fill_value(const word_type w) noexcept {
    assert(is_fill_word(w));
    return (w & word_type(2)) != 0;
  }

  /// Extract the fill length from a fill word. The actual length is:
  /// (returned value + 1) * (word_bitlength - 1)
  static inline constexpr std::size_t
  extract_fill_repetitions(const word_type w) noexcept {
    assert(is_fill_word(w));
    return w >> (2 + dirty_position_bit_cnt);
  }

  /// Returns true if the given fill word piggybacks a dirty bit position.
  static inline constexpr u1
  has_dirty_bit_position(const word_type w) noexcept {
    assert(is_fill_word(w));
    const word_type mask = (1u << dirty_position_bit_cnt) - 1;
    return ((w >> 2) & mask) != 0;
  }

  /// Extract the position of the dirty bit from a fill word.
  static inline constexpr std::size_t
  extract_dirty_bit_position(const word_type w) noexcept {
    assert(is_fill_word(w));
    const word_type mask = (1u << dirty_position_bit_cnt) - 1;
    const auto raw_val = (w >> 2) & mask;
    assert(raw_val != 0);
    return raw_val - 1;
  }

  /// Increments the counter in the given fill word by the given value.
  static inline constexpr word_type
  increment_fill_repetitions(const word_type w, const std::size_t i) noexcept {
    assert(is_fill_word(w));
    const auto val = extract_fill_repetitions(w) + i;
    assert(val <= max_fill_repetitions);
    const word_type mask = (1u << (2 + dirty_position_bit_cnt)) - 1;
    return (w & mask) | (val << (2 + dirty_position_bit_cnt));
  }

  /// Creates a fill word.
  static inline constexpr word_type
  make_fill_word(u1 val, const std::size_t len) noexcept {
    assert(len <= max_fill_repetitions);
    word_type w = val ? word_type(3) : word_type(1);
    w |= len << (2 + dirty_position_bit_cnt);
    return w;
  }

  /// Creates a fill word that contains a dirty bit.
  static inline constexpr word_type
  make_fill_word(u1 val, const std::size_t len,
      const std::size_t dirty_bit_pos) noexcept {
    assert(len <= max_fill_repetitions);
    word_type w = val ? word_type(3) : word_type(1);
    w |= (dirty_bit_pos + 1) << 2;
    w |= len << (2 + dirty_position_bit_cnt);
    return w;
  }

  /// Returns a literal word where all bits are 0.
  static inline constexpr word_type
  make_empty_literal_word() noexcept {
    return word_type(0);
  }

  /// Returns a literal word where all bits are set to 1.
  static inline constexpr word_type
  make_all_ones_literal_word() noexcept {
    return ~word_type(0) << 1;
  }

  /// Tests a bit within a literal word.
  static inline constexpr u1
  test_bit_in_literal_word(const word_type w, const std::size_t i) noexcept {
    assert(i <= payload_bit_cnt);
    return ((w >> (1 + i)) & 1) != 0;
  }

  /// Sets a bit within a literal word.
  static inline constexpr word_type
  set_bit_in_literal_word(const word_type w, const std::size_t i, const u1 val)
      noexcept {
    assert(i <= payload_bit_cnt);
    const u1 cur_val = test_bit_in_literal_word(w, i);
    if (cur_val == val) {
      return w;
    }
    // Toggle the bit.
    return w ^ (word_type(1u << (i + 1)));
  }

  /// Returns true, iff the payload of the given literal word consists of a
  /// single 1-bit.
  static inline constexpr u1
  contains_single_bit(const word_type w) {
    assert(is_literal_word(w));
    return dtl::bits::pop_count(w) == 1;
  }
  //===--------------------------------------------------------------------===//

  /// Try to merge the last two words.
  inline void
  compress() {
    if (data_.size() < 2) {
      return;
    }

    const word_type w_last = data_.back();
    const word_type w_prev = data_[data_.size() - 2];

    assert(is_literal_word(w_last));

    u1 last_word_is_all_zeros = (w_last == make_empty_literal_word());
    u1 last_word_is_all_ones = (w_last == make_all_ones_literal_word());

    if (!last_word_is_all_zeros && !last_word_is_all_ones) {
      return;
    }

    u1 prev_word_is_zero_fill =
        is_fill_word(w_prev) && (extract_fill_value(w_prev) == false);
    u1 prev_word_is_one_fill =
        is_fill_word(w_prev) && (extract_fill_value(w_prev) == true);

    if ((last_word_is_all_zeros && prev_word_is_zero_fill)
        || (last_word_is_all_ones && prev_word_is_one_fill)) {
      data_.pop_back();
      data_.back() = increment_fill_repetitions(data_.back(), 1);
      return;
    }

    if ((last_word_is_all_zeros && prev_word_is_one_fill)
        && (last_word_is_all_ones && prev_word_is_zero_fill)) {
      return;
    }

    assert(is_literal_word(w_prev));

    word_type payload_prev = extract_payload(w_prev);
    if (last_word_is_all_ones) {
      payload_prev = (~payload_prev) & ((1u << payload_bit_cnt) - 1);
    }
    if (payload_prev == 0 || dtl::bits::pop_count(payload_prev) == 1) {
      // Merge the last two literal words into a fill word.
      data_.pop_back();
      u1 fill_value = last_word_is_all_ones;
      if (dtl::bits::pop_count(payload_prev) == 1) {
        const auto dirty_bit_pos = dtl::bits::tz_count(payload_prev);
        data_.back() = make_fill_word(fill_value, 1, dirty_bit_pos);
      }
      else {
        data_.back() = make_fill_word(fill_value, 1);
      }
    }
  }

  /// Append a 1-bit to the encoded bitmap at the given position.
  inline void
  append(const std::size_t i) {
    assert(i >= max_);
    if (data_.empty()) {
      const auto zero_block_cnt = i / payload_bit_cnt;
      if (zero_block_cnt == 1) {
        // Append an all-0 literal word.
        data_.push_back(make_empty_literal_word());
      }
      else if (zero_block_cnt > 1) {
        // Append a 0-fill.
        data_.push_back(make_fill_word(0, zero_block_cnt - 1));
      }
      // Append the literal word that contains the (very first) 1-bit.
      auto lit = make_empty_literal_word();
      const auto bit_pos = i % payload_bit_cnt;
      lit = set_bit_in_literal_word(lit, bit_pos, 1);
      data_.push_back(lit);
      max_ = i;
      return;
    }

    // The bit position to set, relative to the current block.
    auto bit_pos = (max_ % payload_bit_cnt) + i - max_;

    // Add 0-bits between the last set bit and the target bit that is to be
    // appended, if necessary.
    if (bit_pos >= payload_bit_cnt) {
      // The target bit is beyond the current block.
      const auto zero_block_cnt = (bit_pos / payload_bit_cnt) - 1;
      if (zero_block_cnt > 0) {
        // Append a 0-fill.
        if (contains_single_bit(data_.back())) {
          // Convert the last word into a fill word that piggybacks that bit.
          const auto dirty_bit_pos = dtl::bits::tz_count(
              extract_payload(data_.back()));
          data_.back() = make_fill_word(0, zero_block_cnt, dirty_bit_pos);
        }
        else {
          if (zero_block_cnt == 1) {
            data_.push_back(make_empty_literal_word()); // all-0 literal
          }
          else {
            data_.push_back(make_fill_word(0, zero_block_cnt - 1)); // 0-fill
          };
        }
      }

      // Append an all-0 literal word (the active word).
      data_.push_back(make_empty_literal_word());

      // Adjust the relative position of the bit to append.
      bit_pos = bit_pos % payload_bit_cnt;
    }

    assert(is_literal_word(data_.back()));
    data_.back() = set_bit_in_literal_word(data_.back(), bit_pos, 1);
    max_ = i;

    compress();
  }

  /// Append a 0-run of the given length to the compressed bitmap.
  inline void
  append_zero_run(std::size_t len) {
    encoded_bitmap_length_ += len;
  }

  /// Append a 1-run of the given length to the compressed bitmap.
  inline void
  append_one_run(const std::size_t len) {
    for (std::size_t i = 0; i < len; ++i) {
      append(encoded_bitmap_length_ + i);
    }
    encoded_bitmap_length_ += len;
  }

public:
  concise() = default;

  explicit concise(const boost::dynamic_bitset<$u32>& in) {
    // Obtain a 1-run iterator for the input bitmap.
    dtl::plain_bitmap_iter<boost::dynamic_bitset<$u32>> it(in);
    // Append runs to the compressed bitmap.
    std::size_t i = 0;
    if (it.pos() > 0) {
      // The bitmap starts with a 0-run.
      auto len = it.pos() - i;
      append_zero_run(len);
      i += len;
    }
    while (!it.end()) {
      auto len = it.length();
      append_one_run(len);
      i += len;
      it.next();
      if (!it.end()) {
        len = it.pos() - i;
        append_zero_run(len);
        i += len;
      }
    }
    if (i < in.size()) {
      append_zero_run(in.size() - i);
    }
    // Try to reduce the memory consumption.
    shrink();
  }

  ~concise() = default;
  concise(const concise& other) = default;
  concise(concise&& other) noexcept = default;
  concise& operator=(const concise& other) = default;
  concise& operator=(concise&& other) noexcept = default;

  /// Return the size in bytes.
  std::size_t __forceinline__
  size_in_byte() const {
    return data_.size() * sizeof(word_type) /* size of the compressed bitmap */
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
    std::size_t i = 0;
    for (std::size_t word_idx = 0; word_idx < data_.size(); ++word_idx) {
      auto& w = data_[word_idx];
      if (is_fill_word(w)) {
        auto val = extract_fill_value(w);
        auto len = extract_fill_repetitions(w) * payload_bit_cnt;
        ret.set(i, i + len, val);
        i += len;
      }
      else {
        const std::size_t cnt = word_idx == data_.size() - 1
            ? payload_bit_cnt - remaining_bit_cnt_in_last_literal_word_
            : payload_bit_cnt;
        for (std::size_t k = 0; k < cnt; ++k) {
          ret.set(i + k, dtl::bits::bit_test(w, k + 1));
        }
        i += cnt;
      }
    }
    if (i < encoded_bitmap_length_) {
      ret.clear(i, encoded_bitmap_length_);
    }
    return std::move(ret);
  }

  static std::string
  name() {
    return "concise" + std::to_string(word_bitlength);
  }

  /// Returns the value of the bit at the position pos.
  u1 __forceinline__
  test(const std::size_t pos) const {
    std::size_t i = 0;
    std::size_t word_idx = 0;
    // Find the corresponding word.
    for (; word_idx < data_.size(); ++word_idx) {
      auto& w = data_[word_idx];
      if (is_literal_word(w)) {
        if (pos >= i + payload_bit_cnt) {
          i += payload_bit_cnt;
          continue;
        }
        else {
          return dtl::bits::bit_test(w, pos - i + 1); // TODO optimize
        }
      }
      else {
        auto fill_len = (extract_fill_repetitions(w) + 1) * payload_bit_cnt;
        if (pos >= i + fill_len) {
          i += fill_len;
          continue;
        }
        else {
          $u1 val = extract_fill_value(w);
          if (has_dirty_bit_position(w)
              && (pos % payload_bit_cnt == extract_dirty_bit_position(w))) {
              val = !val;
          }
          return val;
        }
      }
    }
    return false;
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
    const concise& outer_;

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

    void __forceinline__
    increment() {
      const auto word_cnt = outer_.data_.size();
      word_type w;
      while (word_idx_ < word_cnt) {
        w = outer_.data_[word_idx_];
        if (is_fill_word(w)) {
          if (!has_dirty_bit_position(w)) {
            if (extract_fill_value(w) == false) {
              pos_ += (extract_fill_repetitions(w) + 1) * payload_bit_cnt;
            }
            else {
              length_ = (extract_fill_repetitions(w) + 1) * payload_bit_cnt;
              in_word_idx_ = 0;
              ++word_idx_;
              break;
            }
          }
          else {
            // Fill with dirty bit.
            const auto dirty_bit_pos = extract_dirty_bit_position(w);
            const auto dirty_bit_val = !extract_fill_value(w);

            if (dirty_bit_val) { // 000001000000
              if (in_word_idx_ == 0) {
                // Skip forward to the dirty bit and produce a 1-run of length 1.
                pos_ += dirty_bit_pos;
                length_ = 1;
                in_word_idx_ = 2; // refers to everything past the dirty bit
                break;
              }
              else {
                assert(in_word_idx_ == 2);
                pos_ += (extract_fill_repetitions(w) + 1) * payload_bit_cnt
                    - (dirty_bit_pos + 1);
                in_word_idx_ = 0;
              }
            }
            else {
              if (dirty_bit_pos == 0) { // Handle special case 011111111
                pos_ += 1; // Skip over the dirty 0-bit.
                length_ = (extract_fill_repetitions(w) + 1) * payload_bit_cnt - 1;
                ++word_idx_;
                in_word_idx_ = 0;
                break;
              }
              else { // 11011111
                if (in_word_idx_ == 0) {
                  length_ += dirty_bit_pos;
                  in_word_idx_ = 1;
                  break;
                }
                else if (in_word_idx_ == 1) {
                  pos_ += 1; // Skip over the dirty 0-bit.
                  length_ = (extract_fill_repetitions(w) + 1) * payload_bit_cnt
                      - (dirty_bit_pos + 1);
                  ++word_idx_;
                  in_word_idx_ = 0;
                  break;
                }
              }
            }
          }
        }
        else { // Literal
          if (in_word_idx_ < payload_bit_cnt) { // TODO decode the entire literal word at once.
            const word_type payload = extract_payload(w) >> in_word_idx_;
            if (payload == 0) {
              pos_ += payload_bit_cnt - in_word_idx_;
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
              in_word_idx_ += e;
              break;
            }
          }
        }
        in_word_idx_ = 0;
        ++word_idx_;
      }
    }

  public:
    explicit __forceinline__
    iter(const concise& outer)
        : outer_(outer),
          word_idx_(0),
          in_word_idx_(0),
          pos_(0), length_(0) {
      increment();
    }

    /// Forward the iterator to the next 1-run.
    void __forceinline__
    next() {
      pos_ += length_;
      length_ = 0;
      increment();
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
        + ",\"size\":" + std::to_string(size_in_byte())
        + ",\"word_size\":" + std::to_string(sizeof(word_type))
        + "}";
  }

  // For debugging purposes.
  void
  print(std::ostream& os) const {
    os << "\nword idx | word type | content" << std::endl;
    for (std::size_t i = 0; i < data_.size(); ++i) {
      os << std::setw(8) << i << " | ";
      auto& w = data_[i];
      if (is_fill_word(w)) {
        os << "fill      | ";
        os << "'" << std::bitset<word_bitlength>(w);
        os << "' => ";
        os << extract_fill_repetitions(w) << " x ";
        os << (extract_fill_value(w) ? "'1'" : "'0'");
        if (has_dirty_bit_position(w)) {
          os << " (dirty bit at " << extract_dirty_bit_position(w) << ")";
        }
        os << std::endl;
      }
      else {
        os << "literal   | ";
        os << "'" << std::bitset<word_bitlength>(w);
        os << "'";
        os << std::endl;
      }
    }
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
