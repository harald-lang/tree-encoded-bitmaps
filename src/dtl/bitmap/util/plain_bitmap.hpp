#pragma once
//===----------------------------------------------------------------------===//
#include "bitmap_fun.hpp"
#include "bitmap_writer.hpp"
#include "buffer.hpp"

#include <dtl/dtl.hpp>

#include <cassert>
#include <cstddef>
#include <memory>
#include <type_traits>
#include <vector>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// A plain uncompressed fixed size bitmap.
template<
    typename _word_type,
    typename _alloc = std::allocator<_word_type>>
class plain_bitmap {
  using word_type = typename std::remove_cv<_word_type>::type;
  static_assert(std::is_integral<word_type>::value,
      "The word type must be an integral type.");
  static_assert(!std::is_same<word_type, bool>::value,
      "The word type must not be a boolean.");

  static constexpr auto word_bitlength = sizeof(_word_type) * 8;

  /// Returns the number of words required to store a bitmap of length n.
  static constexpr auto word_cnt(std::size_t n) {
    return (n + word_bitlength - 1) / word_bitlength;
  }

  /// The number of bits.
  std::size_t n_;
  /// The actual bitmap.
  dtl::buffer<word_type, _alloc> bitmap_;
  /// The allocator.
  _alloc allocator_; // TODO remove

public:
  /// Construct a bitmap of length n.
  explicit plain_bitmap(std::size_t n, const _alloc& alloc = _alloc())
      : n_(n), bitmap_(word_cnt(n)), allocator_(alloc) {}

  /// Construct a bitmap of length n.
  explicit plain_bitmap(std::size_t n, u1 init, const _alloc& alloc = _alloc())
      : n_(n), bitmap_(word_cnt(n), init), allocator_(alloc) {
  }

  /// Return the size of the bitmap.
  std::size_t
  size() const noexcept {
    return n_;
  }

  /// Test the i-th bit.
  u1 __forceinline__
  operator[](std::size_t i) const {
    assert((i / word_bitlength) < bitmap_.size());
    assert(i <= n_);
    return bitmap_fun<word_type>::test(bitmap_.data(), i);
  }

  /// Test the i-th bit.
  u1 __forceinline__
  test(std::size_t i) const {
    assert((i / word_bitlength) < bitmap_.size());
    assert(i <= n_);
    return bitmap_fun<word_type>::test(bitmap_.data(), i);
  }

  /// Set the i-th bit.
  void __forceinline__
  set(std::size_t i, u1 val) {
    assert((i / word_bitlength) < bitmap_.size());
    assert(i <= n_);
    bitmap_fun<word_type>::set(bitmap_.data(), i, val);
  }

  /// Set the i-th bit to 1.
  void __forceinline__
  set(std::size_t i) {
    assert((i / word_bitlength) < bitmap_.size());
    assert(i <= n_);
    bitmap_fun<word_type>::set(bitmap_.data(), i);
  }

  /// Set the bits in [b,e) to 1.
  void __forceinline__
  set(std::size_t b, std::size_t e) {
    assert((b / word_bitlength) < bitmap_.size());
    assert(b <= n_);
    assert(e <= n_);
    assert(b <= e);
    bitmap_fun<word_type>::set(bitmap_.data(), b, e);
  }

  /// Set the i-th bit to 0.
  void __forceinline__
  clear(std::size_t i) {
    assert((i / word_bitlength) < bitmap_.size());
    assert(i <= n_);
    bitmap_fun<word_type>::clear(bitmap_.data(), i);
  }

  /// Set the bits in [b,e) to 0.
  void __forceinline__
  clear(std::size_t b, std::size_t e) {
    assert((b / word_bitlength) < bitmap_.size());
    assert(b <= n_);
    assert(e <= n_);
    assert(b <= e);
    bitmap_fun<word_type>::clear(bitmap_.data(), b, e);
  }

  /// Returns a pointer to the first word of the bitmap.
  __forceinline__ word_type*
  data() noexcept {
    return bitmap_.data();
  }

  __forceinline__ const word_type*
  data() const noexcept {
    return bitmap_.data();
  }

  /// Returns a pointer to the first word of the bitmap.
  __forceinline__ word_type*
  data_begin() noexcept {
    return bitmap_.data();
  }

  /// Returns a pointer to one past the last word of the bitmap.
  __forceinline__ word_type*
  data_end() noexcept {
    return bitmap_.data() + bitmap_.size();
  }

  /// Fetch up to size(word_type)*8 consecutive bits from the range [b,e).
  /// Note: The behavior is undefined when the range spans more bits than fit in
  /// a single word.
  word_type __forceinline__
  fetch_bits(std::size_t b, std::size_t e) const {
    assert((b / word_bitlength) < bitmap_.size());
    assert(b <= n_);
    assert(e <= n_);
    assert(b <= e);
    assert((e - b) <= (sizeof(word_type) * 8));
    return bitmap_fun<word_type>::fetch_bits(bitmap_.data(), b, e);
  }

  void __forceinline__
  store_bits(std::size_t b, std::size_t e, word_type bits_to_store) {
    assert((b / word_bitlength) < bitmap_.size());
    assert(b <= n_);
    assert(e <= n_);
    assert(b <= e);
    assert((e - b) <= (sizeof(word_type) * 8));
    return bitmap_fun<word_type>::store_bits(bitmap_.data(), b, e, bits_to_store);
  }

  /// Bitwise AND
  plain_bitmap __forceinline__
  operator&(const plain_bitmap& other) const {
    assert(size() == other.size());
    plain_bitmap ret(size(), false);
    for (std::size_t w = 0; w < bitmap_.size(); ++w) {
      ret.data()[w] = data()[w] & other.data()[w]; // TODO SIMDize
    }
    return ret;
  }

  /// Bitwise OR
  plain_bitmap __forceinline__
  operator|(const plain_bitmap& other) const {
    assert(size() == other.size());
    plain_bitmap ret(size(), false);
    for (std::size_t w = 0; w < bitmap_.size(); ++w) {
      ret.data()[w] = data()[w] | other.data()[w]; // TODO SIMDize
    }
    return ret;
  }

  /// Bitwise XOR
  plain_bitmap __forceinline__
  operator^(const plain_bitmap& other) const {
    assert(size() == other.size());
    plain_bitmap ret(size(), false);
    for (std::size_t w = 0; w < bitmap_.size(); ++w) {
      ret.data()[w] = data()[w] ^ other.data()[w]; // TODO SIMDize
    }
    return ret;
  }

  /// Bitwise NOT
  plain_bitmap __forceinline__
  operator~() const {
    plain_bitmap ret(size(), false);
    for (std::size_t w = 0; w < bitmap_.size(); ++w) {
      ret.data()[w] = ~(data()[w]); // TODO SIMDize
    }
    return ret;
  }

  /// Returns the position of the first set bit. If no bits are set, the length
  /// of the bitmap is returned.
  std::size_t __attribute__((noinline))
  find_first() const {
    const auto ret_val = bitmap_fun<word_type>::find_first(
        bitmap_.data(), bitmap_.data() + bitmap_.size());
    return std::min(ret_val, n_);
  }

  /// Returns the position of the first set bit within the range [b,e). If no
  /// bits are set, e is returned.
  std::size_t __attribute__((noinline))
  find_first(std::size_t b, std::size_t e) const {
    assert(b <= n_);
    assert(e <= n_);
    assert(b <= e);
    const auto ret_val = bitmap_fun<word_type>::find_first(bitmap_.data(), b, e);
    assert(ret_val >= b);
    assert(ret_val <= e);
    assert(ret_val == bitmap_fun<word_type>::find_first_dense(bitmap_.data(), b, e));
    return ret_val;
  }

  /// Returns the position of the last set bit. If no bits are set, the length
  /// of the bitmap is returned.
  std::size_t __attribute__((noinline))
  find_last() const {
    const auto ret_val = bitmap_fun<word_type>::find_last(
        bitmap_.data(), bitmap_.data() + bitmap_.size());
    return std::min(ret_val, n_);
  }

  /// Finds the position of the next set bit within the range (b,e). If no set
  /// bit is found, the length of the bitmap is returned.
  inline std::size_t
  find_next(std::size_t b) const noexcept {
    if (b >= (n_ - 1)) {
      return n_;
    }
    const auto ret_val = bitmap_fun<word_type>::find_next(
        bitmap_.data(), bitmap_.data() + bitmap_.size(), b);
    return std::min(ret_val, n_);
  }

  /// Counts the number of set bits within the range [b,e).
  std::size_t __attribute__((noinline))
  count(std::size_t b, std::size_t e) const noexcept {
    const auto ret_val = bitmap_fun<word_type>::count(bitmap_.data(), b, e);
#ifndef NDEBUG
    std::size_t cntr = 0;
    for (std::size_t i = b; i <  e; ++i) {
      cntr += test(i);
    }
    assert(cntr == ret_val);
#endif
    return ret_val;
  }

  void
  print(std::ostream& os) const noexcept {
    for (std::size_t i = 0; i < n_; ++i) {
      os << (test(i) ? "1" : "0");
    }
  }

  bitmap_writer<word_type>
  writer(std::size_t start_idx) {
    return bitmap_writer<word_type>(bitmap_.data(), start_idx);
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
