#pragma once
//===----------------------------------------------------------------------===//
#include "bitmap_fun.hpp"

#include <dtl/dtl.hpp>

#include <cassert>
#include <cstddef>
#include <type_traits>
#include <vector>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// A plain uncompressed fixed size bitmap.
template<
    typename _word_type,
    typename _alloc = std::allocator<_word_type>
>
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
  std::vector<word_type, _alloc> bitmap_;
  /// The allocator.
  _alloc allocator_;

public:

  /// Construct a bitmap of length n.
  explicit
  plain_bitmap(std::size_t n, const _alloc& alloc = _alloc())
      : n_(n), bitmap_(word_cnt(n)), allocator_(alloc) {}

  /// Test the i-th bit.
  u1 __forceinline__
  operator[](std::size_t i) const {
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
    assert(b < e);
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
    assert(b < e);
    bitmap_fun<word_type>::clear(bitmap_.data(), b, e);
  }

  __forceinline__ word_type*
  data() {
    return bitmap_.data();
  }

  /// Fetch up to size(word_type)*8 consecutive bits from the range [b,e).
  /// Note: The behavior is undefined when the range spans more bits than fit in
  /// a single word.
  word_type __forceinline__
  fetch_bits(std::size_t b, std::size_t e) const {
    assert((b / word_bitlength) < bitmap_.size());
    assert(b <= n_);
    assert(e <= n_);
    assert(b < e);
    assert((e - b) <= (sizeof(word_type) * 8));
    return bitmap_fun<word_type>::fetch_bits(bitmap_.data(), b, e);
  }

};
//===----------------------------------------------------------------------===//
} // namespace dtl
