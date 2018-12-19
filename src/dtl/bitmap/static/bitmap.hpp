#pragma once

#include <cstddef>
#include <bitset>

#include <dtl/dtl.hpp>


namespace dtl {

//===----------------------------------------------------------------------===//
/// -UN-compressed representation of a bitmap of length N.
/// Wraps a standard bitset. (used for comparisons)
template<std::size_t N>
struct bitmap {

  std::bitset<N> bitmap_;

  bitmap() = default;

  explicit
  bitmap(const std::bitset<N>& in) {
    bitmap_ = in;
  }

  ~bitmap() = default;

  bitmap(const bitmap& other) = default;

  bitmap(bitmap&& other) noexcept = default;

  bitmap&
  operator=(const bitmap& other) = default;

  bitmap&
  operator=(bitmap&& other) noexcept = default;

  /// Return the size in bytes.
  std::size_t
  size_in_byte() {
    return ((N + 7) / 8) + 4 /* length */;
  }

  /// Conversion to an std::bitset.
  std::bitset<N>
  to_bitset() {
    std::bitset<N> ret(bitmap_);
    return ret;
  }

  /// Bitwise AND
  bitmap
  operator&(const bitmap& other) const {
    bitmap ret(*this);
    ret.bitmap_ &= other.bitmap_;
    return ret;
  }

  /// Bitwise AND (range encoding)
  bitmap
  and_re(const bitmap& other) const {
    return *this & other; // nothing special here. fall back to AND
  }

  /// Bitwise XOR
  bitmap
  operator^(const bitmap& other) const {
    bitmap ret;
    ret.bitmap_ = bitmap_;
    ret.bitmap_ ^= other.bitmap_;
    return ret;
  }

  /// Bitwise XOR (range encoding)
  bitmap
  xor_re(const bitmap& other) const {
    return *this ^ other; // nothing special here. fall back to XOR
  }

  /// Computes (a XOR b) & this
  /// Note: this, a and b must be different instances. Otherwise, the behavior is undefined.
  bitmap&
  fused_xor_and(const bitmap& a, const bitmap& b) {
    auto x = a ^ b;
    bitmap_ &= x.bitmap_;
    return *this;
  }

  void
  print(std::ostream& os) const {
    os << "n/a";
  }

  static std::string
  name() {
    return "bitmap";
  }

  /// Returns the value of the bit at the position pos.
  u1
  test(const std::size_t pos) const {
    return bitmap_.test(pos);
  }


};
//===----------------------------------------------------------------------===//


} // namespace dtl
