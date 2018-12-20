#pragma once

#include <cstddef>
#include <bitset>

#include <dtl/dtl.hpp>


namespace dtl {

//===----------------------------------------------------------------------===//
/// -UN-compressed representation of a bitmap of length N.
/// Wraps a standard bitset. (used for comparisons)
template<std::size_t N>
struct uncompressed_bitmap {

  std::bitset<N> bitmap_;

  uncompressed_bitmap() = default;

  explicit
  uncompressed_bitmap(const std::bitset<N>& in) {
    bitmap_ = in;
  }

  ~uncompressed_bitmap() = default;

  uncompressed_bitmap(const uncompressed_bitmap& other) = default;

  uncompressed_bitmap(uncompressed_bitmap&& other) noexcept = default;

  uncompressed_bitmap&
  operator=(const uncompressed_bitmap& other) = default;

  uncompressed_bitmap&
  operator=(uncompressed_bitmap&& other) noexcept = default;

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
  uncompressed_bitmap
  operator&(const uncompressed_bitmap& other) const {
    uncompressed_bitmap ret(*this);
    ret.bitmap_ &= other.bitmap_;
    return ret;
  }

  /// Bitwise AND (range encoding)
  uncompressed_bitmap
  and_re(const uncompressed_bitmap& other) const {
    return *this & other; // nothing special here. fall back to AND
  }

  /// Bitwise XOR
  uncompressed_bitmap
  operator^(const uncompressed_bitmap& other) const {
    uncompressed_bitmap ret;
    ret.bitmap_ = bitmap_;
    ret.bitmap_ ^= other.bitmap_;
    return ret;
  }

  /// Bitwise XOR (range encoding)
  uncompressed_bitmap
  xor_re(const uncompressed_bitmap& other) const {
    return *this ^ other; // nothing special here. fall back to XOR
  }

  /// Computes (a XOR b) & this
  /// Note: this, a and b must be different instances. Otherwise, the behavior is undefined.
  uncompressed_bitmap&
  fused_xor_and(const uncompressed_bitmap& a, const uncompressed_bitmap& b) {
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
