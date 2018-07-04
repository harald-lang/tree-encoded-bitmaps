#pragma once

#include <cstddef>

#include <dtl/dtl.hpp>

#include <roaring/../../cpp/roaring.hh>
#include <roaring/roaring.h>

namespace dtl {

//===----------------------------------------------------------------------===//
/// Compressed representation of a bitmap of length N.
/// Wraps a Roaring bitmap.
template<std::size_t N>
struct roaring_bitmap {

  Roaring bitmap;

  roaring_bitmap() = default;

  explicit
  roaring_bitmap(const std::bitset<N>& in) {
    for (std::size_t i = 0; i < N; i++) {
      if (in[i]) { bitmap.add(i); };
    }
    bitmap.runOptimize();
  }

  ~roaring_bitmap() = default;

  roaring_bitmap(const roaring_bitmap& other) = default;

  roaring_bitmap(roaring_bitmap&& other) noexcept = default;

  roaring_bitmap&
  operator=(const roaring_bitmap& other) = default;

  roaring_bitmap&
  operator=(roaring_bitmap&& other) noexcept = default;

  /// Return the size in bytes.
  std::size_t
  size_in_byte() {
    return bitmap.getSizeInBytes(false);
  }

  /// Conversion to an std::bitset.
  std::bitset<N>
  to_bitset() {
    std::bitset<N> ret;
    for(Roaring::const_iterator i = bitmap.begin() ; i != bitmap.end() ; i++) {
      ret[*i] = true;
    }
    return ret;
  }

  /// Bitwise AND
  roaring_bitmap
  operator&(const roaring_bitmap& other) const {
    roaring_bitmap ret(*this);
    ret.bitmap &= other.bitmap;
    return ret;
  }

  /// Bitwise XOR
  roaring_bitmap
  operator^(const roaring_bitmap& other) const {
    roaring_bitmap ret;
    ret.bitmap = bitmap;
    ret.bitmap ^= other.bitmap;
    return ret;
  }

  /// Computes (a XOR b) & this
  /// Note: this, a and b must be different instances. Otherwise, the behavior is undefined.
  roaring_bitmap&
  fused_xor_and(const roaring_bitmap& a, const roaring_bitmap& b) {
    auto x = a ^ b;
    bitmap &= x.bitmap;
    return *this;
  }

  void
  print(std::ostream& os) const {
    os << "n/a";
  }


};
//===----------------------------------------------------------------------===//


} // namespace dtl
