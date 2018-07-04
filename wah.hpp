#pragma once

#include <cstddef>

#include <dtl/dtl.hpp>

#include <fastbit/bitvector.h>
#include <fastbit/bitvector64.h>
#include <fastbit/fileManager.h>

namespace dtl {

namespace internal {

static const auto& filemanager_instance = ibis::fileManager::instance(); // initialize file manager, which is responsible for memory management in IBIS.

//===----------------------------------------------------------------------===//
/// WAH compressed representation of a bitmap of length N using either
/// 32- or 64-bit words.
template<std::size_t N, typename bitvector_t = ibis::bitvector>
struct wah {

  bitvector_t bv;

  wah() = default;

  explicit
  wah(const std::bitset<N>& in) {
    for (std::size_t i = 0; i < N; i++) {
      bv.setBit(i, in[i]);
    }
    bv.compress();
  }

  ~wah() = default;

  wah(const wah& other) = default;

  wah(wah&& other) noexcept = default;

  wah&
  operator=(const wah& other) = default;

  wah&
  operator=(wah&& other) noexcept = default;

  /// Return the size in bytes.
  std::size_t
  size_in_byte() {
    return bv.bytes();
  }

  /// Conversion to an std::bitset.
  std::bitset<N>
  to_bitset() {
    std::bitset<N> ret;
    typename bitvector_t::pit pit(bv);
    while (*pit < N) {
      ret[*pit] = true;
      pit.next();
    }
    return ret;
  }

//  /// Bitwise AND
//  wah
//  operator&(const wah& other) const {
//    wah ret;
//    ret.bv.copy(bv);
//    ret.bv &= other.bv;
//    return ret;
//  }

  /// Bitwise XOR
  wah
  operator^(const wah& other) const {
    wah ret;
    ret.bv = bv ^ other.bv;
    return ret;
  }

  /// Computes (a XOR b) & this
  /// Note: this, a and b must be different instances. Otherwise, the behavior is undefined.
  wah&
  fused_xor_and(const wah& a, const wah& b) {
    auto x = a ^ b;
    bv &= x.bv;
    return *this;
  }

  void
  print(std::ostream& os) const {
    os << "n/a";
  }

};
//===----------------------------------------------------------------------===//

} // namespace internal

//===----------------------------------------------------------------------===//
/// WAH compressed representation of a bitmap of length N using 32-bit words.
template<std::size_t N>
using wah32 = internal::wah<N, ibis::bitvector>;
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
/// WAH compressed representation of a bitmap of length N using 64-bit words.
template<std::size_t N>
using wah64 = internal::wah<N, ibis::bitvector64>;
//===----------------------------------------------------------------------===//

} // namespace dtl
