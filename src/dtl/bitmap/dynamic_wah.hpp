#pragma once

#include <cstddef>

#include <dtl/dtl.hpp>

#include <fastbit/bitvector.h>
#include <fastbit/bitvector64.h>
#include <fastbit/fileManager.h>
#include <boost/dynamic_bitset.hpp>

namespace dtl {

namespace internal {

static const auto& filemanager_instance = ibis::fileManager::instance(); // initialize file manager, which is responsible for memory management in IBIS.

//===----------------------------------------------------------------------===//
/// WAH compressed representation of a bitmap of length N using either
/// 32- or 64-bit words.
template<typename bitvector_t = ibis::bitvector>
struct dynamic_wah {

  bitvector_t bv;
  std::size_t size_;

  dynamic_wah() = default;

  explicit
  dynamic_wah(const boost::dynamic_bitset<$u32>& in)
      : size_(in.size()) {
    for (std::size_t i = 0; i < size_; i++) {
      bv.setBit(i, in[i]);
    }
    bv.compress();
  }

  ~dynamic_wah() = default;

  dynamic_wah(const dynamic_wah& other) = default;

  dynamic_wah(dynamic_wah&& other) noexcept = default;

  dynamic_wah&
  operator=(const dynamic_wah& other) = default;

  dynamic_wah&
  operator=(dynamic_wah&& other) noexcept = default;

  /// Return the size in bytes.
  std::size_t
  size_in_byte() const {
    return bv.bytes();
  }

  /// Conversion to an std::bitset.
  boost::dynamic_bitset<$u32>
  to_bitset() {
    boost::dynamic_bitset<$u32> ret(size_);
    typename bitvector_t::pit pit(bv);
    while (*pit < size_) {
      ret[*pit] = true;
      pit.next();
    }
    return ret;
  }

  /// Bitwise AND
  dynamic_wah
  operator&(const dynamic_wah& other) const {
    dynamic_wah ret(*this);
    ret.bv &= other.bv;
    return ret;
  }

  /// Bitwise AND (range encoding)
  dynamic_wah
  and_re(const dynamic_wah& other) const {
    return *this & other; // nothing special here. fall back to AND
  }

  /// Bitwise XOR
  dynamic_wah
  operator^(const dynamic_wah& other) const {
    dynamic_wah ret(*this);
    ret.bv ^= other.bv;
    return ret;
  }

  /// Bitwise AND (range encoding)
  dynamic_wah
  xor_re(const dynamic_wah& other) const {
    return *this ^ other; // nothing special here. fall back to XOR
  }

  /// Computes (a XOR b) & this
  /// Note: this, a and b must be different instances. Otherwise, the behavior
  /// is undefined.
  dynamic_wah&
  fused_xor_and(const dynamic_wah& a, const dynamic_wah& b) {
    auto x = a ^ b;
    bv &= x.bv;
    return *this;
  }

  void
  print(std::ostream& os) const {
    os << "n/a";
  }

  static std::string
  name() {
    return "dynamic_wah" + std::to_string(sizeof(typename bitvector_t::word_t) * 8);
  }

  /// Returns the value of the bit at the position pos.
  u1
  test(const std::size_t pos) const {
    return bv.getBit(pos);
  }


};
//===----------------------------------------------------------------------===//

} // namespace internal

//===----------------------------------------------------------------------===//
/// WAH compressed representation of a bitmap of length N using 32-bit words.
using dynamic_wah32 = internal::dynamic_wah<ibis::bitvector>;
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
/// WAH compressed representation of a bitmap of length N using 64-bit words.
using dynamic_wah64 = internal::dynamic_wah<ibis::bitvector64>;
//===----------------------------------------------------------------------===//

} // namespace dtl
