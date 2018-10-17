#pragma once

#include <cstddef>

#include <dtl/dtl.hpp>
#include <dtl/math.hpp>

#include <boost/dynamic_bitset.hpp>

namespace dtl {

//===----------------------------------------------------------------------===//
/// -UN-compressed representation of a bitmap of variable length.
/// Wraps a boost dynamic_bitset. (used for comparisons)
template<typename _block_type = u64>
struct dynamic_bitmap {

  const boost::dynamic_bitset<_block_type> bitmap_;

  // the number of bits
  u64 n_;


  dynamic_bitmap() = default;

  explicit
  dynamic_bitmap(const boost::dynamic_bitset<_block_type>& in)
    : bitmap_(in), n_(in.size()) {

    if (!dtl::is_power_of_two(n_)) {
      throw std::invalid_argument("The length of the bitmask must be a power of two.");
    }

  }

  ~dynamic_bitmap() = default;

  dynamic_bitmap(const dynamic_bitmap& other) = default;

  dynamic_bitmap(dynamic_bitmap&& other) noexcept = default;

  dynamic_bitmap&
  operator=(const dynamic_bitmap& other) = default;

  dynamic_bitmap&
  operator=(dynamic_bitmap&& other) noexcept = default;

  /// Return the size in bytes.
  std::size_t
  size_in_byte() {
    return ((n_ + 7) / 8) + 4 /* length */;
  }

  /// Conversion to an boost::dynamic_bitset.
  boost::dynamic_bitset<_block_type>
  to_bitset() {
    boost::dynamic_bitset<_block_type> ret(bitmap_);
    return ret;
  }

  /// Bitwise AND
  dynamic_bitmap
  operator&(const dynamic_bitmap& other) const {
    dynamic_bitmap ret(*this);
    ret.bitmap_ &= other.bitmap_;
    return ret;
  }

  /// Bitwise AND (range encoding)
  dynamic_bitmap
  and_re(const dynamic_bitmap& other) const {
    return *this & other; // nothing special here. fall back to AND
  }

  /// Bitwise XOR
  dynamic_bitmap
  operator^(const dynamic_bitmap& other) const {
    dynamic_bitmap ret;
    ret.bitmap_ = bitmap_;
    ret.bitmap_ ^= other.bitmap_;
    return ret;
  }

  /// Bitwise XOR (range encoding)
  dynamic_bitmap
  xor_re(const dynamic_bitmap& other) const {
    return *this ^ other; // nothing special here. fall back to XOR
  }

  /// Computes (a XOR b) & this
  /// Note: this, a and b must be different instances. Otherwise, the behavior is undefined.
  dynamic_bitmap&
  fused_xor_and(const dynamic_bitmap& a, const dynamic_bitmap& b) {
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
    return "dynamic_bitmap";
  }

  /// Returns the value of the bit at the position pos.
  u1
  test(const std::size_t pos) const {
    return bitmap_.test(pos);
  }


  //===----------------------------------------------------------------------===//
  /// 1-fill iterator, with skip support.
  class iter {

    const dynamic_bitmap& outer_;

    //===----------------------------------------------------------------------===//
    // Iterator state
    //===----------------------------------------------------------------------===//
    /// points to the beginning of a 1-fill
    $u64 pos_;
    /// the length of the current 1-fill
    $u64 length_;
    //===----------------------------------------------------------------------===//

  public:

    void
    next() {
      pos_ += length_;
      length_ = 0;
      pos_ = outer_.bitmap_.find_next(pos_);
      // determine the length of the current 1fill
      if (pos_ < outer_.n_) {
        length_ = 1;
        while (outer_.bitmap_.test(pos_ + length_)) { // TODO optimize
          length_++;
        }
      }
      else {
        pos_ = outer_.n_;
        length_ = 0;
      }
    }

    explicit
    iter(const dynamic_bitmap& outer)
        : outer_(outer), pos_(outer_.bitmap_.find_first()), length_(0) {
      // determine the length of the current 1fill
      if (pos_ < outer_.n_) {
        length_ = 1;
        while (pos_ + length_ < outer_.n_ && outer_.bitmap_.test(pos_ + length_)) {
          length_++;
        }
      }
      else {
        pos_ = outer_.n_;
        length_ = 0;
      }
    }

    void
    skip_to(const std::size_t to_pos) {
      pos_ = to_pos;
      length_ = 0;
      if (!outer_.bitmap_.test(pos_)) {
        pos_ = outer_.bitmap_.find_next(pos_);
      }
      if (pos_ < outer_.n_) {
        length_ = 1;
        while (pos_ + length_ < outer_.n_ && outer_.bitmap_.test(pos_ + length_)) {
          length_++;
        }
      }
      else {
        pos_ = outer_.n_;
        length_ = 0;
      }
    }

    u1
    end() const noexcept {
      return length_ == 0;
    }

    u64
    pos() const noexcept {
      return pos_;
    }

    u64
    length() const noexcept {
      return length_;
    }

  };
  //===----------------------------------------------------------------------===//

  iter
  it() const {
    return iter(*this);
  }


};
//===----------------------------------------------------------------------===//


} // namespace dtl
