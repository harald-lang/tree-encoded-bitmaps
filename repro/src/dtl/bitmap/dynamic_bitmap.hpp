#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/dtl.hpp>
#include <dtl/math.hpp>

#include <boost/dynamic_bitset.hpp>

#include <cstddef>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// -UN-compressed representation of a bitmap of variable length.
/// Wraps a boost dynamic_bitset. (used for comparisons)
template<typename _block_type = u64>
struct dynamic_bitmap {
  boost::dynamic_bitset<_block_type> bitmap_;

  // The number of bits.
  u64 n_;

  dynamic_bitmap() = default;

  /// Constructs a bitmap from a bitmap, thereby the set bits are copied
  /// one by one, similarly to the other implementations.
  explicit dynamic_bitmap(const boost::dynamic_bitset<_block_type>& in)
      : bitmap_(in.size(), 0), n_(in.size()) {
    auto i = in.find_first();
    while (i != boost::dynamic_bitset<$u32>::npos) {
      bitmap_[i] = true;
      i = in.find_next(i);
    }
  }

  ~dynamic_bitmap() = default;
  dynamic_bitmap(const dynamic_bitmap& other) = default;
  dynamic_bitmap(dynamic_bitmap&& other) noexcept = default;
  dynamic_bitmap& operator=(const dynamic_bitmap& other) = default;
  dynamic_bitmap& operator=(dynamic_bitmap&& other) noexcept = default;

  /// Return the size in bytes.
  std::size_t
  size_in_byte() const {
    const auto bits_per_block = sizeof(_block_type) * 8;
    return ((n_ + bits_per_block - 1) / bits_per_block) * sizeof(_block_type) /* bitmap */
        + 4 /* length */;
  }

  /// Returns the length of the bitmap.
  inline std::size_t
  size() const {
    return n_;
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
  dynamic_bitmap __forceinline__
  operator^(const dynamic_bitmap& other) const {
    dynamic_bitmap ret(bitmap_);
    ret.bitmap_ ^= other.bitmap_;
    return ret;
  }

  /// Bitwise XOR (range encoding)
  dynamic_bitmap __forceinline__
  xor_re(const dynamic_bitmap& other) const {
    return *this ^ other; // nothing special here. fall back to XOR
  }

  /// Computes (a XOR b) & this
  /// Note: this, a and b must be different instances. Otherwise, the behavior
  /// is undefined.
  __forceinline__ dynamic_bitmap&
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
    return "bitmap";
  }

  /// Returns the value of the bit at the position pos.
  u1 __forceinline__
  test(const std::size_t pos) const {
    return bitmap_.test(pos);
  }

  //===--------------------------------------------------------------------===//
  /// 1-fill iterator, with skip support.
  class iter {
    const dynamic_bitmap& outer_;

    //===------------------------------------------------------------------===//
    // Iterator state
    //===------------------------------------------------------------------===//
    /// points to the beginning of a 1-fill
    $u64 pos_;
    /// the length of the current 1-fill
    $u64 length_;
    //===------------------------------------------------------------------===//

  public:
    void __forceinline__
    next() {
      pos_ += length_;
      length_ = 0;
      pos_ = outer_.bitmap_.find_next(pos_);
      // determine the length of the current 1fill
      if (pos_ < outer_.n_) {
        length_ = 1;
        while (pos_ + length_ < outer_.n_
            && outer_.bitmap_.test(pos_ + length_)) {
          length_++;
        }
      }
      else {
        pos_ = outer_.n_;
        length_ = 0;
      }
    }

    explicit __forceinline__
    iter(const dynamic_bitmap& outer)
        : outer_(outer), pos_(outer_.bitmap_.find_first()), length_(0) {
      // determine the length of the current 1fill
      if (pos_ < outer_.n_) {
        length_ = 1;
        while (pos_ + length_ < outer_.n_
            && outer_.bitmap_.test(pos_ + length_)) {
          length_++;
        }
      }
      else {
        pos_ = outer_.n_;
        length_ = 0;
      }
    }

    void __forceinline__
    skip_to(const std::size_t to_pos) {
      if (to_pos >= outer_.n_) {
        pos_ = outer_.n_;
        length_ = 0;
        return;
      }
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
  //===--------------------------------------------------------------------===//

  using skip_iter_type = iter;
  using scan_iter_type = iter;

  skip_iter_type __forceinline__
  it() const {
    return skip_iter_type(*this);
  }

  scan_iter_type __forceinline__
  scan_it() const {
    return scan_iter_type(*this);
  }

  /// Returns the name of the instance including the most important parameters
  /// in JSON.
  std::string
  info() const {
    return "{\"name\":\"" + name() + "\""
        + ",\"n\":" + std::to_string(n_)
        + ",\"size\":" + std::to_string(size_in_byte())
        + "}";
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
