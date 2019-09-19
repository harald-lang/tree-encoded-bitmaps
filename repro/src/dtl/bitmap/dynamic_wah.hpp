#pragma once
//===----------------------------------------------------------------------===//
#include <fastbit/bitvector.h>
#include <fastbit/bitvector64.h>
#include <fastbit/fileManager.h>

#include <dtl/dtl.hpp>

#include <boost/dynamic_bitset.hpp>

#include <cstddef>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
namespace internal {
//===----------------------------------------------------------------------===//
// Initialize file manager, which is responsible for memory management in IBIS.
static const auto& filemanager_instance = ibis::fileManager::instance();
//===----------------------------------------------------------------------===//
/// WAH compressed representation of a bitmap of length N using either
/// 32- or 64-bit words.
template<typename bitvector_t = ibis::bitvector>
struct dynamic_wah {
  bitvector_t bv;
  std::size_t size_;

  dynamic_wah() = default;

  explicit dynamic_wah(const boost::dynamic_bitset<$u32>& in)
      : size_(in.size()) {
    auto i = in.find_first();
    while (i != boost::dynamic_bitset<$u32>::npos) {
      bv.setBit(i, in[i]);
      i = in.find_next(i);
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
    return bv.bytes() /* size of the compressed bitmap */
        + sizeof(size_); /* bit-length of the original bitmap */
  }

  /// Returns the size of the bitmap.
  std::size_t
  size() const {
    return size_;
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
    const auto word_size = sizeof(typename bitvector_t::word_t) * 8;
    if (word_size == 4) {
      return "wah";
    }
    else {
      return "wah" + std::to_string(word_size);
    }
  }

  /// Returns the value of the bit at the position pos.
  u1 test(const std::size_t pos) const {
    return bv.getBit(pos);
  }

  //===--------------------------------------------------------------------===//
  /// 1-fill iterator, with skip support.
  class iter {
    const dynamic_wah& outer_;

    /// Position iterator.
    typename bitvector_t::pit pit;

    //===------------------------------------------------------------------===//
    // Iterator state
    //===------------------------------------------------------------------===//
    /// points to the beginning of a 1-fill
    $u64 pos_;
    /// the length of the current 1-fill
    $u64 length_;
    //===------------------------------------------------------------------===//

  public:
    explicit __forceinline__
    iter(const dynamic_wah& outer)
        : outer_(outer),
          pit(outer_.bv),
          pos_(*pit), length_(0) {
      // determine the length of the current 1fill
      if (*pit < outer_.size_) {
        length_ = 1;
        pit.next();
        while (*pit < outer_.size_
            && *pit == pos_ + length_) {
          ++length_;
          pit.next();
        }
      }
      else {
        pos_ = outer_.size_;
        length_ = 0;
      }
    }

    void __forceinline__
    next() {
      // Determine the length of the current 1-fill.
      if (*pit < outer_.size_) {
        pos_ = *pit;
        length_ = 1;
        pit.next();
        while (*pit < outer_.size_
            && *pit == pos_ + length_) {
          ++length_;
          pit.next();
        }
      }
      else {
        pos_ = outer_.size_;
        length_ = 0;
      }
    }

    void __forceinline__
    skip_to(const std::size_t to_pos) {
      assert(pos_ <= to_pos);
      if (to_pos >= outer_.size_) {
        pos_ = outer_.size_;
        length_ = 0;
        return;
      }
      if (to_pos < (pos_ + length_)) {
        length_ -= to_pos - pos_;
        pos_ = to_pos;
        return;
      }
      while (*pit < to_pos) {
        pit.next();
      }
      // Determine the length of the current 1-fill.
      if (*pit < outer_.size_) {
        pos_ = *pit;
        length_ = 1;
        pit.next();
        while (*pit < outer_.size_
            && *pit == pos_ + length_) {
          ++length_;
          pit.next();
        }
      }
      else {
        pos_ = outer_.size_;
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

  iter __forceinline__
  it() const {
    return iter(*this);
  }

  iter __forceinline__
  scan_it() const {
    return iter(*this);
  }

  /// Returns the name of the instance including the most important parameters
  /// in JSON.
  std::string
  info() const {
    return "{\"name\":\"" + name() + "\""
        + ",\"n\":" + std::to_string(size_)
        + ",\"size\":" + std::to_string(size_in_byte())
        + ",\"word_size\":" + std::to_string(sizeof(typename bitvector_t::word_t))
        + "}";
  }
};
//===----------------------------------------------------------------------===//
} // namespace internal
//===----------------------------------------------------------------------===//
/// WAH compressed representation of a bitmap of length N using 32-bit words.
using dynamic_wah32 = internal::dynamic_wah<ibis::bitvector>;
/// WAH compressed representation of a bitmap of length N using 64-bit words.
using dynamic_wah64 = internal::dynamic_wah<ibis::bitvector64>;
//===----------------------------------------------------------------------===//
} // namespace dtl
