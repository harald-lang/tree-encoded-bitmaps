#pragma once
//===----------------------------------------------------------------------===//
#include "bitmap_fun.hpp"

#include <dtl/dtl.hpp>
#include <dtl/math.hpp>

#include <boost/dynamic_bitset.hpp>

#include <vector>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
template<typename _word_type = $u32>
struct rank1_naive {

  using word_type = typename std::remove_cv<_word_type>::type;

  using size_type = $u32;

  ~rank1_naive() = default;

  std::size_t bitmap_length_ = 0;

  void
  init(const boost::dynamic_bitset<word_type>& bitmap) {
    bitmap_length_ = bitmap.size();
  }

  size_type __forceinline__
  operator()(u64 idx, const word_type* bitmap_ptr) const {
    assert(idx < bitmap_length_);
    size_type ret_val = 0;
    for (std::size_t i = 0; i <= idx; ++i) {
      ret_val += bitmap_fun<word_type>::test(bitmap_ptr, i);
    }
    return ret_val;
  }

  u64
  size_in_bytes() const {
    return 0;
  }

  static u64
  estimate_size_in_bytes(u64 bitmap_size) {
    return 0;
  }

  /// Returns the important properties in JSON.
  std::string
  info() const noexcept {
    return "{\"name\":" + std::string("\"naive\"")
        + ",\"size\":" + std::to_string(size_in_bytes())
        + "}";
  }

};
//===----------------------------------------------------------------------===//
} // namespace dtl
