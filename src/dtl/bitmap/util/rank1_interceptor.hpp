#pragma once

#include <vector>

#include <dtl/dtl.hpp>
#include <dtl/math.hpp>
#include <boost/dynamic_bitset.hpp>

namespace dtl {

//===----------------------------------------------------------------------===//
// Rank wrapper, used for profiling and debugging.
//===----------------------------------------------------------------------===//

static thread_local $u64 last_idx = ~0ull;

template<typename _wrapped_type>
struct rank1_interceptor {

  using wrapped_type = _wrapped_type;
  using word_type = typename _wrapped_type::word_type;
  using size_type = typename _wrapped_type::size_type;

  wrapped_type wrapped_instance_;

  ~rank1_interceptor() = default;

  void
  init(const boost::dynamic_bitset<word_type>& bitmap) {
    last_idx = ~0ull;
    wrapped_instance_.init(bitmap);
  }

  size_type __forceinline__
  operator()(u64 idx, const word_type* bitmap_ptr) const {
    std::cout << "rank " << idx << std::endl;
    if (idx == last_idx) {
      std::cout << "breakpoint - idx=" << idx << std::endl;
    }
    last_idx = idx;
    return wrapped_instance_(idx, bitmap_ptr);
  }

  u64
  size_in_bytes() const {
    return wrapped_instance_.size_in_bytes();
  }

  static u64
  estimate_size_in_bytes(u64 bitmap_size) {
    return wrapped_type::estimate_size_in_bytes(bitmap_size);
  }

  /// Returns the important properties in JSON.
  std::string
  info() const noexcept {
    return "{\"name\":" + std::string("\"Interceptor\"") + "}";
  }

};
//===----------------------------------------------------------------------===//
} // namespace dtl

