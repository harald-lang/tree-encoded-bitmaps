#pragma once

#include <dtl/dtl.hpp>
#include <boost/dynamic_bitset.hpp>

namespace dtl {

//===----------------------------------------------------------------------===//
template<typename bitset_t>
boost::dynamic_bitset<$u32>
to_dynamic_bitset(const bitset_t& b) {
  boost::dynamic_bitset<$u32> ret(b.size(), false);
  for (std::size_t i = 0; i < b.size(); i++) {
    ret[i] = b[i];
  }
  return ret;
}
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
/// Reconstruct a plain bitmap using the range iterator of the type under test.
template<typename T>
boost::dynamic_bitset<$u32>
to_bitmap_using_iterator(const T& encoded_bitmap) {
  boost::dynamic_bitset<$u32> bm(encoded_bitmap.size());
  auto it = encoded_bitmap.it();
  while (!it.end()) {
    for (std::size_t i = it.pos(); i < it.pos() + it.length(); ++i) {
      assert(bm[i] == false);
      bm[i] = true;
    }
    it.next();
  }
  return bm;
}
//===----------------------------------------------------------------------===//


} // namespace dtl