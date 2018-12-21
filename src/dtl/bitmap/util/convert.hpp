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

} // namespace dtl