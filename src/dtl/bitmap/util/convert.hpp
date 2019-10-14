#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/bitmap/util/bitmap_fun.hpp>
#include <dtl/dtl.hpp>

#include <boost/dynamic_bitset.hpp>
//===----------------------------------------------------------------------===//
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
/// Reconstruct a plain bitmap using the run iterator of the type under test.
template<typename T>
boost::dynamic_bitset<$u32>
to_bitmap_using_iterator(const T& encoded_bitmap) {
  boost::dynamic_bitset<$u32> bm(encoded_bitmap.size());
  auto it = encoded_bitmap.scan_it();
  while (!it.end()) {
#ifndef BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS
    for (std::size_t i = it.pos(); i < it.pos() + it.length(); ++i) {
      assert(bm[i] == false);
      bm[i] = true;
    }
#else
    // HACK: This gives access to the private members of the boost::dynamic_bitset.
    const std::size_t b = it.pos();
    const std::size_t e = it.length() + b;
    dtl::bitmap_fun<$u32>::set(bm.m_bits.data(), b, e);
#endif
    it.next();
  }
  return bm;
}
//===----------------------------------------------------------------------===//
template<typename It>
boost::dynamic_bitset<$u32>
to_bitmap_from_iterator(It& it, u64 n) {
  boost::dynamic_bitset<$u32> bm(n);
  while (!it.end()) {
#ifndef BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS
    for (std::size_t i = it.pos(); i < it.pos() + it.length(); ++i) {
      assert(bm[i] == false);
      bm[i] = true;
    }
#else
    // HACK: This gives access to the private members of the boost::dynamic_bitset.
    const std::size_t b = it.pos();
    const std::size_t e = it.length() + b;
    dtl::bitmap_fun<$u32>::set(bm.m_bits.data(), b, e);
#endif
    it.next();
  }
  return bm;
}
//===----------------------------------------------------------------------===//
/// Convert the given bitmap into a position list.
template<typename T>
static std::vector<$u32>
to_position_list(const T& bitmap) {
  std::vector<$u32> positions;
  auto it = bitmap.scan_it();
  while (!it.end()) {
    for (std::size_t i = it.pos(); i < it.pos() + it.length(); ++i) {
      positions.push_back(i);
    }
    it.next();
  }
  return positions;
}
static std::vector<$u32>
to_position_list(const boost::dynamic_bitset<$u32>& bitmap) {
  std::vector<$u32> positions;
  auto i = bitmap.find_first();
  while (i != boost::dynamic_bitset<$u32>::npos) {
    positions.push_back(i);
    i = bitmap.find_next(i);
  }
  return positions;
}
//===----------------------------------------------------------------------===//
} // namespace dtl