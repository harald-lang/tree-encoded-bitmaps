#pragma once

#include <bitset>
#include <cstddef>

#include <boost/dynamic_bitset.hpp>

//===----------------------------------------------------------------------===//
// Utility functions for bitsets.
//===----------------------------------------------------------------------===//

/// Count the number of 1-fills in the given bitset.
template<typename bitset_t>
std::size_t
count_1fills(const bitset_t& b) {
  if (b.size() == 0) return 0;
  bool last_bit = b[0];
  std::size_t cntr = last_bit;
  for (std::size_t i = 1; i < b.size(); i++) {
    const bool current_bit = b[i];
    cntr += current_bit & !last_bit;
    last_bit = current_bit;
  }
  return cntr;
}

/// Count the number of 0-fills in the given bitset.
template<typename bitset_t>
std::size_t
count_0fills(const bitset_t& b) {
  if (b.size() == 0) return 0;
  bool last_bit = b[0];
  std::size_t cntr = !last_bit;
  for (std::size_t i = 1; i < b.size(); i++) {
    const bool current_bit = b[i];
    cntr += !current_bit & last_bit;
    last_bit = current_bit;
  }
  return cntr;
}

template<typename bitset_t>
std::vector<$u1>
to_vector(const bitset_t& b) {
  std::vector<$u1> ret(b.size(), false);
  for (std::size_t i = 0; i < b.size(); i++) {
    ret[i] = b[i];
  }
  return ret;
}


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
