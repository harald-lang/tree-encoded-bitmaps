#pragma once

#include <dtl/bitmap.hpp>
#include "two_state_markov_process.hpp"

namespace dtl {

//===----------------------------------------------------------------------===//
inline dtl::bitmap
gen_random_bitmap_markov(u64 n, $f64 f, $f64 d) {

  // init bitset
  f64 f_min = d >= 1.0 ? n : d/(1-d);
  f64 f_actual = std::max(f, f_min);
  two_state_markov_process mp(f_actual, d);
  dtl::bitmap bs(n);
  for ($u64 i = 0; i < n; i++) {
    bs[i] = mp.next();
  }
  return bs;
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
/// Count the number of 1-fills in the given bitmap.
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
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
/// Count the number of 1-fills in the given bitmap.
template<typename bitset_t>
f64
determine_bit_density(const bitset_t& b) {
  if (b.size() == 0) return 0.0;
  return (b.count() * 1.0) / b.size();
}
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
/// Count the number of 1-fills in the given bitmap.
template<typename bitset_t>
f64
determine_clustering_factor(const bitset_t& b) {
  if (b.size() == 0) return 0.0;
  if (b.count() == 0) return 0.0;
  return (b.count() * 1.0) / count_1fills(b);
}
//===----------------------------------------------------------------------===//

} // namespace dtl
