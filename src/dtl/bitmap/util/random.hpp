#pragma once

#include <dtl/bitmap.hpp>
#include "two_state_markov_process.hpp"

namespace dtl {

//===----------------------------------------------------------------------===//
inline dtl::bitmap
gen_random_bitmap(u64 n, $f64 f, $f64 d) {

  // init bitset
  f64 f_min = d >= 1.0 ? n : d/(1-d);
  f64 f_actual = std::max(f, f_min);
  two_state_markov_process mp(f_actual, d);
  dtl::bitmap bs(n);
  for ($u64 i = 0; i < n; i++) {
    bs[i] = mp.next();
  }

  $f64 d_actual = (bs.count() * 1.0) / n;
  if (std::abs(d - d_actual) > 1
      || std::abs(f - f_actual) > 0.25) {
    throw std::invalid_argument(
        "Failed to construct a random bitmap with the given parameters.");
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

} // namespace dtl
