#include "test_utils.hpp"

//===----------------------------------------------------------------------===//
dtl::bitmap
gen_bitmap(u64 n, $f64 f, $f64 d) {

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
    throw std::invalid_argument("Failed to construct a random bitmap with the given parameters.");
  }
  return bs;
}
//===----------------------------------------------------------------------===//

