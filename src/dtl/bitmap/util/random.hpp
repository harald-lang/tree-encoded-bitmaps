#pragma once

#include <random>
#include <dtl/bitmap.hpp>
#include "two_state_markov_process.hpp"

namespace dtl {

//===----------------------------------------------------------------------===//
static inline dtl::bitmap
gen_random_bitmap_markov(u64 n, $f64 f, $f64 d) {

  // init bitset
  f64 f_min = d >= 1.0 ? n : d/(1-d);
  f64 f_actual = std::max(f, f_min);
  two_state_markov_process mp(f_actual, d);
  dtl::bitmap bs(n);
  for ($u64 i = 0; i < n; ++i) {
    bs[i] = mp.next();
  }
  return bs;
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
static inline dtl::bitmap
gen_random_bitmap_uniform(u64 n, $f64 d) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<std::size_t> dis(0, RAND_MAX);

  // Adapted from https://stackoverflow.com/questions/2509679/how-to-generate-a-
  //              random-integer-number-from-within-a-range
  auto random_at_most = [&](long max) {
    unsigned long
    // max <= RAND_MAX < ULONG_MAX, so this is okay.
        num_bins = (unsigned long) max + 1,
        num_rand = (unsigned long) RAND_MAX + 1,
        bin_size = num_rand / num_bins,
        defect   = num_rand % num_bins;

    long x;
    do {
      x = dis(gen);
    }
      // This is carefully written not to overflow
    while (num_rand - defect <= (unsigned long)x);

    // Truncated division is intentional
    return x / bin_size;
  };

  dtl::bitmap bs(n);
  if (d <= 0) {
    return bs;
  }
  if (d >= 1.0) {
    bs.flip();
    return bs;
  }

  const auto m = static_cast<std::size_t>(n * d);
  for (std::size_t i = 0; i < m; ++i) {
    bs[i] = true;
  }
  // Fisher-Yates shuffle
  for (std::size_t i = n - 1; i > 0; --i) {
    const auto j = random_at_most(i - 1);
    u1 t = bs[i];
    bs[i] = bs[j];
    bs[j] = t;
  }
  return bs;
}
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
/// Count the number of 1-fills in the given bitmap.
template<typename bitset_t>
static std::size_t
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
static f64
determine_bit_density(const bitset_t& b) {
  if (b.size() == 0) return 0.0;
  return (b.count() * 1.0) / b.size();
}
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
/// Determine the average length of the 1-fills in the given bitmap.
template<typename bitset_t>
static f64
determine_clustering_factor(const bitset_t& b) {
  if (b.size() == 0) return 0.0;
  if (b.count() == 0) return 0.0;
  return (b.count() * 1.0) / count_1fills(b);
}
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
/// Count the number of runs of the of the given value.
static std::size_t
count_value_runs(const std::vector<$u32>& b, u32 val) {
  if (b.size() == 0) return 0;
  auto last_val = b[0];
  std::size_t cntr = b[0] == val;
  for (std::size_t i = 1; i < b.size(); i++) {
    const auto current_val = b[i];
    cntr += last_val != val
      && current_val == val;
    last_val = current_val;
  }
  assert(cntr > 0);
  return cntr;
}
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
/// Count the occurrences of the given value.
static std::size_t
count_value(const std::vector<$u32>& b, u32 val) {
  if (b.size() == 0) return 0;
  std::size_t cntr = b[0] == val;
  for (std::size_t i = 1; i < b.size(); i++) {
    cntr += b[i] == val;
  }
  return cntr;
}
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
/// Determine the average run-length of the given value.
static f64
determine_clustering_factor(const std::vector<$u32>& b, u32 val) {
  if (b.size() == 0) return 0.0;
  const auto cnt = count_value(b, val);
  if (cnt == 0) return 0.0;
  return (cnt * 1.0) / count_value_runs(b, val);
}
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
/// Determine the average run-length.
static f64
determine_clustering_factor(const std::vector<$u32>& b) {
  if (b.size() == 0) return 0;
  auto last_val = b[0];
  std::size_t run_cntr = 1;
  for (std::size_t i = 1; i < b.size(); i++) {
    const auto current_val = b[i];
    run_cntr += last_val != current_val;
    last_val = current_val;
  }
  assert(run_cntr > 0);
  return (b.size() * 1.0) / run_cntr;
}
//===----------------------------------------------------------------------===//

} // namespace dtl
