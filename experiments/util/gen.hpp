#pragma once

#include <iostream>
#include <sstream>
#include <dtl/dtl.hpp>
#include <dtl/bitmap/util/random.hpp>
#include <dtl/bitmap/util/two_state_markov_process.hpp>

//===----------------------------------------------------------------------===//
static u1
markov_parameters_are_valid(u64 n, $f64 f, $f64 d) {
  return d >= 0
      && d <= 1.0
      && f >= d/(1-d)
      && f <= d*n;
}

/// Generate a random bitmap with the given parameters.
/// The actual f and d are at most 2% off.
///
/// Throws an exception if the parameters are invalid or a random bitmap could
/// not be constructed after several retries.
static dtl::bitmap
gen_random_bitmap_markov(u64 n, $f64 f, $f64 d) {
  if (!markov_parameters_are_valid(n, f, d)) {
    throw std::invalid_argument("Invalid parameters for the Markov process.");
  }
  two_state_markov_process mp(f, d);
  dtl::bitmap bs(n);

  // Error bounds.
  f64 e = 0.05;
  f64 f_min = f - f * e;
  f64 f_max = f + f * e;
  f64 d_min = d - d * e;
  f64 d_max = d + d * e;

  // Retry if error is too large.
  for ($u64 r = 0; r < 10000; ++r) {
    bs.reset();
    // Populate the bitmap.
    for ($u64 i = 0; i < n; i++) {
      bs[i] = mp.next();
    }
    const auto f_actual = dtl::determine_clustering_factor(bs);
    const auto d_actual = dtl::determine_bit_density(bs);

    if (f_actual >= f_min
        && f_actual <= f_max
        && d_actual >= d_min
        && d_actual <= d_max) {
      return bs;
    }
  }
  std::stringstream err;
  err << "Failed: n=" << n
      << ", f=" << f
      << ", d=" << d
      << std::endl;
  std::cerr << err.str();
  throw std::invalid_argument(
      "Failed to construct a random bitmap with the given parameters after "
      "several retries.");
}
//===----------------------------------------------------------------------===//
/// Generate a uniformly populated random bitmap with the given density.
static dtl::bitmap
gen_random_bitmap_uniform(u64 n, $f64 d) {
  if (d < 0 || d > 1.0) {
    throw std::invalid_argument("Invalid bit density.");
  }

  dtl::bitmap bs = dtl::gen_random_bitmap_uniform(n, d);

  return bs;
}
//===----------------------------------------------------------------------===//
