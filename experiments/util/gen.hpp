#pragma once
//===----------------------------------------------------------------------===//
#include "bitmap_db.hpp"
#include "params.hpp"

#include <dtl/bitmap/util/markov_process.hpp>
#include <dtl/bitmap/util/random.hpp>
#include <dtl/bitmap/util/two_state_markov_process.hpp>
#include <dtl/dtl.hpp>

#include <iostream>
#include <sstream>
#include <vector>
//===----------------------------------------------------------------------===//
static u1
markov_parameters_are_valid(u64 n, $f64 f, $f64 d) {
  return d >= 0
      && d <= 1.0
      && f >= d / (1 - d)
      && f <= d * n;
}
//===----------------------------------------------------------------------===//
/// Generate a random bitmap with the given parameters.
/// The actual f and d are at most 2% off.
///
/// Throws an exception if the parameters are invalid or a random bitmap could
/// not be constructed after several retries.
static dtl::bitmap
gen_random_bitmap_markov(u64 n, $f64 f, $f64 d) {
  if (!markov_parameters_are_valid(n, f, d)) {
    throw std::invalid_argument(
        "Invalid parameters for the Markov process: n="
        + std::to_string(n) + ", f=" + std::to_string(f) + ", d="
        + std::to_string(d));
  }
  two_state_markov_process mp(f, d);
  dtl::bitmap bs(n);

  // Error bounds.
  f64 e = 0.02;
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
// TODO max error
static dtl::bitmap
gen_random_bitmap_uniform(u64 n, $f64 d) {
  if (d < 0 || d > 1.0) {
    throw std::invalid_argument("Invalid bit density.");
  }

  dtl::bitmap bs = dtl::gen_random_bitmap_uniform(n, d);

  return bs;
}
//===----------------------------------------------------------------------===//
/// Generate a random integer sequence of length n with integers in [0, c).
/// The parameter f refers to the clustering factor, which is the average length
/// of runs with identical values. If f is 1, the result is uniformly
/// distributed. - The actual f is at most 2% off.
///
/// Throws an exception if the parameters are invalid or a random bitmap could
/// not be constructed after several retries.
static std::vector<$u32>
gen_random_integer_sequence_markov(u64 n, u32 c, $f64 f) {
  //  if (!markov_parameters_are_valid(n, f, d)) { // TODO
  //    throw std::invalid_argument("Invalid parameters for the Markov process.");
  //  }
  assert(c < n);
  assert(c * f < n);

  markov_process mp(c, f);
  std::vector<$u32> seq(n);

  // Error bounds.
  f64 e = 0.02;
  f64 f_min = f - f * e;
  f64 f_max = f + f * e;

  // Retry if error is too large.
  for ($u64 r = 0; r < 10000; ++r) {
    std::fill(seq.begin(), seq.end(), 0);
    // Populate the integer vector.
    for ($u64 i = 0; i < n; i++) {
      seq[i] = mp.next();
    }
    const auto f_actual = dtl::determine_clustering_factor(seq);
    std::vector<$u32> cpy = seq;
    std::sort(cpy.begin(), cpy.end());
    const auto unique_cnt = std::unique(cpy.begin(), cpy.end()) - cpy.begin();
    if (f_actual >= f_min
        && f_actual <= f_max
        && unique_cnt == c) {
      return seq;
    }
  }
  std::stringstream err;
  err << "Failed: n=" << n
      << ", c=" << c
      << ", f=" << f
      << std::endl;
  std::cerr << err.str();
  throw std::invalid_argument(
      "Failed to construct a random integer sequence with the given parameters "
      "after several retries.");
}
//===----------------------------------------------------------------------===//
/// Generates multiple (Markov) random bitmaps in parallel and stores them in
/// the database. When the return values is > 0, then some of the bitmap
/// couldn't be generated.
std::size_t
gen(std::vector<params_markov>& params,
    bitmap_db& db // the database, where the generated bitmaps are stored
);
//===----------------------------------------------------------------------===//
std::size_t
gen(std::vector<params_uniform>& params,
    bitmap_db& db // the database, where the generated bitmaps are stored
);
//===----------------------------------------------------------------------===//
