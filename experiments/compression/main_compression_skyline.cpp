#include <iostream>

#include <dtl/dtl.hpp>

#include "common.hpp"

//===----------------------------------------------------------------------===//
// Experiment: Comparison of compression ratios.
//===----------------------------------------------------------------------===//


$i32 main() {

//  //===----------------------------------------------------------------------===//
//  // Bit density
//  u64 d_lo = dtl::next_power_of_two(dtl::env<$u64>::get("D_LO", 1ull << 10));
//  u64 d_lo_log2 = dtl::log_2(d_lo);
//  u64 d_hi = dtl::next_power_of_two(dtl::env<$u64>::get("D_HI", 1ull << 28));
//  u64 d_hi_log2 = dtl::log_2(d_hi);
//  //===----------------------------------------------------------------------===//
//
//  //===----------------------------------------------------------------------===//
//  // Clustering factor
//  u64 f_lo = dtl::next_power_of_two(dtl::env<$u64>::get("F_LO", 1ull));
//  u64 f_lo_log2 = dtl::log_2(d_lo);
//  u64 f_hi = dtl::next_power_of_two(dtl::env<$u64>::get("F_HI", 1ull << 28));
//  u64 f_hi_log2 = dtl::log_2(d_hi);
//  //===----------------------------------------------------------------------===//
//
//  // all valid n's
//  const std::set<$u64> n_s = [&]() {
//    std::set<$u64> n_s;
//
//    for ($u64 d_log2 = d_lo_log2; d_log2 <= d_hi_log2; d_log2++) {
//      const std::vector<$f64> exp {
//          d_log2 +  0 * 0.0625,
//          d_log2 +  1 * 0.0625,
//          d_log2 +  2 * 0.0625,
//          d_log2 +  3 * 0.0625,
//          d_log2 +  4 * 0.0625,
//          d_log2 +  5 * 0.0625,
//          d_log2 +  6 * 0.0625,
//          d_log2 +  7 * 0.0625,
//          d_log2 +  8 * 0.0625,
//          d_log2 +  9 * 0.0625,
//          d_log2 + 10 * 0.0625,
//          d_log2 + 11 * 0.0625,
//          d_log2 + 12 * 0.0625,
//          d_log2 + 13 * 0.0625,
//          d_log2 + 14 * 0.0625,
//          d_log2 + 15 * 0.0625,
//          d_log2 + 16 * 0.0625,
//      };
//
//      for (auto e : exp) {
//        u64 n = std::pow(2, e);
//        if ((n * b_lo) > m_hi) continue; // make sure to not exceed the max filter size
//        n_s.insert(n);
//      }
//    }
//    return n_s;
//  }();
//
//  std::cout << "---------------------" << std::endl;
//  std::cout << "d's:" << std::endl;
//  for (auto d : n_s) {
//    std::cout << d << " " << std::endl;
//  }
//  std::cout << "---------------------" << std::endl;
//


  std::cerr << "run_id=" << RUN_ID << std::endl;
  std::vector<$f64> clustering_factors;
  std::cerr << "f:" << std::endl;
  for ($f64 f = 1; f <= N; f *= 2) {
    clustering_factors.push_back(f);
    std::cerr << f << std::endl;
  }

  std::vector<$f64> bit_densities;
  std::cerr << "d:" << std::endl;
  for ($f64 d = 1; d <= 10000; d *= 1.25) {
    bit_densities.push_back(d/10000);
    std::cerr << (d/10000) << std::endl;
  }

  std::vector<config> configs;
  for (auto f: clustering_factors) {
    for (auto d: bit_densities) {

      if (f > N*d) continue;

      config c;
      c.clustering_factor = f;
      c.density = d;

      c.bitmap_type = bitmap_t::bitmap;
      configs.push_back(c);
      c.bitmap_type = bitmap_t::roaring;
      configs.push_back(c);
      c.bitmap_type = bitmap_t::tree_mask_lo;
      configs.push_back(c);
      c.bitmap_type = bitmap_t::tree_mask_po;
      configs.push_back(c);
      c.bitmap_type = bitmap_t::wah;
      configs.push_back(c);

    }
  }

  {
    // # of runs
    std::vector<config> c(configs.begin(), configs.end());
    for ($u64 r = 1; r < RUNS; r++) {
      c.insert(c.end(), configs.begin(), configs.end());
    }
    std::swap(c, configs);
  }

  {
    // Shuffle the configurations to better predict the overall runtime of the benchmark.
    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(configs.begin(), configs.end(), gen);
  }

  std::function<void(const config&, std::ostream&)> fn = [](const config c, std::ostream& os) -> void {
    run(c, os);
  };
  dispatch(configs, fn);

}
