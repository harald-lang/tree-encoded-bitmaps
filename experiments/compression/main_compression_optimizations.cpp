#include "common.hpp"
#include "experiments/util/gen.hpp"

#include <dtl/bitmap/position_list.hpp>
#include <dtl/bitmap/range_list.hpp>
#include <dtl/bitmap/teb_legacy.hpp>
#include <dtl/bitmap/teb_wrapper.hpp>
#include <dtl/bitmap/util/convert.hpp>
#include <dtl/bitmap/util/random.hpp>
#include <dtl/dtl.hpp>

#include <functional>
#include <iostream>
#include <sstream>
#include <vector>
//===----------------------------------------------------------------------===//
// Experiment: Comparison of compression ratios for the different TEB
//             optimization levels.
//             The results are required to plot the Figures 2 and 5.
//===----------------------------------------------------------------------===//
struct config_t {
  $u64 n;
  $f64 bit_density;
  $f64 clustering_factor;

  void
  print(std::ostream& os) const {
    os << "config[n=" << n
       << ",d=" << bit_density
       << ",f=" << clustering_factor
       << "]";
  }
};
//===----------------------------------------------------------------------===//
void run(const config_t& config, std::ostream& os) {
  $f64 size_wah = 0;
  $f64 size_roaring = 0;
  $f64 size_teb_o0 = 0;
  $f64 size_teb_o1 = 0;
  $f64 size_teb_o2 = 0;
  $f64 size_teb_o3 = 0;

  // # of runs
  for ($u64 r = 0; r < RUNS; r++) {
    auto bm = gen_random_bitmap_markov(config.n,
        config.clustering_factor,
        config.bit_density);

    dtl::dynamic_wah32 wah(bm);
    size_wah += wah.size_in_bytes();
    dtl::dynamic_roaring_bitmap roaring(bm);
    size_roaring += roaring.size_in_bytes();
    dtl::teb<0> teb_o0(bm);
    size_teb_o0 += teb_o0.size_in_bytes();
    //    dtl::teb<1> teb_o1(bm); // deprecated. in the paper we only distinguish between -o0 and -o3
    //    size_teb_o1 += teb_o1.size_in_bytes();
    //    dtl::teb<2> teb_o2(bm);
    //    size_teb_o2 += teb_o2.size_in_bytes();
    dtl::teb<3> teb_o3(bm);
    size_teb_o3 += teb_o3.size_in_bytes();

    // Validation
    {
      if (bm != dtl::to_bitmap_using_iterator(wah)) {
        std::cerr << "Validation failed. (WAH)" << std::endl;
        std::exit(1);
      }
      if (bm != dtl::to_bitmap_using_iterator(roaring)) {
        std::cerr << "Validation failed. (Roaring)" << std::endl;
        std::exit(1);
      }
      if (bm != dtl::to_bitmap_using_iterator(teb_o0)) {
        std::cerr << "Validation failed. (TEBo0)" << std::endl;
        std::exit(1);
      }
      //      if (bm != dtl::to_bitmap_using_iterator(teb_o1)) {
      //        std::cerr << "Validation failed. (TEBo1)" << std::endl;
      //        std::exit(1);
      //      }
      //      if (bm != dtl::to_bitmap_using_iterator(teb_o2)) {
      //        std::cerr << "Validation failed. (TEBo2)" << std::endl;
      //        std::exit(1);
      //      }
      if (bm != dtl::to_bitmap_using_iterator(teb_o3)) {
        std::cerr << "Validation failed. (TEBo3)" << std::endl;
        std::exit(1);
      }
    }
  }

  std::stringstream out;
  out << config.n / 8 / 1024.0
      << "," << config.bit_density
      << "," << config.clustering_factor
      << "," << (size_teb_o0 / RUNS) / 1024.0
      << "," << (size_teb_o1 / RUNS) / 1024.0
      << "," << (size_teb_o2 / RUNS) / 1024.0
      << "," << (size_teb_o3 / RUNS) / 1024.0
      << "," << (size_roaring / RUNS) / 1024.0
      << "," << (size_wah / RUNS) / 1024.0
      << std::endl;
  os << out.str();
}
//===----------------------------------------------------------------------===//
$i32 main() {
  std::vector<$u64> n_s { 1ull << 20 };
  std::vector<$f64> clustering_factors { 8 };

  std::vector<$f64> bit_densities { 0.01 };
  for ($f64 d = 5; d <= 85; d += 5) {
    bit_densities.push_back(d / 100);
  }

  std::vector<config_t> configs;
  for (auto n : n_s) {
    for (auto d : bit_densities) {
      for (auto f : clustering_factors) {
        if (f > n * d) {
          std::stringstream err;
          err << "Skipping n=" << n << ", d=" << d << ", f=" << f << "."
              << std::endl;
          std::cerr << err.str();
          continue;
        }
        config_t c;
        c.n = n;
        c.bit_density = d;
        c.clustering_factor = f;
        configs.push_back(c);
      }
    }
  }

  std::function<void(const config_t&, std::ostream&)> fn =
      [](const config_t c, std::ostream& os) -> void {
    try {
      run(c, os);
    }
    catch (...) {
      std::stringstream err;
      err << "Failed to run " << c << "." << std::endl;
      std::cerr << err.str();
    }
  };
  dispatch<config_t>(configs, fn);
}
//===----------------------------------------------------------------------===//
