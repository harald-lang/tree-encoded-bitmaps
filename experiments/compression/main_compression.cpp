#include <iostream>

#include <dtl/dtl.hpp>

#include "common.hpp"

//===----------------------------------------------------------------------===//
// Experiment: Comparison of compression ratios for varying densities and
// fixed clustering factors.
//===----------------------------------------------------------------------===//

$i32 main() {

  std::cerr << "run_id=" << RUN_ID << std::endl;
  std::vector<$f64> clustering_factors { 8, 64, 128 };

  std::vector<$f64> bit_densities { 0.01 };
  std::cerr << "d:" << std::endl;
  for ($f64 d = 5; d <= 100; d += 5) {
    bit_densities.push_back(d/100);
    std::cerr << (d/100) << std::endl;
  }

  std::vector<config> configs;
  for (auto f: clustering_factors) {
    for (auto d: bit_densities) {

      if (f > N * d) continue;

      config c;
      c.clustering_factor = f;
      c.density = d;

      c.bitmap_type = bitmap_t::bitmap;
      configs.push_back(c);
      c.bitmap_type = bitmap_t::roaring;
      configs.push_back(c);
      c.bitmap_type = bitmap_t::teb;
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
