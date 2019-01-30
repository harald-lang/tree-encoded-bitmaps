#include <iostream>

#include <dtl/dtl.hpp>

#include "common.hpp"

//===----------------------------------------------------------------------===//
// Experiment: Comparison of compression ratios.
//===----------------------------------------------------------------------===//


$i32 main() {

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
      c.bitmap_type = bitmap_t::wah;
      configs.push_back(c);
      c.bitmap_type = bitmap_t::teb;
      configs.push_back(c);
      c.bitmap_type = bitmap_t::position_list;
      configs.push_back(c);
      c.bitmap_type = bitmap_t::partitioned_position_list_u8;
      configs.push_back(c);
      c.bitmap_type = bitmap_t::partitioned_position_list_u16;
      configs.push_back(c);
      c.bitmap_type = bitmap_t::range_list;
      configs.push_back(c);
      c.bitmap_type = bitmap_t::partitioned_range_list_u8;
      configs.push_back(c);
      c.bitmap_type = bitmap_t::partitioned_range_list_u16;
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
