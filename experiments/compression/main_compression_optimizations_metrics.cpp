#include "common.hpp"
#include "experiments/util/gen.hpp"

#include <dtl/bitmap/partitioned_position_list.hpp>
#include <dtl/bitmap/partitioned_range_list.hpp>
#include <dtl/bitmap/position_list.hpp>
#include <dtl/bitmap/range_list.hpp>
#include <dtl/bitmap/teb.hpp>
#include <dtl/bitmap/teb_wrapper.hpp>
#include <dtl/bitmap/util/convert.hpp>
#include <dtl/bitmap/util/random.hpp>
#include <dtl/dtl.hpp>

#include <functional>
#include <iostream>
#include <sstream>
#include <vector>
//===----------------------------------------------------------------------===//
// Experiment: Comparison of "basic TEBs" and "space optimized TEBs".
//             The results are required to plot the Figures 6 and 7.
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
void
run(const config_t& config, std::ostream& os) {
  const auto num_runs = 10; // RUNS;
  $f64 size_teb_o0 = 0;
//  $f64 size_teb_o1 = 0; // deprecated. in the paper we only distinguish between -o0 and -o3
//  $f64 size_teb_o2 = 0;
  $f64 size_teb_o3 = 0;

  $u64 tree_bit_cnt_o0 = 0;
  $u64 label_bit_cnt_o0 = 0;
  $u64 tree_bit_cnt_o3 = 0;
  $u64 tree_leading_1bit_cnt_o3 = 0;
  $u64 tree_trailing_0bit_cnt_o3 = 0;
  $u64 label_bit_cnt_o3 = 0;
  $u64 label_leading_0bit_cnt = 0;
  $u64 label_trailing_0bit_cnt = 0;

  // # of runs
  for ($u64 r = 0; r < num_runs; r++) {

    auto bm = dtl::gen_random_bitmap_markov(config.n,
                                            config.clustering_factor,
                                            config.bit_density);

    dtl::teb<0> teb_o0(bm);
    size_teb_o0 += teb_o0.size_in_byte();
    tree_bit_cnt_o0 += teb_o0._get_tree_bit_cnt();
    label_bit_cnt_o0 += teb_o0._get_label_bit_cnt();

    dtl::teb<3> teb_o3(bm);
    size_teb_o3 += teb_o3.size_in_byte();
    tree_bit_cnt_o3 += teb_o3._get_tree_bit_cnt();
    tree_leading_1bit_cnt_o3 += teb_o3._get_tree_leading_1bit_cnt();
    tree_trailing_0bit_cnt_o3 += teb_o3._get_tree_trailing_0bit_cnt();
    label_bit_cnt_o3 += teb_o3._get_label_bit_cnt();
    label_leading_0bit_cnt += teb_o3._get_label_leading_0bit_cnt();
    label_trailing_0bit_cnt += teb_o3._get_label_trailing_0bit_cnt();
  }

  std::stringstream out;
  out << config.n / 8 / 1024.0
      << "," << config.bit_density
      << "," << config.clustering_factor
      << "," << (size_teb_o0 / num_runs) / 1024.0
      << "," << (1.0 * tree_bit_cnt_o0) / num_runs
      << "," << (1.0 * label_bit_cnt_o0) / num_runs
      << "," << (size_teb_o3 / num_runs) / 1024.0
      << "," << (1.0 * tree_bit_cnt_o3) / num_runs
      << "," << (1.0 * tree_leading_1bit_cnt_o3) / num_runs
      << "," << (1.0 * tree_trailing_0bit_cnt_o3) / num_runs
      << "," << (1.0 * label_bit_cnt_o3) / num_runs
      << "," << (1.0 * label_leading_0bit_cnt) / num_runs
      << "," << (1.0 * label_trailing_0bit_cnt) / num_runs
      << std::endl;
  os << out.str();
}
//===----------------------------------------------------------------------===//
$i32 main() {

  std::vector<$u64> n_s {
      1ull << 20 };
  std::vector<$f64> clustering_factors { 8 };

  std::vector<$f64> bit_densities { 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01 };
  for ($f64 d = 10; d <= 80; d += 10) {
    bit_densities.push_back(d/100);
  }

  std::vector<config_t> configs;
  for (auto n: n_s) {
    for (auto d: bit_densities) {
      for (auto f: clustering_factors) {
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
