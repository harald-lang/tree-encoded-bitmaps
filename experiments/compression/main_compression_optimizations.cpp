#include <iostream>

#include <dtl/dtl.hpp>
#include <dtl/bitmap/position_list.hpp>
#include <dtl/bitmap/partitioned_position_list.hpp>
#include <dtl/bitmap/range_list.hpp>
#include <dtl/bitmap/partitioned_range_list.hpp>
#include <dtl/bitmap/util/random.hpp>
#include "common.hpp"

//===----------------------------------------------------------------------===//
// Experiment: Comparison of compression ratios for the different TEB
//             optimization levels.
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

void
run(const config_t& config, std::ostream& os) {
  $f64 size_wah = 0;
  $f64 size_roaring = 0;
  $f64 size_teb_o0 = 0;
  $f64 size_teb_o1 = 0;
  $f64 size_teb_o2 = 0;
  $f64 size_pl = 0;
  $f64 size_ppl_u8 = 0;
  $f64 size_ppl_u16 = 0;
  $f64 size_rl = 0;
  $f64 size_prl_u8 = 0;
  $f64 size_prl_u16 = 0;

  // # of runs
  for ($u64 r = 0; r < RUNS; r++) {

    auto bm = dtl::gen_random_bitmap_markov(config.n,
                                            config.clustering_factor,
                                            config.bit_density);

    dtl::dynamic_wah32 wah(bm);
    size_wah += wah.size_in_byte();
    dtl::dynamic_roaring_bitmap roaring(bm);
    size_roaring += roaring.size_in_byte();
    dtl::teb teb_o0(bm, 0);
    size_teb_o0 += teb_o0.size_in_byte();
    dtl::teb teb_o1(bm, 1);
    size_teb_o1 += teb_o1.size_in_byte();
    dtl::teb teb_o2(bm, 2);
    teb_o2.run_optimize();
    size_teb_o2 += teb_o2.size_in_byte();
    dtl::position_list<$u32> pl(bm);
    size_pl += pl.size_in_byte();
    dtl::partitioned_position_list<$u32, $u8> ppl_u8(bm);
    size_ppl_u8 += ppl_u8.size_in_byte();
    dtl::partitioned_position_list<$u32, $u16> ppl_u16(bm);
    size_ppl_u16 += ppl_u16.size_in_byte();
    dtl::range_list<$u32> rl(bm);
    size_rl += rl.size_in_byte();
    dtl::partitioned_range_list<$u32, $u8> prl_u8(bm);
    size_prl_u8 += prl_u8.size_in_byte();
    dtl::partitioned_range_list<$u32, $u16> prl_u16(bm);
    size_prl_u16 += prl_u16.size_in_byte();

    // Validation
    {
      if (bm != wah.to_bitset()) {
        std::cerr << "Validation failed. (WAH)" << std::endl;
        std::exit(1);
      }
      if (bm != roaring.to_bitset()) {
        std::cerr << "Validation failed. (Roaring)" << std::endl;
        std::exit(1);
      }
      if (bm != teb_o0.to_bitset()) {
        std::cerr << "Validation failed. (TEBo0)" << std::endl;
        std::exit(1);
      }
      if (bm != teb_o1.to_bitset()) {
        std::cerr << "Validation failed. (TEBo1)" << std::endl;
        std::exit(1);
      }
      if (bm != teb_o2.to_bitset()) {
        std::cerr << "Validation failed. (TEBo2)" << std::endl;
        std::cerr << "expected: " << bm << std::endl;
        std::cerr << "actual:   " << teb_o2.to_bitset() << std::endl;
        std::exit(1);
      }
      if (bm != pl.to_bitset()) {
        std::cerr << "Validation failed. (PL)" << std::endl;
        std::exit(1);
      }
      if (bm != ppl_u8.to_bitset()) {
        std::cerr << "Validation failed. (PPLu8)" << std::endl;
        std::exit(1);
      }
      if (bm != ppl_u16.to_bitset()) {
        std::cerr << "Validation failed. (PPLu16)" << std::endl;
        std::exit(1);
      }
      if (bm != rl.to_bitset()) {
        std::cerr << "Validation failed. (RL)" << std::endl;
        std::exit(1);
      }
      if (bm != prl_u8.to_bitset()) {
        std::cerr << "Validation failed. (PRLu8)" << std::endl;
        std::exit(1);
      }
      if (bm != prl_u16.to_bitset()) {
        std::cerr << "Validation failed. (PRLu16)" << std::endl;
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
      << "," << (size_roaring / RUNS) / 1024.0
      << "," << (size_wah / RUNS) / 1024.0
      << "," << (size_pl / RUNS) / 1024.0
      << "," << (size_ppl_u8 / RUNS) / 1024.0
      << "," << (size_ppl_u16 / RUNS) / 1024.0
      << "," << (size_rl / RUNS) / 1024.0
      << "," << (size_prl_u8 / RUNS) / 1024.0
      << "," << (size_prl_u16 / RUNS) / 1024.0
      << std::endl;
  os << out.str();
}

$i32 main() {

//  u64 N = 1ull << 20;
  std::vector<$u64> n_s {
      1ull << 10, 1ull << 12, 1ull << 14, 1ull << 16, 1ull << 18, 1ull << 20 };
//  std::vector<$f64> clustering_factors { 8 };
  std::vector<$f64> clustering_factors { 1, 2, 4, 8, 16, 32, 64, 128 };
  std::vector<$f64> bit_densities { 0.01 };
  for ($f64 d = 5; d <= 85; d += 5) {
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

  std::random_shuffle(configs.begin(), configs.end());

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
