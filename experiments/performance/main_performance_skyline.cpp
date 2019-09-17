#include <iostream>
#include <dtl/dtl.hpp>
#include <experiments/util/bitmap_db.hpp>
#include <experiments/util/gen.hpp>
#include <navin/experiments/util/prep_data.hpp>

#include "common.hpp"
//===----------------------------------------------------------------------===//
// Experiment: Exhaustive performance measurements, used to compute a skyline
//             of the fastest compression schemes.
//===----------------------------------------------------------------------===//
$i32 main() {

  // Prepare benchmark settings.
  u64 n_min = 1ull << 20;
  u64 n_max = 1ull << 20;

  // Use the same setting as with the compression skyline.
  std::cerr << "run_id=" << RUN_ID << std::endl;
  std::cerr << "build_id=" << BUILD_ID << std::endl;
  std::vector<$f64> clustering_factors;
  for ($f64 f = 1; f <= n_max; f *= 2) {
    clustering_factors.push_back(f);
  }

  std::vector<$f64> bit_densities;
  for ($f64 d = 1; d <= 10000; d *= 1.25) {
    bit_densities.push_back(d/10000);
  }

  std::vector<$u64> n_values;
  for ($u64 n = n_min; n <= n_max; n <<= 1) {
    n_values.push_back(n);
  }

  if (GEN_DATA) {
    std::vector<params_markov> params;
    for (auto f: clustering_factors) {
      for (auto d: bit_densities) {
        for (auto n: n_values) {
          if (!markov_parameters_are_valid(n, f, d)) continue;
          params_markov p;
          p.n = n;
          p.clustering_factor = f;
          p.density = d;
          params.push_back(p);
        }
      }
    }
    prep_data(params, RUNS, db);
    std::exit(0);
  }
  else {
    if (db.empty()) {
      std::cerr << "Bitmap database is empty. Use GEN_DATA=1 to populate the "
          "database." << std::endl;
      std::exit(1);
    }
  }

  // The implementations under test.
  std::vector<bitmap_t> bitmap_types;
  for (auto bitmap_type = static_cast<int>(bitmap_t::_first);
       bitmap_type <= static_cast<int>(bitmap_t::_last);
       ++bitmap_type) {
    bitmap_types.push_back(static_cast<bitmap_t>(bitmap_type));
  }

  std::vector<config> configs;
  for (auto f: clustering_factors) {
    for (auto d: bit_densities) {
      for (auto n: n_values) {
        if (!markov_parameters_are_valid(n, f, d)) continue;

        config c;
        c.n = n;
        c.clustering_factor = f;
        c.density = d;

        auto bitmap_ids = db.find_bitmaps(n, f, d);
        if (bitmap_ids.empty()) {
          continue;
        }
        if (bitmap_ids.size() < RUNS) {
          std::cerr << "There are only " << bitmap_ids.size() << " prepared "
              << "bitmaps for the parameters n=" << n << ", f=" << f
              << ", d=" << d << ", but " << RUNS << " are required."
              << std::endl;
        }
        for (auto bitmap_id : bitmap_ids) {
          c.bitmap_id = bitmap_id;
          for (auto b : bitmap_types) {
            c.bitmap_type = b;
            configs.push_back(c);
          }
        }
      }
    }
  }

  {
    // Shuffle the configurations to better predict the overall runtime of the
    // benchmark.
    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(configs.begin(), configs.end(), gen);
  }

  // Run the actual benchmark.
  run(configs);
}
//===----------------------------------------------------------------------===//
