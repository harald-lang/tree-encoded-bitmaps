#include "common.hpp"
#include "common_bitwise.hpp"
#include "experiments/util/bitmap_db.hpp"
#include "experiments/util/gen.hpp"
#include "experiments/util/prep_data.hpp"

#include <dtl/dtl.hpp>

#include <iostream>
#include <random>
#include <set>
#include <vector>
//===----------------------------------------------------------------------===//
// Experiment: Measure the intersection time of two bitmaps.
//             Setting 1: d1=0.01, f1=8, d2=VARYING, f2=4
//             Setting 2: d1=0.01, f1=8, d2=0.25, f2=VARYING
//===----------------------------------------------------------------------===//
$i32 main() {
  std::cerr << "run_id=" << RUN_ID << std::endl;
  std::cerr << "build_id=" << BUILD_ID << std::endl;

  // Prepare benchmark settings.
  std::vector<config_pair> configs;

  u64 n_min = 1ull << 20;
  u64 n_max = 1ull << 20;

  // Loop over n.
  for ($u64 n = n_min; n <= n_max; n <<= 1) {
    //===------------------------------------------------------------------===//
    // Fix d1, f1, f2, and vary d2.
    {
      config_pair c;
      c.n = n;
      c.density1 = 0.01;
      c.clustering_factor1 = 8.0;
      c.clustering_factor2 = 4.0;

      if (!markov_parameters_are_valid(c.n, c.clustering_factor1, c.density1)) {
        std::cerr
            << "The characteristics of the outer bitmap are invalid: "
            << c.first()
            << std::endl;
        std::exit(1);
      }

      c.density2 = 0.01;
      if (!markov_parameters_are_valid(c.n, c.clustering_factor2, c.density2)) {
        continue;
      }
      configs.push_back(c);
      for ($f64 d = 5; d <= 100; d += 5) {
        c.density2 = d / 100;
        if (!markov_parameters_are_valid(c.n, c.clustering_factor2, c.density2)) {
          continue;
        }
        configs.push_back(c);
      }
    }
    //===------------------------------------------------------------------===//

    //===------------------------------------------------------------------===//
    // Fix d1, f1, d2, and vary f2.
    {
      config_pair c;
      c.n = n;
      c.density1 = 0.01;
      c.clustering_factor1 = 8.0;
      c.density2 = 0.25;
      if (!markov_parameters_are_valid(c.n, c.clustering_factor1, c.density1)) {
        std::cerr
            << "The characteristics of the outer bitmap are invalid: "
            << c.first()
            << std::endl;
        std::exit(1);
      }

      c.clustering_factor2 = 1;
      if (!markov_parameters_are_valid(c.n, c.clustering_factor2, c.density2)) {
        continue;
      }
      for ($u64 f = 1; f <= 32; f += 1) {
        c.clustering_factor2 = f;
        if (!markov_parameters_are_valid(c.n, c.clustering_factor2, c.density2)) {
          continue;
        }
        configs.push_back(c);
      }
    }
    //===------------------------------------------------------------------===//
  }

  if (GEN_DATA) {
    std::vector<params_markov> params;
    {
      std::set<config> required_bitmaps;
      for (auto pair : configs) {
        required_bitmaps.insert(pair.first());
        required_bitmaps.insert(pair.second());
      }
      for (auto& c : required_bitmaps) {
        params_markov p;
        p.n = c.n;
        p.clustering_factor = c.clustering_factor;
        p.density = c.density;
        params.push_back(p);
      }
    }
    prep_data(params, RUNS, db);
    std::exit(0);
  }

  // The implementations under test.
  std::vector<bitmap_t> bitmap_types;
  for (auto bitmap_type = static_cast<int>(bitmap_t::_first);
       bitmap_type <= static_cast<int>(bitmap_t::_last);
       ++bitmap_type) {
    bitmap_types.push_back(static_cast<bitmap_t>(bitmap_type));
  }

  std::vector<config_pair> benchmark_configs;

  for (auto c : configs) {
    auto bitmap_ids1 = db.find_bitmaps(c.n, c.clustering_factor1, c.density1);
    auto bitmap_ids2 = db.find_bitmaps(c.n, c.clustering_factor2, c.density2);
    if (bitmap_ids1.size() < RUNS) {
      std::cerr << "There are only " << bitmap_ids1.size() << " prepared "
                << "bitmaps for the parameters n=" << c.n << ", f="
                << c.clustering_factor1 << ", d=" << c.density1 << ", but " << RUNS
                << " are required."
                << std::endl;
      continue;
    }
    if (bitmap_ids2.size() < RUNS) {
      std::cerr << "There are only " << bitmap_ids2.size() << " prepared "
                << "bitmaps for the parameters n=" << c.n << ", f="
                << c.clustering_factor2 << ", d=" << c.density2 << ", but " << RUNS
                << " are required."
                << std::endl;
      continue;
    }

    for (std::size_t i = 0; i < RUNS; ++i) {
      c.bitmap_id1 = bitmap_ids1[i];
      c.bitmap_id2 = bitmap_ids2[i];
      for (auto b : bitmap_types) {
        c.bitmap_type = b;
        benchmark_configs.push_back(c);
      }
    }
  }

  std::cerr << "Prepared " << benchmark_configs.size() << " benchmark(s)."
            << std::endl;

  {
    // Shuffle the configurations to better predict the overall runtime of the
    // benchmark.
    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(configs.begin(), configs.end(), gen);
  }

  // Run the actual benchmark.
  run_intersect(benchmark_configs);
}
//===----------------------------------------------------------------------===//
