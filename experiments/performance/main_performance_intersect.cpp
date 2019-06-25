#include <iostream>
#include <set>

#include <dtl/dtl.hpp>
#include <experiments/util/bitmap_db.hpp>
#include <experiments/util/gen.hpp>
#include "common.hpp"
#include "common_bitwise.hpp"

//===----------------------------------------------------------------------===//
// Experiment: TBD
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
/// Data generation.
void gen_data(const std::set<config>& configs) {
  // Prepare the random bitmaps.
  std::cout << "Preparing the data set." << std::endl;

  std::vector<config> missing_bitmaps;
  for (auto c : configs) {
    auto ids = db.find_bitmaps(c.n, c.clustering_factor, c.density);
    if (ids.size() < RUNS) {
      for (std::size_t i = ids.size(); i < RUNS; ++i) {
        missing_bitmaps.push_back(c);
      }
    }
  }

  if (!missing_bitmaps.empty()) {
    std::cout << "Generating " << missing_bitmaps.size()
              << " random bitmaps." << std::endl;

    {
      // Shuffle the configurations to better predict the overall runtime of the
      // benchmark.
      std::random_device rd;
      std::mt19937 gen(rd());
      std::shuffle(missing_bitmaps.begin(), missing_bitmaps.end(), gen);
    }

    std::atomic<std::size_t> failure_cntr {0};
    std::function<void(const config&, std::ostream&)> fn =
        [&](const config c, std::ostream& os) -> void {
          try {
            const auto b = gen_random_bitmap_markov(
                c.n, c.clustering_factor, c.density);
            const auto id =
                db.store_bitmap(c.n, c.clustering_factor, c.density, b);
            // Validation.
            const auto loaded = db.load_bitmap(id);
            if (b != loaded) {
              // Fatal!
              std::cerr << "Validation failed" << std::endl;
              std::exit(1);
            }
          }
          catch (std::exception& ex) {
            ++failure_cntr;
          }
        };
    dispatch(missing_bitmaps, fn);

    if (failure_cntr > 0) {
      std::cerr << "Failed to generate all required bitmaps. "
                << failure_cntr << " bitmaps are still missing."
                << std::endl;
    }
  }

  std::size_t pass = 2;
  while (true) {
    std::cout << "Preparing the data set. (pass " << pass << ")" << std::endl;
    std::vector<config> incomplete_bitmaps;
    for (auto c : configs) {
      auto ids = db.find_bitmaps(c.n, c.clustering_factor, c.density);
      if (ids.size() > 0 && ids.size() < RUNS) {
        for (std::size_t i = ids.size(); i < RUNS; ++i) {
          incomplete_bitmaps.push_back(c);
        }
      }
    }

    std::cout << incomplete_bitmaps.size() << " remaining." << std::endl;
    if (!incomplete_bitmaps.empty()) {
      std::cout << "Generating " << incomplete_bitmaps.size()
                << " random bitmaps. (pass " << pass << ")" << std::endl;

      {
        // Shuffle the configurations to better predict the overall runtime of the
        // benchmark.
        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(incomplete_bitmaps.begin(), incomplete_bitmaps.end(), gen);
      }

      std::atomic<std::size_t> failure_cntr {0};
      std::function<void(const config&, std::ostream&)> fn =
          [&](const config c, std::ostream& os) -> void {
            try {
              const auto b = gen_random_bitmap_markov(
                  c.n, c.clustering_factor, c.density);
              const auto id = db.store_bitmap(
                  c.n, c.clustering_factor, c.density, b);
              // Validation.
              const auto loaded = db.load_bitmap(id);
              if (b != loaded) {
                // Fatal!
                std::cerr << "Validation failed" << std::endl;
                std::exit(1);
              }
            }
            catch (std::exception& ex) {
              ++failure_cntr;
            }
          };
      dispatch(incomplete_bitmaps, fn);
      if (failure_cntr == 0) {
        break;
      }

    }
    else {
      break;
    }
    pass++;
  }
  std::cerr << "Done generating random bitmaps after "
            << pass << " passes." << std::endl;
}
//===----------------------------------------------------------------------===//
/// Data generation. (for all valid combinations of n, f, and d)
void gen_data(const std::vector<$u64>& n_values,
              const std::vector<$f64>& clustering_factors,
              const std::vector<$f64>& bit_densities) {
  // Prepare the random bitmaps.
  std::set<config> bitmaps;
  for (auto f: clustering_factors) {
    for (auto d: bit_densities) {
      for (auto n: n_values) {

        if (!markov_parameters_are_valid(n, f, d)) continue;

        auto ids = db.find_bitmaps(n, f, d);
        if (ids.size() < RUNS) {
          config c;
          c.n = n;
          c.clustering_factor = f;
          c.density = d;
          for (std::size_t i = ids.size(); i < RUNS; ++i) {
            bitmaps.insert(c);
          }
        }

      }
    }
  }
  gen_data(bitmaps);
}
//===----------------------------------------------------------------------===//
$i32 main() {
  std::cerr << "run_id=" << RUN_ID << std::endl;
  std::cerr << "build_id=" << BUILD_ID << std::endl;

  // Prepare benchmark settings.
  std::vector<config_pair> configs;

//  u64 n_min = 1ull << 16;
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
    std::set<config> required_bitmaps;
    for (auto pair : configs) {
      required_bitmaps.insert(pair.first());
      required_bitmaps.insert(pair.second());
    }
//    for (auto& c : required_bitmaps) {
//      std::cerr << c << std::endl;
//    }
    gen_data(required_bitmaps);
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
