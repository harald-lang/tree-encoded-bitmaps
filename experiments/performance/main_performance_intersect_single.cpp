#include <iostream>
#include <vector>

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
void gen_data(const std::vector<$u64>& n_values,
              const std::vector<$f64>& clustering_factors,
              const std::vector<$f64>& bit_densities) {
  // Prepare the random bitmaps.
  std::cout << "Preparing the data set." << std::endl;
  std::vector<config> missing_bitmaps;
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
            missing_bitmaps.push_back(c);
          }
        }

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
    for (auto f: clustering_factors) {
      for (auto d: bit_densities) {
        for (auto n: n_values) {

          if (!markov_parameters_are_valid(n, f, d)) continue;

          auto ids = db.find_bitmaps(n, f, d);
          if (ids.size() > 0 && ids.size() < RUNS) {
            config c;
            c.n = n;
            c.clustering_factor = f;
            c.density = d;
            for (std::size_t i = ids.size(); i < RUNS; ++i) {
              incomplete_bitmaps.push_back(c);
            }
          }

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
$i32 main() {

  // Prepare benchmark settings.
//  u64 n_min = 1ull << 10;
//  u64 n_min = 1ull << 14;
//  u64 n_max = 1ull << 19;
//  u64 n_min = 1ull << 20;
//  u64 n_max = 1ull << 20;

  u64 n = 1ull << dtl::env<$u64>::get("N_LOG2", 0);

  const auto f1 = dtl::env<$f64>::get("F1", 0.0);
  const auto f2 = dtl::env<$f64>::get("F2", 0.0);
  const auto d1 = dtl::env<$f64>::get("D1", 0.0);
  const auto d2 = dtl::env<$f64>::get("D2", 0.0);
  if (n == 1 || f1 == 0.0 || f2 == 0.0 || d1 == 0.0 || d2 == 0.0) {
    std::cerr << "Invalid arguments." << std::endl;
    std::exit(1);
  }

  // Use the same setting as with the compression skyline.
  std::cerr << "run_id=" << RUN_ID << std::endl;
  std::cerr << "build_id=" << BUILD_ID << std::endl;
  std::vector<$f64> clustering_factors;
  clustering_factors.push_back(f1);
  clustering_factors.push_back(f2);

  std::vector<$f64> bit_densities;
  bit_densities.push_back(d1);
  bit_densities.push_back(d2);

  std::vector<$u64> n_values;
  n_values.push_back(n);

//  if (GEN_DATA) {
    gen_data(n_values, clustering_factors, bit_densities);
//    std::exit(0);
//  }
//  else {
//    if (db.empty()) {
//      std::cerr << "Bitmap database is empty. Use GEN_DATA=1 to populate the "
//          "database." << std::endl;
//      std::exit(1);
//    }
//  }

  // The implementations under test.
  std::vector<bitmap_t> bitmap_types;
  for (auto bitmap_type = static_cast<int>(bitmap_t::_first);
       bitmap_type <= static_cast<int>(bitmap_t::_last);
       ++bitmap_type) {
    bitmap_types.push_back(static_cast<bitmap_t>(bitmap_type));
  }

  std::vector<config_pair> configs;

  config_pair c;
  c.n = n;
  c.clustering_factor1 = f1;
  c.clustering_factor2 = f2;
  c.density1 = d1;
  c.density2 = d2;

  auto bitmap_ids1 = db.find_bitmaps(n, f1, d1);
  auto bitmap_ids2 = db.find_bitmaps(n, f2, d2);
  if (bitmap_ids1.size() < RUNS) {
    std::cerr << "There are only " << bitmap_ids1.size() << " prepared "
        << "bitmaps for the parameters n=" << n << ", f=" << f1
        << ", d=" << d1 << ", but " << RUNS << " are required."
        << std::endl;
  }
  if (bitmap_ids2.size() < RUNS) {
    std::cerr << "There are only " << bitmap_ids2.size() << " prepared "
        << "bitmaps for the parameters n=" << n << ", f=" << f2
        << ", d=" << d2 << ", but " << RUNS << " are required."
        << std::endl;
  }
  for (auto bitmap_id1 : bitmap_ids1) {
    c.bitmap_id1 = bitmap_id1;
    for (auto bitmap_id2 : bitmap_ids2) {
      c.bitmap_id2 = bitmap_id2;
      for (auto b : bitmap_types) {
        c.bitmap_type = b;
        configs.push_back(c);
      }
      break;
    }
    break;
  }

//  {
//    // Shuffle the configurations to better predict the overall runtime of the
//    // benchmark.
//    std::random_device rd;
//    std::mt19937 gen(rd());
//    std::shuffle(configs.begin(), configs.end(), gen);
//  }

  // Run the actual benchmark.
  run_intersect(configs);
}
//===----------------------------------------------------------------------===//
