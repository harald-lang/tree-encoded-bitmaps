#include <iostream>
#include <dtl/dtl.hpp>
#include <experiments/util/bitmap_db.hpp>
#include <experiments/util/gen.hpp>
#include "common.hpp"

//===----------------------------------------------------------------------===//
// Experiment: Comparison of compression ratios of (uniform) random bitmaps.
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
/// Data generation.
void gen_data_uniform(const std::vector<$u64>& n_values,
                      const std::vector<$f64>& bit_densities) {
  // Prepare the random bitmaps.
  std::cout << "Preparing the data set." << std::endl;
  std::vector<config> missing_bitmaps;
  for (auto d: bit_densities) {
    for (auto n: n_values) {

      if (d < 0 || d > 1.0) continue;

      auto ids = db.find_bitmaps(n, d);
      if (ids.size() < RUNS) {
        config c;
        c.n = n;
        c.clustering_factor = 0.0;
        c.density = d;
        for (std::size_t i = ids.size(); i < RUNS; ++i) {
          missing_bitmaps.push_back(c);
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
            const auto b = gen_random_bitmap_uniform(
                c.n, c.density);
            const auto id =
                db.store_bitmap(c.n, c.density, b);
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
  std::cerr << "Done generating random bitmaps." << std::endl;
}
//===----------------------------------------------------------------------===//
$i32 main() {

  // Prepare benchmark settings.
//  u64 n_min = 1ull << 10;
  u64 n_min = 1ull << 20;
  u64 n_max = 1ull << 20;

  std::cerr << "run_id=" << RUN_ID << std::endl;

  std::vector<$f64> bit_densities;
  // log scale
  for ($f64 d = 1; d <= 10000; d *= 1.25) {
    bit_densities.push_back(1.0 - d/10000);
  }
  // linear scale
  bit_densities.push_back(0.01);
  for ($f64 d = 5; d <= 100; d += 5) {
//    bit_densities.push_back(d/100);
  }

  std::vector<$u64> n_values;
  for ($u64 n = n_min; n <= n_max; n <<= 1) {
    n_values.push_back(n);
  }

  if (GEN_DATA) {
    gen_data_uniform(n_values, bit_densities);
    std::exit(0);
  }
  else {
    if (db.empty()) {
      std::cerr << "Bitmap database is empty. Use GEN_DATA=1 to populate the "
          "database." << std::endl;
      std::exit(1);
    }
  }

  std::vector<config> configs;
  for (auto d: bit_densities) {
    for (auto n: n_values) {

      if (d < 0 || d > 1.0) continue;

      config c;
      c.n = n;
      c.clustering_factor = 0.0;
      c.density = d;

      auto bitmap_ids = db.find_bitmaps(n, d);
      if (bitmap_ids.empty()) {
        continue;
      }
      if (bitmap_ids.size() < RUNS) {
        std::cerr << "There are only " << bitmap_ids.size() << " prepared "
            << "bitmaps for the parameters n=" << n
            << ", d=" << d << ", but " << RUNS << " are required."
            << std::endl;
      }
      for (auto bitmap_id : bitmap_ids) {
        c.bitmap_id = bitmap_id;
        for (auto bitmap_type = static_cast<int>(bitmap_t::_first);
             bitmap_type <= static_cast<int>(bitmap_t::_last);
             ++bitmap_type) {
          c.bitmap_type = static_cast<bitmap_t>(bitmap_type);
          configs.push_back(c);
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
