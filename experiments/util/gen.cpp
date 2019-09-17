#include "gen.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include <dtl/dtl.hpp>

#include "bitmap_db.hpp"
#include "params.hpp"
#include "threading.hpp"
//===----------------------------------------------------------------------===//
std::size_t
gen(
    std::vector<params_markov>& params,
    bitmap_db& db // the database, where the generated bitmaps are stored
) {
  // Shuffle the configurations to better predict the overall runtime of the
  // benchmark.
  std::random_device rd;
  std::mt19937 gen(rd());
  std::shuffle(params.begin(), params.end(), gen);

  std::atomic<std::size_t> failure_cntr {0};
  std::function<void(const params_markov&, std::ostream&)> fn =
      [&](const params_markov c, std::ostream& os) -> void {
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
  dispatch(params, fn);
  return failure_cntr;
}
//===----------------------------------------------------------------------===//
std::size_t
gen(
    std::vector<params_uniform>& params,
    bitmap_db& db // the database, where the generated bitmaps are stored
) {
  // Shuffle the configurations to better predict the overall runtime of the
  // benchmark.
  std::random_device rd;
  std::mt19937 gen(rd());
  std::shuffle(params.begin(), params.end(), gen);

  std::atomic<std::size_t> failure_cntr {0};
  std::function<void(const params_uniform&, std::ostream&)> fn =
      [&](const params_uniform c, std::ostream& os) -> void {
        try {
          const auto b = gen_random_bitmap_uniform(c.n, c.density);
          const auto id = db.store_bitmap(c.n, c.density, b);
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
  dispatch(params, fn);
  return failure_cntr;
}
//===----------------------------------------------------------------------===//
