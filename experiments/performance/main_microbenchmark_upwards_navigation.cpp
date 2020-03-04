#include "experiments/util/bitmap_db.hpp"
#include "experiments/util/gen.hpp"
#include "experiments/util/threading.hpp"
#include "thirdparty/perfevent/PerfEvent.hpp"

#include <dtl/bitmap/teb_legacy.hpp>
#include <dtl/bitmap/util/convert.hpp>
#include <dtl/dtl.hpp>
#include <dtl/env.hpp>

#include <boost/algorithm/string.hpp>

#include <chrono>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>
//===----------------------------------------------------------------------===//
// Micro-Experiment: Determine the costs of upward navigational steps.
//===----------------------------------------------------------------------===//
/// The database file where the bitmaps are stored.
static const std::string DB_FILE =
    dtl::env<std::string>::get("DB_FILE", "./random_bitmaps.sqlite3");
/// The database instance where the bitmaps are stored.
static bitmap_db db(DB_FILE);
/// The clustering factor.
f64 F = 1.0;
/// The bit densities we test. Note that with higher density, the number of
/// perfect levels increases and therefore the number of navigational steps
/// decreases. If the number of steps is less than 5, the results are
/// (significantly) distorted due to other overheads.
std::vector<$f64> bit_densities = { 0.0001, 0.001, 0.01, 0.1 };
/// The bitmap length.
u64 N = 1ull << 20;
//===----------------------------------------------------------------------===//
// Helper
auto now_nanos = []() {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(
      std::chrono::steady_clock::now().time_since_epoch())
      .count();
};
//===----------------------------------------------------------------------===//
void run(u64 n, f64 f, f64 d, i64 bitmap_id) {
  const auto bid = bitmap_id;
  const auto plain_bitmap = db.load_bitmap(bid);

  // We need to use an unoptimized TEB, otherwise the stack would be too small
  // to get reliable results.
  using T = dtl::teb<0>;

  T enc_bitmap(plain_bitmap);

  std::vector<std::size_t> probe_positions;
  std::size_t up_step_sum = 0;
  {
    auto it = enc_bitmap.it();
    while (!it.end()) {
      if (probe_positions.size() > 1024) break;

      auto up_steps = it.bench_nav_upwards_get_stack_size();
      if (up_steps >= 5) {
        probe_positions.push_back(it.pos());
        up_step_sum += up_steps;
      }

      it.next();
    }
  }
  if (up_step_sum < 1024) {
    std::cerr << "The given bitmap is not suitable for this experiment." << std::endl;
    std::cerr << enc_bitmap.info() << std::endl;
    std::cerr << " # of upward steps = " << up_step_sum
              << ", # of probe positions = " << probe_positions.size()
              << std::endl;
    return;
  }

  auto it = enc_bitmap.it();
  auto tsc_cntr = 0ull;
  for (std::size_t i = 0; i < probe_positions.size(); ++i) {
    std::size_t to_pos = probe_positions[i];
    it.nav_from_root_to(to_pos);
    const auto tsc_begin = _rdtsc();
    it.bench_nav_upwards(1, 0);
    const auto tsc_end = _rdtsc();
    tsc_cntr += tsc_end - tsc_begin;
  }
  std::cout << n << "," << f << "," << d
            << "," << (tsc_cntr * 1.0 / up_step_sum)
            << "," << enc_bitmap.info()
            << std::endl;
}
//===----------------------------------------------------------------------===//
$i32 main() {
  std::vector<$i64> bitmap_ids;

  for (auto d : bit_densities) {
    const auto ids = db.find_bitmaps(N, F, d);
    if (ids.empty()) {
      const auto bitmap = gen_random_bitmap_markov(N, F, d);
      const auto id = db.store_bitmap(N, F, d, bitmap);
      bitmap_ids.push_back(id);
    }
  }

  std::cerr << "n,f,d,cycles,info" << std::endl;
  for (auto d : bit_densities) {
    const auto ids = db.find_bitmaps(N, F, d);
    if (!ids.empty()) {
      run(N, F, d, ids[0]);
    }
  }
}
//===----------------------------------------------------------------------===//
