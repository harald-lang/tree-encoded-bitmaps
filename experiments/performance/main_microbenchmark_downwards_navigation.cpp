#include "experiments/util/bitmap_db.hpp"
#include "experiments/util/gen.hpp"
#include "experiments/util/threading.hpp"
#include "thirdparty/perfevent/PerfEvent.hpp"

#include <dtl/bitmap/teb.hpp>
#include <dtl/bitmap/teb_wrapper.hpp>
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
// Micro-Experiment: Determine the costs of downward navigational steps.
//===----------------------------------------------------------------------===//
/// The database file where the bitmaps are stored.
static const std::string DB_FILE =
    dtl::env<std::string>::get("DB_FILE", "./random_bitmaps.sqlite3");
/// The database instance where the bitmaps are stored.
bitmap_db db(DB_FILE);
/// The clustering factor needs to be 1, to ensure that we navigate downwards
/// the tree until the very last level.
f64 F = 1.0;
/// The bit densities we test. Note that with higher density, the number of
/// perfect levels increases and therefore the number of downward steps
/// decreases. If the number of downward steps is less than 5, the results
/// are (significantly) distorted due to the overhead in the
/// 'nav_from_root_to(pos)' function.
std::vector<$f64> bit_densities = {0.0001, 0.001, 0.01, 0.1};
/// The bitmap length.
u64 N = 1ull << 20;
/// Each measurement is repeated until the time below is elapsed.
static $u64 RUN_DURATION_NANOS = 1000e6; // 1s
//===----------------------------------------------------------------------===//
// Helper
auto now_nanos = []() {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(
      std::chrono::steady_clock::now().time_since_epoch()).count();
};
//===----------------------------------------------------------------------===//
void run(u64 n, f64 f, f64 d, i64 bitmap_id) {
  const auto bid = bitmap_id;
  const auto plain_bitmap = db.load_bitmap(bid);

  using T = dtl::teb_wrapper;
//  using T = dtl::teb<3>;

  T enc_bitmap(plain_bitmap);

  std::vector<std::size_t> probe_positions;
  std::size_t down_step_sum = 0;
  {
    auto it = enc_bitmap.it();
    while (!it.end()) {
      if (probe_positions.size() > 1024) break;
      probe_positions.push_back(it.pos());
      auto down_steps = dtl::teb<>::determine_tree_height(n) - (it.perfect_levels() - 1);
      down_step_sum += down_steps;
      it.next();
    }
  }
  if (down_step_sum == 0) return;

  const auto nanos_end = now_nanos() + RUN_DURATION_NANOS;

  auto it = enc_bitmap.it();
  std::size_t repeat_cnt = 0;
  PerfEvent e;
  e.startCounters();
  std::size_t chksum = 0;
  while (now_nanos() < nanos_end) {
    ++repeat_cnt;
    for (std::size_t i = 0; i < probe_positions.size(); ++i) {
      std::size_t to_pos = probe_positions[i];
      it.nav_from_root_to(to_pos);
      chksum += it.pos();
    }
  }
  e.stopCounters();
  down_step_sum *= repeat_cnt;

  std::string type_info = enc_bitmap.info();
  boost::replace_all(type_info, "\"", "\"\""); // Escape JSON for CSV output.

  std::cout << n << "," << f << "," << d
      << "," << (e.getCounter("cycles") / down_step_sum)
      << "," << (e.getCounter("instructions") / down_step_sum)
      << "," << (e.getCounter("branch-misses") / down_step_sum)
      << "," << e.getIPC()
      << "," << type_info
      << "," << chksum
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

  // CSV header
  std::cerr << "n,f,d,cycles,instructions,branch-misses,ipc,info,dontcare" << std::endl;
  for (auto d : bit_densities) {
    const auto ids = db.find_bitmaps(N, F, d);
    if (!ids.empty()) {
      run(N, F, d, ids[0]);
    }
  }

}
//===----------------------------------------------------------------------===//
