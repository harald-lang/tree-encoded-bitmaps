#include "common.hpp"
#include "experiments/util/bitmap_db.hpp"
#include "experiments/util/bitmap_types.hpp"
#include "experiments/util/differential_bitmap_types.hpp"
#include "experiments/util/gen.hpp"
#include "experiments/util/params.hpp"
#include "experiments/util/prep_data.hpp"
#include "experiments/util/prep_updates.hpp"

#include <dtl/dtl.hpp>

#include <dtl/bitmap/bitwise_operations.hpp>
#include <dtl/bitmap/diff/diff.hpp>
#include <dtl/bitmap/diff/merge.hpp>
#include <dtl/bitmap/diff/merge_teb.hpp>
#include <iostream>
#include <random>
#include <vector>
//===----------------------------------------------------------------------===//
// Experiment: Measures the update performance for varying number of merge
//             thresholds.
//===----------------------------------------------------------------------===//
struct task {
//  diff_bitmap_t bitmap_type; // under test
//  diff_merge_t merge_type;
  $u64 bitmap_id_a;
  $u64 bitmap_id_b;
  $f64 update_threshold;
  $f64 merge_threshold;
};
//===----------------------------------------------------------------------===//
struct benchmark_results {
  std::size_t n;
  std::string type_name;
  std::string type_info_before_updates;
  std::string type_info_after_updates;
  std::string merge_type_name;
  std::size_t encoded_size_before_updates;
  std::size_t encoded_size_after_updates;
  std::size_t update_cnt;
  std::size_t merge_cnt;
  $u64 runtime_nanos = 0;
  $u64 runtime_cycles = 0;
  $u64 runtime_update_nanos = 0;
  $u64 runtime_update_cycles = 0;
  $u64 runtime_merge_nanos = 0;
  $u64 runtime_merge_cycles = 0;
  // Don't care.
  std::size_t checksum = 0;
};
//===----------------------------------------------------------------------===//
template<typename T, typename M>
benchmark_results __attribute__((noinline))
run_benchmark(
    const dtl::bitmap b,
    const std::vector<update_entry>& updates,
    $u64 pending_updates_threshold) {

  const std::size_t total_update_cnt = updates.size();
  pending_updates_threshold = std::min(pending_updates_threshold, updates.size());

  const auto duration_nanos = 1250e6;
  const std::size_t MIN_REPS = 5;

  benchmark_results res;
  res.n = b.size();
  {
    T enc(b);
    res.encoded_size_before_updates = enc.size_in_byte();
    res.type_name = enc.name();
    res.type_info_before_updates = enc.info();
    res.merge_type_name = M::name();
  }

  // The actual measurement.
  res.runtime_nanos = 0;
  res.runtime_cycles = 0;
  res.runtime_update_nanos = 0;
  res.runtime_update_cycles = 0;
  res.runtime_merge_nanos = 0;
  res.runtime_merge_cycles = 0;

  // The repetition loop.
  const auto nanos_begin = now_nanos();
  const auto nanos_term = nanos_begin + duration_nanos;
  std::size_t rep_cntr = 0;
  std::size_t checksum = 0;
  while (now_nanos() < nanos_term || rep_cntr < MIN_REPS) {
    ++rep_cntr;

    // Encode the bitmap.
    T enc(b);

    std::size_t applied_update_cnt = 0;
    std::size_t applied_merge_cnt = 0;
    while (applied_update_cnt < total_update_cnt) {

      // The number of updates before the next merge happens.
      const std::size_t batch_size =
          std::min(pending_updates_threshold, total_update_cnt - applied_update_cnt);

      const std::size_t batch_begin = applied_update_cnt;
      const std::size_t batch_end = applied_update_cnt + batch_size;

      const auto nanos_begin_updates = now_nanos();
      const auto cycles_begin_updates = _rdtsc();

      // Apply the updates.
      for (std::size_t i = batch_begin; i < batch_end; ++i) {
        const auto& u = updates[i];
        enc.set(u.pos, u.value);
      }

      const auto cycles_end_updates = _rdtsc();
      const auto nanos_end_updates = now_nanos();

      res.runtime_update_nanos += nanos_end_updates - nanos_begin_updates;
      res.runtime_update_cycles += cycles_end_updates - cycles_begin_updates;
      applied_update_cnt += batch_size;

      if (batch_size == pending_updates_threshold) {
        // Merge the updates.
        const auto nanos_begin_merge = now_nanos();
        const auto cycles_begin_merge = _rdtsc();

        enc.template merge<M>();

        const auto cycles_end_merge = _rdtsc();
        const auto nanos_end_merge = now_nanos();

        res.runtime_merge_nanos += nanos_end_merge - nanos_begin_merge;
        res.runtime_merge_cycles += cycles_end_merge - cycles_begin_merge;
        ++applied_merge_cnt;
      }
    }
    res.encoded_size_after_updates = enc.size_in_byte();
    checksum += res.encoded_size_after_updates;
    res.type_info_after_updates = enc.info();
    res.update_cnt = applied_update_cnt;
    res.merge_cnt = applied_merge_cnt;
  }
  res.checksum = checksum;

  res.runtime_update_nanos /= rep_cntr;
  res.runtime_update_cycles /= rep_cntr;
  res.runtime_merge_nanos /= rep_cntr;
  res.runtime_merge_cycles /= rep_cntr;

  res.runtime_nanos = res.runtime_update_nanos + res.runtime_merge_nanos;
  res.runtime_cycles = res.runtime_update_cycles + res.runtime_merge_cycles;

  return res;
}
//===----------------------------------------------------------------------===//
template<typename T, typename M>
void
do_measurement(task& t) {
  // The bitmap that is used as starting point.
  const auto bm_a = db.load_bitmap(t.bitmap_id_a);
  const auto bm_a_density = dtl::determine_bit_density(bm_a);
  const auto bm_a_clustering_factor = dtl::determine_clustering_factor(bm_a);
  // The second bitmap determines the updates to perform.
  const auto bm_b = db.load_bitmap(t.bitmap_id_b);
  const auto bm_b_density = dtl::determine_bit_density(bm_b);
  const auto bm_b_clustering_factor = dtl::determine_clustering_factor(bm_b);

  assert(bm_a.size() == bm_b.size());

  // Prepare the updates. The updates are clustered, so that the characteristics
  // of the bitmap remains roughly the same during modifications.
  auto range_updates = prepare_range_updates(bm_a, bm_b);

  // Shuffle the update order.
  std::mt19937 gen(42); // for reproducible results, same order among all benchmarks
  std::shuffle(range_updates.begin(), range_updates.end(), gen);

  // The maximum number of updates to perform.
  const std::size_t max_updates = t.update_threshold;

  std::vector<update_entry> updates;
  updates.reserve(bm_b.count());
  for (auto& range : range_updates) {
    for (std::size_t i = range.pos; i < range.pos + range.length; ++i) {
      if (updates.size() < max_updates) {
        updates.emplace_back(i, range.value);
      }
    }
  }

  auto res = run_benchmark<T, M>(bm_a, updates, t.merge_threshold);

  std::string type_info_before_updates = res.type_info_before_updates;
  boost::replace_all(type_info_before_updates, "\"", "\"\""); // Escape JSON for CSV output.
  std::string type_info_after_updates = res.type_info_after_updates;
  boost::replace_all(type_info_after_updates, "\"", "\"\""); // Escape JSON for CSV output.

  std::cout << RUN_ID
      << ",\"" << BUILD_ID << "\""
      << "," << res.n
      << "," << "\"" << res.type_name << "\""
      << "," << "\"" << res.merge_type_name << "\""
      << "," << t.bitmap_id_a
      << "," << t.bitmap_id_b
      << "," << (res.runtime_nanos / res.update_cnt)
      << "," << (res.runtime_nanos / res.merge_cnt)
      << "," << std::min(std::size_t(t.update_threshold), updates.size())
      << "," << t.merge_threshold
      << "," << res.update_cnt
      << "," << res.merge_cnt
      << "," << res.runtime_nanos
//      << "," << res.runtime_cycles
      << "," << res.runtime_update_nanos
//      << "," << res.runtime_update_cycles
      << "," << res.runtime_merge_nanos
//      << "," << res.runtime_merge_cycles
      << "," << dtl::determine_bit_density(bm_a)
      << "," << dtl::determine_clustering_factor(bm_a)
      << "," << dtl::determine_bit_density(bm_b)
      << "," << dtl::determine_clustering_factor(bm_b)
      << "," << res.encoded_size_before_updates
      << "," << res.encoded_size_after_updates
      << "," << "\"" << type_info_before_updates << "\""
      << "," << "\"" << type_info_after_updates << "\""
      << "," << res.checksum
      << std::endl;
}
//===----------------------------------------------------------------------===//
template<diff_bitmap_t T, diff_merge_t M>
void
do_measurement(task& t) {
  for (std::size_t mt : {
      1000, 2500, 5000, 7500, 10000, 12500, 15000, 17500, 20000, 22500, 25000,
      27500, 30000, 32500, 35000, 37500, 40000, 42500, 45000, 47500, 50000}) {
    if (t.update_threshold < mt) break;
    t.merge_threshold = mt;
    do_measurement<
        typename diff_type_of<T>::type,
        typename diff_merge_type_of<T, M>::type>(t);
  }
}
//===----------------------------------------------------------------------===//
$i32 main() {
  // Prepare benchmark settings.
  u64 n_min = 1ull << 20;
  u64 n_max = 1ull << 20;

  // Use the same setting as with the compression skyline.
  std::cerr << "run_id=" << RUN_ID << std::endl;
  std::cerr << "build_id=" << BUILD_ID << std::endl;
  std::vector<$f64> clustering_factors;
  clustering_factors.push_back(8.0);

  std::vector<$f64> bit_densities;
//  bit_densities.push_back(0.01);
  bit_densities.push_back(0.1);
//  bit_densities.push_back(0.25);

  std::vector<$u64> n_values;
  for ($u64 n = n_min; n <= n_max; n <<= 1) {
    n_values.push_back(n);
  }

  if (GEN_DATA) {
    std::vector<params_markov> params;
    for (auto f : clustering_factors) {
      for (auto d : bit_densities) {
        for (auto n : n_values) {
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
                   "database."
                << std::endl;
      std::exit(1);
    }
  }

  for (auto f : clustering_factors) {
    for (auto d : bit_densities) {
      for (auto n : n_values) {
        if (!markov_parameters_are_valid(n, f, d)) continue;

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
        if (bitmap_ids.size() >= 2) {
          task t;
          t.bitmap_id_a = bitmap_ids[0];
          t.bitmap_id_b = bitmap_ids[1];
          const auto bm_b = db.load_bitmap(t.bitmap_id_b);
          t.update_threshold = std::min(std::size_t(100000), bm_b.count());
          std::cout<< t.update_threshold << std::endl;
//          do_measurement<diff_bitmap_t::part_teb,        diff_merge_t::no>(t); // DON'T! VERY SLOW
//          do_measurement<diff_bitmap_t::part_wah,        diff_merge_t::no>(t); // DON'T! VERY SLOW
          do_measurement<diff_bitmap_t::part_diff_teb,   diff_merge_t::naive>(t);
          do_measurement<diff_bitmap_t::part_diff_teb,   diff_merge_t::naive_iter>(t);
          do_measurement<diff_bitmap_t::part_diff_teb,   diff_merge_t::tree>(t);

          do_measurement<diff_bitmap_t::part_diff_wah,   diff_merge_t::naive>(t);
          do_measurement<diff_bitmap_t::part_diff_wah,   diff_merge_t::naive_iter>(t);

          do_measurement<diff_bitmap_t::roaring_roaring, diff_merge_t::naive>(t);
          do_measurement<diff_bitmap_t::roaring_roaring, diff_merge_t::naive_iter>(t);
          do_measurement<diff_bitmap_t::roaring_roaring, diff_merge_t::inplace>(t);

          do_measurement<diff_bitmap_t::teb_roaring,     diff_merge_t::naive>(t);
          do_measurement<diff_bitmap_t::teb_roaring,     diff_merge_t::naive_iter>(t);
          do_measurement<diff_bitmap_t::teb_roaring,     diff_merge_t::tree>(t);

          // DON'T use WAH (without fence pointers or partitioning).
          // Recall, every point update cause a point lookup, with is in O(n)
          // WAH.
//          do_measurement<diff_bitmap_t::wah_roaring,     diff_merge_t::naive>(t);
//          do_measurement<diff_bitmap_t::wah_roaring,     diff_merge_t::naive_iter>(t);

          do_measurement<diff_bitmap_t::part_wah_roaring,  diff_merge_t::naive>(t);
          do_measurement<diff_bitmap_t::part_wah_roaring,  diff_merge_t::naive_iter>(t);

          // DON'T use WAH as diff structure! SLOW
//          do_measurement<diff_bitmap_t::roaring_wah,     diff_merge_t::naive>(t);
//          do_measurement<diff_bitmap_t::teb_wah,         diff_merge_t::naive>(t);
//          do_measurement<diff_bitmap_t::wah_wah,         diff_merge_t::naive>(t);
        }
      }
    }
  }
}
//===----------------------------------------------------------------------===//
