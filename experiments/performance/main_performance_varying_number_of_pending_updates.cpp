#include "common.hpp"
#include "experiments/util/bitmap_db.hpp"
#include "experiments/util/bitmap_types.hpp"
#include "experiments/util/gen.hpp"
#include "experiments/util/params.hpp"
#include "experiments/util/prep_data.hpp"
#include "experiments/util/prep_updates.hpp"

#include <dtl/dtl.hpp>

#include <dtl/bitmap/bitwise_operations.hpp>
#include <dtl/bitmap/diff/diff.hpp>
#include <dtl/bitmap/diff/merge.hpp>
#include <dtl/bitmap/diff/merge_teb.hpp>
#include <dtl/bitmap/part/part.hpp>
#include <dtl/bitmap/part/part_upforward.hpp>
#include <iostream>
#include <random>
#include <vector>
//===----------------------------------------------------------------------===//
// Experiment: Measures the read performance for varying number of pending
//             updates.
//===----------------------------------------------------------------------===//
struct benchmark_results {
  std::size_t n;
  std::string type_name;
  std::string type_info;
  std::size_t encoded_size;
  $u64 runtime_nanos_scan_it = 0;
  $u64 runtime_cycles_scan_it = 0;
  $u64 runtime_nanos_skip_it = 0;
  $u64 runtime_cycles_skip_it = 0;
  // Don't care.
  std::size_t checksum = 0;
};
//===----------------------------------------------------------------------===//
template<typename T>
benchmark_results __attribute__((noinline))
run_benchmark(const T& enc_bs, const dtl::bitmap expected, std::ostream& os) {
//  const auto duration_nanos = RUN_DURATION_NANOS;
  const auto duration_nanos = 1250e6;
  const std::size_t MIN_REPS = 10;

  benchmark_results res;
  res.n = expected.size();
  res.encoded_size = enc_bs.size_in_bytes();
  res.type_name = enc_bs.name();
  res.type_info = enc_bs.info();

  // Validation code.
  {
    const auto dec_bs = dtl::to_bitmap_using_iterator(enc_bs);
    if (expected != dec_bs) {
      std::cerr << "Validation failed. " << std::endl;
      std::cerr << "Failed to decode the bitmap " << enc_bs.info()
          << "." << std::endl;
      std::exit(1);
    }
  }

  // Warm up run
  res.checksum = 0;
  {
    std::size_t pos_sink = 0;
    std::size_t length_sink = 0;
    auto it = enc_bs.scan_it();
    while (!it.end()) {
      pos_sink += it.pos();
      length_sink += it.length();
      it.next();
    }

    res.checksum += pos_sink + length_sink;

    // Validation. (Fail fast)
    if ((length_sink != expected.count())) {
      std::cerr << "Validation failed (during warm-up run): " << std::endl;
      std::cerr << "Expected length to be " << expected.count()
          << " but got " << length_sink << std::endl;
      std::exit(1);
    }
  }

  // The actual measurement.
  res.runtime_nanos_scan_it = 0;
  res.runtime_cycles_scan_it = 0;
  res.runtime_nanos_skip_it = 0;
  res.runtime_cycles_skip_it = 0;

  {
    std::size_t pos_sink = 0;
    std::size_t length_sink = 0;
    const auto nanos_begin = now_nanos();
    const auto nanos_term = nanos_begin + duration_nanos;
    const auto tsc_begin = _rdtsc();
    std::size_t rep_cntr = 0;
    while (now_nanos() < nanos_term || rep_cntr < MIN_REPS) {
      ++rep_cntr;
      auto it = enc_bs.scan_it();
      while (!it.end()) {
        pos_sink += it.pos();
        length_sink += it.length();
        it.next();
      }
    }
    const auto tsc_end = _rdtsc();
    const auto nanos_end = now_nanos();

    if ((length_sink / rep_cntr) != expected.count()) {
      std::cerr << "Validation failed: " << std::endl;
      std::cerr << "Expected length to be " << expected.count()
          << " but got " << (length_sink / rep_cntr) << std::endl;
      std::exit(1);
    }

    res.runtime_nanos_scan_it = (nanos_end - nanos_begin) / rep_cntr;
    res.runtime_cycles_scan_it = (tsc_end - tsc_begin) / rep_cntr;
    res.checksum += pos_sink + length_sink;
  }

  using scan_iter_t = decltype(std::declval<T>().scan_it());
  using skip_iter_t = decltype(std::declval<T>().it());

  if (std::is_same<scan_iter_t, skip_iter_t>::value) {
    // Scan and skip iterator are of the same type, thus replicate the results.
    res.runtime_nanos_skip_it = res.runtime_nanos_scan_it;
    res.runtime_cycles_skip_it = res.runtime_cycles_scan_it;
  }
  else {
    // Re-run the experiment with the other iterator.
    std::size_t pos_sink = 0;
    std::size_t length_sink = 0;
    const auto nanos_begin = now_nanos();
    const auto tsc_begin = _rdtsc();
    std::size_t rep_cntr = 0;
    while (now_nanos() - nanos_begin < duration_nanos
        || rep_cntr < MIN_REPS) {
      ++rep_cntr;
      auto it = enc_bs.it();
      while (!it.end()) {
        pos_sink += it.pos();
        length_sink += it.length();
        it.next();
      }
    }
    const auto tsc_end = _rdtsc();
    const auto nanos_end = now_nanos();

    if ((length_sink / rep_cntr) != expected.count()) {
      std::cerr << "Validation failed: " << std::endl;
      std::cerr << "Expected length to be " << expected.count()
          << " but got " << (length_sink / rep_cntr) << std::endl;
      std::exit(1);
    }

    res.runtime_nanos_skip_it = (nanos_end - nanos_begin) / rep_cntr;
    res.runtime_cycles_skip_it = (tsc_end - tsc_begin) / rep_cntr;
    res.checksum += pos_sink + length_sink;
  }

  return res;
}
//===----------------------------------------------------------------------===//
//template<typename B, typename D>
template<typename T, typename Tnodiff>
void
do_measurement(uint64_t bid_a, uint64_t bid_b) {
  // The bitmap that is used as starting point.
  const auto bm_a = db.load_bitmap(bid_a);
  const auto bm_a_density = dtl::determine_bit_density(bm_a);
  const auto bm_a_clustering_factor = dtl::determine_clustering_factor(bm_a);
  // The second bitmap determines the updates to perform.
  const auto bm_b = db.load_bitmap(bid_b);
  const auto bm_b_density = dtl::determine_bit_density(bm_b);
  const auto bm_b_clustering_factor = dtl::determine_clustering_factor(bm_b);

  assert(bm_a.size() == bm_b.size());

  // Prepare the updates. The updates are clustered, so that the characteristics
  // of the bitmap remains roughly the same during modifications.
  auto range_updates = prepare_range_updates(bm_a, bm_b);

  // Shuffle the update order.
  std::mt19937 gen(42); // for reproducible results
  std::shuffle(range_updates.begin(), range_updates.end(), gen);

  // The maximum number of pending updates.
  const std::size_t max_updates = 20000;

  std::vector<update_entry> updates;
  updates.reserve(bm_b.count());
  for (auto& range : range_updates) {
    for (std::size_t i = range.pos; i < range.pos + range.length; ++i) {
      if (updates.size() < max_updates) {
        updates.emplace_back(i, range.value);
      }
    }
  }

  // The encoded bitmap.
  T b(bm_a);
  // A plain bitmap used for validation.
  auto bm_expected = bm_a;

  // The number of data points to collect.
  const std::size_t data_points_cnt = 1;

  const std::size_t update_cnt_per_iteration =
      updates.size() / data_points_cnt;


  std::size_t pending_update_cnt = 0;

  // Baseline measurement.
  {
    // - diff, without pending updates
    auto res = run_benchmark(b, bm_expected, std::cout);

    // - without a diff structure
    auto baseline_res = res;
    {
      Tnodiff no_diff(bm_expected);
      baseline_res = run_benchmark(no_diff, bm_expected, std::cout);
    }

    std::string type_info = res.type_info;
    boost::replace_all(type_info, "\"", "\"\""); // Escape JSON for CSV output.

    std::cout << RUN_ID
        << ",\"" << BUILD_ID << "\""
        << "," << res.n
        << "," << "\"" << b.name() << "\""
        << "," << bid_a
        << "," << bid_b
        << "," << pending_update_cnt
        << "," << res.runtime_nanos_skip_it
        << "," << res.runtime_nanos_scan_it
        << "," << baseline_res.runtime_nanos_skip_it
        << "," << baseline_res.runtime_nanos_scan_it
        << "," << dtl::determine_bit_density(bm_expected)
        << "," << dtl::determine_clustering_factor(bm_expected)
        << "," << res.encoded_size
        << "," << baseline_res.encoded_size
        << "," << "\"" << type_info << "\""
        << "," << res.checksum
        << std::endl;
  }

  for (std::size_t i = 0; i < updates.size(); ++i) {
    const auto& u = updates[i];
    b.set(u.pos, u.value);
    bm_expected[u.pos] = u.value;
    ++pending_update_cnt;

    if (pending_update_cnt % update_cnt_per_iteration != 0) {
      continue;
    }

    auto res = run_benchmark(b, bm_expected, std::cout);
    auto baseline_res = res;
    {
      // Repeat the measurement without a differential data structure.
      Tnodiff no_diff(bm_expected);
      baseline_res = run_benchmark(no_diff, bm_expected, std::cout);
    }

    std::string type_info = res.type_info;
    boost::replace_all(type_info, "\"", "\"\""); // Escape JSON for CSV output.

    std::cout << RUN_ID
        << ",\"" << BUILD_ID << "\""
        << "," << res.n
        << "," << "\"" << b.name() << "\""
        << "," << bid_a
        << "," << bid_b
        << "," << pending_update_cnt
        << "," << res.runtime_nanos_skip_it
        << "," << res.runtime_nanos_scan_it
        << "," << baseline_res.runtime_nanos_skip_it
        << "," << baseline_res.runtime_nanos_scan_it
        << "," << dtl::determine_bit_density(bm_expected)
        << "," << dtl::determine_clustering_factor(bm_expected)
        << "," << res.encoded_size
        << "," << baseline_res.encoded_size
        << "," << "\"" << type_info << "\""
        << "," << res.checksum
        << std::endl;
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
          const auto bid_a = bitmap_ids[0];
          const auto bid_b = bitmap_ids[1];
          auto bm_b = db.load_bitmap(bid_b);
          std::cerr << "popcount(bm_b)=" << bm_b.count() << std::endl;
          // The bitmap type to use as differential data structure.
          using D = dtl::dynamic_roaring_bitmap;
          // The partition size to use.
          static constexpr std::size_t P = 1ull << 16;
          // Aliases
          using T = dtl::teb_wrapper;
          using R = dtl::dynamic_roaring_bitmap;
          using W = dtl::dynamic_wah32;
          using pW = dtl::part<dtl::dynamic_wah32, P>;

          do_measurement<dtl::diff<R,  D>,                        R >(bid_a, bid_b);
          do_measurement<dtl::diff<T,  D>,                        T >(bid_a, bid_b);
          do_measurement<dtl::diff<W,  D>,                        W >(bid_a, bid_b);
          do_measurement<dtl::diff<pW, D>,                        pW>(bid_a, bid_b);
          do_measurement<dtl::part_upforward<dtl::diff<R, D>, P>, dtl::part<R,P>>(bid_a, bid_b);
          do_measurement<dtl::part_upforward<dtl::diff<T, D>, P>, dtl::part<T,P>>(bid_a, bid_b);
          do_measurement<dtl::part_upforward<dtl::diff<W, D>, P>, dtl::part<W,P>>(bid_a, bid_b);
        }
      }
    }
  }
}
//===----------------------------------------------------------------------===//
