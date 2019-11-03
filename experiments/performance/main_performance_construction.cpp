#include "common.hpp"
#include "experiments/util/bitmap_db.hpp"
#include "experiments/util/gen.hpp"
#include "experiments/util/prep_data.hpp"

#include <dtl/bitmap/part/part.hpp>
#include <dtl/dtl.hpp>

#include <iostream>
//===----------------------------------------------------------------------===//
// Experiment: Measure the construction time (compression time) for various
//             bit densities.
//===----------------------------------------------------------------------===//
/// Data generation.
void gen_data(const std::vector<$u64>& n_values,
    const std::vector<$f64>& clustering_factors,
    const std::vector<$f64>& bit_densities) {
  // Prepare the random bitmaps.
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
//===----------------------------------------------------------------------===//
template<typename T>
void __attribute__((noinline))
run_construction_benchmark(const config& c, std::ostream& os) {
  const auto duration_nanos = RUN_DURATION_NANOS;
#ifndef NDEBUG
  const std::size_t MIN_REPS = 1;
#else
  const std::size_t MIN_REPS = 10;
#endif

  // Load the bitmap from DB.
  const auto bs = db.load_bitmap(c.bitmap_id);
  // Encode the bitmap.
  const T enc_bs(bs);
  // Validation code.
  {
    const auto dec_bs = dtl::to_bitmap_using_iterator(enc_bs);
    if (bs != dec_bs) {
      std::cerr << "Validation failed: " << c << std::endl;
      std::cerr << "Failed to decode the bitmap " << enc_bs.info()
                << "." << std::endl;
      std::exit(1);
    }
  }

  // The actual measurement.
  $u64 runtime_nanos = 0;
  $u64 runtime_cycles = 0;
  std::size_t checksum = 0;

  {
    const auto nanos_begin = now_nanos();
    const auto nanos_term = nanos_begin + duration_nanos;
    const auto tsc_begin = _rdtsc();
    std::size_t rep_cntr = 0;
    while (now_nanos() < nanos_term || rep_cntr < MIN_REPS) {
      ++rep_cntr;
      // Encode the bitmap.
      const T enc(bs);
      checksum += enc.size();
    }
    const auto tsc_end = _rdtsc();
    const auto nanos_end = now_nanos();

    runtime_nanos = (nanos_end - nanos_begin) / rep_cntr;
    runtime_cycles = (tsc_end - tsc_begin) / rep_cntr;
  }

  std::string type_info = enc_bs.info();
  boost::replace_all(type_info, "\"", "\"\""); // Escape JSON for CSV output.

  os << RUN_ID
     << ",\"" << BUILD_ID << "\""
     << "," << c.n
     << "," << "\"" << T::name() << "\""
     << "," << runtime_nanos
     << "," << runtime_cycles
     << "," << c.density
     << "," << dtl::determine_bit_density(bs)
     << "," << c.clustering_factor
     << "," << dtl::determine_clustering_factor(bs)
     << "," << c.bitmap_id
     << "," << enc_bs.size_in_byte()
     << ","
     << "\"" << type_info << "\""
     << "," << checksum
     << std::endl;
}
//===----------------------------------------------------------------------===//
$i32 main() {
  // Pin main thread to CPU 0.
  dtl::this_thread::set_cpu_affinity(0);

  // Prepare benchmark settings.
  u64 n_min = 1ull << 20;
  u64 n_max = 1ull << 20;

  // Use the same setting as with the compression skyline.
  std::cerr << "run_id=" << RUN_ID << std::endl;
  std::cerr << "build_id=" << BUILD_ID << std::endl;
  std::vector<$f64> clustering_factors;
  clustering_factors.push_back(8);

  std::vector<$f64> bit_densities;
  bit_densities.push_back(0.01);
  for ($f64 d = 5; d <= 85; d += 5) {
    bit_densities.push_back(d / 100);
  }
  std::vector<$u64> n_values;
  for ($u64 n = n_min; n <= n_max; n <<= 1) {
    n_values.push_back(n);
  }

  if (GEN_DATA) {
    gen_data(n_values, clustering_factors, bit_densities);
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

  std::vector<config> configs;
  for (auto f : clustering_factors) {
    for (auto d : bit_densities) {
      for (auto n : n_values) {
        std::cerr << "n=" << n << ", f=" << f << ", d=" << d << std::endl;
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
          c.bitmap_type = bitmap_t::bitmap; // Don't care.
          configs.push_back(c);
        }
      }
    }
  }

  // Run the actual benchmark.
  std::function<void(const config&, std::ostream&)> fn =
      [](const config c, std::ostream& os) -> void {
        run_construction_benchmark<dtl::dynamic_wah32>(c, os);
        run_construction_benchmark<dtl::dynamic_roaring_bitmap>(c, os);
        run_construction_benchmark<dtl::teb_wrapper>(c, os);
        run_construction_benchmark<dtl::part<dtl::teb_wrapper, 1ull << 16>>(c, os);
  };
  const auto thread_cnt = 1; // run performance measurements single-threaded
  dispatch(configs, fn, thread_cnt);
}
//===----------------------------------------------------------------------===//
