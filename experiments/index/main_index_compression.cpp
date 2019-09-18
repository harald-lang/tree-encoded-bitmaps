#include "experiments/util/bitmap_types.hpp"
#include "experiments/util/gen.hpp"
#include "experiments/util/seq_db.hpp"
#include "experiments/util/threading.hpp"
#include "version.h"

#include <dtl/bitmap.hpp>
#include <dtl/dtl.hpp>
#include <dtl/env.hpp>

#include <atomic>
#include <chrono>
#include <functional>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <string>
//===----------------------------------------------------------------------===//
// Experiment: This experiment constructs and compresses bitmap indexes from
//             randomly generated integer sequences. Thereby varying the
//             parameters n, c, and f, where c refers to the number of distinct
//             values.
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
// The number of independent runs.
static constexpr u64 RUNS = 10;
// Identifies the benchmark run.
static const i64 RUN_ID = std::chrono::duration_cast<std::chrono::seconds>(
    std::chrono::system_clock::now().time_since_epoch()).count();
// The data set.
static const std::string DB_FILE =
    dtl::env<std::string>::get("DB_FILE", "./random_sequences.sqlite3");
static seq_db db(DB_FILE);

static u64 GEN_DATA = dtl::env<$u64>::get("GEN_DATA", 0);

//===----------------------------------------------------------------------===//
// An instance of config contains all information to run a single measurement.
struct config {
  /// Sequence ID (in the database).
  $u64 seq_id;
  /// The bitmap type.
  bitmap_t bitmap_type;
  /// Total number of elements.
  $u64 n;
  /// Attribute cardinality (# of distinct values).
  $u32 c;
  /// Clustering factor.
  $f64 f;
};
//===----------------------------------------------------------------------===//
// Forward declaration.
void bmi_benchmark(const config& conf, std::ostream& os);
//===----------------------------------------------------------------------===//
/// Data generation.
void gen_data(const std::vector<$u64>& n_values,
              const std::vector<$u32>& cardinalities,
              const std::vector<$f64>& clustering_factors) {

  // Prepare the random integer sequences.
  std::cout << "Preparing the data set." << std::endl;
  std::vector<config> missing_config;
  for (auto f: clustering_factors) {
    for (auto c: cardinalities) {
      for (auto n: n_values) {

        if (c * f >= n / 4) continue;

        auto ids = db.find(n, c, f);
        if (ids.size() < RUNS) {
          config conf;
          conf.n = n;
          conf.c = c;
          conf.f = f;
          for (std::size_t i = ids.size(); i < RUNS; ++i) {
            missing_config.push_back(conf);
          }
        }

      }
    }
  }

  if (!missing_config.empty()) {
    std::cout << "Generating " << missing_config.size()
              << " random integer sequences." << std::endl;

    {
      // Shuffle the configurations to better predict the overall runtime of the
      // benchmark.
      std::random_device rd;
      std::mt19937 gen(rd());
      std::shuffle(missing_config.begin(), missing_config.end(), gen);
    }

    std::atomic<std::size_t> failure_cntr {0};
    std::function<void(const config&, std::ostream&)> fn =
        [&](const config conf, std::ostream& os) -> void {
          try {
            const auto seq = gen_random_integer_sequence_markov(
                conf.n, conf.c, conf.f);
            const auto id =
                db.put(conf.n, conf.c, conf.f, seq);
            // Validation.
            const auto loaded = db.get(id);
            for (std::size_t i = 0; i < seq.size(); ++i) {
              if (seq[i] != loaded[i]) {
                // Fatal!
                std::cerr << "Validation failed" << std::endl;
                std::exit(1);
              }
            }
          }
          catch (std::exception& ex) {
            ++failure_cntr;
          }
        };
    dispatch(missing_config, fn);

    if (failure_cntr > 0) {
      std::cerr << "Failed to generate all required bitmaps. "
                << failure_cntr << " bitmaps are still missing."
                << std::endl;
    }
  }

  std::size_t pass = 2;
  while (true) {
    std::cout << "Preparing the data set. (pass " << pass << ")" << std::endl;
    std::vector<config> incomplete;
    for (auto f: clustering_factors) {
      for (auto c: cardinalities) {
        for (auto n: n_values) {

          if (c * f >= n / 4) continue;

          auto ids = db.find(n, c, f);
          if (!ids.empty() && ids.size() < RUNS) {
            config conf;
            conf.n = n;
            conf.c = c;
            conf.f = f;
            for (std::size_t i = ids.size(); i < RUNS; ++i) {
              incomplete.push_back(conf);
            }
          }

        }
      }
    }
    std::cout << incomplete.size() << " remaining." << std::endl;
    if (!incomplete.empty()) {
      std::cout << "Generating " << incomplete.size()
                << " random integer sequences. (pass " << pass << ")"
                << std::endl;

      {
        // Shuffle the configurations to better predict the overall runtime of the
        // benchmark.
        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(incomplete.begin(), incomplete.end(), gen);
      }

      std::atomic<std::size_t> failure_cntr {0};
      std::function<void(const config&, std::ostream&)> fn =
          [&](const config conf, std::ostream& os) -> void {
            try {
              const auto seq = gen_random_integer_sequence_markov(
                  conf.n, conf.c, conf.f);
              const auto id = db.put(
                  conf.n, conf.c, conf.f, seq);
              // Validation.
              const auto loaded = db.get(id);
              for (std::size_t i = 0; i < seq.size(); ++i) {
                if (seq[i] != loaded[i]) {
                  // Fatal!
                  std::cerr << "Validation failed" << std::endl;
                  std::exit(1);
                }
              }
            }
            catch (std::exception& ex) {
              ++failure_cntr;
            }
          };
      dispatch(incomplete, fn);
      if (failure_cntr == 0) {
        break;
      }

    }
    else {
      break;
    }
    pass++;
  }
  std::cerr << "Done generating random integer sequences after "
            << pass << " passes." << std::endl;
}
//===----------------------------------------------------------------------===//
$i32 main() {

  // Prepare benchmark settings.
  u64 n_min = 1ull << 10;
  u64 n_max = 1ull << 20;

  std::cerr << "run_id=" << RUN_ID << std::endl;
  std::vector<$f64> clustering_factors;
  for ($f64 f = 1; f <= 128; f *= 2) {
    clustering_factors.push_back(f);
  }

  std::vector<$u32> cardinalities;
  for ($u32 c = 8; c <= 8192; c *= 2) {
    cardinalities.push_back(c);
  }

  std::vector<$u64> n_values;
  for ($u64 n = n_min; n <= n_max; n <<= 1) {
    n_values.push_back(n);
  }

  if (GEN_DATA) {
    gen_data(n_values, cardinalities, clustering_factors);
    std::exit(0);
  }
  else {
    if (db.empty()) {
      std::cerr << "Integer sequence database is empty. Use GEN_DATA=1 to "
          "populate the database." << std::endl;
      std::exit(1);
    }
  }

  std::vector<config> configs;
  for (auto f: clustering_factors) {
    for (auto c: cardinalities) {
      for (auto n: n_values) {

        if (c * f >= n / 4) continue;

        config conf;
        conf.n = n;
        conf.c = c;
        conf.f = f;

        auto seq_ids = db.find(n, c, f);
        if (seq_ids.empty()) {
          continue;
        }
        if (seq_ids.size() < RUNS) {
          std::cerr << "There are only " << seq_ids.size() << " prepared "
              << "sequences for the parameters n=" << n << ", c=" << c
              << ", f=" << f << ", but " << RUNS << " are required."
              << std::endl;
        }
        for (auto seq_id : seq_ids) {
          conf.seq_id = seq_id;
          for (auto bitmap_type = static_cast<int>(bitmap_t::_first);
               bitmap_type <= static_cast<int>(bitmap_t::_last);
               ++bitmap_type) {
            conf.bitmap_type = static_cast<bitmap_t>(bitmap_type);
            configs.push_back(conf);
          }
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
  std::function<void(const config&, std::ostream&)> fn =
      [&](const config conf, std::ostream& os) -> void {
        const auto int_seq = db.get(conf.seq_id);
        bmi_benchmark(conf, os);
      };
  dispatch(configs, fn);
}

template<typename T>
void __attribute__((noinline))
run(const config& conf, std::ostream& os) {
  // Load the integer values.
  const auto values = db.get(conf.seq_id);
  if (values.size() != conf.n) {
    std::cerr << "Validation failed. Number of values does not match."
              << std::endl;
    std::exit(1);
  }

  // Construct the bitmap index.
  std::vector<dtl::bitmap*> bmi(conf.c, nullptr);
  for (std::size_t i = 0; i < conf.c; ++i) {
    bmi[i] = new dtl::bitmap(conf.n);
  }
  for (std::size_t i = 0; i < conf.n; ++i) {
    const auto v = values[i];
    (*bmi[v])[i] = true;
  }

  std::size_t valid_cntr = 0;
  for (std::size_t i = 0; i < conf.c; ++i) {
    valid_cntr += bmi[i]->count();
  }
  if (valid_cntr != conf.n) {
    std::cerr << "Validation failed. Number of set bits does not match."
              << std::endl;
    std::exit(1);
  }

  // Compress the bitmap index.
  std::size_t compressed_size = 0;
  for (auto bm : bmi) {
    T bmc(*bm);
    auto s = bmc.size_in_byte();
    compressed_size += s;
  }

  // Range encode the bitmap index.
  for (std::size_t i = 1; i < bmi.size(); ++i) {
    *bmi[i] = *bmi[i] | *bmi[i-1];
  }
  if (bmi.back()->count() != conf.n) {
    std::cout << "Validation failed. All bits in the last bitmap of a "
        "range-encoded bitmap index are supposed to be set." << std::endl;
    std::exit(1);
  }
  delete bmi[conf.c - 1];
  bmi.pop_back();

  // Compress the range-encoded bitmap index.
  std::size_t re_compressed_size = 0;
  for (auto bm : bmi) {
    T bmc(*bm);
    auto s = bmc.size_in_byte();
    re_compressed_size += s;
  }

  // Release memory.
  for (auto bm : bmi) {
    delete bm;
  }

  // Print results.
  os << RUN_ID
     << ",\"" << BUILD_ID << "\""
     << "," << conf.n
     << "," << conf.c
     << "," << conf.f
     << "," << conf.seq_id
     << "," << T::name()
     << "," << compressed_size
     << "," << re_compressed_size
     << std::endl;
}
//===----------------------------------------------------------------------===//
template<>
void __attribute__((noinline))
run<type_of<bitmap_t::bitmap>::type>(const config& conf, std::ostream& os) {
  // Print results.
  using T = type_of<bitmap_t::bitmap>::type;
  std::size_t size_in_bytes = (conf.n * conf.c) / 8;
  os << RUN_ID
     << ",\"" << BUILD_ID << "\""
     << "," << conf.n
     << "," << conf.c
     << "," << conf.f
     << "," << conf.seq_id
     << "," << T::name()
     << "," << size_in_bytes
     << "," << size_in_bytes
     << std::endl;
}
//===----------------------------------------------------------------------===//
void bmi_benchmark(const config& conf, std::ostream& os) {
  switch (conf.bitmap_type) {
    case bitmap_t::bitmap:
      run<type_of<bitmap_t::bitmap>::type>(conf, os);
      break;
    case bitmap_t::roaring:
      run<type_of<bitmap_t::roaring>::type>(conf, os);
      break;
    case bitmap_t::teb:
      run<type_of<bitmap_t::teb>::type>(conf, os);
      break;
    case bitmap_t::teb_scan:
//      run<type_of<bitmap_t::teb_scan>::type>(conf, os); // TODO remove
      break;
    case bitmap_t::wah:
      run<type_of<bitmap_t::wah>::type>(conf, os);
      break;
    case bitmap_t::position_list:
      run<type_of<bitmap_t::position_list>::type>(conf, os);
      break;
    case bitmap_t::partitioned_position_list_u8:
      run<type_of<bitmap_t::partitioned_position_list_u8>::type>(conf, os);
      break;
    case bitmap_t::partitioned_position_list_u16:
      run<type_of<bitmap_t::partitioned_position_list_u16>::type>(conf, os);
      break;
    case bitmap_t::range_list:
      run<type_of<bitmap_t::range_list>::type>(conf, os);
      break;
    case bitmap_t::partitioned_range_list_u8:
      run<type_of<bitmap_t::partitioned_range_list_u8>::type>(conf, os);
      break;
    case bitmap_t::partitioned_range_list_u16:
      run<type_of<bitmap_t::partitioned_range_list_u16>::type>(conf, os);
      break;
  }
}
//===----------------------------------------------------------------------===//
//$i32 main_off() {
//  config conf;
//  conf.n = 1ull << 10;
//  conf.c = 10;
//  conf.f = 4;
//
//  const auto seq = gen_random_integer_sequence_markov(
//      conf.n, conf.c, conf.f);
//  for (auto val : seq) {
//    std::cout << val << ", ";
//  }
//  std::cout << std::endl;
//  std::cout << "f_actual=" << dtl::determine_clustering_factor(seq) << std::endl;
//  const auto id =
//      db.put(conf.n, conf.c, conf.f, seq);
//  std::cout << "ID=" << id << std::endl;
//  // Validation.
//  const auto loaded = db.get(id);
//  for (auto val : loaded) {
//    std::cout << val << ", ";
//  }
//  std::cout << std::endl;
//  for (std::size_t i = 0; i < seq.size(); ++i) {
//    if (seq[i] != loaded[i]) {
//      // Fatal!
//      std::cerr << "Validation failed" << std::endl;
//      std::exit(1);
//    }
//  }
//  std::cout << "count=" << db.count() << std::endl;
//}
//===----------------------------------------------------------------------===//
