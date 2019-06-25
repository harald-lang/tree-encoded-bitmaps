#pragma once

#include <bitset>
#include <chrono>
#include <cmath>
#include <limits>
#include <utility>
#include <iostream>

#include <boost/algorithm/string.hpp>

#include <dtl/dtl.hpp>
#include <dtl/env.hpp>
#include <dtl/thread.hpp>

#include <dtl/bitmap/dynamic_bitmap.hpp>
#include <dtl/bitmap/dynamic_roaring_bitmap.hpp>
#include <dtl/bitmap/teb.hpp>
#include <dtl/bitmap/teb_scan.hpp>
#include <dtl/bitmap/dynamic_wah.hpp>
#include <dtl/bitmap/util/two_state_markov_process.hpp>
#include <dtl/bitmap/util/random.hpp>
#include <dtl/bitmap/util/convert.hpp>
#include <dtl/bitmap/position_list.hpp>
#include <dtl/bitmap/partitioned_position_list.hpp>
#include <dtl/bitmap/partitioned_range_list.hpp>

#include <dtl/bitmap/range_list.hpp>
#include <experiments/util/config.hpp>
#include <experiments/util/threading.hpp>
#include <experiments/util/bitmap_types.hpp>

#include "version.h"

//===----------------------------------------------------------------------===//
// The number of independent runs.
static $u64 RUNS = 10;
static constexpr u64 N = 1u << 20;

static const i64 RUN_ID = std::chrono::duration_cast<std::chrono::seconds>(
    std::chrono::system_clock::now().time_since_epoch()).count();

static const std::string DB_FILE =
    dtl::env<std::string>::get("DB_FILE", "./random_bitmaps.sqlite3");
static u64 GEN_DATA = dtl::env<$u64>::get("GEN_DATA", 0);

static bitmap_db db(DB_FILE);

//static $u64 RUN_DURATION_NANOS = 250e6; // run for at least 250ms
static $u64 RUN_DURATION_NANOS = 1250e6; // run for at least 250ms
//===----------------------------------------------------------------------===//
u64
now_nanos() {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(
      std::chrono::steady_clock::now().time_since_epoch()).count();
}
//===----------------------------------------------------------------------===//
template<typename T>
void __attribute__ ((noinline))
run(const config& c, std::ostream& os) {
  const auto duration_nanos = RUN_DURATION_NANOS;
  const std::size_t MIN_REPS = 10;
  // Load the bitmap from DB.
  const auto bs = db.load_bitmap(c.bitmap_id);
  const auto bs_count = bs.count();
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

  // Warm up run
  std::size_t checksum = 0;
  {
    std::size_t pos_sink = 0;
    std::size_t length_sink = 0;
    auto it = enc_bs.scan_it();
    while (!it.end()) {
      pos_sink += it.pos();
      length_sink += it.length();
      it.next();
    }

    checksum += pos_sink + length_sink;

    // Validation. (Fail fast)
    if ((length_sink != bs_count)) {
      std::cerr << "Validation failed (during warm-up run): " << c << std::endl;
      std::cerr << "Expected length to be " << bs.count()
                << " but got " << length_sink << std::endl;
      std::exit(1);
    }
  }


  // The actual measurement.
  $u64 runtime_nanos_scan_it = 0;
  $u64 runtime_cycles_scan_it = 0;
  $u64 runtime_nanos_skip_it = 0;
  $u64 runtime_cycles_skip_it = 0;

  {
    std::size_t pos_sink = 0;
    std::size_t length_sink = 0;
    const auto nanos_begin = now_nanos();
    const auto tsc_begin = _rdtsc();
    std::size_t rep_cntr = 0;
    while (now_nanos() - nanos_begin < duration_nanos
        || rep_cntr < MIN_REPS) {
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

    if ((length_sink / rep_cntr) != bs_count) {
      std::cerr << "Validation failed: " << c << std::endl;
      std::cerr << "Expected length to be " << bs_count
          << " but got " << (length_sink / rep_cntr) << std::endl;
      std::exit(1);
    }

    runtime_nanos_scan_it = (nanos_end - nanos_begin) / rep_cntr;
    runtime_cycles_scan_it = (tsc_end - tsc_begin) / rep_cntr;
    checksum += pos_sink + length_sink;
  }

  using scan_iter_t = decltype(std::declval<T>().scan_it());
  using skip_iter_t = decltype(std::declval<T>().it());

  if (std::is_same<scan_iter_t, skip_iter_t>::value) {
    // Scan and skip iterator are of the same type, thus replicate the results.
    runtime_nanos_skip_it = runtime_nanos_scan_it;
    runtime_cycles_skip_it = runtime_cycles_scan_it;
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

    if ((length_sink / rep_cntr) != bs_count) {
      std::cerr << "Validation failed: " << c << std::endl;
      std::cerr << "Expected length to be " << bs_count
          << " but got " << (length_sink / rep_cntr) << std::endl;
      std::exit(1);
    }

    runtime_nanos_skip_it = (nanos_end - nanos_begin) / rep_cntr;
    runtime_cycles_skip_it = (tsc_end - tsc_begin) / rep_cntr;
    checksum += pos_sink + length_sink;
  }

  std::string type_info = enc_bs.info();
  boost::replace_all(type_info, "\"", "\"\""); // Escape JSON for CSV output.

  os << RUN_ID
     << ",\"" << BUILD_ID << "\""
     << "," << c.n
     << "," << T::name()
     << "," << runtime_nanos_skip_it
     << "," << runtime_cycles_skip_it
     << "," << runtime_nanos_scan_it
     << "," << runtime_cycles_scan_it
     << "," << c.density
     << "," << dtl::determine_bit_density(bs)
     << "," << c.clustering_factor
     << "," << dtl::determine_clustering_factor(bs)
     << "," << c.bitmap_id
     << "," << enc_bs.size_in_byte()
     << "," << "\"" << type_info << "\""
     << "," << checksum
     << std::endl;
}
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
void run(config c, std::ostream& os) {
  switch (c.bitmap_type) {
    case bitmap_t::bitmap:
      run<dtl::dynamic_bitmap<$u32>>(c, os);
      break;
    case bitmap_t::roaring:
      run<dtl::dynamic_roaring_bitmap>(c, os);
      break;
    case bitmap_t::teb:
      run<dtl::teb<>>(c, os);
      break;
//    case bitmap_t::teb_scan:
//      run<dtl::teb_scan<>>(c, os);
//      break;
    case bitmap_t::wah:
      run<dtl::dynamic_wah32>(c, os);
      break;

//    case bitmap_t::position_list:
//      run<dtl::position_list<$u32>>(c, os);
//      break;
//    case bitmap_t::partitioned_position_list_u8:
//      run<dtl::partitioned_position_list<$u32, $u8>>(c, os);
//      break;
//    case bitmap_t::partitioned_position_list_u16:
//      run<dtl::partitioned_position_list<$u32, $u16>>(c, os);
//      break;
//    case bitmap_t::range_list:
//      run<dtl::range_list<$u32>>(c, os);
//      break;
//    case bitmap_t::partitioned_range_list_u8:
//      run<dtl::partitioned_range_list<$u32, $u8>>(c, os);
//      break;
//    case bitmap_t::partitioned_range_list_u16:
//      run<dtl::partitioned_range_list<$u32, $u16>>(c, os);
//      break;
  }
}
//===----------------------------------------------------------------------===//
void run(const std::vector<config>& configs) {
  std::function<void(const config&, std::ostream&)> fn =
      [](const config c, std::ostream& os) -> void {
        run(c, os);
      };
  const auto thread_cnt = 1; // run performance measurements single-threaded
  dispatch(configs, fn, thread_cnt);
}
//===----------------------------------------------------------------------===//
