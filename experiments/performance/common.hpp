#pragma once

#include <bitset>
#include <chrono>
#include <cmath>
#include <limits>
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
#include <dtl/bitmap/position_list.hpp>
#include <dtl/bitmap/partitioned_position_list.hpp>
#include <dtl/bitmap/partitioned_range_list.hpp>
#include <dtl/bitmap/range_list.hpp>

#include <experiments/util/threading.hpp>
#include <dtl/bitmap/teb_scan.hpp>

#include "version.h"

// The number of independent runs.
static constexpr u64 RUNS = 10;
static constexpr u64 N = 1u << 20;

static const i64 RUN_ID = std::chrono::duration_cast<std::chrono::seconds>(
    std::chrono::system_clock::now().time_since_epoch()).count();

static const std::string DB_FILE =
    dtl::env<std::string>::get("DB_FILE", "./random_bitmaps.sqlite3");
static u64 GEN_DATA = dtl::env<$u64>::get("GEN_DATA", 0);

static bitmap_db db(DB_FILE);

static $u64 RUN_DURATION_NANOS = 250e6; // run for at least 250ms

u64
now_nanos() {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(
      std::chrono::steady_clock::now().time_since_epoch()).count();
}


//===----------------------------------------------------------------------===//
enum class bitmap_t {
  bitmap,
  roaring,
  teb,
  teb_scan,
  wah,
  position_list,
  partitioned_position_list_u8,
  partitioned_position_list_u16,
  range_list,
  partitioned_range_list_u8,
  partitioned_range_list_u16,
  _first = bitmap,
  _last = partitioned_range_list_u16
};
std::vector<std::string> bitmap_names {
    "bitmap",
    "roaring",
    "teb",
    "teb_scan",
    "wah",
    "position_list",
    "partitioned_position_list_u8",
    "partitioned_position_list_u16",
    "range_list",
    "partitioned_range_list_u8",
    "partitioned_range_list_u16",
};
std::ostream& operator<<(std::ostream& out, const bitmap_t& b) {
  const auto i = static_cast<int>(b);
  assert(i < bitmap_names.size());
  out << bitmap_names[i];
  return out;
}
//===----------------------------------------------------------------------===//
struct config {
  bitmap_t bitmap_type;
  $u64 n;
  $f64 density;
  $f64 clustering_factor;
  $i64 bitmap_id; // in DB

  void
  print(std::ostream& os) const noexcept {
    os << "config[bitmap_type=" << bitmap_type
       << ",n=" << n
       << ",d=" << density
       << ",f=" << clustering_factor
       << ",bitmap_id=" << bitmap_id
       << "]";
  }
};
//===----------------------------------------------------------------------===//
template<typename T>
void __attribute__ ((noinline))
run(const config& c, std::ostream& os) {
  const auto duration_nanos = RUN_DURATION_NANOS;
  const std::size_t MIN_REPS = 10;
  // Load the bitmap from DB.
  auto bs = db.load_bitmap(c.bitmap_id);
  const auto bs_count = bs.count();
  // Encode the bitmap.
  T enc_bs(bs);

  std::size_t pos_sink = 0;
  std::size_t length_sink = 0;

  // Warm up run
  {
    auto it = enc_bs.it();
    while (!it.end()) {
      pos_sink += it.pos();
      length_sink += it.length();
      it.next();
    }
  }

  // Validation. (Fail fast)
  if (!(length_sink == bs_count)) {
    std::cerr << "Validation failed: " << c << std::endl;
    std::cerr << "Expected length to be " << bs.count()
              << " but got " << length_sink << std::endl;
    std::exit(1);
  }

  // The actual measurement.
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


  if (!(length_sink / (rep_cntr + 1) == bs_count)) {
    std::cerr << "Validation failed: " << c << std::endl;
    std::cerr << "Expected length to be " << bs_count
              << " but got " << (length_sink / (rep_cntr + 1)) << std::endl;
    std::exit(1);
  }

  const std::size_t checksum = pos_sink + length_sink;

  std::string type_info = enc_bs.info();
  boost::replace_all(type_info, "\"", "\"\""); // Escape JSON for CSV output.

  os << RUN_ID
     << ",\"" << BUILD_ID << "\""
     << "," << c.n
     << "," << T::name()
     << "," << (nanos_end - nanos_begin) / rep_cntr
     << "," << (tsc_end - tsc_begin) / rep_cntr
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
    case bitmap_t::teb_scan:
      run<dtl::teb_scan<>>(c, os);
      break;
    case bitmap_t::wah:
      run<dtl::dynamic_wah32>(c, os);
      break;
    case bitmap_t::position_list:
      run<dtl::position_list<$u32>>(c, os);
      break;
    case bitmap_t::partitioned_position_list_u8:
      run<dtl::partitioned_position_list<$u32, $u8>>(c, os);
      break;
    case bitmap_t::partitioned_position_list_u16:
      run<dtl::partitioned_position_list<$u32, $u16>>(c, os);
      break;
    case bitmap_t::range_list:
      run<dtl::range_list<$u32>>(c, os);
      break;
    case bitmap_t::partitioned_range_list_u8:
      run<dtl::partitioned_range_list<$u32, $u8>>(c, os);
      break;
    case bitmap_t::partitioned_range_list_u16:
      run<dtl::partitioned_range_list<$u32, $u16>>(c, os);
      break;
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
