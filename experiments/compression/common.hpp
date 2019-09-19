#pragma once
//===----------------------------------------------------------------------===//
#include "experiments/util/bitmap_db.hpp"
#include "experiments/util/threading.hpp"

#include <dtl/bitmap/dynamic_bitmap.hpp>
#include <dtl/bitmap/dynamic_roaring_bitmap.hpp>
#include <dtl/bitmap/dynamic_wah.hpp>
#include <dtl/bitmap/partitioned_position_list.hpp>
#include <dtl/bitmap/partitioned_range_list.hpp>
#include <dtl/bitmap/position_list.hpp>
#include <dtl/bitmap/range_list.hpp>
#include <dtl/bitmap/teb.hpp>
#include <dtl/bitmap/teb_wrapper.hpp>
#include <dtl/bitmap/util/convert.hpp>
#include <dtl/bitmap/util/random.hpp>
#include <dtl/bitmap/util/two_state_markov_process.hpp>
#include <dtl/dtl.hpp>
#include <dtl/env.hpp>
#include <dtl/thread.hpp>

#include <boost/algorithm/string.hpp>

#include <bitset>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
//===----------------------------------------------------------------------===//
/// The number of independent runs.
static constexpr u64 RUNS = 10;
/// The (default) size of a randomly generated bitmap.
static constexpr u64 N = 1u << 20;
/// A timestamp that identifies the current experiment.
static const i64 RUN_ID =
    std::chrono::duration_cast<std::chrono::seconds>(
        std::chrono::system_clock::now().time_since_epoch())
        .count();
/// The database file where the bitmaps are stored.
static const std::string DB_FILE =
    dtl::env<std::string>::get("DB_FILE", "./random_bitmaps.sqlite3");
/// 0 = run experiment, 1 = generate data for the experiment
static u64 GEN_DATA = dtl::env<$u64>::get("GEN_DATA", 0);
/// The database instance where the bitmaps are stored.
static bitmap_db db(DB_FILE);
//===----------------------------------------------------------------------===//
/// All implementations that should be considered in the experiment.
enum class bitmap_t {
  bitmap,
  roaring,
  teb,
  teb_scan,
  teb_wrapper,
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
//===----------------------------------------------------------------------===//
/// Refers to a single independent task that is executed during the experiment.
struct config {
  bitmap_t bitmap_type;
  $u64 n;
  $f64 density;
  $f64 clustering_factor;
  $i64 bitmap_id; // in DB
};
//===----------------------------------------------------------------------===//
template<typename T>
static void
run(const config& c, std::ostream& os) {
  // Load the bitmap from DB.
  auto bs = db.load_bitmap(c.bitmap_id);
  // Encode the bitmap.
  T enc_bs(bs);

  const auto size_in_bytes = enc_bs.size_in_byte();
  std::string type_info = enc_bs.info();

  // Validation. - Reconstruct the original bitmap.
  auto dec_bs = dtl::to_bitmap_using_iterator(enc_bs);
  if (bs != dec_bs) {
    std::cerr << "Validation failed for "
              << type_info << "."
              << "\nBitmap ID: " << c.bitmap_id
              << std::endl;
    std::exit(1);
  }

  boost::replace_all(type_info, "\"", "\"\""); // Escape JSON for CSV output.

  os << RUN_ID
     << "," << c.n
     << "," << T::name()
     << "," << c.density
     << "," << dtl::determine_bit_density(bs)
     << "," << c.clustering_factor
     << "," << dtl::determine_clustering_factor(bs)
     << "," << size_in_bytes
     << ","
     << "\"" << type_info << "\""
     << std::endl;
}
//===----------------------------------------------------------------------===//
static void
run(config c, std::ostream& os) {
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
    case bitmap_t::teb_wrapper:
      run<dtl::teb_wrapper>(c, os);
      break;
    case bitmap_t::wah:
      run<dtl::dynamic_wah32>(c, os);
      break;
      // clang-format off
//    case bitmap_t::teb_scan: /* deprecated*/
//      run<dtl::teb_scan<>>(c, os);
//      break;
    // EXPERIMENTAL
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
      // clang-format on
  }
}
//===----------------------------------------------------------------------===//
static void
run(const std::vector<config>& configs) {
  std::function<void(const config&, std::ostream&)> fn =
      [](const config c, std::ostream& os) -> void {
    run(c, os);
  };
  dispatch(configs, fn);
}
//===----------------------------------------------------------------------===//
