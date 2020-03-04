#pragma once
//===----------------------------------------------------------------------===//
#include "experiments/util/bitmap_db.hpp"
#include "experiments/util/bitmap_types.hpp"
#include "experiments/util/threading.hpp"
#include "version.h"

#include <dtl/bitmap/dynamic_bitmap.hpp>
#include <dtl/bitmap/dynamic_roaring_bitmap.hpp>
#include <dtl/bitmap/dynamic_wah.hpp>
#include <dtl/bitmap/position_list.hpp>
#include <dtl/bitmap/range_list.hpp>
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
#include <utility>
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

  const auto size_in_bytes = enc_bs.size_in_bytes();
  std::string type_info = enc_bs.info();
  std::string type_name = T::name();

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
  boost::replace_all(type_name, "\"", "\"\""); // Escape JSON for CSV output.

  os << RUN_ID
     << "," << c.n
     << ","
     << "\"" << type_name << "\""
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
#define __GENERATE_CASE(name)                  \
  case bitmap_t::name:                         \
    run<type_of<bitmap_t::name>::type>(c, os); \
    break;

    __GENERATE_CASE(bitmap)

    __GENERATE_CASE(roaring)

    __GENERATE_CASE(teb_wrapper)
    __GENERATE_CASE(partitioned_teb_wrapper)

    __GENERATE_CASE(wah)
    __GENERATE_CASE(partitioned_wah)

    __GENERATE_CASE(position_list)
    __GENERATE_CASE(partitioned_position_list_u8)
    __GENERATE_CASE(partitioned_position_list_u16)

    __GENERATE_CASE(range_list)
    __GENERATE_CASE(partitioned_range_list_u8)
    __GENERATE_CASE(partitioned_range_list_u16)

    __GENERATE_CASE(uah8)
    __GENERATE_CASE(uah8_skip)
    __GENERATE_CASE(uah16)
    __GENERATE_CASE(uah16_skip)
    __GENERATE_CASE(uah32)
    __GENERATE_CASE(uah32_skip)
    __GENERATE_CASE(uah64)
    __GENERATE_CASE(uah64_skip)
    __GENERATE_CASE(partitioned_uah8)
    __GENERATE_CASE(partitioned_uah8_skip)
    __GENERATE_CASE(partitioned_uah16)
    __GENERATE_CASE(partitioned_uah16_skip)
    __GENERATE_CASE(partitioned_uah32)
    __GENERATE_CASE(partitioned_uah32_skip)
    __GENERATE_CASE(partitioned_uah64)
    __GENERATE_CASE(partitioned_uah64_skip)

    __GENERATE_CASE(xah8)
    __GENERATE_CASE(xah8_skip)
    __GENERATE_CASE(xah16)
    __GENERATE_CASE(xah16_skip)
    __GENERATE_CASE(xah32)
    __GENERATE_CASE(xah32_skip)
    __GENERATE_CASE(xah64)
    __GENERATE_CASE(xah64_skip)
    __GENERATE_CASE(partitioned_xah8)
    __GENERATE_CASE(partitioned_xah8_skip)
    __GENERATE_CASE(partitioned_xah16)
    __GENERATE_CASE(partitioned_xah16_skip)
    __GENERATE_CASE(partitioned_xah32)
    __GENERATE_CASE(partitioned_xah32_skip)
    __GENERATE_CASE(partitioned_xah64)
    __GENERATE_CASE(partitioned_xah64_skip)

    __GENERATE_CASE(bah)
    __GENERATE_CASE(partitioned_bah)
#undef __GENERATE_CASE
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
