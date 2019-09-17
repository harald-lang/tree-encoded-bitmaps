#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/bitmap/bitwise_operations.hpp>
#include "common.hpp"
//===----------------------------------------------------------------------===//
template<typename T>
void __attribute__ ((noinline))
run_intersect(const config_pair& c, std::ostream& os) {
  const auto duration_nanos = RUN_DURATION_NANOS;
  const std::size_t MIN_REPS = 10;
  // Load the bitmap from DB.
  const auto bs1 = db.load_bitmap(c.bitmap_id1);
  const auto bs2 = db.load_bitmap(c.bitmap_id2);
  // Encode the bitmap.
  const T enc_bs1(bs1);
  const T enc_bs2(bs2);

  // Validation code.
  {
    const auto expected_bs = bs1 & bs2;
    {
      auto it1 = enc_bs1.scan_it();
      auto it2 = enc_bs2.it();
      auto result_it = dtl::bitwise_and_it(it1, it2);
      const auto result_bs = dtl::to_bitmap_from_iterator(result_it, bs1.size());
      if (result_bs != expected_bs) {
        std::cerr << "Validation failed: " << c << std::endl;
        std::exit(1);
      }
    }
    {
      auto it1 = enc_bs1.it();
      auto it2 = enc_bs2.it();
      auto result_it = dtl::bitwise_and_it(it1, it2);
      const auto result_bs = dtl::to_bitmap_from_iterator(result_it, bs1.size());
      if (result_bs != expected_bs) {
        std::cerr << "Validation failed: " << c << std::endl;
        std::exit(1);
      }
    }
  }

  // Warm up run
  std::size_t checksum = 0;
  {
    std::size_t pos_sink = 0;
    std::size_t length_sink = 0;
    auto it1 = enc_bs1.scan_it();
    auto it2 = enc_bs2.it();
    auto result_it = dtl::bitwise_and_it(it1, it2);
    while (!result_it.end()) {
      pos_sink += result_it.pos();
      length_sink += result_it.length();
      result_it.next();
    }

    checksum += pos_sink + length_sink;
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
      auto it1 = enc_bs1.scan_it();
      auto it2 = enc_bs2.it();
      auto result_it = dtl::bitwise_and_it(it1, it2);
      while (!result_it.end()) {
        pos_sink += result_it.pos();
        length_sink += result_it.length();
        result_it.next();
      }
    }
    const auto tsc_end = _rdtsc();
    const auto nanos_end = now_nanos();

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
      auto it1 = enc_bs1.it();
      auto it2 = enc_bs2.it();
      auto result_it = dtl::bitwise_and_it(it1, it2);
      while (!result_it.end()) {
        pos_sink += result_it.pos();
        length_sink += result_it.length();
        result_it.next();
      }
    }
    const auto tsc_end = _rdtsc();
    const auto nanos_end = now_nanos();

    runtime_nanos_skip_it = (nanos_end - nanos_begin) / rep_cntr;
    runtime_cycles_skip_it = (tsc_end - tsc_begin) / rep_cntr;
    checksum += pos_sink + length_sink;
  }

  std::string type_info1 = enc_bs1.info();
  boost::replace_all(type_info1, "\"", "\"\""); // Escape JSON for CSV output.
  std::string type_info2 = enc_bs2.info();
  boost::replace_all(type_info2, "\"", "\"\""); // Escape JSON for CSV output.

  os << RUN_ID
     << ",\"" << BUILD_ID << "\""
     << "," << c.n
     << "," << T::name()
     << "," << runtime_nanos_skip_it
     << "," << runtime_cycles_skip_it
     << "," << runtime_nanos_scan_it
     << "," << runtime_cycles_scan_it
     << "," << c.density1
     << "," << c.density2
     << "," << dtl::determine_bit_density(bs1)
     << "," << dtl::determine_bit_density(bs2)
     << "," << c.clustering_factor1
     << "," << c.clustering_factor2
     << "," << dtl::determine_clustering_factor(bs1)
     << "," << dtl::determine_clustering_factor(bs2)
     << "," << c.bitmap_id1
     << "," << c.bitmap_id2
     << "," << enc_bs1.size_in_byte()
     << "," << enc_bs2.size_in_byte()
     << "," << "\"" << type_info1 << "\""
     << "," << "\"" << type_info2 << "\""
     << "," << checksum
     << std::endl;
}
//===----------------------------------------------------------------------===//
void run_intersect(config_pair c, std::ostream& os) {
  switch (c.bitmap_type) {
    case bitmap_t::bitmap:
      run_intersect<dtl::dynamic_bitmap<$u32>>(c, os);
      break;
    case bitmap_t::roaring:
      run_intersect<dtl::dynamic_roaring_bitmap>(c, os);
      break;
    case bitmap_t::teb:
      run_intersect<dtl::teb<>>(c, os);
      break;
//    case bitmap_t::teb_scan:
//      run_intersect<dtl::teb_scan<>>(c, os);
//      break;
    case bitmap_t::wah:
      run_intersect<dtl::dynamic_wah32>(c, os);
      break;
    case bitmap_t::teb_wrapper:
      run_intersect<dtl::teb_wrapper>(c, os);
      break;
    // EXPERIMENTAL
//    case bitmap_t::position_list:
//      run_intersect<dtl::position_list<$u32>>(c, os);
//      break;
//    case bitmap_t::partitioned_position_list_u8:
//      run_intersect<dtl::partitioned_position_list<$u32, $u8>>(c, os);
//      break;
//    case bitmap_t::partitioned_position_list_u16:
//      run_intersect<dtl::partitioned_position_list<$u32, $u16>>(c, os);
//      break;
//    case bitmap_t::range_list:
//      run_intersect<dtl::range_list<$u32>>(c, os);
//      break;
//    case bitmap_t::partitioned_range_list_u8:
//      run_intersect<dtl::partitioned_range_list<$u32, $u8>>(c, os);
//      break;
//    case bitmap_t::partitioned_range_list_u16:
//      run_intersect<dtl::partitioned_range_list<$u32, $u16>>(c, os);
//      break;
  }
}
//===----------------------------------------------------------------------===//
void run_intersect(const std::vector<config_pair>& configs) {
  std::function<void(const config_pair&, std::ostream&)> fn =
      [](const config_pair c, std::ostream& os) -> void {
        run_intersect(c, os);
      };
  const auto thread_cnt = 1; // run performance measurements single-threaded
  dispatch(configs, fn, thread_cnt);
}
//===----------------------------------------------------------------------===//
