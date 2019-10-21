#include "experiments/util/bitmap_db.hpp"
#include "experiments/util/threading.hpp"

#include <dtl/bitmap/dynamic_bitmap.hpp>
#include <dtl/bitmap/dynamic_roaring_bitmap.hpp>
#include <dtl/bitmap/dynamic_wah.hpp>
#include <dtl/bitmap/teb_wrapper.hpp>
#include <dtl/bitmap/util/convert.hpp>
#include <dtl/bitmap/util/random.hpp>
#include <dtl/dtl.hpp>
#include <dtl/env.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include <atomic>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <limits>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
//===----------------------------------------------------------------------===//
namespace fs = boost::filesystem;
//===----------------------------------------------------------------------===//
void run(const std::string& dir, std::ostream& result_out) {
//  result_out
//      << "//===----------------------------------------------------------------------===//"
//      << "\ndata set: '" << dir << "'"
//      << std::endl;

  std::size_t min_val = std::numeric_limits<std::size_t>::max();
  std::size_t max_val = 0;
  std::size_t file_cnt = 0;

  std::vector<std::string> filenames;
  fs::directory_iterator end_it;
  for (fs::directory_iterator dir_it(dir); dir_it != end_it; ++dir_it) {
    const auto file = dir_it->path().string();
    filenames.push_back(file);
    ++file_cnt;
    std::cout << "reading file: " << file << std::endl;
    std::ifstream is(file);
    std::string token;
    while (std::getline(is, token, ',')) {
      const auto val = std::stoull(token);
      if (val < min_val) min_val = val;
      if (val > max_val) max_val = val;
    }
    is.close();
  }
  std::sort(filenames.begin(), filenames.end());

  const auto n = max_val + 1;
  const auto n_pow2 = dtl::next_power_of_two(n);
  const auto c = file_cnt;

//  result_out << "min_val: " << min_val
//             << " max_val: " << max_val
//             << " bitmap_cnt: " << file_cnt
//             << " n: " << n
//             << std::endl;

  $f64 min_d = std::numeric_limits<f64>::max();
  $f64 max_d = std::numeric_limits<f64>::min();
  $f64 sum_d = 0.0;

  $f64 min_f = std::numeric_limits<f64>::max();
  $f64 max_f = std::numeric_limits<f64>::min();
  $f64 sum_f = 0.0;

  std::atomic<std::size_t> bytes_roaring { 0 };
  std::atomic<std::size_t> bytes_wah { 0 };
  std::atomic<std::size_t> bytes_wah64 { 0 };
  std::atomic<std::size_t> bytes_teb { 0 };
  std::atomic<std::size_t> bytes_teb_rank { 0 };

  std::vector<dtl::bitmap> bitmaps;
  std::vector<dtl::bitmap> bitmaps_pow2;

  std::size_t total_bit_cnt = 0;
  for (auto& file : filenames) {
    std::ifstream is(file);
    std::string token;

    dtl::bitmap bm(n);
    dtl::bitmap bm_pow2(n_pow2);

    while (std::getline(is, token, ',')) {
      const auto val = std::stoull(token);
      bm[val] = true;
      bm_pow2[val] = true;
    }
    is.close();

    bitmaps.push_back(bm);
    bitmaps_pow2.push_back(bm_pow2);

    const auto d = dtl::determine_bit_density(bm);
    min_d = std::min(d, min_d);
    max_d = std::max(d, max_d);
    sum_d += d;

    const auto f = dtl::determine_clustering_factor(bm);
    min_f = std::min(f, min_f);
    max_f = std::max(f, max_f);
    sum_f += f;

    total_bit_cnt += bm.count();
  }
//  result_out << "cardinality: " << total_bit_cnt << std::endl;
//  result_out << "d (min/avg/max): "
//             << min_d << " / "
//             << (sum_d / file_cnt) << " / "
//             << max_d << std::endl;
//  result_out << "f (min/avg/max): "
//             << min_f << " / "
//             << (sum_f / file_cnt) << " / "
//             << max_f << std::endl;

  auto thread_fn = [&](const std::size_t bid, std::ostream& os) {
    auto& bm = bitmaps[bid];
    auto& bm_pow2 = bitmaps_pow2[bid];

    std::size_t r = 0;
    std::size_t t = 0;
    std::size_t t_rank = 0;
    std::size_t w = 0;
    std::size_t w64 = 0;
    {
      dtl::dynamic_roaring_bitmap roaring(bm);
      r = roaring.size_in_byte();
    }
    {
//      dtl::dynamic_wah32 wah(bm);
//      w = wah.size_in_byte();
    }
    {
//      dtl::dynamic_wah64 wah(bm);
//      w64 = wah.size_in_byte();
    }
    {
      dtl::teb<> teb(bm_pow2);
      t = teb.size_in_byte();
      t_rank = teb.rank_.size_in_bytes();
      const auto dec = dtl::to_bitmap_using_iterator(teb);
      if ((bm_pow2 & dec) != bm_pow2) {
        std::cerr << "Validation failed." << std::endl;
        std::exit(1);
      }
    }
    // clang-format off
// TEB (lossy compressed)
//    {
//      const auto fpr = 0.0001;
//      dtl::teb_wrapper teb(bm_pow2, fpr);
//      t = teb.size_in_byte();
//      const auto dec = dtl::to_bitmap_using_iterator(teb);
//      if ((bm_pow2 & dec) != bm_pow2) {
//        std::cerr << "Validation failed." << std::endl;
//        std::exit(1);
//      }
//      u64 max_fp_cnt = static_cast<u64>(teb.size() * fpr);
//      u64 fp_cnt = (bm_pow2 ^ dec).count();
//      std::cout << "fp_cnt=" << fp_cnt << std::endl;
//      if (fp_cnt > max_fp_cnt) {
//        std::cerr << "Validation failed. Max FP count exceeded." << std::endl;
//        std::exit(1);
//      }
//      os << teb.info() << std::endl;
//    }
    // clang-format on

    bytes_roaring += r;
    bytes_wah += w;
    bytes_wah64 += w64;
    bytes_teb += t;
    bytes_teb_rank += t_rank;
  };

  dispatch(0, bitmaps.size(), thread_fn);

//  result_out << "roaring: " << std::setw(15) << bytes_roaring << " bytes, " << std::setw(15) << std::setprecision(4) << ((bytes_roaring * 8.0) / total_bit_cnt) << " bits/int" << std::endl;
//  result_out << "teb:     " << std::setw(15) << bytes_teb << " bytes, " << std::setw(15) << std::setprecision(4) << ((bytes_teb * 8.0) / total_bit_cnt) << " bits/int" << std::endl;
//  result_out << "wah:     " << std::setw(15) << bytes_wah << " bytes, " << std::setw(15) << std::setprecision(4) << ((bytes_wah * 8.0) / total_bit_cnt) << " bits/int" << std::endl;
//  result_out << "wah64:   " << std::setw(15) << bytes_wah64 << " bytes, " << std::setw(15) << std::setprecision(4) << ((bytes_wah64 * 8.0) / total_bit_cnt) << " bits/int" << std::endl;

    result_out
      << "\"" << dir << "\","
      << total_bit_cnt << ","
      << bytes_teb << "," << bytes_teb_rank << ","
      << std::setprecision(4) << ((bytes_teb * 8.0) / total_bit_cnt)
      << std::endl;

//    result_out
//      << "\"" << dir << "\","
//      << total_bit_cnt << ","
//      << bytes_roaring << ","
//      << std::setprecision(4) << ((bytes_roaring * 8.0) / total_bit_cnt)
//      << std::endl;
}
//===----------------------------------------------------------------------===//
$i32 main() {
  const auto basedir = dtl::env<std::string>::get("DIR",
      "./real-roaring-dataset");
  std::vector<std::string> dirs {
    basedir + "/census1881/",
    basedir + "/census1881_srt/",
    basedir + "/census-income/",
    basedir + "/census-income_srt/",
    basedir + "/weather_sept_85/",
    basedir + "/weather_sept_85_srt/",
    basedir + "/wikileaks-noquotes/",
    basedir + "/wikileaks-noquotes_srt/",
  };

  std::stringstream results;
  for (auto& dir : dirs) {
    run(dir, results);
  }
  std::cout << "Results" << std::endl;
  std::cout << results.str() << std::endl;
}
//===----------------------------------------------------------------------===//
