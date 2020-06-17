#include "experiments/util/bitmap_db.hpp"
#include "experiments/util/threading.hpp"

#include <dtl/bitmap/bah.hpp>
#include <dtl/bitmap/bbc.hpp>
#include <dtl/bitmap/concise.hpp>
#include <dtl/bitmap/dynamic_bitmap.hpp>
#include <dtl/bitmap/dynamic_roaring_bitmap.hpp>
#include <dtl/bitmap/dynamic_wah.hpp>
#include <dtl/bitmap/teb_wrapper.hpp>
#include <dtl/bitmap/uah.hpp>
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
template<typename T>
std::size_t
get_compressed_size(dtl::bitmap& bm) {
  T enc(bm);
  std::size_t size = enc.size_in_bytes();
  // Validation
  const auto dec = dtl::to_bitmap_using_iterator(enc);
  if ((bm & dec) != bm) {
    std::cerr << "Validation failed for type " << enc.name() << std::endl;
    std::exit(1);
  }
  return size;
}
// Disable validation for WAH64, since the iterator on that type is not
// implemented.
template<>
std::size_t
get_compressed_size<dtl::dynamic_wah64>(dtl::bitmap& bm) {
  dtl::dynamic_wah64 enc(bm);
  std::size_t size = enc.size_in_bytes();
  return size;
}
//===----------------------------------------------------------------------===//
void run(const std::string& dir, std::ostream& result_out,
    std::ostream& result_csv_out) {
  result_out
      << "//===----------------------------------------------------------------------===//"
      << "\ndata set: '" << dir << "'"
      << std::endl;

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

  result_out << "min_val: " << min_val
             << " max_val: " << max_val
             << " bitmap_cnt: " << file_cnt
             << " n: " << n
             << std::endl;

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
  std::atomic<std::size_t> bytes_uah8 { 0 };
  std::atomic<std::size_t> bytes_uah16 { 0 };
  std::atomic<std::size_t> bytes_uah32 { 0 };
  std::atomic<std::size_t> bytes_bah { 0 };
  std::atomic<std::size_t> bytes_concise { 0 };
  std::atomic<std::size_t> bytes_bbc { 0 };

  std::vector<dtl::bitmap> bitmaps;
  std::vector<dtl::bitmap> bitmaps_pow2;

  std::size_t total_bit_cnt = 0;
  for (auto& file : filenames) {
    std::ifstream is(file);
    std::string token;

    dtl::bitmap bm(n);

    while (std::getline(is, token, ',')) {
      const auto val = std::stoull(token);
      bm[val] = true;
    }
    is.close();

    bitmaps.push_back(bm);

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
  result_out << "cardinality: " << total_bit_cnt << std::endl;
  result_out << "d (min/avg/max): "
             << min_d << " / "
             << (sum_d / file_cnt) << " / "
             << max_d << std::endl;
  result_out << "f (min/avg/max): "
             << min_f << " / "
             << (sum_f / file_cnt) << " / "
             << max_f << std::endl;

  auto thread_fn = [&](const std::size_t bid, std::ostream& os) {
    auto& bm = bitmaps[bid];
    auto& bm_pow2 = bitmaps_pow2[bid];

    std::size_t r = 0;
    std::size_t t = 0;
    std::size_t w = 0;
    std::size_t w64 = 0;
    std::size_t uah8 = 0;
    std::size_t uah16 = 0;
    std::size_t uah32 = 0;
    std::size_t bah = 0;
    std::size_t concise = 0;
    std::size_t bbc = 0;
    {
      using T = dtl::dynamic_roaring_bitmap;
      r = get_compressed_size<T>(bm);
    }
    {
      using T = dtl::dynamic_wah32;
      w = get_compressed_size<T>(bm);
    }
    {
      using T = dtl::dynamic_wah64;
      w64 = get_compressed_size<T>(bm);
    }
    {
      using T = dtl::teb_wrapper;
      t = get_compressed_size<T>(bm);
    }
    {
      using T = dtl::uah8;
      uah8 = get_compressed_size<T>(bm);
    }
    {
      using T = dtl::uah16;
      uah16 = get_compressed_size<T>(bm);
    }
    {
      using T = dtl::uah32;
      uah32 = get_compressed_size<T>(bm);
    }
    {
      using T = dtl::bah;
      bah = get_compressed_size<T>(bm);
    }
    {
      using T = dtl::concise;
      concise = get_compressed_size<T>(bm);
    }
    {
      using T = dtl::bbc;
      bbc = get_compressed_size<T>(bm);
//      dtl::bbc x(bm);
//      bbc = x.size_in_bytes();
//      const auto dec = dtl::to_bitmap_using_iterator(x);
//      if (bm != dec) {
//        std::cerr << "BBC: Validation failed." << std::endl;
//        std::cout << "n=" << bm.size() << std::endl;
//        std::cout << "expected: ";
//        auto it_a = dtl::plain_bitmap_iter<dtl::bitmap>(bm);
//        while (!it_a.end()) {
//          std::cout << "[" << it_a.pos() << "," << (it_a.pos() + it_a.length())
//            << "),";
//          it_a.next();
//        }
//        std::cout << std::endl;
//        std::cout << "but got:  ";
//        auto it_b = x.it();
//        while (!it_b.end()) {
//          std::cout << "[" << it_b.pos() << "," << (it_b.pos() + it_b.length())
//              << "),";
//          it_b.next();
//        }
//        std::cout << std::endl;
//        std::cout << x << std::endl;
//        std::exit(1);
//      }
    }
    // clang-format off
// TEB (lossy compressed)
//    {
//      const auto fpr = 0.0001;
//      dtl::teb_wrapper teb(bm_pow2, fpr);
//      t = teb.size_in_bytes();
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
    bytes_uah8 += uah8;
    bytes_uah16 += uah16;
    bytes_uah32 += uah32;
    bytes_bah += bah;
    bytes_concise += concise;
    bytes_bbc += bbc;
  };

  dispatch(0, bitmaps.size(), thread_fn);

  result_out << "roaring: " << std::setw(15) << bytes_roaring << " bytes, " << std::setw(15) << std::setprecision(4) << ((bytes_roaring * 8.0) / total_bit_cnt) << " bits/int" << std::endl;
  result_out << "teb:     " << std::setw(15) << bytes_teb << " bytes, " << std::setw(15) << std::setprecision(4) << ((bytes_teb * 8.0) / total_bit_cnt) << " bits/int" << std::endl;
  result_out << "wah:     " << std::setw(15) << bytes_wah << " bytes, " << std::setw(15) << std::setprecision(4) << ((bytes_wah * 8.0) / total_bit_cnt) << " bits/int" << std::endl;
  result_out << "wah64:   " << std::setw(15) << bytes_wah64 << " bytes, " << std::setw(15) << std::setprecision(4) << ((bytes_wah64 * 8.0) / total_bit_cnt) << " bits/int" << std::endl;
  result_out << "uah8 :   " << std::setw(15) << bytes_uah8 << " bytes, " << std::setw(15) << std::setprecision(4) << ((bytes_uah8 * 8.0) / total_bit_cnt) << " bits/int" << std::endl;
  result_out << "uah16 :  " << std::setw(15) << bytes_uah16 << " bytes, " << std::setw(15) << std::setprecision(4) << ((bytes_uah16 * 8.0) / total_bit_cnt) << " bits/int" << std::endl;
  result_out << "uah32 :  " << std::setw(15) << bytes_uah32 << " bytes, " << std::setw(15) << std::setprecision(4) << ((bytes_uah32 * 8.0) / total_bit_cnt) << " bits/int" << std::endl;
  result_out << "bah :    " << std::setw(15) << bytes_bah << " bytes, " << std::setw(15) << std::setprecision(4) << ((bytes_bah * 8.0) / total_bit_cnt) << " bits/int" << std::endl;
  result_out << "concise: " << std::setw(15) << bytes_concise << " bytes, " << std::setw(15) << std::setprecision(4) << ((bytes_concise * 8.0) / total_bit_cnt) << " bits/int" << std::endl;
  result_out << "bbc:     " << std::setw(15) << bytes_bbc << " bytes, " << std::setw(15) << std::setprecision(4) << ((bytes_bbc * 8.0) / total_bit_cnt) << " bits/int" << std::endl;

  result_csv_out << "\"" << dir << "\"";
  result_csv_out << "," << std::setprecision(4) << ((bytes_roaring * 8.0) / total_bit_cnt);
  result_csv_out << "," << std::setprecision(4) << ((bytes_teb * 8.0) / total_bit_cnt);
  result_csv_out << "," << std::setprecision(4) << ((bytes_wah * 8.0) / total_bit_cnt);
  result_csv_out << "," << std::setprecision(4) << ((bytes_bah * 8.0) / total_bit_cnt);
  result_csv_out << "," << std::setprecision(4) << ((bytes_concise * 8.0) / total_bit_cnt);
  result_csv_out << "," << std::setprecision(4) << ((bytes_bbc * 8.0) / total_bit_cnt);
  result_csv_out << std::endl;
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
  std::stringstream results_csv;
  for (auto& dir : dirs) {
    run(dir, results, results_csv);
  }
  std::cout << "Results:" << std::endl;
  std::cout << results.str() << std::endl;
  std::cout << std::endl;
  std::cout << "Results [bits per value]:" << std::endl;
  std::cout << "=========================" << std::endl;
  std::cout << std::endl;
  std::cout << "data_set,roaring,teb,wah,bah,concise,bbc" << std::endl;
  std::cout << results_csv.str() << std::endl;
}
//===----------------------------------------------------------------------===//
