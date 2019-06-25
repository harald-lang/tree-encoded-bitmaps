#include <atomic>
#include <vector>

#include <dtl/dtl.hpp>
#include <dtl/env.hpp>
#include <dtl/bitmap/util/convert.hpp>
#include <dtl/bitmap/teb.hpp>
#include <dtl/bitmap/dynamic_bitmap.hpp>
#include <dtl/bitmap/dynamic_roaring_bitmap.hpp>
#include <dtl/bitmap/dynamic_wah.hpp>

#include <boost/algorithm/string.hpp>
#include <dtl/bitmap/util/random.hpp>
#include <fstream>

#include "../util/bitmap_db.hpp"
#include "../util/threading.hpp"

#include "boost/filesystem.hpp"

//===----------------------------------------------------------------------===//
namespace fs = boost::filesystem;
//===----------------------------------------------------------------------===//
void run(const std::string& dir, const std::string& name,
    std::ostream& result_out) {

  std::size_t min_val = std::numeric_limits<std::size_t>::max();
  std::size_t max_val = 0;
  std::size_t file_cnt = 0;

  std::vector<std::string> filenames;
  fs::directory_iterator end_it;
  for (fs::directory_iterator dir_it(dir); dir_it != end_it; ++dir_it ) {
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

  std::cout << "min=" << min_val
      << ", max_val=" << max_val
      << ", file_cnt=" << file_cnt
      << std::endl;

  const auto n = max_val + 1;
  const auto n_pow2 = dtl::next_power_of_two(n);
  const auto c = file_cnt;

  std::atomic<std::size_t> bytes_teb1 { 0 };
  std::atomic<std::size_t> bytes_teb2 { 0 };
  std::atomic<std::size_t> bytes_teb3 { 0 };
  std::atomic<std::size_t> bytes_teb4 { 0 };
  std::atomic<std::size_t> bytes_teb5 { 0 };

  std::vector<dtl::bitmap> bitmaps_pow2;

  std::size_t total_bit_cnt = 0;
  for (auto& file : filenames) {
//    std::cout << "reading file: " << file << std::endl;
    std::ifstream is(file);
    std::string token;

    dtl::bitmap bm_pow2(n_pow2);

    while (std::getline(is, token, ',')) {
      const auto val = std::stoull(token);
      bm_pow2[val] = true;
    }
    is.close();

    bitmaps_pow2.push_back(bm_pow2);
    total_bit_cnt += bm_pow2.count();
  }
  std::cout << "total bit cnt: " << total_bit_cnt << std::endl;

  auto thread_fn = [&](const std::size_t bid, std::ostream& os) {
    auto& bm_pow2 = bitmaps_pow2[bid];

    auto compressed_size = [&](f64 fpr) {
//      const auto fpr = (bm_pow2.count() * p) / n;
      dtl::teb<> teb(bm_pow2, fpr);
      const auto dec = dtl::to_bitmap_using_iterator(teb);
      if ((bm_pow2 & dec) != bm_pow2) {
        std::cerr << "Validation failed. False negative encountered."
            << std::endl;
        std::exit(1);
      }
      u64 max_fp_cnt = static_cast<u64>(teb.size() * fpr);
      u64 fp_cnt = (bm_pow2 ^ dec).count();
      std::cout << "fp_cnt=" << fp_cnt << std::endl;
      if (fp_cnt > max_fp_cnt) {
        std::cerr << "Validation failed. Max false positive count exceeded."
            << std::endl;
        std::exit(1);
      }
      return teb.size_in_byte();
    };

    bytes_teb1 += compressed_size(0.0);
    bytes_teb2 += compressed_size(0.0001);
    bytes_teb3 += compressed_size(0.001);
    bytes_teb4 += compressed_size(0.01);
    bytes_teb5 += compressed_size(0.1);
  };

  dispatch(0, bitmaps_pow2.size(), thread_fn);

  auto savings_percent = [&](const std::size_t lossy_compressed_size) {
    const auto compressed_size = bytes_teb1.load();
    return 100.0 - ((100.0 / compressed_size) * lossy_compressed_size);
  };

  result_out << "% " << name
      << " | " << std::setw(10) << bytes_teb1
      << " | " << std::setw(10) << bytes_teb2
      << " | " << std::setw(10) << bytes_teb3
      << " | " << std::setw(10) << bytes_teb4
      << " | " << std::setw(10) << bytes_teb5
      << std::endl;
  result_out << std::setw(30) << name
      << " & " << std::setw(5) << std::fixed << std::setprecision(1)
        << savings_percent(bytes_teb2) << "\\,\\%"
      << " & " << std::setw(5) << std::fixed << std::setprecision(1)
        << savings_percent(bytes_teb3) << "\\,\\%"
      << " & " << std::setw(5) << std::fixed << std::setprecision(1)
        << savings_percent(bytes_teb4) << "\\,\\%"
      << " & " << std::setw(5) << std::fixed << std::setprecision(1)
        << savings_percent(bytes_teb5) << "\\,\\%"
      << " \\\\" << std::endl;
}
//===----------------------------------------------------------------------===//
$i32 main() {
  const auto basedir = dtl::env<std::string>::get("DIR",
      "/home/hl/git/storage/RoaringBitmap/real-roaring-dataset/src/main/resources/real-roaring-dataset");
  std::vector<std::string> dirs {
    basedir + "/census-income/",
    basedir + "/census-income_srt/",
    basedir + "/census1881/",
    basedir + "/census1881_srt/",
    basedir + "/weather_sept_85/",
    basedir + "/weather_sept_85_srt/",
    basedir + "/wikileaks-noquotes/",
    basedir + "/wikileaks-noquotes_srt/",
  };

  std::vector<std::string> names {
    "{\\sc \\small Census Income}",
    "{\\sc \\small Census Income}\\tiny{ (sorted)}",
    "{\\sc \\small Census 1881}",
    "{\\sc \\small Census 1881}\\tiny{ (sorted)}",
    "{\\sc \\small Weather}",
    "{\\sc \\small Weather}\\tiny{ (sorted)}",
    "{\\sc \\small WikiLeaks}",
    "{\\sc \\small WikiLeaks}\\tiny{ (sorted)}",
  };

  std::stringstream results;
  for (std::size_t i = 0; i < dirs.size(); ++i) {
    run(dirs[i], names[i], results);
  }
  std::cout << results.str() << std::endl;
}
//===----------------------------------------------------------------------===//
