#include <atomic>
#include <regex>
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
i32 extract_bitmap_id_from_filename(const std::string file) {
  std::regex exp("csv(\\d+)\\.txt$",
      std::regex_constants::ECMAScript | std::regex_constants::icase);
  auto match_it =
      std::sregex_iterator(file.begin(), file.end(), exp);
  if (match_it == std::sregex_iterator()) {
    std::cerr << "Unexpected file name. Skipping " << file << std::endl;
    return -1;
  }
  std::smatch match = *match_it;
  std::string match_str = match.str(1);

  const auto bitmap_id = std::stoi(match_str);
  return bitmap_id;
}
//===----------------------------------------------------------------------===//
void run(const std::string& dir, const std::string& name,
    std::ostream& result_out) {

  std::size_t min_val = std::numeric_limits<std::size_t>::max();
  std::size_t max_val = 0;
  std::size_t file_cnt = 0;

  // Examine the bitmaps in the given directory.
  std::vector<std::string> filenames;
  fs::directory_iterator end_it;
  for (fs::directory_iterator dir_it(dir); dir_it != end_it; ++dir_it ) {
    const auto file = dir_it->path().string();
    filenames.push_back(file);
    ++file_cnt;
    std::cerr << "reading file: " << file << std::endl;
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

  std::cerr << "min=" << min_val
      << ", max_val=" << max_val
      << ", file_cnt=" << file_cnt
      << std::endl;

  const auto n = max_val + 1;
  const auto n_pow2 = dtl::next_power_of_two(n);
  const auto c = file_cnt;

  // Read the raw bitmaps.
  std::vector<dtl::bitmap> bitmaps_pow2(file_cnt);
  std::size_t total_bit_cnt = 0;
  for (auto& file : filenames) {
    const auto bitmap_id = extract_bitmap_id_from_filename(file);
    if (bitmap_id < 0) continue;

    std::ifstream is(file);
    std::string token;

    dtl::bitmap bm_pow2(n_pow2);

    while (std::getline(is, token, ',')) {
      const auto val = std::stoull(token);
      bm_pow2[val] = true;
    }
    is.close();

    bitmaps_pow2[bitmap_id] = bm_pow2;
    total_bit_cnt += bm_pow2.count();
  }
  std::cerr << "total bit cnt: " << total_bit_cnt << std::endl;

  for (std::size_t i = 0; i < bitmaps_pow2.size(); ++i) {
    if (bitmaps_pow2[i].count() == 0) {
      std::cerr << "Error. Found an empty bitmap. Exiting." << std::endl;
      std::exit(1);
    }
  }


  std::atomic<std::size_t> bytes_teb1 { 0 };
  std::atomic<std::size_t> bytes_teb2 { 0 };
  std::atomic<std::size_t> bytes_teb3 { 0 };
  std::atomic<std::size_t> bytes_teb4 { 0 };
  std::atomic<std::size_t> bytes_teb5 { 0 };

  auto thread_fn = [&](const std::size_t bid, std::ostream& os) {
    auto& bm_pow2 = bitmaps_pow2[bid];

    auto compress = [&](f64 fpr) {
//      const auto fpr = (bm_pow2.count() * p) / n;
      dtl::teb<> teb(bm_pow2, fpr);

      // Validation.
      const auto dec = dtl::to_bitmap_using_iterator(teb);
      if ((bm_pow2 & dec) != bm_pow2) {
        std::cerr << "Validation failed. False negative encountered."
            << std::endl;
        std::exit(1);
      }
      u64 max_fp_cnt = static_cast<u64>(teb.size() * fpr);
      u64 fp_cnt = (bm_pow2 ^ dec).count();
      if (fp_cnt > max_fp_cnt) {
        std::cerr << "Validation failed. Max false positive count exceeded."
            << std::endl;
        std::exit(1);
      }
      return teb;
    };

    for (auto fpr_rec : {0, 10000, 1000, 100, 10}) {
      const auto fpr = fpr_rec == 0 ? 0.0 : 1.0 / fpr_rec;
      const auto teb = compress(fpr);

      const auto dec = dtl::to_bitmap_using_iterator(teb);
      const auto fp_cnt = (bm_pow2 ^ dec).count();
      const auto fpr_actual = (fp_cnt * 1.0) / n;
      const auto d = dtl::determine_bit_density(dec);
      const auto f = dtl::determine_clustering_factor(dec);

      std::string type_info = teb.info();
      boost::replace_all(type_info, "\"", "\"\""); // Escape JSON for CSV output.

      os << name
          << "," << bid
          << "," << n
          << "," << dec.count()
          << "," << fpr_rec
          << "," << fpr
          << "," << fpr_actual
          << "," << fp_cnt
          << "," << d
          << "," << f
          << "," << teb.size_in_byte()
          << ",\"" << type_info << "\""
          << std::endl;
    }

  };

  dispatch(0, bitmaps_pow2.size(), thread_fn);

  auto savings_percent = [&](const std::size_t lossy_compressed_size) {
    const auto compressed_size = bytes_teb1.load();
    return 100.0 - ((100.0 / compressed_size) * lossy_compressed_size);
  };

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
    "income",
    "income-srt",
    "census1881",
    "census1881-srt",
    "weather",
    "weather-srt",
    "wikileaks",
    "wikileaks-srt",
  };

  std::stringstream results;
  for (std::size_t i = 0; i < dirs.size(); ++i) {
    run(dirs[i], names[i], results);
  }
//  std::cout << results.str() << std::endl;
}
//===----------------------------------------------------------------------===//
