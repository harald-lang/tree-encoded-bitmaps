#include <dtl/dtl.hpp>
#include <dtl/env.hpp>
#include <dtl/bitmap/util/convert.hpp>
#include <dtl/bitmap/teb.hpp>

#include "../util/bitmap_db.hpp"
#include "../util/threading.hpp"

//===----------------------------------------------------------------------===//
const std::string DB_FILE =
    dtl::env<std::string>::get("DB_FILE", "./random_bitmaps.sqlite3");
bitmap_db db(DB_FILE);

constexpr auto opt_level = 2;

std::vector<$f64> fprs = {
    0.0,
    0.0001,
    0.00025,
    0.0005,
    0.00075,
    0.001,
    0.0025,
    0.005,
    0.0075,
    0.01,
    0.025,
    0.05,
    0.075,
    0.1
};
//===----------------------------------------------------------------------===//
void run(u64 bitmap_id, std::ostream& os) {
  const auto bm = db.load_bitmap(bitmap_id);
  const auto n = bm.size();
//    if (n < 1ull << 20) continue;
  const auto c = bm.count();

  dtl::teb<opt_level> teb_lossless(bm);

  for (auto fpr : fprs) {
    dtl::teb<opt_level> teb_lossy(bm, fpr);
    const auto dec_lossy = dtl::to_bitmap_using_iterator(teb_lossy);
    const auto fp_cnt = dec_lossy.count() - c;
    os << bitmap_id
        << "," << fpr
        << "," << teb_lossless.size_in_byte()
        << "," << teb_lossy.size_in_byte()
        << "," << fp_cnt
        << "," << (teb_lossy.size_in_byte() * 1.0 / teb_lossless.size_in_byte())
        << std::endl;

    // Validation
    {
      if ((bm & dec_lossy) != bm) {
        std::cerr << "Validation failed. Aborting." << std::endl;
        std::exit(1);
      }
    }
  }
}
//===----------------------------------------------------------------------===//
$i32 main() {
  auto bitmap_ids = db.ids();
  std::random_shuffle(bitmap_ids.begin(), bitmap_ids.end());
  std::function<void(i64&, std::ostream&)> fn =
      [](i64 bitmap_id, std::ostream& os) -> void {
        run(bitmap_id, os);
      };
  const auto thread_cnt = std::thread::hardware_concurrency() / 2;
  dispatch(bitmap_ids, fn, thread_cnt);
}
//===----------------------------------------------------------------------===//
