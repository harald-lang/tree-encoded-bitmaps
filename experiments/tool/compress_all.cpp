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

#include "../util/bitmap_db.hpp"
#include "../util/threading.hpp"

//===----------------------------------------------------------------------===//
const std::string DB_FILE =
    dtl::env<std::string>::get("DB_FILE", "./random_bitmaps.sqlite3");
bitmap_db db(DB_FILE);

std::vector<$f64> fprs = {
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
void prepare_db() {
  auto* sqlite_db = db.get_db();
  char* err_msg = nullptr;
  auto rc = sqlite3_exec(sqlite_db,
      "drop table if exists compression_info;"
      "create table compression_info ("
      "  bitmap_id  integer not null,"
      "  bitmap_type varchar(128) not null, "
      "  fpr double not null, "
      "  fpc bigint not null, "
      "  compressed_size bigint not null, "
      "  bitmap_info json not null);",
      nullptr, nullptr, &err_msg);
}
//===----------------------------------------------------------------------===//
template<typename T>
void run(u64 bitmap_id, std::ostream& os) {
  const auto bm = db.load_bitmap(bitmap_id);

  T enc_bm(bm);
  const auto dec_bm = dtl::to_bitmap_using_iterator(enc_bm);

  std::string type_info = enc_bm.info();
  boost::replace_all(type_info, "\"", "\"\""); // Escape JSON for CSV output.

  os << bitmap_id
      << "," << enc_bm.name()
      << "," << 0.0
      << "," << 0
      << "," << enc_bm.size_in_byte()
      << ",\"" << type_info << "\""
      << std::endl;

  // Validation
  {
    if (dec_bm != bm) {
      std::cerr << "Validation failed. Aborting." << std::endl;
      std::cerr << bitmap_id
          << "," << enc_bm.name()
          << "," << 0.0
          << "," << 0
          << "," << enc_bm.size_in_byte()
          << "," << enc_bm.info()
          << std::endl;
      std::exit(1);
    }
  }


  std::stringstream sql;
  sql << "insert into compression_info values ("
      << std::to_string(bitmap_id)
      << ", '" << enc_bm.name() << "'"
      << ",0.0"
      << ",0"
      << "," + std::to_string(enc_bm.size_in_byte())
      << ",'" + enc_bm.info() + "');";
  auto* sqlite_db = db.get_db();
  char* err_msg = nullptr;
  auto rc = sqlite3_exec(sqlite_db,
      sql.str().c_str(),
      nullptr, nullptr, &err_msg);
  if (rc) {
    std::stringstream err;
    err << "Can't write to database: "
        << sqlite3_errmsg(sqlite_db)
        << std::endl;
    std::cerr << err.str();
    std::exit(1);
  }

}
//===----------------------------------------------------------------------===//
template<typename T>
void run_lossy(u64 bitmap_id, f64 fpr, std::ostream& os) {
  const auto bm = db.load_bitmap(bitmap_id);
  const auto c = bm.count();

  T enc_bm(bm, fpr);
  const auto dec_bm = dtl::to_bitmap_using_iterator(enc_bm);
  const auto fp_cnt = dec_bm.count() - c;

  std::string type_info = enc_bm.info();
  boost::replace_all(type_info, "\"", "\"\""); // Escape JSON for CSV output.

  os << bitmap_id
      << "," << enc_bm.name()
      << "," << fpr
      << "," << fp_cnt
      << "," << enc_bm.size_in_byte()
      << ",\"" << type_info << "\""
      << std::endl;

  // Validation
  {
    if ((bm & dec_bm) != bm) {
      std::cerr << "Validation failed. Aborting." << std::endl;
      std::cerr << bitmap_id
          << "," << enc_bm.name()
          << "," << fpr
          << "," << fp_cnt
          << "," << enc_bm.size_in_byte()
          << "," << enc_bm.info()
          << std::endl;
      std::exit(1);
    }
  }
}
//===----------------------------------------------------------------------===//
void run(u64 bitmap_id, std::ostream& os) {
  run<dtl::dynamic_bitmap<$u32>>(bitmap_id, os);
  run<dtl::dynamic_roaring_bitmap>(bitmap_id, os);
  run<dtl::dynamic_wah32>(bitmap_id, os);
  run<dtl::teb<0>>(bitmap_id, os);
  run<dtl::teb<1>>(bitmap_id, os);
  run<dtl::teb<2>>(bitmap_id, os);
  run<dtl::teb<3>>(bitmap_id, os);
//  for (auto fpr : fprs) {
//    run_lossy<dtl::teb<0>>(bitmap_id, fpr, os);
//    run_lossy<dtl::teb<1>>(bitmap_id, fpr, os);
//    run_lossy<dtl::teb<2>>(bitmap_id, fpr, os);
//    run_lossy<dtl::teb<3>>(bitmap_id, fpr, os);
//  }
}
//===----------------------------------------------------------------------===//
$i32 main() {
  prepare_db();
  std::vector<$i64> all_bitmap_ids = db.ids();
  std::vector<$i64> bitmap_ids;
  std::copy_if (all_bitmap_ids.begin(), all_bitmap_ids.end(),
      std::back_inserter(bitmap_ids), [&](auto bid) {
         return db.load_bitmap(bid).size() == (1ull << 20);
      }
  );
  std::random_shuffle(bitmap_ids.begin(), bitmap_ids.end());
  std::function<void(i64&, std::ostream&)> fn =
      [](i64 bitmap_id, std::ostream& os) -> void {
        run(bitmap_id, os);
      };
  const auto thread_cnt = dtl::env<$u64>::get("THREAD_CNT",
      std::thread::hardware_concurrency() / 2);
  dispatch(bitmap_ids, fn, thread_cnt);
}
//===----------------------------------------------------------------------===//
