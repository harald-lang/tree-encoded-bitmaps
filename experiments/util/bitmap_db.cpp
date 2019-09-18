#include "bitmap_db.hpp"

#include <dtl/bitmap.hpp>
#include <dtl/bitmap/util/base64.hpp>
#include <dtl/bitmap/util/random.hpp>
#include <dtl/dtl.hpp>

#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
//===----------------------------------------------------------------------===//
bitmap_db::bitmap_db(const std::string& file)
    : file_(file),
      mutex_(),
      db_(nullptr),
      insert_stmt_(nullptr),
      select_by_id_stmt_(nullptr),
      select_ids_stmt_(nullptr),
      delete_by_id_stmt_(nullptr),
      count_stmt_(nullptr),
      select_all_ids_stmt_(nullptr) {
  open();
  init();
}
//===----------------------------------------------------------------------===//
bitmap_db::~bitmap_db() {
  close();
}
//===----------------------------------------------------------------------===//
void
bitmap_db::open() {
  auto rc = sqlite3_open(file_.c_str(), &db_);
  if (rc) {
    std::stringstream err;
    err << "Can't open database: " << file_ << ". Error: "
        << sqlite3_errmsg(db_)
        << std::endl;
    sqlite3_close(db_);
    throw std::invalid_argument(err.str());
  }
}
//===----------------------------------------------------------------------===//
void
bitmap_db::init() {
  auto rc = sqlite3_open(file_.c_str(), &db_);
  if (rc) {
    std::stringstream err;
    err << "Can't open database: " << file_ << ". Error: "
        << sqlite3_errmsg(db_)
        << std::endl;
    sqlite3_close(db_);
    throw std::invalid_argument(err.str());
  }

  const std::string sql_create_table =
      "create table if not exists bitmaps (\n"
          "  id       integer primary key,\n"
          "  n        bigint not null,\n"
          "  d        decimal(7,6) not null,\n"
          "  d_actual decimal(7,6) not null,\n"
          "  f        decimal(7,6) not null,\n"
          "  f_actual decimal(7,6) not null,\n"
          "  bitmap   blob\n"
          ");";

  char* err_msg = nullptr;
  rc = sqlite3_exec(db_, sql_create_table.c_str(),
                    nullptr, nullptr, &err_msg);
  if (rc) {
    std::stringstream err;
    err << "Can't create table 'bitmaps': "
        << sqlite3_errmsg(db_)
        << std::endl;
    close();
    throw std::runtime_error(err.str());
  }

  //===--------------------------------------------------------------------===//
  // Prepare the SQL statements.
  {
    const std::string sql_stmt =
        "insert into bitmaps (n, d, d_actual, f, f_actual, bitmap)\n"
            "  values (:n, :d, :d_actual, :f, :f_actual, :bitmap )";
    rc = sqlite3_prepare_v2(db_, sql_stmt.c_str(), -1, &insert_stmt_, nullptr);
    if (rc) {
      std::stringstream err;
      err << "Can't prepare SQL statement."
          << sqlite3_errmsg(db_)
          << std::endl;
      throw std::runtime_error(err.str());
    }
  }
  {
    const std::string sql_stmt =
        "select * from bitmaps\n"
        " where id = :id";
    rc = sqlite3_prepare_v2(db_, sql_stmt.c_str(), -1, &select_by_id_stmt_,
        nullptr);
    if (rc) {
      std::stringstream err;
      err << "Can't prepare SQL statement."
          << sqlite3_errmsg(db_)
          << std::endl;
      throw std::runtime_error(err.str());
    }
  }
  {
    const std::string sql_stmt =
        "select id from bitmaps\n"
        " where n = :n and f = :f and d = :d"
        " order by id";
    rc = sqlite3_prepare_v2(db_, sql_stmt.c_str(), -1, &select_ids_stmt_,
        nullptr);
    if (rc) {
      std::stringstream err;
      err << "Can't prepare SQL statement."
          << sqlite3_errmsg(db_)
          << std::endl;
      throw std::runtime_error(err.str());
    }
  }
  {
    const std::string sql_stmt =
        "delete from bitmaps\n"
        " where id = :id";
    rc = sqlite3_prepare_v2(db_, sql_stmt.c_str(), -1, &delete_by_id_stmt_,
        nullptr);
    if (rc) {
      std::stringstream err;
      err << "Can't prepare SQL statement."
          << sqlite3_errmsg(db_)
          << std::endl;
      throw std::runtime_error(err.str());
    }
  }
  {
    const std::string sql_stmt = "select count(*) from bitmaps";
    rc = sqlite3_prepare_v2(db_, sql_stmt.c_str(), -1, &count_stmt_,
        nullptr);
    if (rc) {
      std::stringstream err;
      err << "Can't prepare SQL statement."
          << sqlite3_errmsg(db_)
          << std::endl;
      throw std::runtime_error(err.str());
    }
  }
  {
    const std::string sql_stmt = "select id from bitmaps order by id";
    rc = sqlite3_prepare_v2(db_, sql_stmt.c_str(), -1, &select_all_ids_stmt_,
        nullptr);
    if (rc) {
      std::stringstream err;
      err << "Can't prepare SQL statement."
          << sqlite3_errmsg(db_)
          << std::endl;
      throw std::runtime_error(err.str());
    }
  }
  //===--------------------------------------------------------------------===//
}
//===----------------------------------------------------------------------===//
void
bitmap_db::close() {
  if (db_ != nullptr) {

    sqlite3_finalize(insert_stmt_);
    sqlite3_finalize(select_by_id_stmt_);
    sqlite3_finalize(select_ids_stmt_);
    sqlite3_finalize(delete_by_id_stmt_);
    sqlite3_finalize(count_stmt_);
    sqlite3_finalize(select_all_ids_stmt_);

    db_ = nullptr;
    auto* db = db_;
    auto rc = sqlite3_close(db);
    if (rc) {
      std::stringstream err;
      err << "Can't close database. Error: "
          << sqlite3_errmsg(db)
          << std::endl;
      throw std::runtime_error(err.str());
    }
  }
}
//===----------------------------------------------------------------------===//
i64
bitmap_db::store_bitmap(u64 n, f64 f, f64 d, const dtl::bitmap& b) {

  // Determine the actual bit density and clustering factor.
  const auto d_actual = dtl::determine_bit_density(b);
  const auto f_actual = dtl::determine_clustering_factor(b);
  // Base64-encode the bitmap.
  const auto encoded_bitmap = dtl::base64_encode_bitmap(b);

  std::lock_guard<std::mutex> lock(mutex_);
  // Critical section.
  sqlite3_reset(insert_stmt_);
  sqlite3_bind_int64(insert_stmt_,
      sqlite3_bind_parameter_index(insert_stmt_, ":n"), n);
  sqlite3_bind_double(insert_stmt_,
      sqlite3_bind_parameter_index(insert_stmt_, ":d"), d);
  sqlite3_bind_double(insert_stmt_,
      sqlite3_bind_parameter_index(insert_stmt_, ":d_actual"), d_actual);
  sqlite3_bind_double(insert_stmt_,
      sqlite3_bind_parameter_index(insert_stmt_, ":f"), f);
  sqlite3_bind_double(insert_stmt_,
      sqlite3_bind_parameter_index(insert_stmt_, ":f_actual"), f_actual);
  sqlite3_bind_blob64(insert_stmt_,
      sqlite3_bind_parameter_index(insert_stmt_, ":bitmap"),
      encoded_bitmap.data(), encoded_bitmap.size(), SQLITE_TRANSIENT);

  $i32 rc;
  while(true) {
    rc = sqlite3_step(insert_stmt_);
    if (!(rc == SQLITE_LOCKED || rc == SQLITE_BUSY)) {
      break;
    }
    std::this_thread::sleep_for(std::chrono::microseconds(1));
  }
  if (rc != SQLITE_DONE) {
    std::stringstream err;
    err << "Can't insert data. Error: " << rc << " - "
        << sqlite3_errmsg(db_)
        << std::endl;
    throw std::runtime_error(err.str());
  }
  return sqlite3_last_insert_rowid(db_);
}
//===----------------------------------------------------------------------===//
dtl::bitmap
bitmap_db::load_bitmap(i64 id) {
  std::lock_guard<std::mutex> lock(mutex_);
  // Critical section.
  sqlite3_reset(select_by_id_stmt_);
  sqlite3_bind_int64(select_by_id_stmt_,
      sqlite3_bind_parameter_index(select_by_id_stmt_, ":id"), id);

  $i32 rc;
  while((rc = sqlite3_step(select_by_id_stmt_)) == SQLITE_ROW) {
    auto* blob = reinterpret_cast<u8*>(
        sqlite3_column_blob(select_by_id_stmt_, 6));
    const auto bytes = sqlite3_column_bytes(select_by_id_stmt_, 6);
    return dtl::base64_decode_bitmap(dtl::data_view<u8>{blob, blob + bytes});
  }
  std::stringstream err;
  err << "Can't fetch data. Error: " << rc << " - "
      << sqlite3_errmsg(db_)
      << std::endl;
  throw std::runtime_error(err.str());
}
//===----------------------------------------------------------------------===//
std::vector<$i64>
bitmap_db::find_bitmaps(u64 n, f64 f, f64 d) {
  std::lock_guard<std::mutex> lock(mutex_);
  // Critical section.
  sqlite3_reset(select_ids_stmt_);
  sqlite3_bind_int64(select_ids_stmt_,
      sqlite3_bind_parameter_index(select_ids_stmt_, ":n"), n);
  sqlite3_bind_double(select_ids_stmt_,
      sqlite3_bind_parameter_index(select_ids_stmt_, ":d"), d);
  sqlite3_bind_double(select_ids_stmt_,
      sqlite3_bind_parameter_index(select_ids_stmt_, ":f"), f);

  std::vector<$i64> ret;
  $i32 rc;
  while((rc = sqlite3_step(select_ids_stmt_)) == SQLITE_ROW) {
    ret.push_back(sqlite3_column_int64(select_ids_stmt_, 0));
  }
  if (rc != SQLITE_DONE) {
    std::stringstream err;
    err << "Can't fetch data. Error: " << rc << " - "
        << sqlite3_errmsg(db_)
        << std::endl;
    throw std::runtime_error(err.str());
  }
  return ret;
}
//===----------------------------------------------------------------------===//
void
bitmap_db::delete_bitmap(i64 id) {
  std::lock_guard<std::mutex> lock(mutex_);
  // Critical section.
  sqlite3_reset(delete_by_id_stmt_);
  sqlite3_bind_int64(delete_by_id_stmt_,
      sqlite3_bind_parameter_index(delete_by_id_stmt_, ":id"), id);

  $i32 rc;
  while((rc = sqlite3_step(delete_by_id_stmt_)) == SQLITE_ROW) {
  }
  if (rc != SQLITE_DONE) {
    std::stringstream err;
    err << "Can't fetch data. Error: " << rc << " - "
        << sqlite3_errmsg(db_)
        << std::endl;
    throw std::runtime_error(err.str());
  }
}
//===----------------------------------------------------------------------===//
std::size_t
bitmap_db::count() {
  std::lock_guard<std::mutex> lock(mutex_);
  // Critical section.
  sqlite3_reset(count_stmt_);

  std::size_t ret = 0;
  $i32 rc;
  while((rc = sqlite3_step(count_stmt_)) == SQLITE_ROW) {
    ret = sqlite3_column_int64(count_stmt_, 0);
  }
  if (rc != SQLITE_DONE) {
    std::stringstream err;
    err << "Can't fetch data. Error: " << rc << " - "
        << sqlite3_errmsg(db_)
        << std::endl;
    throw std::runtime_error(err.str());
  }
  return ret;
}
//===----------------------------------------------------------------------===//
u1
bitmap_db::empty() {
  return count() == 0;
}
//===----------------------------------------------------------------------===//
std::vector<$i64>
bitmap_db::ids() {
  std::lock_guard<std::mutex> lock(mutex_);
  // Critical section.
  sqlite3_reset(select_all_ids_stmt_);

  std::vector<$i64> ret;
  $i32 rc;
  while((rc = sqlite3_step(select_all_ids_stmt_)) == SQLITE_ROW) {
    ret.push_back(sqlite3_column_int64(select_all_ids_stmt_, 0));
  }
  if (rc != SQLITE_DONE) {
    std::stringstream err;
    err << "Can't fetch data. Error: " << rc << " - "
        << sqlite3_errmsg(db_)
        << std::endl;
    throw std::runtime_error(err.str());
  }
  return ret;
}
//===----------------------------------------------------------------------===//
