#include "seq_db.hpp"

#include <dtl/bitmap.hpp>
#include <dtl/bitmap/util/base64.hpp>
#include <dtl/bitmap/util/random.hpp>
#include <dtl/dtl.hpp>
#include <dtl/math.hpp>

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
//===----------------------------------------------------------------------===//
seq_db::seq_db(const std::string& file)
    : file_(file),
      mutex_(),
      db_(nullptr),
      insert_def_stmt_(nullptr),
      select_by_id_stmt_(nullptr),
      select_ids_stmt_(nullptr),
      delete_seq_values_by_id_stmt_(nullptr),
      count_stmt_(nullptr),
      select_all_ids_stmt_(nullptr) {
  open();
  init();
}
//===----------------------------------------------------------------------===//
seq_db::~seq_db() {
  close();
}
//===----------------------------------------------------------------------===//
void
seq_db::open() {
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
seq_db::init() {
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
      "create table if not exists seq_def (\n"
      "  id integer primary key,\n"
      "  n bigint not null,\n"
      "  c int not null,\n"
      "  f decimal(7,6) not null,\n"
      "  f_actual decimal(7,6) not null\n"
      ");\n"
      "create table if not exists seq_data (\n"
      "  id int not null,\n"
//      "  value_seq_no bigint not null,\n"
//      "  value int not null,\n"
      "  data blob,\n"
      "  foreign key (id) references seq_def(id)\n"
      ");";

  char* err_msg = nullptr;
  rc = sqlite3_exec(db_, sql_create_table.c_str(),
                    nullptr, nullptr, &err_msg);
  if (rc) {
    std::stringstream err;
    err << "Can't create table: "
        << sqlite3_errmsg(db_)
        << std::endl;
    close();
    throw std::runtime_error(err.str());
  }

  //===--------------------------------------------------------------------===//
  // Prepare the SQL statements.
  {
    const std::string sql_stmt =
        "insert into seq_def (n, c, f, f_actual)\n"
         "  values (:n, :c, :f, :f_actual )";
    rc = sqlite3_prepare_v2(db_, sql_stmt.c_str(), -1, &insert_def_stmt_, nullptr);
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
        "insert into seq_data (id, data)\n"
         "  values (:id, :data )";
//        "insert into seq_data (id, value_seq_no, value)\n"
//         "  values (:id, :value_seq_no, :value )";
    rc = sqlite3_prepare_v2(db_, sql_stmt.c_str(), -1, &insert_value_stmt_, nullptr);
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
        "select * from seq_data\n"
        " where id = :id\n";
//        " order by value_seq_no";
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
        "select id from seq_def\n"
        " where n = :n and c = :c and f = :f";
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
        "delete from seq_data\n"
        " where id = :id";
    rc = sqlite3_prepare_v2(db_, sql_stmt.c_str(), -1, &delete_seq_values_by_id_stmt_,
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
        "delete from seq_def\n"
        " where id = :id";
    rc = sqlite3_prepare_v2(db_, sql_stmt.c_str(), -1, &delete_seq_def_by_id_stmt_,
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
    const std::string sql_stmt = "select count(*) from seq_def";
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
    const std::string sql_stmt = "select id from seq_def order by id";
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
seq_db::close() {
  if (db_ != nullptr) {

    sqlite3_finalize(insert_def_stmt_);
    sqlite3_finalize(insert_value_stmt_);
    sqlite3_finalize(select_by_id_stmt_);
    sqlite3_finalize(select_ids_stmt_);
    sqlite3_finalize(delete_seq_values_by_id_stmt_);
    sqlite3_finalize(delete_seq_def_by_id_stmt_);
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
seq_db::put(u64 n, u32 c, f64 f, const seq_t& s) {

  // Determine the actual bit density and clustering factor.
  const auto f_actual = dtl::determine_clustering_factor(s);

  std::lock_guard<std::mutex> lock(mutex_);
  sqlite3_exec(db_, "BEGIN", nullptr, nullptr, nullptr);

  // Critical section.
  sqlite3_reset(insert_def_stmt_);
  sqlite3_bind_int64(insert_def_stmt_,
      sqlite3_bind_parameter_index(insert_def_stmt_, ":n"), n);
  sqlite3_bind_double(insert_def_stmt_,
      sqlite3_bind_parameter_index(insert_def_stmt_, ":c"), c);
  sqlite3_bind_double(insert_def_stmt_,
      sqlite3_bind_parameter_index(insert_def_stmt_, ":f"), f);
  sqlite3_bind_double(insert_def_stmt_,
      sqlite3_bind_parameter_index(insert_def_stmt_, ":f_actual"), f_actual);
  $i32 rc;
  while(true) {
    rc = sqlite3_step(insert_def_stmt_);
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
    sqlite3_exec(db_, "ROLLBACK", nullptr, nullptr, nullptr);
    throw std::runtime_error(err.str());
  }
  i64 seq_id = sqlite3_last_insert_rowid(db_);

  sqlite3_reset(insert_value_stmt_);
  sqlite3_bind_int64(insert_value_stmt_,
      sqlite3_bind_parameter_index(insert_value_stmt_, ":id"), seq_id);
//    sqlite3_bind_int64(insert_value_stmt_,
//        sqlite3_bind_parameter_index(insert_value_stmt_, ":value_seq_no"), i);
//    sqlite3_bind_int64(insert_value_stmt_,
//        sqlite3_bind_parameter_index(insert_value_stmt_, ":value"), s[i]);
  sqlite3_bind_blob64(insert_value_stmt_,
      sqlite3_bind_parameter_index(insert_value_stmt_, ":data"),
      reinterpret_cast<const void*>(s.data()), s.size() * 4, SQLITE_TRANSIENT);

  while(true) {
    rc = sqlite3_step(insert_value_stmt_);
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
    sqlite3_exec(db_, "ROLLBACK", nullptr, nullptr, nullptr);
    throw std::runtime_error(err.str());
  }
  sqlite3_exec(db_, "COMMIT", nullptr, nullptr, nullptr);
  return seq_id;
}
//===----------------------------------------------------------------------===//
seq_db::seq_t
seq_db::get(i64 id) {
  std::lock_guard<std::mutex> lock(mutex_);
  // Critical section.
  sqlite3_reset(select_by_id_stmt_);
  sqlite3_bind_int64(select_by_id_stmt_,
      sqlite3_bind_parameter_index(select_by_id_stmt_, ":id"), id);

  $i32 rc;
  seq_t ret;
  while((rc = sqlite3_step(select_by_id_stmt_)) == SQLITE_ROW) {
    auto* blob = reinterpret_cast<u8*>(
        sqlite3_column_blob(select_by_id_stmt_, 1));
    const auto bytes = sqlite3_column_bytes(select_by_id_stmt_, 1);
    dtl::data_view<u32> ints {
        reinterpret_cast<u32*>(blob),
        reinterpret_cast<u32*>(blob + bytes)
    };
    seq_t ret(ints.begin(), ints.end());
    return ret;
  }
  std::stringstream err;
  err << "Can't fetch data. Error: " << rc << " - "
      << sqlite3_errmsg(db_)
      << std::endl;
  throw std::runtime_error(err.str());
}
//===----------------------------------------------------------------------===//
std::vector<$i64>
seq_db::find(u64 n, u32 c, f64 f) {
  std::lock_guard<std::mutex> lock(mutex_);
  // Critical section.
  sqlite3_reset(select_ids_stmt_);
  sqlite3_bind_int64(select_ids_stmt_,
      sqlite3_bind_parameter_index(select_ids_stmt_, ":n"), n);
  sqlite3_bind_double(select_ids_stmt_,
      sqlite3_bind_parameter_index(select_ids_stmt_, ":c"), c);
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
seq_db::remove(i64 id) {
  std::lock_guard<std::mutex> lock(mutex_);
  // Critical section.
  sqlite3_reset(delete_seq_values_by_id_stmt_);
  sqlite3_bind_int64(delete_seq_values_by_id_stmt_,
      sqlite3_bind_parameter_index(delete_seq_values_by_id_stmt_, ":id"), id);

  $i32 rc;
  while((rc = sqlite3_step(delete_seq_values_by_id_stmt_)) == SQLITE_ROW) {
  }
  if (rc != SQLITE_DONE) {
    std::stringstream err;
    err << "Can't fetch data. Error: " << rc << " - "
        << sqlite3_errmsg(db_)
        << std::endl;
    throw std::runtime_error(err.str());
  }

  sqlite3_reset(delete_seq_def_by_id_stmt_);
  sqlite3_bind_int64(delete_seq_def_by_id_stmt_,
      sqlite3_bind_parameter_index(delete_seq_def_by_id_stmt_, ":id"), id);

  while((rc = sqlite3_step(delete_seq_def_by_id_stmt_)) == SQLITE_ROW) {
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
seq_db::count() {
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
seq_db::empty() {
  return count() == 0;
}
//===----------------------------------------------------------------------===//
std::vector<$i64>
seq_db::ids() {
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
