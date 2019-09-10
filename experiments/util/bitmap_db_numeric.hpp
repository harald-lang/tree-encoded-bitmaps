#pragma once

#include <mutex>
#include <string>
#include <dtl/dtl.hpp>
#include <dtl/bitmap.hpp>
#include <sqlite/sqlite3.h>
#include "types.hpp"

class bitmap_db {

  const std::string file_;
  std::mutex mutex_;
  sqlite3* db_;
  sqlite3_stmt* insert_stmt_;
  sqlite3_stmt* select_by_id_stmt_;
  sqlite3_stmt* select_ids_stmt_;
  sqlite3_stmt* delete_by_id_stmt_;
  sqlite3_stmt* count_stmt_;
  sqlite3_stmt* select_all_ids_stmt_;

public:

  using num = numeric<6>;
  explicit bitmap_db(const std::string& file);
  virtual ~bitmap_db();

  /// Puts the given bitmap in the DB and return the ID (thread-safe).
  i64 store_bitmap(u64 n, num f, num d, const dtl::bitmap& b);
  i64 store_bitmap(u64 n, num d, const dtl::bitmap& b) {
    return store_bitmap(n, 0.0, d, b);
  };
  /// Loads the bitmap with the given ID (thread-safe).
  dtl::bitmap load_bitmap(i64 id);
  /// Returns the IDs of the bitmaps that match n, f, and d (thread-safe).
  std::vector<$i64> find_bitmaps(u64 n, num f, num d);
  std::vector<$i64> find_bitmaps(u64 n, num d) {
    return find_bitmaps(n, num(0), d);
  };
  /// Deletes the bitmap with the given ID (thread-safe).
  void delete_bitmap(i64 id);
  /// Returns the number of bitmaps in the database.
  std::size_t count();
  /// Returns true if the database is empty.
  u1 empty();
  /// Returns all bitmap IDs.
  std::vector<$i64> ids();

private:

  /// Opens the database.
  void open();
  /// Initialize the database schema and prepare the SQL statements.
  void init();
  /// Close the database. Called by the constructor.
  void close();

};