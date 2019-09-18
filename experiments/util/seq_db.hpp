#pragma once
//===----------------------------------------------------------------------===//
#include "sqlite/sqlite3.h"

#include <dtl/bitmap.hpp>
#include <dtl/dtl.hpp>

#include <mutex>
#include <string>
#include <vector>
//===----------------------------------------------------------------------===//
class seq_db {

  const std::string file_;
  std::mutex mutex_;
  sqlite3* db_;
  sqlite3_stmt* insert_def_stmt_;
  sqlite3_stmt* insert_value_stmt_;
  sqlite3_stmt* select_by_id_stmt_;
  sqlite3_stmt* select_ids_stmt_;
  sqlite3_stmt* delete_seq_values_by_id_stmt_;
  sqlite3_stmt* delete_seq_def_by_id_stmt_;
  sqlite3_stmt* count_stmt_;
  sqlite3_stmt* select_all_ids_stmt_;

public:

  using seq_t = std::vector<$u32>;

  explicit seq_db(const std::string& file);
  virtual ~seq_db();

  /// Puts the given sequence in the DB and return the ID (thread-safe).
  i64 put(u64 n, u32 c, f64 f, const seq_t& b);
  /// Loads the sequence with the given ID (thread-safe).
  seq_t get(i64 id);
  /// Returns the IDs of the sequences that match n, f, and d (thread-safe).
  std::vector<$i64> find(u64 n, u32 c, f64 f);
  /// Deletes the sequence with the given ID (thread-safe).
  void remove(i64 id);
  /// Returns the number of sequences in the database.
  std::size_t count();
  /// Returns true if the database is empty.
  u1 empty();
  /// Returns all sequence IDs.
  std::vector<$i64> ids();

private:

  /// Opens the database.
  void open();
  /// Initialize the database schema and prepare the SQL statements.
  void init();
  /// Close the database. Called by the destructor.
  void close();

};
//===----------------------------------------------------------------------===//
