#pragma once
//===----------------------------------------------------------------------===//
#include <utility>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
// Helper to obtain either a skip or a scan iterator.
//===----------------------------------------------------------------------===//
/// The two kinds of run iterators.
enum class run_iterator_type {
  SKIP,
  SCAN
};
//===----------------------------------------------------------------------===//
/// Base template.
template<typename bitmap_type, run_iterator_type it_type>
struct obtain_run_iterator {};
//===----------------------------------------------------------------------===//
/// Helper to obtain either the (C++) type or an instance.
template<typename bitmap_type>
struct obtain_run_iterator<bitmap_type, run_iterator_type::SKIP> {
  /// The C++ type.
  using type = typename bitmap_type::skip_iter_type;

  /// Create and return an instance.
  static type
  from(bitmap_type& b) noexcept {
    return std::move(b.it());
  }
};
//===----------------------------------------------------------------------===//
/// Helper to obtain either the (C++) type or an instance.
template<typename bitmap_type>
struct obtain_run_iterator<bitmap_type, run_iterator_type::SCAN> {
  /// The C++ type.
  using type = typename bitmap_type::scan_iter_type;

  /// Create and return an instance.
  static type
  from(bitmap_type& b) noexcept {
    return std::move(b.scan_it());
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
