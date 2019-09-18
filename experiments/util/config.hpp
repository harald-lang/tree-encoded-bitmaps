#pragma once
//===----------------------------------------------------------------------===//
#include "bitmap_types.hpp"

#include <dtl/dtl.hpp>

#include <ostream>
//===----------------------------------------------------------------------===//
struct config {
  bitmap_t bitmap_type;
  $u64 n;
  $f64 density;
  $f64 clustering_factor;
  $i64 bitmap_id; // in DB

  bool
  operator<(const config& o) const {
    const auto t = static_cast<$u32>(bitmap_type);
    const auto ot = static_cast<$u32>(o.bitmap_type);
    return t < ot
        || (t == ot && n  < o.n)
        || (t == ot && n == o.n && density  < o.density)
        || (t == ot && n == o.n && density == o.density && clustering_factor  < o.clustering_factor)
        ;
  }

  void
  print(std::ostream& os) const noexcept {
    os << "config[bitmap_type=" << bitmap_type
        << ",n=" << n
        << ",d=" << density
        << ",f=" << clustering_factor
        << ",bitmap_id=" << bitmap_id
        << "]";
  }
};
//===----------------------------------------------------------------------===//
struct config_pair {
  bitmap_t bitmap_type = bitmap_t::bitmap;
  $u64 n = 0;
  $f64 density1 = 0;
  $f64 density2 = 0;
  $f64 clustering_factor1 = 0;
  $f64 clustering_factor2 = 0;
  $i64 bitmap_id1 = 0; // in DB
  $i64 bitmap_id2 = 0; // in DB

  bool
  operator<(const config_pair& o) const {
    const auto t = static_cast<uint32_t>(bitmap_type);
    const auto ot = static_cast<uint32_t>(o.bitmap_type);
    return t < ot
        || (t == ot && n  < o.n)
        || (t == ot && n == o.n && density1  < o.density1)
        || (t == ot && n == o.n && density1 == o.density1 && clustering_factor1  < o.clustering_factor1)
        || (t == ot && n == o.n && density1 == o.density1 && clustering_factor1 == o.clustering_factor1 && density2  < o.density2)
        || (t == ot && n == o.n && density1 == o.density1 && clustering_factor1 == o.clustering_factor1 && density2 == o.density2 && clustering_factor2 < o.clustering_factor2);
  }

  void
  print(std::ostream& os) const noexcept {
    os << "config_pair[bitmap_type=" << bitmap_type
        << ",n=" << n
        << ",d1=" << density1
        << ",d2=" << density2
        << ",f1=" << clustering_factor1
        << ",f2=" << clustering_factor2
        << ",bitmap_id1=" << bitmap_id1
        << ",bitmap_id2=" << bitmap_id2
        << "]";
  }

  config
  first() {
    config ret_val;
    ret_val.n = n;
    ret_val.density = density1;
    ret_val.clustering_factor = clustering_factor1;
    ret_val.bitmap_type = bitmap_type;
    ret_val.bitmap_id = bitmap_id1;
    return ret_val;
  }

  config
  second() {
    config ret_val;
    ret_val.n = n;
    ret_val.density = density2;
    ret_val.clustering_factor = clustering_factor2;
    ret_val.bitmap_type = bitmap_type;
    ret_val.bitmap_id = bitmap_id2;
    return ret_val;
  }

};
//===----------------------------------------------------------------------===//
