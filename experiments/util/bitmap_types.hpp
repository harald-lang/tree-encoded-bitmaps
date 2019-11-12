#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/bitmap/dynamic_bitmap.hpp>
#include <dtl/bitmap/dynamic_roaring_bitmap.hpp>
#include <dtl/bitmap/dynamic_wah.hpp>
#include <dtl/bitmap/part/part.hpp>
#include <dtl/bitmap/position_list.hpp>
#include <dtl/bitmap/range_list.hpp>
#include <dtl/bitmap/teb.hpp>
#include <dtl/bitmap/teb_wrapper.hpp>
#include <dtl/dtl.hpp>

#include <dtl/bitmap/uah.hpp>
#include <dtl/bitmap/uah_skip.hpp>
#include <dtl/bitmap/xah_skip.hpp>
#include <iostream>
#include <vector>
//===----------------------------------------------------------------------===//
enum class bitmap_t {
  bitmap,
  roaring,
  teb, /* deprecated */
  teb_scan, /* deprecated */
  teb_wrapper,
  wah,
  position_list,
  partitioned_position_list_u8,
  partitioned_position_list_u16,
  range_list,
  partitioned_range_list_u8,
  partitioned_range_list_u16,

  uah8,
  uah8_skip,
  uah16,
  uah16_skip,
  uah32,
  uah32_skip,
  uah64,
  uah64_skip,

  partitioned_uah8,
  partitioned_uah8_skip,
  partitioned_uah16,
  partitioned_uah16_skip,
  partitioned_uah32,
  partitioned_uah32_skip,
  partitioned_uah64,
  partitioned_uah64_skip,

  xah8,
  xah8_skip,
  xah16,
  xah16_skip,
  xah32,
  xah32_skip,
  xah64,
  xah64_skip,

  partitioned_xah8,
  partitioned_xah8_skip,
  partitioned_xah16,
  partitioned_xah16_skip,
  partitioned_xah32,
  partitioned_xah32_skip,
  partitioned_xah64,
  partitioned_xah64_skip,

  partitioned_teb_wrapper,
  partitioned_wah,

  _first = bitmap,
  _last = partitioned_wah
};
static const std::vector<bitmap_t>
    bitmap_t_list = []() {
      std::vector<bitmap_t> l;
      for (auto bitmap_type = static_cast<int>(bitmap_t::_first);
           bitmap_type <= static_cast<int>(bitmap_t::_last);
           ++bitmap_type) {
        l.push_back(static_cast<bitmap_t>(bitmap_type));
      }
      return l;
    }();
//===----------------------------------------------------------------------===//
std::vector<std::string> bitmap_names { // TODO remove
  "bitmap",
  "roaring",
  "teb",
  "teb_scan", /* deprecated */
  "teb_wrapper",
  "wah",
  "position_list",
  "partitioned_position_list_u8",
  "partitioned_position_list_u16",
  "range_list",
  "partitioned_range_list_u8",
  "partitioned_range_list_u16",
  "uah8",
  "uah8_skip",
  "uah16",
  "uah16_skip",
  "uah32",
  "uah32_skip",
  "uah64",
  "uah64_skip",
  "partitioned_uah8",
  "partitioned_uah8_skip",
  "partitioned_uah16",
  "partitioned_uah16_skip",
  "partitioned_uah32",
  "partitioned_uah32_skip",
  "partitioned_uah64",
  "partitioned_uah64_skip",
  "xah8",
  "xah8_skip",
  "xah16",
  "xah16_skip",
  "xah32",
  "xah32_skip",
  "xah64",
  "xah64_skip",
  "partitioned_xah8",
  "partitioned_xah8_skip",
  "partitioned_xah16",
  "partitioned_xah16_skip",
  "partitioned_xah32",
  "partitioned_xah32_skip",
  "partitioned_xah64",
  "partitioned_xah64_skip",
};
std::ostream& operator<<(std::ostream& out, const bitmap_t& b) {
  const auto i = static_cast<int>(b);
  assert(i < bitmap_names.size());
  out << bitmap_names[i];
  return out;
}
//===----------------------------------------------------------------------===//
template<bitmap_t B>
struct type_of {
  using type = void;
};
template<>
struct type_of<bitmap_t::bitmap> {
  using type = dtl::dynamic_bitmap<$u32>;
};
template<>
struct type_of<bitmap_t::roaring> {
  using type = dtl::dynamic_roaring_bitmap;
};
//template<>
//struct type_of<bitmap_t::teb> { /* deprecated */
//  using type = dtl::teb<>;
//};
//template<> struct type_of<bitmap_t::teb_scan> { /* deprecated */
//  using type = dtl::teb_scan<>; };
template<>
struct type_of<bitmap_t::teb_wrapper> {
  using type = dtl::teb_wrapper;
};
template<>
struct type_of<bitmap_t::partitioned_teb_wrapper> {
  using type = dtl::part<dtl::teb_wrapper, 1ull << 16>;
};
template<>
struct type_of<bitmap_t::wah> {
  using type = dtl::dynamic_wah32;
};
template<>
struct type_of<bitmap_t::partitioned_wah> {
  using type = dtl::part<dtl::dynamic_wah32, 1ull << 16>;
};
template<>
struct type_of<bitmap_t::position_list> {
  using type = dtl::position_list<$u32>;
};
template<>
struct type_of<bitmap_t::partitioned_position_list_u8> {
  using type = dtl::part<dtl::position_list<$u8>, 1ull << 8>;
};
template<>
struct type_of<bitmap_t::partitioned_position_list_u16> {
  using type = dtl::part<dtl::position_list<$u16>, 1ull << 16>;
};
template<>
struct type_of<bitmap_t::range_list> {
  using type = dtl::range_list<$u32>;
};
template<>
struct type_of<bitmap_t::partitioned_range_list_u8> {
  using type = dtl::part<dtl::range_list<$u8>, 1ull << 8>;
};
template<>
struct type_of<bitmap_t::partitioned_range_list_u16> {
  using type = dtl::part<dtl::range_list<$u16>, 1ull << 16>;
};

template<>
struct type_of<bitmap_t::uah8> {
  using type = dtl::uah8;
};
template<>
struct type_of<bitmap_t::uah8_skip> {
  using type = dtl::uah8_skip;
};
template<>
struct type_of<bitmap_t::uah16> {
  using type = dtl::uah16;
};
template<>
struct type_of<bitmap_t::uah16_skip> {
  using type = dtl::uah16_skip;
};
template<>
struct type_of<bitmap_t::uah32> {
  using type = dtl::uah32;
};
template<>
struct type_of<bitmap_t::uah32_skip> {
  using type = dtl::uah32_skip;
};
template<>
struct type_of<bitmap_t::uah64> {
  using type = dtl::uah64;
};
template<>
struct type_of<bitmap_t::uah64_skip> {
  using type = dtl::uah64_skip;
};

template<>
struct type_of<bitmap_t::partitioned_uah8> {
  using type = dtl::part<dtl::uah8, 1ull << 16>;
};
template<>
struct type_of<bitmap_t::partitioned_uah8_skip> {
  using type = dtl::part<dtl::uah8_skip, 1ull << 16>;
};
template<>
struct type_of<bitmap_t::partitioned_uah16> {
  using type = dtl::part<dtl::uah16, 1ull << 16>;
};
template<>
struct type_of<bitmap_t::partitioned_uah16_skip> {
  using type = dtl::part<dtl::uah16_skip, 1ull << 16>;
};
template<>
struct type_of<bitmap_t::partitioned_uah32> {
  using type = dtl::part<dtl::uah32, 1ull << 16>;
};
template<>
struct type_of<bitmap_t::partitioned_uah32_skip> {
  using type = dtl::part<dtl::uah32_skip, 1ull << 16>;
};
template<>
struct type_of<bitmap_t::partitioned_uah64> {
  using type = dtl::part<dtl::uah64, 1ull << 16>;
};
template<>
struct type_of<bitmap_t::partitioned_uah64_skip> {
  using type = dtl::part<dtl::uah64_skip, 1ull << 16>;
};

template<>
struct type_of<bitmap_t::xah8> {
  using type = dtl::xah8;
};
template<>
struct type_of<bitmap_t::xah8_skip> {
  using type = dtl::xah8_skip;
};
template<>
struct type_of<bitmap_t::xah16> {
  using type = dtl::xah16;
};
template<>
struct type_of<bitmap_t::xah16_skip> {
  using type = dtl::xah16_skip;
};
template<>
struct type_of<bitmap_t::xah32> {
  using type = dtl::xah32;
};
template<>
struct type_of<bitmap_t::xah32_skip> {
  using type = dtl::xah32_skip;
};
template<>
struct type_of<bitmap_t::xah64> {
  using type = dtl::xah64;
};
template<>
struct type_of<bitmap_t::xah64_skip> {
  using type = dtl::xah64_skip;
};

template<>
struct type_of<bitmap_t::partitioned_xah8> {
  using type = dtl::part<dtl::xah8, 1ull << 16>;
};
template<>
struct type_of<bitmap_t::partitioned_xah8_skip> {
  using type = dtl::part<dtl::xah8_skip, 1ull << 16>;
};
template<>
struct type_of<bitmap_t::partitioned_xah16> {
  using type = dtl::part<dtl::xah16, 1ull << 16>;
};
template<>
struct type_of<bitmap_t::partitioned_xah16_skip> {
  using type = dtl::part<dtl::xah16_skip, 1ull << 16>;
};
template<>
struct type_of<bitmap_t::partitioned_xah32> {
  using type = dtl::part<dtl::xah32, 1ull << 16>;
};
template<>
struct type_of<bitmap_t::partitioned_xah32_skip> {
  using type = dtl::part<dtl::xah32_skip, 1ull << 16>;
};
template<>
struct type_of<bitmap_t::partitioned_xah64> {
  using type = dtl::part<dtl::xah64, 1ull << 16>;
};
template<>
struct type_of<bitmap_t::partitioned_xah64_skip> {
  using type = dtl::part<dtl::xah64_skip, 1ull << 16>;
};
//===----------------------------------------------------------------------===//
