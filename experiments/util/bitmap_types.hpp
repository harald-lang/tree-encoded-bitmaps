#pragma once

#include <vector>

#include <dtl/dtl.hpp>
#include <dtl/bitmap/dynamic_bitmap.hpp>
#include <dtl/bitmap/dynamic_roaring_bitmap.hpp>
#include <dtl/bitmap/teb.hpp>
#include <dtl/bitmap/teb_scan.hpp>
#include <dtl/bitmap/dynamic_wah.hpp>
#include <dtl/bitmap/position_list.hpp>
#include <dtl/bitmap/partitioned_position_list.hpp>
#include <dtl/bitmap/partitioned_range_list.hpp>
#include <dtl/bitmap/range_list.hpp>

//===----------------------------------------------------------------------===//
enum class bitmap_t {
  bitmap,
  roaring,
  teb,
  teb_scan,
  wah,
  position_list,
  partitioned_position_list_u8,
  partitioned_position_list_u16,
  range_list,
  partitioned_range_list_u8,
  partitioned_range_list_u16,

  _first = bitmap,
  _last = partitioned_range_list_u16
};
static const std::vector<bitmap_t>
    bitmap_t_list = [&]() {
        std::vector<bitmap_t> l;
        for (auto bitmap_type = static_cast<int>(bitmap_t::_first);
             bitmap_type <= static_cast<int>(bitmap_t::_last);
             ++bitmap_type) {
          l.push_back(static_cast<bitmap_t>(bitmap_type));
        }
        return l;
      }();
//===----------------------------------------------------------------------===//
std::vector<std::string> bitmap_names {
    "bitmap",
    "roaring",
    "teb",
    "teb_scan",
    "wah",
    "position_list",
    "partitioned_position_list_u8",
    "partitioned_position_list_u16",
    "range_list",
    "partitioned_range_list_u8",
    "partitioned_range_list_u16",
};
std::ostream& operator<<(std::ostream& out, const bitmap_t& b) {
  const auto i = static_cast<int>(b);
  assert(i < bitmap_names.size());
  out << bitmap_names[i];
  return out;
}
//===----------------------------------------------------------------------===//
template<bitmap_t B> struct type_of {
  using type = void; };
template<> struct type_of<bitmap_t::bitmap> {
  using type = dtl::dynamic_bitmap<$u32>; };
template<> struct type_of<bitmap_t::roaring> {
  using type = dtl::dynamic_roaring_bitmap; };
template<> struct type_of<bitmap_t::teb> {
  using type = dtl::teb<>; };
template<> struct type_of<bitmap_t::teb_scan> {
  using type = dtl::teb_scan<>; };
template<> struct type_of<bitmap_t::wah> {
  using type = dtl::dynamic_wah32; };
template<> struct type_of<bitmap_t::position_list> {
  using type = dtl::position_list<$u32>; };
template<> struct type_of<bitmap_t::partitioned_position_list_u8> {
  using type = dtl::partitioned_position_list<$u32, $u8>; };
template<> struct type_of<bitmap_t::partitioned_position_list_u16> {
  using type = dtl::partitioned_position_list<$u32, $u16>; };
template<> struct type_of<bitmap_t::range_list> {
  using type = dtl::range_list<$u32>; };
template<> struct type_of<bitmap_t::partitioned_range_list_u8> {
  using type = dtl::partitioned_range_list<$u32, $u8>; };
template<> struct type_of<bitmap_t::partitioned_range_list_u16> {
  using type = dtl::partitioned_range_list<$u32, $u16>; };
//===----------------------------------------------------------------------===//
