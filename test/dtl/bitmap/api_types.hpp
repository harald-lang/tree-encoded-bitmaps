#pragma once
//===----------------------------------------------------------------------===//
#include "gtest/gtest.h"

#include <dtl/bitmap/dynamic_bitmap.hpp>
#include <dtl/bitmap/dynamic_roaring_bitmap.hpp>
#include <dtl/bitmap/dynamic_wah.hpp>
#include <dtl/bitmap/partitioned_position_list.hpp>
#include <dtl/bitmap/partitioned_range_list.hpp>
#include <dtl/bitmap/position_list.hpp>
#include <dtl/bitmap/range_list.hpp>
#include <dtl/bitmap/teb.hpp>
#include <dtl/bitmap/teb_wrapper.hpp>
//===----------------------------------------------------------------------===//
// Types under test.
using teb_v2 = dtl::teb_wrapper;
using teb_o0 = dtl::teb<0>;
using teb_o1 = dtl::teb<1>;
using teb_o2 = dtl::teb<2>;
using teb_o3 = dtl::teb<3>;
// Competitors
using plain_bitmap_32 = dtl::dynamic_bitmap<$u32>;
using roaring_bitmap = dtl::dynamic_roaring_bitmap;
using wah = dtl::dynamic_wah32;
// EXPERIMENTAL
using position_list_32 = dtl::position_list<$u32>;
using partitioned_position_list_32_u8 = dtl::partitioned_position_list<$u32, $u8>;
using partitioned_position_list_32_u16 = dtl::partitioned_position_list<$u32, $u16>;
using range_list_32 = dtl::range_list<$u32>;
using partitioned_range_list_32_u8 = dtl::partitioned_range_list<$u32, $u8>;
using partitioned_range_list_32_u16 = dtl::partitioned_range_list<$u32, $u16>;
//===----------------------------------------------------------------------===//
// The types for which we want to run the API tests.
using types_under_test = ::testing::Types<
    // Tree-encoded Bitmap
    teb_v2,

    // A Tree-encoded Bitmap implementation which allows to enable/disable the
    // space optimizations described in the paper.
    teb_o0,
    teb_o1,
    teb_o2,
    teb_o3,

    // Competitors
    plain_bitmap_32,
    roaring_bitmap,
    wah

    // EXPERIMENTAL
//    position_list_32,
//    partitioned_position_list_32_u8,
//    partitioned_position_list_32_u16,
//    range_list_32,
//    partitioned_range_list_32_u8,
//    partitioned_range_list_32_u16
>;
//===----------------------------------------------------------------------===//
