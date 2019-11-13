#pragma once
//===----------------------------------------------------------------------===//
#include "gtest/gtest.h"

#include <dtl/bitmap/bah.hpp>
#include <dtl/bitmap/diff/diff.hpp>
#include <dtl/bitmap/diff/merge.hpp>
#include <dtl/bitmap/dynamic_bitmap.hpp>
#include <dtl/bitmap/dynamic_roaring_bitmap.hpp>
#include <dtl/bitmap/dynamic_wah.hpp>
#include <dtl/bitmap/part/part.hpp>
#include <dtl/bitmap/part/part_run.hpp>
#include <dtl/bitmap/part/part_updirect.hpp>
#include <dtl/bitmap/part/part_upforward.hpp>
#include <dtl/bitmap/position_list.hpp>
#include <dtl/bitmap/range_list.hpp>
#include <dtl/bitmap/teb.hpp>
#include <dtl/bitmap/teb_wrapper.hpp>
#include <dtl/bitmap/uah.hpp>
#include <dtl/bitmap/uah_skip.hpp>
#include <dtl/bitmap/xah.hpp>
#include <dtl/bitmap/xah_skip.hpp>
//===----------------------------------------------------------------------===//
// Types under test.
using teb_v2 = dtl::teb_wrapper;
using teb_o0 = dtl::teb<0>; // DEPRECATED
using teb_o1 = dtl::teb<1>; // DEPRECATED
using teb_o2 = dtl::teb<2>; // DEPRECATED
using teb_o3 = dtl::teb<3>; // DEPRECATED

// Competitors
using plain_bitmap_32 = dtl::dynamic_bitmap<$u32>;
using roaring_bitmap = dtl::dynamic_roaring_bitmap;
using wah = dtl::dynamic_wah32;

// Differential (EXPERIMENTAL)
using diff_teb_roaring = dtl::diff<teb_v2, roaring_bitmap>;
using diff_teb_wah = dtl::diff<teb_v2, wah>;
using diff_roaring_roaring = dtl::diff<roaring_bitmap, roaring_bitmap>;
using diff_roaring_wah = dtl::diff<roaring_bitmap, wah>;
using diff_wah_roaring = dtl::diff<wah, roaring_bitmap>;
using diff_wah_wah = dtl::diff<wah, wah>;

// Partitioned (EXPERIMENTAL)
// - Updates are handled directly (by re-compressing affected partitions)..
using part_8_teb = dtl::part_updirect<teb_v2, 1ull << 8>;
using part_16_teb = dtl::part_updirect<teb_v2, 1ull << 16>;
using part_8_wah = dtl::part_updirect<wah, 1ull << 8>;
using part_16_wah = dtl::part_updirect<wah, 1ull << 16>;
// - Updates are forwarded to the internal type.
using part_8_upfwd_wah = dtl::part_upforward<wah, 1ull << 8>;
using part_16_upfwd_wah = dtl::part_upforward<wah, 1ull << 16>;
using part_8_upfwd_diff_teb = dtl::part_upforward<diff_teb_roaring, 1ull << 8>;
using part_16_upfwd_diff_teb = dtl::part_upforward<diff_teb_roaring, 1ull << 16>;

// Partitioned, with single-value optimization (EXPERIMENTAL)
using part_run_8_teb = dtl::part_run<teb_v2, 1ull << 8>;
using part_run_16_teb = dtl::part_run<teb_v2, 1ull << 16>;
using part_run_8_wah = dtl::part_run<wah, 1ull << 8>;
using part_run_16_wah = dtl::part_run<wah, 1ull << 16>;

// Partitioned + Differential (EXPERIMENTAL, using Roaring for diffs)
using diff_part_8_teb = dtl::diff<part_8_teb, roaring_bitmap>;
using part_8_diff_teb = dtl::part<diff_teb_roaring, 1ull << 8>;

// EXPERIMENTAL
using position_list_32 = dtl::position_list<$u32>;
using part_position_list_8 = dtl::part_upforward<dtl::position_list<$u8>, 1ull << 8>;
using part_position_list_16 = dtl::part_upforward<dtl::position_list<$u16>, 1ull << 16>;

using range_list_32 = dtl::range_list<$u32>;
using part_range_list_8 = dtl::part<dtl::range_list<$u8>, 1ull << 8>;
using part_range_list_16 = dtl::part<dtl::range_list<$u16>, 1ull << 16>;
//===----------------------------------------------------------------------===//
// The types for which we want to run the API tests.
using types_under_test = ::testing::Types<
    // Tree-encoded Bitmap
    teb_v2,

    // A Tree-encoded Bitmap implementation which allows to enable/disable the
    // space optimizations described in the paper.
    // teb_o0,
    // teb_o1,
    // teb_o2,
    teb_o3,

    // Competitors
    plain_bitmap_32,
    roaring_bitmap,
    wah,

    // Differential (EXPERIMENTAL)
    diff_teb_roaring,
    diff_teb_wah,
    diff_roaring_roaring,
    diff_roaring_wah,
    diff_wah_roaring,
    diff_wah_wah,

    // Partitioned (EXPERIMENTAL)
    part_8_teb,
    part_run_8_teb,
    part_8_wah,
    part_run_8_wah,

    // Partitioned + Differential (EXPERIMENTAL)
    diff_part_8_teb,
    part_8_diff_teb,
    // EXPERIMENTAL
    dtl::uah8,
    dtl::uah16,
    dtl::uah32,
    dtl::uah64,
    dtl::uah_skip<u8, 2>, // Skip distance is intentionally chosen small, as the bitmaps in the test are also rather small.
    dtl::uah_skip<u16, 2>,
    dtl::uah_skip<u32, 2>,
    dtl::uah_skip<u64, 2>,
    dtl::xah8,
    dtl::xah16,
    dtl::xah32,
    dtl::xah64,
    dtl::xah_skip<u8, 2>, // Skip distance is intentionally chosen small, as the bitmaps in the test are also rather small.
    dtl::xah_skip<u16, 2>,
    dtl::xah_skip<u32, 2>,
    dtl::xah_skip<u64, 2>,
    position_list_32,
    part_position_list_8,
    part_position_list_16,
    range_list_32,
    part_range_list_8,
    part_range_list_16,

    dtl::bah,
    dtl::part<dtl::bah, 1ull << 16>
    >;
//===----------------------------------------------------------------------===//
