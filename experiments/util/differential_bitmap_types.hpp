#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/bitmap/diff/diff.hpp>
#include <dtl/bitmap/diff/merge.hpp>
#include <dtl/bitmap/diff/merge_teb.hpp>
#include <dtl/bitmap/dynamic_roaring_bitmap.hpp>
#include <dtl/bitmap/dynamic_wah.hpp>
#include <dtl/bitmap/part/part_updirect.hpp>
#include <dtl/bitmap/part/part_upforward.hpp>
#include <dtl/bitmap/teb_wrapper.hpp>
//===----------------------------------------------------------------------===//
using partitioned_teb = dtl::part_updirect<dtl::teb_wrapper, 1ull << 16>;
using partitioned_teb_diff_structure = dtl::dynamic_roaring_bitmap;
using partitioned_differential_teb =
dtl::part_upforward<dtl::diff<dtl::teb_wrapper,partitioned_teb_diff_structure>,
    1ull << 16>;

using partitioned_wah = dtl::part<dtl::dynamic_wah32, 1ull << 16>;
using partitioned_wah_diff_structure = dtl::dynamic_roaring_bitmap;
//using partitioned_wah_diff_structure = dtl::dynamic_wah32; // WAH disqualifies as diff structure, in general
using partitioned_differential_wah =
dtl::part_upforward<dtl::diff<dtl::dynamic_wah32,partitioned_wah_diff_structure>,
    1ull << 16>;
//===----------------------------------------------------------------------===//
#define __DIFF_NAME_ROARING_WAH roaring_wah
#define __DIFF_TYPE_ROARING_WAH dtl::diff<dtl::dynamic_roaring_bitmap, dtl::dynamic_wah32>

#define __DIFF_NAME_ROARING_ROARING roaring_roaring
#define __DIFF_TYPE_ROARING_ROARING dtl::diff<dtl::dynamic_roaring_bitmap, dtl::dynamic_roaring_bitmap>

#define __DIFF_NAME_TEB_WAH teb_wah
#define __DIFF_TYPE_TEB_WAH dtl::diff<dtl::teb_wrapper, dtl::dynamic_wah32>

#define __DIFF_NAME_TEB_ROARING teb_roaring
#define __DIFF_TYPE_TEB_ROARING dtl::diff<dtl::teb_wrapper, dtl::dynamic_roaring_bitmap>

#define __DIFF_NAME_WAH_WAH wah_wah
#define __DIFF_TYPE_WAH_WAH dtl::diff<dtl::dynamic_wah32, dtl::dynamic_wah32>

#define __DIFF_NAME_WAH_ROARING wah_roaring
#define __DIFF_TYPE_WAH_ROARING dtl::diff<dtl::dynamic_wah32, dtl::dynamic_roaring_bitmap>

#define __DIFF_NAME_PART_WAH_ROARING part_wah_roaring
#define __DIFF_TYPE_PART_WAH_ROARING dtl::diff<partitioned_wah, dtl::dynamic_roaring_bitmap>
//===----------------------------------------------------------------------===//
// FIXME misnormer ... updatable bitmap types
enum class diff_bitmap_t {
  __DIFF_NAME_ROARING_WAH,
  __DIFF_NAME_ROARING_ROARING,
  __DIFF_NAME_TEB_WAH,
  __DIFF_NAME_TEB_ROARING,
  __DIFF_NAME_WAH_WAH,
  __DIFF_NAME_WAH_ROARING,
  __DIFF_NAME_PART_WAH_ROARING,
  roaring,
  part_teb, // a partitioned TEB which supports direct updates
  part_diff_teb, // a partitioned TEB where each partition has a diff structure
  part_wah, // a partitioned WAH which supports direct updates
  part_diff_wah, // a partitioned WAH where each partition has a diff structure
};
//===----------------------------------------------------------------------===//
template<diff_bitmap_t B>
struct diff_type_of {
  using type = void;
};

template<>
struct diff_type_of<diff_bitmap_t::roaring> {
  using type = dtl::dynamic_roaring_bitmap;
};

template<>
struct diff_type_of<diff_bitmap_t::part_teb> {
  using type = partitioned_teb;
};

template<>
struct diff_type_of<diff_bitmap_t::part_diff_teb> {
  using type = partitioned_differential_teb;
};

template<>
struct diff_type_of<diff_bitmap_t::part_wah> {
  using type = partitioned_wah;
};

template<>
struct diff_type_of<diff_bitmap_t::part_diff_wah> {
  using type = partitioned_differential_wah;
};

#define __GENERATE(NAME, TYPE)               \
  template<>                                 \
  struct diff_type_of<diff_bitmap_t::NAME> { \
    using type = TYPE;                       \
  };

__GENERATE(__DIFF_NAME_ROARING_WAH, __DIFF_TYPE_ROARING_WAH)
__GENERATE(__DIFF_NAME_ROARING_ROARING, __DIFF_TYPE_ROARING_ROARING)
__GENERATE(__DIFF_NAME_TEB_WAH, __DIFF_TYPE_TEB_WAH)
__GENERATE(__DIFF_NAME_TEB_ROARING, __DIFF_TYPE_TEB_ROARING)
__GENERATE(__DIFF_NAME_WAH_WAH, __DIFF_TYPE_WAH_WAH)
__GENERATE(__DIFF_NAME_WAH_ROARING, __DIFF_TYPE_WAH_ROARING)
__GENERATE(__DIFF_NAME_PART_WAH_ROARING, __DIFF_TYPE_PART_WAH_ROARING)
#undef __GENERATE
//===----------------------------------------------------------------------===//
enum class diff_merge_t {
  naive,
  naive_iter,
  inplace,
  tree,
  no,
  forward,
};
//===----------------------------------------------------------------------===//
template<diff_bitmap_t B, diff_merge_t M>
struct diff_merge_type_of {
  using type = void;
};

#define __GENERATE(BNAME, MNAME, MTYPE)                                  \
  template<>                                                             \
  struct diff_merge_type_of<diff_bitmap_t::BNAME, diff_merge_t::MNAME> { \
    using type = MTYPE<                                                  \
        diff_type_of<diff_bitmap_t::BNAME>::type::bitmap_type,           \
        diff_type_of<diff_bitmap_t::BNAME>::type::diff_type>;            \
  };

__GENERATE(__DIFF_NAME_ROARING_ROARING, naive, dtl::merge_naive)
__GENERATE(__DIFF_NAME_ROARING_ROARING, naive_iter, dtl::merge_naive_iter)
__GENERATE(__DIFF_NAME_ROARING_ROARING, inplace, dtl::merge_inplace)

__GENERATE(__DIFF_NAME_ROARING_WAH, naive, dtl::merge_naive)
__GENERATE(__DIFF_NAME_ROARING_WAH, naive_iter, dtl::merge_naive_iter)

__GENERATE(__DIFF_NAME_TEB_ROARING, naive, dtl::merge_naive)
__GENERATE(__DIFF_NAME_TEB_ROARING, naive_iter, dtl::merge_naive_iter)
__GENERATE(__DIFF_NAME_TEB_ROARING, tree, dtl::merge_tree)

__GENERATE(__DIFF_NAME_TEB_WAH, naive, dtl::merge_naive)
__GENERATE(__DIFF_NAME_TEB_WAH, naive_iter, dtl::merge_naive_iter)
__GENERATE(__DIFF_NAME_TEB_WAH, tree, dtl::merge_tree)

__GENERATE(__DIFF_NAME_WAH_ROARING, naive, dtl::merge_naive)
__GENERATE(__DIFF_NAME_WAH_ROARING, naive_iter, dtl::merge_naive_iter)

__GENERATE(__DIFF_NAME_PART_WAH_ROARING, naive, dtl::merge_naive)
__GENERATE(__DIFF_NAME_PART_WAH_ROARING, naive_iter, dtl::merge_naive_iter)

__GENERATE(__DIFF_NAME_WAH_WAH, naive, dtl::merge_naive)
__GENERATE(__DIFF_NAME_WAH_WAH, naive_iter, dtl::merge_naive_iter)
__GENERATE(__DIFF_NAME_WAH_WAH, inplace, dtl::merge_inplace)

#undef __GENERATE

// Partitioned TEB.
template<>
struct diff_merge_type_of<diff_bitmap_t::part_teb, diff_merge_t::no> {
  using type = dtl::merge_not<partitioned_teb, void>; // no merge
};

// Partitioned WAH.
template<>
struct diff_merge_type_of<diff_bitmap_t::part_wah, diff_merge_t::no> {
  using type = dtl::merge_not<partitioned_wah, void>; // no merge
};

// Partitioned/differential TEB.
template<>
struct diff_merge_type_of<diff_bitmap_t::part_diff_teb, diff_merge_t::naive> {
  // Merge is handled by the bitmap itself.
  using type = dtl::merge_naive<dtl::teb_wrapper,partitioned_teb_diff_structure>;
};
template<>
struct diff_merge_type_of<diff_bitmap_t::part_diff_teb, diff_merge_t::naive_iter> {
  // Merge is handled by the bitmap itself.
  using type = dtl::merge_naive_iter<dtl::teb_wrapper,partitioned_teb_diff_structure>;
};
template<>
struct diff_merge_type_of<diff_bitmap_t::part_diff_teb, diff_merge_t::tree> {
  // Merge is handled by the bitmap itself.
  using type = dtl::merge_tree<dtl::teb_wrapper,partitioned_teb_diff_structure>;
};

// Partitioned/differential WAH.
template<>
struct diff_merge_type_of<diff_bitmap_t::part_diff_wah, diff_merge_t::naive> {
  // Merge is handled by the bitmap itself.
  using type = dtl::merge_naive<dtl::dynamic_wah32,partitioned_wah_diff_structure>;
};
template<>
struct diff_merge_type_of<diff_bitmap_t::part_diff_wah, diff_merge_t::naive_iter> {
  // Merge is handled by the bitmap itself.
  using type = dtl::merge_naive_iter<dtl::dynamic_wah32,partitioned_wah_diff_structure>;
};
//===----------------------------------------------------------------------===//
#undef __DIFF_NAME_ROARING_WAH
#undef __DIFF_NAME_ROARING_ROARING
#undef __DIFF_NAME_TEB_WAH
#undef __DIFF_NAME_TEB_ROARING
#undef __DIFF_NAME_WAH_WAH
#undef __DIFF_NAME_WAH_ROARING

#undef __DIFF_TYPE_ROARING_WAH
#undef __DIFF_TYPE_ROARING_ROARING
#undef __DIFF_TYPE_TEB_WAH
#undef __DIFF_TYPE_TEB_ROARING
#undef __DIFF_TYPE_WAH_WAH
#undef __DIFF_TYPE_WAH_ROARING
//===----------------------------------------------------------------------===//
