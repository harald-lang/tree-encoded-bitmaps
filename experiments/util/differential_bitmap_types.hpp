#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/bitmap/diff/diff.hpp>
#include <dtl/bitmap/dynamic_roaring_bitmap.hpp>
#include <dtl/bitmap/dynamic_wah.hpp>
#include <dtl/bitmap/teb_wrapper.hpp>
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
//===----------------------------------------------------------------------===//

enum class diff_bitmap_t {
  __DIFF_NAME_ROARING_WAH,
  __DIFF_NAME_ROARING_ROARING,
  __DIFF_NAME_TEB_WAH,
  __DIFF_NAME_TEB_ROARING,
  __DIFF_NAME_WAH_WAH,
  __DIFF_NAME_WAH_ROARING,
};
//===----------------------------------------------------------------------===//
template<diff_bitmap_t B>
struct diff_type_of {
  using type = void;
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

#undef __GENERATE

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
