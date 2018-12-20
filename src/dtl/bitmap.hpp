#pragma once

#include <boost/dynamic_bitset.hpp>
#include <dtl/dtl.hpp>

#include "bitmap/dynamic_bitmap.hpp"
#include "bitmap/dynamic_partitioned_tree_mask.hpp"
#include "bitmap/dynamic_roaring_bitmap.hpp"
#include "bitmap/dynamic_tree_mask_lo.hpp"
#include "bitmap/dynamic_wah.hpp"
#include "bitmap/static/bitmap.hpp"
#include "bitmap/static/partitioned_tree_mask.hpp"
#include "bitmap/static/roaring_bitmap.hpp"
#include "bitmap/static/tree_mask_lo.hpp"
#include "bitmap/static/tree_mask_po.hpp"
#include "bitmap/static/wah.hpp"

//===----------------------------------------------------------------------===//
namespace dtl {
using bitmap = boost::dynamic_bitset<$u32>;
} // namespace dtl
//===----------------------------------------------------------------------===//
