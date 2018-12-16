#pragma once

#include <boost/dynamic_bitset.hpp>

#include <dtl/dtl.hpp>

#include "two_state_markov_process.hpp"


//===----------------------------------------------------------------------===//
namespace dtl {
using bitmap = boost::dynamic_bitset<$u32>;
} // namespace dtl
//===----------------------------------------------------------------------===//


//===----------------------------------------------------------------------===//
dtl::bitmap
gen_bitmap(u64 n, $f64 f, $f64 d);
//===----------------------------------------------------------------------===//

