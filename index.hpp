#pragma once

#include "bitmap/dynamic_roaring_bitmap.hpp"
#include "bitmap/static/tree_mask_po.hpp"
#include "bitmap/static/wah.hpp"


// API:
// 1) construct from bitset
// 2) size_in_bytes() in bytes
// 3) binary op: a XOR b (range encoding, first predicate)
// 3.1) binary op: a OR b (without range encoding, first predicate)
// 4) ternary op: a & (b XOR c) (range encoding, second predicate)
// 5) to position list (maybe batchwise)
// 6) to bitset (validation code)
// 7) equality (RE, reference to an equal mask instead of storing it redundantly)

// constexpr u1 is_lossy, is_precise, ....



// set() set all bits
// reset() clear all bits
// == (compare neighboring entries)