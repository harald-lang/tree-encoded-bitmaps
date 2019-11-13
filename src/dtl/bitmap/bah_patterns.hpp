#pragma once
//===----------------------------------------------------------------------===//
#include "bah_types.hpp"

#include <dtl/dtl.hpp>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
// Pre-defined lookup tables for the one- and two-byte encodable patterns, as
// described in 'BAH: A Bitmap Index Compression Algorithm for Fast Data
// Retrieval' by Li et al.
//===----------------------------------------------------------------------===//
/// Special value that signals, that no code for the given bit pattern exists.
static constexpr bah_word_t bah_code_not_found = bah_word_t(~0ull);
//===----------------------------------------------------------------------===//
/// Determines the code for the given word. If no code exists,
/// 'bah_code_not_found' is returned.
/// Time complexity: Logarithmic in the LuT size.
bah_word_t
get_oep_code(bah_word_t word); // TODO should be in O(1)
//===----------------------------------------------------------------------===//
/// Returns the word (or pattern) associated with the given code.
bah_word_t
get_oep(bah_word_t code);
//===----------------------------------------------------------------------===//
/// Determines the code for the given word. If no code exists,
/// 'bah_code_not_found' is returned.
/// Time complexity: Logarithmic in the LuT size.
bah_word_t
get_tep_code(bah_word_t word); // TODO should be in O(1)
//===----------------------------------------------------------------------===//
/// Returns the word (or pattern) associated with the given code.
bah_word_t
get_tep(bah_word_t code);
//===----------------------------------------------------------------------===//
} // namespace dtl
