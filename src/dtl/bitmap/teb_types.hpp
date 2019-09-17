#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/dtl.hpp>
#include <dtl/bitmap/util/rank1_surf.hpp>
//===----------------------------------------------------------------------===//
#if !defined(__teb_inline__)
#if defined(NDEBUG)
// Release build.
#define __teb_inline__ inline __attribute__((always_inline))
//#define __teb_inline__ __attribute__((noinline))
#else
#define __teb_inline__
#endif
#endif
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// The size type.
using teb_size_type = $u32;
/// The storage type. The size of a TEB is a multiple of sizeof(teb_word_type).
using teb_word_type = $u64;
/// Support data structure for rank1 operations on the tree structure.
using teb_rank_type = dtl::rank1_surf<teb_word_type, true>;
//===----------------------------------------------------------------------===//
#pragma pack(push, 1)
/// The header of a TEB.
struct teb_header {
  /// The number of bits in the bitmap.
  teb_size_type n = 0;
  /// The length of the encoded tree.
  teb_size_type tree_bit_cnt = 0;
  /// The number of implicit inner nodes in the tree structure.
  teb_size_type implicit_inner_node_cnt = 0;
  /// The number of label bits.
  teb_size_type label_bit_cnt = 0;
  /// The number of implicit leading 0-labels.
  teb_size_type implicit_leading_label_cnt = 0;
  /// The number of perfect levels.
  $u8 perfect_level_cnt = 0; // FIXME redundant, as is can be computed from the number of implicit inner nodes
  /// The height of the encoded (pruned) tree.
  $u8 encoded_tree_height = 0;
  /// True if the TEB contains level offsets at the very end.
  $u1 has_level_offsets = false;
  /// Padding.
  $u8 padding = 0;
};
#pragma pack(pop)
//===----------------------------------------------------------------------===//
static_assert(sizeof(teb_header) == 24,
    "A TEB header is supposed to be 24 bytes in size.");
static_assert(sizeof(teb_header) % sizeof(teb_word_type) == 0,
    "A TEB header is supposed to be a multiple of the word size.");
//===----------------------------------------------------------------------===//
/// Used to represent empty trees and labels.
static const teb_word_type teb_null_word = 0;
//===----------------------------------------------------------------------===//
}; // namespace dtl
