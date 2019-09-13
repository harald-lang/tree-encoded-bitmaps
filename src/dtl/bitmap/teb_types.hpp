#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/dtl.hpp>
#include <dtl/bitmap/util/rank1_surf.hpp>
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
  teb_size_type n;
  /// The length of the encoded tree.
  teb_size_type tree_bit_cnt;
  /// The number of implicit inner nodes in the tree structure.
  teb_size_type implicit_inner_node_cnt;
  /// The number of label bits.
  teb_size_type label_bit_cnt;
  /// The number of implicit leading 0-labels.
  teb_size_type implicit_leading_label_cnt;
  /// The number of perfect levels.
  $u8 perfect_level_cnt;
  /// The height of the encoded (pruned) tree.
  $u8 encoded_tree_height;
  /// True if the TEB contains level offsets at the very end.
  $u1 has_level_offsets;
  /// Padding.
  $u8 padding = 0;
};
#pragma pack(pop)
static_assert(sizeof(teb_header) == 24,
    "A TEB header is supposed to be 24 bytes in size.");
static_assert(sizeof(teb_header) % sizeof(teb_word_type) == 0,
    "A TEB header is supposed to be a multiple of the word size.");
//===----------------------------------------------------------------------===//
}; // namespace dtl
