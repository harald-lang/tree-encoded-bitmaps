#pragma once
//===----------------------------------------------------------------------===//
#include "teb_types.hpp"

#include <dtl/dtl.hpp>
#include <dtl/bits.hpp>
#include <dtl/math.hpp>
//===----------------------------------------------------------------------===//
namespace dtl {
namespace teb_util {
//===--------------------------------------------------------------------===//
// Helper functions.
//===--------------------------------------------------------------------===//
/// Computes the number of perfect levels in the tree structure based on the
/// number of implicit inner nodes. The number of perfect levels is at least
/// one, because there is at least one node in the tree structure.
static inline u64
determine_perfect_tree_levels(u64 implicit_inner_node_cnt) noexcept {
  return dtl::log_2(implicit_inner_node_cnt + 1) + 1;
}
//===--------------------------------------------------------------------===//
/// Computes the height of the tree based on n. - Note that a tree, consisting
/// of a single (root) node has a height of 0.
static inline u64
determine_tree_height(u64 n) noexcept {
  return dtl::log_2(n);
}
//===----------------------------------------------------------------------===//
/// Determines the path to the common ancestor of the two nodes specified
/// by there paths.
static inline u64
determine_common_ancestor_path(u64 src_path, u64 dst_path) noexcept {
  //TODO should use positions instead of paths (at least for the second argument)
  assert(src_path != dst_path);
  const auto t0 = dtl::bits::lz_count(src_path) + 1;
  const auto a = src_path << t0;
  const auto b = dst_path << (dtl::bits::lz_count(dst_path) + 1);
  assert(a < b);
  const auto src_path_len = sizeof(src_path) * 8 - t0;
  const auto common_prefix_length = dtl::bits::lz_count(a ^ b);
  const auto common_ancestor_path =
      src_path >> (src_path_len - common_prefix_length);
  return common_ancestor_path;
}
//===----------------------------------------------------------------------===//
static inline void
determine_common_ancestor_path(u64 src_path, u64 dst_path,
    $u64& out_common_ancestor_path, $u64& out_common_ancestor_level) noexcept {
  //TODO should use positions instead of paths (at least for the second argument)
  assert(src_path != dst_path);
  const auto t0 = dtl::bits::lz_count(src_path) + 1;
  const auto a = src_path << t0;
  const auto b = dst_path << (dtl::bits::lz_count(dst_path) + 1);
  assert(a < b);
  const auto src_path_len = sizeof(src_path) * 8 - t0;
  const auto common_prefix_length = dtl::bits::lz_count(a ^ b);
  out_common_ancestor_path =
      src_path >> (src_path_len - common_prefix_length);
  out_common_ancestor_level = common_prefix_length;
}
//===----------------------------------------------------------------------===//
static inline void
determine_common_ancestor_path2(u64 src_path, u64 dst_pos, u64 tree_height,
    $u64& out_common_ancestor_path, $u64& out_common_ancestor_level) noexcept {
  const auto t0 = dtl::bits::lz_count(src_path) + 1;
  const auto a = src_path << t0;
  const auto b = dst_pos << (sizeof(dst_pos) * 8 - tree_height);
  assert(a < b);
  const auto src_path_len = sizeof(src_path) * 8 - t0;
  const auto common_prefix_length = dtl::bits::lz_count(a ^ b);
  out_common_ancestor_path =
      src_path >> (src_path_len - common_prefix_length);
  out_common_ancestor_level = common_prefix_length;
}
//===----------------------------------------------------------------------===//
static inline u64
determine_level_of(u64 path) noexcept {
  const auto lz_cnt_path = dtl::bits::lz_count(path);
  const auto level = sizeof(u64) * 8 - 1 - lz_cnt_path;
  return level;
}
//===----------------------------------------------------------------------===//
} // namespace teb_util
} // namespace dtl
//===----------------------------------------------------------------------===//
