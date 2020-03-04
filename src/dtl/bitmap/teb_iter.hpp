#pragma once
//===----------------------------------------------------------------------===//
#include "teb_flat.hpp"
#include "teb_util.hpp"

#include <dtl/dtl.hpp>
#include <dtl/static_stack.hpp>

#include <cassert>
#include <ostream>
#include <string>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// A 1-run iterator for TEBs with efficient skip support.
class teb_iter {
  static constexpr u32 optimization_level_ = 3; // TODO remove

  /// The fundamental type to encode paths within the tree.
  using path_t = $u64;

  struct stack_entry {
    $u64 node_idx;
    path_t path;
    $u64 rank;
    $u64 level;
    $u64 is_inner;
    void
    print(std::ostream& os) const noexcept {
      os << "node[idx=" << node_idx
         << ",rank=" << rank
         << ",level=" << level
         << ",path=" << std::bitset<32>(path)
         << "]";
    }
  };

  /// Reference to the TEB instance.
  const teb_flat& teb_;
  /// The height of the tree structure. // TODO refers to max height
  u64 tree_height_;
  /// The number of perfect levels in the tree structure.
  u64 perfect_levels_;
  // TODO
  u64 partition_shift_;
  /// First node index in the last perfect level.
  u64 top_node_idx_begin_;
  /// Last node index + 1 in the last perfect level.
  u64 top_node_idx_end_;
  /// The current top node.
  $u64 top_node_idx_current_;
  /// The stack contains the tree nodes that need to be visited.
  static_stack<stack_entry, 32> stack_;
  /// Points to the beginning of the current 1-fill.
  $u64 pos_;
  /// The length of the current 1-fill.
  $u64 length_;
  /// The current node.
  $u64 node_idx_;
  /// The path to the current node
  /// (required only for leaf-to-leaf navigation).
  path_t path_;

public:
  /// Constructs an iterator for the given TEB instance. After construction,
  /// the iterator points to the first 1-fill.
  explicit __teb_inline__
  teb_iter(const teb_flat& teb) noexcept
      : teb_(teb),
        tree_height_(dtl::teb_util::determine_tree_height(teb.n_)),
        perfect_levels_(teb.perfect_level_cnt_),
        partition_shift_(tree_height_ - (teb.perfect_level_cnt_ - 1)),
        top_node_idx_begin_((1ull << (teb.perfect_level_cnt_ - 1)) - 1),
        top_node_idx_end_((1ull << teb.perfect_level_cnt_) - 1),
        top_node_idx_current_((1ull << (teb.perfect_level_cnt_ - 1)) - 1),
        pos_(0),
        length_(0),
        node_idx_((1ull << (teb.perfect_level_cnt_ - 1)) - 1),
        path_(path_t(1) << (teb.perfect_level_cnt_ - 1)) {
    // Initialize the stack.
    --top_node_idx_current_;
    next_top_node();
    next();
  }

  teb_iter(teb_iter&&) noexcept = default;

  /// Forwards the iterator to the next top node that is either an inner node
  /// or a leaf node with a 1-label.
  void __teb_inline__
  next_top_node() {
    stack_entry& next_node = stack_.push();
    while (++top_node_idx_current_ < top_node_idx_end_) {
      const auto is_inner = 0ull + teb_.is_inner_node(top_node_idx_current_);
      const auto rank = teb_.rank_inclusive(top_node_idx_current_);
      if (is_inner) {
        next_node.node_idx = top_node_idx_current_;
        next_node.path = path_t(top_node_idx_current_ - top_node_idx_begin_);
        // Set the sentinel bit.
        next_node.path |= path_t(1) << (perfect_levels_ - 1);
        next_node.level = perfect_levels_ - 1;
        next_node.rank = rank;
        next_node.is_inner = is_inner;
        return;
      }
      else {
        // Check label.
        // We don't want to have leaf nodes with a 0-label on the stack.
        const auto label_idx = top_node_idx_current_ - rank;
        const auto label = teb_.get_label_by_idx(label_idx);
        if (label) {
          next_node.node_idx = top_node_idx_current_;
          next_node.path = path_t(top_node_idx_current_ - top_node_idx_begin_);
          // Set the sentinel bit.
          next_node.path |= path_t(1) << (perfect_levels_ - 1);
          next_node.level = perfect_levels_ - 1;
          next_node.rank = rank;
          next_node.is_inner = is_inner;
          return;
        }
      }
    }
  }

  /// Forwards the iterator to the next 1-fill (if any).
  /// Use the functions pos() and length() to get the 1-fill the iterator
  /// is currently pointing to.
  void __teb_inline__
  next_off() noexcept __attribute__((flatten, hot)) {
    while (top_node_idx_current_ < top_node_idx_end_) {
    outer_loop_begin:
      while (!stack_.empty()) {
        // The current node.
        stack_entry node_info = stack_.top();
        stack_.pop();

        $u1 label;
        while (teb_.is_inner_node(node_info.node_idx)) {
        loop_begin:
          label = false;
          // Determine left and right child. - Both exist, because the tree
          // is full binary.
          // Note: If the current node is an inner node, the rank is always
          //       required.
          u64 right_child_idx = 2 * node_info.rank;
          u64 left_child_idx = right_child_idx - 1;

          // Determine whether the children are inner or leaf nodes.
          const auto right_child_is_inner = 0
              + teb_.is_inner_node(right_child_idx);
          const auto left_child_is_inner = 0
              + teb_.is_inner_node(left_child_idx);
          const auto children_are_inner = (left_child_is_inner << 1) | right_child_is_inner; // TODO fetch both bits in one go

          // Compute the rank for one child, and derive the rank of the
          // other one.
          // Rank is required to compute the label index.
          u64 left_child_rank = teb_.rank_inclusive(left_child_idx);
          u64 left_child_label_idx = left_child_idx - left_child_rank
              + left_child_is_inner; // prevent underflow
          u64 right_child_label_idx = left_child_label_idx + 1
              - left_child_is_inner; // adjust index if necessary

          //          // TODO Eagerly fetch the labels.
          //          u1 left_child_label = teb_.get_label_by_idx(left_child_label_idx);
          //          u1 right_child_label = teb_.get_label_by_idx(right_child_label_idx);

          // Switch over the different cases.
          switch (children_are_inner) {
            case 0b00: {
              //===------------------------------------------------------===//
              // Both childs are leaf nodes.
              //===------------------------------------------------------===//

              // Note for UN-optimized TEBs the following holds:
              //   One child has a 1-label the other has
              //   a 0-label, which is guaranteed by the bottom-up pruning.
              //   Thus, right_label == !left_label.
              //
              // However, gradual decompression (optimization level >= 2) may
              // expand nodes and replicate labels. Thus, the sub-tree is no
              // longer guaranteed to be 'compressed', and therefore, both
              // labels need to be inspected.

              // Fetch the labels.
              u1 left_child_label = teb_.get_label_by_idx(left_child_label_idx);
              u1 right_child_label = teb_.get_label_by_idx(right_child_label_idx);

              if (optimization_level_ < 2) {
                // Go to the node which has the 1-label.
                node_info.node_idx = left_child_label ? left_child_idx : right_child_idx;
                node_info.path = (node_info.path << 1) | !left_child_label;
                node_info.level++;
                node_info.rank = left_child_rank;
                goto produce_output;
              }
              else {
                // The price we pay for better compression ratios.
                u64 both_labels = 2 * left_child_label + right_child_label;
                switch (both_labels) {
                  case 0b00: {
                    goto outer_loop_begin;
                  }
                  case 0b01: {
                    node_info.node_idx = right_child_idx;
                    node_info.path = (node_info.path << 1) | 1;
                    node_info.level++;
                    node_info.rank = left_child_rank + 1;
                    goto produce_output;
                  }
                  case 0b10: {
                    node_info.node_idx = left_child_idx;
                    node_info.path <<= 1;
                    node_info.level++;
                    node_info.rank = left_child_rank;
                    goto produce_output;
                  }
                  case 0b11: {
                    goto produce_output;
                  }
                }
              }
            }
            case 0b01: {
              //===------------------------------------------------------===//
              // Left child is a leaf, right child is an inner node.
              //===------------------------------------------------------===//
              // Determine whether the left child produces an output
              // (label = 1).
              u1 left_child_label = teb_.get_label_by_idx(left_child_label_idx);

              // Derive rank of the right child.  The following works, because
              // the left child is leaf (0-bit) and the right is inner (1-bit).
              u64 right_child_rank = left_child_rank + 1;

              if (left_child_label) {
                // Produce the output for the left child iff it has a 1-label,
                // otherwise it can be ignored.
                node_info.node_idx = left_child_idx;
                //                label = true; // TODO remove
                node_info.path <<= 1;
                node_info.rank = left_child_rank;
                node_info.level++;
                // Push the right child on the stack.
                stack_entry& right_child_info = stack_.push();
                right_child_info.node_idx = right_child_idx;
                right_child_info.path = node_info.path | 1;
                right_child_info.level = node_info.level;
                right_child_info.rank = right_child_rank;
                goto produce_output;
              }
              // Else, go to right child, ignoring the left child.
              node_info.node_idx = right_child_idx;
              node_info.path = (node_info.path << 1) | 1;
              node_info.level++;
              node_info.rank = right_child_rank;
              goto loop_begin;
            }
            case 0b10: {
              //===------------------------------------------------------===//
              // Left child is an inner node, right child is a leaf node.
              //===------------------------------------------------------===//

              // Determine whether the right child produces an output
              // (label = 1).
              u1 right_child_label = teb_.get_label_by_idx(right_child_label_idx);
              if (right_child_label) { // FIXME DEP
                // Push the right child on the stack iff it has a 1-label,
                // otherwise the right child is ignored.
                stack_entry& right_child_info = stack_.push();
                right_child_info.node_idx = right_child_idx;
                right_child_info.path = (node_info.path << 1) | 1;
                right_child_info.level = node_info.level + 1;
                // Rank of the right child is equal to the rank of the left
                // child.
                u64 right_child_rank = left_child_rank;
                right_child_info.rank = right_child_rank;
              }
              // Go to left child.
              node_info.node_idx = left_child_idx;
              node_info.path <<= 1;
              node_info.level++;
              node_info.rank = left_child_rank;
              goto loop_begin;
            }
            case 0b11: {
              //===------------------------------------------------------===//
              // Both children are an inner nodes.
              //===------------------------------------------------------===//

              u64 right_child_rank = left_child_rank + 1;
              // Go to left child.
              node_info.node_idx = left_child_idx;
              node_info.path <<= 1;
              node_info.level++;
              node_info.rank = left_child_rank;
              // Push the right child on the stack.
              stack_entry& right_child_info = stack_.push();
              right_child_info.node_idx = right_child_idx;
              right_child_info.path = node_info.path | 1;
              right_child_info.level = node_info.level;
              right_child_info.rank = right_child_rank;
              goto loop_begin;
            }
            default:
              __builtin_unreachable();
          }
        }
        // Reached a leaf node.
        label = teb_.get_label_by_idx(node_info.node_idx - node_info.rank);
        if (label) {
        produce_output:
          // Produce output (a 1-fill).
          //          const auto lz_cnt_path = dtl::bits::lz_count(node_info.path);
          //          const auto level = sizeof(path_t) * 8 - 1 - lz_cnt_path;
          const auto level = node_info.level;
          // Toggle sentinel bit (= highest bit set) and add offset.
          pos_ = (node_info.path ^ (1ull << level)) << (tree_height_ - level);
          // The length of the 1-fill.
          length_ = teb_.n_ >> level;
          node_idx_ = node_info.node_idx;
          path_ = node_info.path;
          return;
        }
      }

      // TODO bypass the stack here
      ++top_node_idx_current_;
      stack_entry& next_node = stack_.push();
      next_node.node_idx = top_node_idx_current_;
      next_node.path = path_t(top_node_idx_current_ - top_node_idx_begin_);
      // Set the sentinel bit.
      next_node.path |= path_t(1) << (perfect_levels_ - 1);
      next_node.level = perfect_levels_ - 1;
      next_node.rank = teb_.rank_inclusive(top_node_idx_current_);
    }
    pos_ = teb_.n_actual_;
    length_ = 0;
  }

  void __teb_inline__
  next() noexcept __attribute__((flatten, hot)) {
    assert(!end());
    while (top_node_idx_current_ < top_node_idx_end_) {
      while (!stack_.empty()) {
        // The current node.
        stack_entry node_info = stack_.top();
        stack_.pop();
        $u1 node_label = node_info.is_inner == 0;

        while (node_info.is_inner) {
          // Determine left and right child. - Both exist, because the tree
          // is full binary.
          u64 right_child_idx = 2 * node_info.rank;
          u64 left_child_idx = right_child_idx - 1;

          // Determine whether the children are inner or leaf nodes.
          const auto right_child_is_inner = 0ull
              + teb_.is_inner_node(right_child_idx);
          const auto left_child_is_inner = 0ull
              + teb_.is_inner_node(left_child_idx);

          // Compute the rank for one child, and derive the rank of the
          // other one.
          u64 left_child_rank = teb_.rank_inclusive(left_child_idx);
          u64 right_child_rank = left_child_rank + right_child_is_inner;

          // Determine the label indices in L.
          u64 left_child_label_idx = left_child_idx - left_child_rank
              + left_child_is_inner; // prevent underflow
          u64 right_child_label_idx = left_child_label_idx + 1
              - left_child_is_inner; // adjust index if necessary

          // Eagerly fetch the labels.
          u1 left_child_label = teb_.get_label_by_idx(left_child_label_idx);
          node_label = left_child_label;
          u1 right_child_label = teb_.get_label_by_idx(right_child_label_idx);

          // Push right child on the stack and go to left child.
          node_info.node_idx = left_child_idx;
          node_info.path <<= 1;
          node_info.level++;
          node_info.rank = left_child_rank;
          node_info.is_inner = left_child_is_inner;
          // Push the right child only if necessary.
          u1 push_right_child = (right_child_is_inner | right_child_label);
          stack_entry& right_child_info = stack_.push();
          right_child_info.node_idx = right_child_idx;
          right_child_info.path = node_info.path | 1;
          right_child_info.level = node_info.level;
          right_child_info.rank = right_child_rank;
          right_child_info.is_inner = right_child_is_inner;
          stack_.cnt_ = push_right_child ? stack_.cnt_ : stack_.cnt_ - 1;
        }
        // Reached a leaf node.
        if (node_label) {
          // Produce output (a 1-fill).
          const auto level = node_info.level;
          // Toggle sentinel bit (= highest bit set) and add offset.
          pos_ = (node_info.path ^ (1ull << level)) << (tree_height_ - level); // TODO there might be no sentinel bit when the stack is used
          // The length of the 1-fill.
          length_ = teb_.n_ >> level;
          node_idx_ = node_info.node_idx;
          path_ = node_info.path;
          return;
        }
      }

      // Push the next top node on the stack (if any).
      next_top_node();
    }
    pos_ = teb_.n_actual_;
    length_ = 0;
  }

  /// Navigate to the desired position, starting from the trees' root node.
  void __teb_inline__
  nav_from_root_to(const std::size_t to_pos) noexcept {
    //===------------------------------------------------------------------===//
    // (Re-)initialize the iterator state.
    stack_.clear();

    // Determine the top-node path.
    $u64 level = perfect_levels_ - 1;
    const auto foo = to_pos >> (tree_height_ - level);
    path_ = (path_t(1) << level) | foo;

    // Determine the top-node idx.
    top_node_idx_current_ = top_node_idx_begin_ + foo;
    //===------------------------------------------------------------------===//

    node_idx_ = top_node_idx_current_;
    nav_downwards(to_pos);
    return;
  }

  /// Navigate downwards the tree to the desired position, starting from the
  /// current node. The destination position must be part of the current sub-
  /// tree, otherwise the behavior is undefined.
  void __teb_inline__
  nav_downwards(const std::size_t to_pos) noexcept {
    $u64 level = dtl::teb_util::determine_level_of(path_);
    auto rank = teb_.rank_inclusive(node_idx_);
    std::size_t i = tree_height_ - level - 1;
    while (true) {
      // First check, if this is already a leaf node.
      if (teb_.is_leaf_node(node_idx_)) {
        // Reached the desired position.
        const auto label_idx = node_idx_ - rank;
        const auto label = teb_.get_label_by_idx(label_idx);
        if (label) {
          // Found the corresponding leaf node.
          // Toggle sentinel bit (= highest bit set) and add offset.
          pos_ = (path_ ^ (1ull << level)) << (tree_height_ - level);
          // The length of the 1-fill.
          length_ = teb_.n_ >> level;
          // Adjust the current position and fill-length.
          length_ -= to_pos - pos_;
          pos_ = to_pos;
          return;
        }
        else {
          // Search forward to the next 1-fill.
          next();
          return;
        }
      }

      // Navigate downwards the tree.
      // 0 -> go to left child, 1 -> go to right child
      u1 direction_bit = dtl::bits::bit_test(to_pos, i);
      i--;
      const auto right_child_idx = 2 * rank;
      const auto left_child_idx = right_child_idx - 1;
      const auto right_child_rank = teb_.rank_inclusive(right_child_idx);
      const auto right_child_is_inner = teb_.is_inner_node(right_child_idx);
      const auto left_child_rank = right_child_rank - right_child_is_inner;
      level++;
      if (!direction_bit) {
        // Push the right child only if necessary.
        if (right_child_is_inner
            || teb_.get_label_by_idx(right_child_idx - right_child_rank)) {
          stack_entry& right_child_info = stack_.push();
          right_child_info.node_idx = right_child_idx;
          right_child_info.path = (path_ << 1) | 1;
          right_child_info.level = level;
          right_child_info.rank = right_child_rank;
          right_child_info.is_inner = 0ull + right_child_is_inner;
        }
        // Go to left child.
        path_ <<= 1;
        node_idx_ = left_child_idx;
        rank = left_child_rank;
      }
      else {
        // Go to right child.
        path_ = (path_ << 1) | 1;
        node_idx_ = right_child_idx;
        rank = right_child_rank;
      }
    }
  }

  /// Forwards the iterator to the given position. The function first
  /// determines the common ancestor node and then decides whether it is
  /// cheaper to navigate starting from the current node or from the root
  /// node.
  void __teb_inline__
  nav_to(const std::size_t to_pos) noexcept {
    assert(to_pos >= pos_ + length_);
    assert(perfect_levels_ > 0);
    // Fast path.  If the skip distance is larger than the range spanned by
    // the current subtree, we immediately start navigating downwards from the
    // root node.  Thus, we do not need to compute the common ancestor node.
    if (pos_ >> partition_shift_ != to_pos >> partition_shift_
        || stack_.empty()) {
      nav_from_root_to(to_pos);
      return;
    }

    // Determine the common ancestor node.  Note that the common ancestor is
    // guaranteed to be in the lower (non-perfect) tree part. Otherwise, we
    // would have taken the fast path above.
    const path_t from_path = path_;
    path_t common_ancestor_path;
    $u64 common_ancestor_level;
    dtl::teb_util::determine_common_ancestor_path2(
        from_path, to_pos, tree_height_,
        common_ancestor_path, common_ancestor_level);

    // Decide whether to start navigating from the current or from the root node.

    // Determine the number of upward steps from the current node to the common
    // ancestor node. - We use a fast approximation. The actual step count is
    // at most 2 steps higher.
    const auto upstep_cnt = stack_.top().level - common_ancestor_level;

    // Determine the number of downward steps when starting from the root node.
    const auto downstep_cnt = common_ancestor_level - (perfect_levels_ - 1);
    assert(downstep_cnt < tree_height_);

    if ((upstep_cnt > tree_height_) // underflow happens if the common ancestor is not on the stack
        || (upstep_cnt << 3) > downstep_cnt) { // cost(downstep) is approx. 9 x cost(upstep), however, we use << 3 instead of * 9
      nav_from_root_to(to_pos);
      return;
    }

    assert(!stack_.empty());

    // Determine the right child of the common ancestor node.
    const auto right_child_of_common_ancestor_path = (common_ancestor_path << 1) | 1ull;
    const auto right_child_of_common_ancestor_level = common_ancestor_level + 1;

    $u64 node_idx = node_idx_;
    auto path = from_path;
    $u64 level = tree_height_;
    while (path != right_child_of_common_ancestor_path) {
      if (stack_.empty()) {
        nav_from_root_to(to_pos);
        return;
      };
      const stack_entry& node = stack_.top();
      node_idx = node.node_idx;
      level = node.level;
      path = node.path;
      // Check if we missed the right child of the common ancestor (could
      // happen because it might not be on the stack).
      //
      // Note: It is no longer guaranteed, that the common ancestor is on
      //   the stack since we push only inner nodes and leaf nodes with 1-labels
      //   on the stack.  Which means that if we cannot find the common ancestor
      //   on the stack, it is either (i) an implicit node or (ii) it is a leaf
      //   node with a 0-label. The first case is handled before we walk the
      //   tree upwards, and in the second case we forward the iterator to
      //   the next 1-fill.
      if (level < right_child_of_common_ancestor_level) {
        next();
        return;
      }
      stack_.pop();
    }
    node_idx_ = node_idx;
    path_ = path;
    nav_downwards(to_pos);
  }

  /// ------------------ FOR BENCHMARKING PURPOSES ONLY. -----------------------
  /// Used to determine the costs for upwards navigation.
  void __attribute__((noinline))
  bench_nav_upwards(
      const path_t right_child_of_common_ancestor_path,
      u64 right_child_of_common_ancestor_level) {
    // clang-format off
// --- slower ---
//    // Walk up the tree to the right child of the common ancestor.
//    path_t path = 1;
//    $u1 found = false;
//    $u1 missed = false;
//    std::size_t i = stack_.size();
//    for (; i > 0; --i) {
//      path = stack_[i - 1].path;
//      const auto level = stack_[i - 1].level;
//      found = (path == right_child_of_common_ancestor_path);
//      missed = (level < right_child_of_common_ancestor_level);
//      if (found | missed) break;
//    }
//
//    if (missed) {
//      node_idx_ = 0;
//      path_ = 1;
//      return;
//    }
//
//    if (found) {
//      stack_.rewind(i - 1);
//      node_idx_ = stack_[i - 1].node_idx;
//      path_ = path;
//      return;
//    }
//
//    node_idx_ = 0;
//    path_ = 1;
//    return;
// --- ------ ---
    // clang-format on

    // Walk up the tree to the right child of the common ancestor.
    $u64 node_idx = node_idx_;
    auto path = path_;
    $u64 level = tree_height_;
    while (path != right_child_of_common_ancestor_path) {
      if (stack_.empty()) {
        node_idx_ = 0;
        path_ = 1;
        return;
      };
      const stack_entry& node = stack_.top();
      node_idx = node.node_idx;
      level = node.level;
      path = node.path;
      if (level < right_child_of_common_ancestor_level) {
        node_idx_ = 0;
        path_ = 1;
        return;
      }
      stack_.pop();
    }
    node_idx_ = node_idx;
    path_ = path;
  }
  auto __teb_inline__
  bench_nav_upwards_get_stack_size() {
    return stack_.size();
  }
  /// --------------------------------------------------------------------------

  /// Fast-forwards the iterator to the given position.
  void __teb_inline__
  skip_to(const std::size_t to_pos) noexcept {
    if (to_pos >= teb_.n_actual_) { // TODO remove condition
      pos_ = teb_.n_actual_;
      length_ = 0;
      return;
    }
    if (to_pos < (pos_ + length_)) {
      length_ -= to_pos - pos_;
      pos_ = to_pos;
      return;
    }
    nav_to(to_pos);
  }

  /// Returns true if the iterator reached the end, false otherwise.
  u1 __forceinline__
  end() const noexcept {
    return pos_ == teb_.n_actual_;
  }

  /// Returns the starting position of the current 1-fill.
  u64 __forceinline__
  pos() const noexcept {
    return pos_;
  }

  /// Returns the length of the current 1-fill.
  u64 __forceinline__
  length() const noexcept {
    return length_;
  }

  /// Returns the path of the current tree node.
  u64 __teb_inline__
  path() const noexcept {
    return path_;
  }

  /// Returns the level of the current tree node.
  u64 __teb_inline__
  level() const noexcept {
    return dtl::teb_util::determine_level_of(path_);
  }

  /// Returns the number of perfect tree levels.
  u64 __teb_inline__
  perfect_levels() const noexcept {
    return perfect_levels_;
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
