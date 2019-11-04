#pragma once
//===----------------------------------------------------------------------===//
#include "plain_bitmap.hpp"

#include <dtl/dtl.hpp>
#include <dtl/math.hpp>

#include <bitset>
#include <cassert>
#include <functional>
#include <queue>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// Represents a binary tree structure with a limited number of nodes. The Tree
/// is stored as an implicit data structure in breadth-first order. Each node
/// has an unique identifier which allow for easy navigation inside the tree.
/// I.e., for a given node i, the left childs identifier is 2i+1, and 2i+2 is
/// the identifier of the right child node; the parent node ID is (i-1)/2.
class binary_tree_structure {
  using bitmap_t = plain_bitmap<$u64>;

public: // TODO make private
  /// The number of leaf nodes.
  $u64 n_;
  /// The max number of nodes.
  $u64 max_node_cnt_;
  /// The tree height.
  $u64 height_;
  /// Indicates whether a node is an inner or a leaf node.
  bitmap_t is_inner_node_;
  /// Indicates whether a node is part of the tree. I.e., its parent is an inner
  /// node.
  bitmap_t is_active_node_;

  /// The node with index i is stored in the bitmap at position i+offset. The +1
  /// offset results in better word alignment, which helps to speed up several
  /// tree algorithms.
  static constexpr std::size_t offset = 1;

public:
  /// Constructs a perfect binary tree (structure) with n leaf nodes and n-1
  /// inner nodes.  Note, n must be a power of two.
  explicit binary_tree_structure(u64 n, u1 init = true)
      : n_(n),
        max_node_cnt_(2 * n_ - 1),
        height_(dtl::log_2(n_)),
        is_inner_node_(max_node_cnt_ + offset),
        is_active_node_(max_node_cnt_ + offset) { // NOLINT

    if (!dtl::is_power_of_two(n_)) {
      throw std::invalid_argument(
          "The number of leaf nodes must be a power of two.");
    }

    // Initialize a perfect binary tree.
    if (init) {
      // The first n-1 nodes are inner nodes.
      is_inner_node_.set(0, (max_node_cnt_ / 2) + offset);
      // All nodes are active in the beginning.
      is_active_node_.set(0, is_active_node_.size()); // TODO decide how to init the bits [0, offset)
    }
  }

  binary_tree_structure(const binary_tree_structure& other) = default;
  binary_tree_structure(binary_tree_structure&& other) noexcept = default;
  binary_tree_structure& operator=(const binary_tree_structure& other) = default;
  binary_tree_structure& operator=(binary_tree_structure&& other) noexcept = default;
  virtual ~binary_tree_structure() = default;

  /// Counts the number of nodes in the given subtree.
  u64 subtree_size(u64 node_idx) const {
    if (!is_inner_node(node_idx)) return 1;
    return 1
        + subtree_size(left_child_of(node_idx))
        + subtree_size(right_child_of(node_idx));
  };

  /// Counts the number of leaf nodes in the given subtree.
  u64 count_leaf_nodes(u64 node_idx) const {
    if (!is_inner_node(node_idx)) return 1;
    return count_leaf_nodes(left_child_of(node_idx))
        + count_leaf_nodes(right_child_of(node_idx));
  };

  /// Returns the ID of the root node.
  static constexpr inline u64
  root() {
    return 0;
  }

  /// Returns the ID of the parent node.
  static constexpr inline u64
  parent_of(u64 node_idx) {
    return (node_idx - 1) / 2;
  }

  /// Returns the ID of left child node.
  static constexpr inline u64
  left_child_of(u64 node_idx) {
    return 2 * node_idx + 1;
  }

  /// Returns the ID of right child node.
  static constexpr inline u64
  right_child_of(u64 node_idx) {
    return 2 * node_idx + 2;
  }

  /// Returns the ID of the sibling node.
  static constexpr inline u64
  sibling_of(u64 node_idx) {
    return node_idx ^ u64(1);
  }

  /// Returns the level of the given node.
  static constexpr inline u64
  level_of(u64 node_idx) {
    return log_2(node_idx + 1);
  }

  /// Returns the number of the last tree level.
  inline u64
  last_level() {
    return dtl::log_2(n_);
  }

  /// Returns the idx of first node at the given level.
  static constexpr inline u64
  first_node_idx_at_level(u64 level) {
    return (1ull << level) - 1;
  }

  /// Returns true, if the given node is active (exists); false otherwise.
  inline u1
  is_active_node(u64 node_idx) const {
    return is_active_node_[node_idx + offset];
  }

  /// Returns true, if the given node is an inner node; false otherwise.
  inline u1
  is_inner_node(u64 node_idx) const {
    return is_inner_node_[node_idx + offset];
  }

  /// Returns true, if the given node is a leaf node; false otherwise.
  inline u1
  is_leaf_node(u64 node_idx) const {
    return !is_inner_node_[node_idx + offset];
  }

  /// Turns the given node into a leaf node.
  inline void
  set_leaf(u64 node_idx) {
    set_leaf_rec(node_idx);
  }

  /// Traverses the tree starting at the given node index. For each visited
  /// node, the function 'fn' is called whereas the first argument is the
  /// current node id.
  inline void
  visit(u64 node_idx, const std::function<void(u64)>& fn) {
    visit_rec(node_idx, fn);
  }

  /// Turns the given node into an inner node.
  inline void
  set_inner(u64 node_idx) {
    assert(node_idx == root() || is_inner_node(parent_of(node_idx)));
    assert((left_child_of(node_idx) + offset) / (sizeof($u64) * 8)
        == (right_child_of(node_idx) + offset) / (sizeof($u64) * 8));
    is_inner_node_.set(node_idx + offset);
    // Set the two child nodes active. The specialized function below sets two
    // bits at once.
    is_active_node_.set2(left_child_of(node_idx) + offset);
  }

  /// Turns the given node into an inner node.
  /// WARNING: This function does NOT update the active node bitmap for
  ///          performance reasons. The caller has to take care, that the
  ///          active flags are set manually.
  inline void
  set_inner_fast(u64 node_idx) {
    assert(node_idx == root() || is_inner_node(parent_of(node_idx)));
    is_inner_node_.set(node_idx + offset);
    assert((left_child_of(node_idx) + offset) / (sizeof($u64) * 8)
        == (right_child_of(node_idx) + offset) / (sizeof($u64) * 8));
  }

  //===--------------------------------------------------------------------===//
  class breadth_first_iterator
      : public std::iterator<
            std::input_iterator_tag, // iterator_category
            $u64, // value_type
            u64, // difference_type
            $u64*, // pointer
            $u64 // reference
            > {
    const binary_tree_structure& tree_;
    $u64 idx_;

  public:
    explicit inline breadth_first_iterator(const binary_tree_structure& tree,
        u64 start_node_idx)
        : tree_(tree),
          idx_(start_node_idx) {
      assert(start_node_idx <= tree_.max_node_cnt_);
      if (start_node_idx < tree_.max_node_cnt_) {
        assert(tree_.is_active_node(start_node_idx));
      }
    }

    inline breadth_first_iterator&
    operator++() {
      assert(idx_ < tree_.max_node_cnt_);
      idx_ = tree_.is_active_node_.find_next(idx_ + tree_.offset) - tree_.offset;
      return *this;
    }

    inline breadth_first_iterator
    operator++(int) {
      breadth_first_iterator ret_val = *this;
      ++(*this);
      return ret_val;
    }

    inline bool
    operator==(const breadth_first_iterator& other) const {
      return idx_ == other.idx_;
    }

    inline bool
    operator!=(const breadth_first_iterator& other) const {
      return !(*this == other);
    }

    inline reference
    operator*() const {
      assert(idx_ <= tree_.max_node_cnt_);
      return idx_;
    }
  };
  //===--------------------------------------------------------------------===//

  /// Returns a breadth-first iterator.
  inline breadth_first_iterator
  breadth_first_begin() const {
    return breadth_first_iterator(*this, root());
  }

  /// Returns a breadth-first iterator.
  inline breadth_first_iterator
  breadth_first_begin(std::size_t start_node_idx) const {
    return breadth_first_iterator(*this, start_node_idx);
  }

  /// Returns a breadth-first iterator that points one past the last node.
  inline breadth_first_iterator
  breadth_first_end() const {
    return breadth_first_iterator(*this, max_node_cnt_);
  }

  //===--------------------------------------------------------------------===//
  // Validation code for debug builds.
  //===--------------------------------------------------------------------===//
  /// Validates that the active flags are set correctly.
  void
  validate_active_nodes() {
    for (std::size_t i = 0; i < this->max_node_cnt_; ++i) {
      u1 is_inner = this->is_inner_node(i);
      if (is_inner) {
        assert(this->is_active_node(i));
        assert(this->is_active_node(this->left_child_of(i)));
        assert(this->is_active_node(this->right_child_of(i)));
        if (i != this->root()) {
          assert(this->is_active_node(this->parent_of(i)));
        }
      }
      else {
        if (i == this->root()) {
          assert(this->is_active_node(i));
        }
        else {
          const auto parent = this->parent_of(i);
          if (this->is_inner_node(parent)) {
            assert(this->is_active_node(parent));
            assert(this->is_active_node(i));
          }
          else {
            assert(this->is_active_node(i) == false);
          }
        }
      }
    }
  }

  /// Sets the active flags.
  void
  fix_active_nodes() {
    auto act = [&](std::size_t i, u1 val) {
      this->is_active_node_.set(i + this->offset, val);
    };
    act(0, true); // root
    for (std::size_t i = 1; i < this->max_node_cnt_; ++i) {
      u1 is_inner = this->is_inner_node(i);
      assert(this->is_inner_node(this->parent_of(i)));
      if (is_inner) {
        act(i, true);
        act(this->left_child_of(i), true);
        act(this->right_child_of(i), true);
      }
      else {
        const auto parent = this->parent_of(i);
        if (this->is_inner_node(parent)) {
          act(i, true);
        }
      }
    }
  }
  //===--------------------------------------------------------------------===//

private:
  /// Mark the given node as a leaf node. The function propagates the call
  /// to the child nodes when those are inner nodes.
  inline void
  set_leaf_rec(u64 node_idx) {
    u1 recurse = is_inner_node_[node_idx + offset];
    is_inner_node_.clear(node_idx + offset);
    if (recurse) {
      const auto left_child_idx = left_child_of(node_idx);
      const auto right_child_idx = right_child_of(node_idx);
      is_active_node_.clear(left_child_idx + offset);
      is_active_node_.clear(right_child_idx + offset);
      set_leaf_rec(left_child_idx);
      set_leaf_rec(right_child_idx);
    }
  }

  inline void
  visit_rec(u64 node_idx, const std::function<void(u64)>& fn) const {
    u1 recurse = is_inner_node_[node_idx + offset];
    fn(node_idx);
    if (recurse) {
      visit_rec(left_child_of(node_idx), fn);
      visit_rec(right_child_of(node_idx), fn);
    }
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
