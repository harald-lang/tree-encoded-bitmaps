#pragma once

#include <bitset>

#include <dtl/dtl.hpp>
#include <dtl/math.hpp>

namespace dtl {

//===----------------------------------------------------------------------===//
/// Represents a binary tree structure with a limited number of nodes. The Tree
/// is stored as an implicit data structure in breadth-first order. Each node
/// has an unique identifier which allow for easy navigation inside the tree.
/// I.e., for a given node i, the left childs identifier is 2i+1, and 2i+2 is
/// the identifier of the right child node; the parent node ID is (i-1)/2.
class binary_tree_structure {

  using bitmap_t = boost::dynamic_bitset<$u32>;

public:

  /// The number of leaf nodes.
  $u64 n_;

  /// The max number of nodes.
  $u64 max_node_cnt_;

  /// The tree height.
  $u64 height_;

  /// Indicates whether a node is an internal or a leaf node.
  bitmap_t is_inner_node_;

public:

  /// Constructs a perfect binary tree (structure) with n leaf nodes and n-1
  /// inner nodes.
  explicit
  binary_tree_structure(u64 n)
      : n_(n),
        max_node_cnt_(2 * n_ - 1),
        height_(dtl::log_2(n_)),
        is_inner_node_(max_node_cnt_) {

    // Initialize a perfect binary tree.
    for ($u64 i = 0; i < max_node_cnt_ / 2; i++) {
      is_inner_node_[i] = true;
    }
  }

  binary_tree_structure(const binary_tree_structure& other) = default;
  binary_tree_structure(binary_tree_structure&& other) noexcept = default;
  binary_tree_structure& operator=(const binary_tree_structure& other) = default;
  binary_tree_structure& operator=(binary_tree_structure&& other) noexcept = default;
  virtual ~binary_tree_structure() = default;

  /// Counts the number of nodes in the given subtree.
  u64
  subtree_size(u64 node_idx) const {
    if (!is_inner_node(node_idx)) return 1;
    return 1
         + subtree_size(left_child_of(node_idx))
         + subtree_size(right_child_of(node_idx));
  };

  /// Counts the number of leaf nodes in the given subtree.
  u64
  count_leaf_nodes(u64 node_idx) const {
    if (!is_inner_node(node_idx)) return 1;
    return count_leaf_nodes(left_child_of(node_idx))
         + count_leaf_nodes(right_child_of(node_idx));
  };

  /// Returns the ID of the root node.
  static inline u64
  root() {
    return 0;
  }

  /// Returns the ID of the parent node.
  static inline u64
  parent_of(u64 node_idx) {
    return (node_idx - 1) / 2;
  }

  /// Returns the ID of left child node.
  static inline u64
  left_child_of(u64 node_idx) {
    return 2 * node_idx + 1;
  }

  /// Returns the ID of right child node.
  static inline u64
  right_child_of(u64 node_idx) {
    return 2 * node_idx + 2;
  }

  /// Returns the level of the given node.
  static inline u64
  level_of(u64 node_idx) {
    return log_2(node_idx + 1);
  }

  /// Returns true, if the given node is an inner node; false otherwise.
  inline u1
  is_inner_node(u64 node_idx) const {
    return is_inner_node_[node_idx];
  }

  /// Returns true, if the given node is a leaf node; false otherwise.
  inline u1
  is_leaf_node(u64 node_idx) const {
    return ! is_inner_node_[node_idx];
  }

  /// Turns the given node into a leaf node.
  inline void
  set_leaf(u64 node_idx) {
    set_leaf_rec(node_idx);
  }

  /// Turns the given node into an inner node.
  inline void
  set_inner(u64 node_idx) {
    is_inner_node_[node_idx] = true;
  }

  //===--------------------------------------------------------------------===//
  class breadth_first_iterator : public std::iterator<
      std::input_iterator_tag,  // iterator_category
      u64,                      // value_type
      u64,                      // difference_type
      const u64*,               // pointer
      u64                       // reference
    >{
    const binary_tree_structure& tree_;
    std::queue<$u64> fifo_;

  public:

    explicit
    breadth_first_iterator(const binary_tree_structure& tree,
        u64 start_node_idx)
        : tree_(tree) {
      fifo_.push(start_node_idx);
    }

    inline breadth_first_iterator&
    operator++() {
      if (!fifo_.empty()) {
        const auto current_node_idx = fifo_.front();
        if (tree_.is_inner_node(current_node_idx)) {
          fifo_.push(tree_.left_child_of(current_node_idx));
          fifo_.push(tree_.right_child_of(current_node_idx));
        }
        fifo_.pop();
      }
      return *this;
    }

    inline breadth_first_iterator
    operator++(int) {
      breadth_first_iterator ret_val = *this;
      ++(*this);
      return ret_val;
    }

    inline bool
    operator==(breadth_first_iterator other) const {
      return **this == *other;
    }

    inline bool
    operator!=(breadth_first_iterator other) const {
      return !(*this == other);
    }

    inline reference
    operator*() const {
      if (fifo_.empty()) {
        return tree_.max_node_cnt_;
      }
      return fifo_.front();
    }
  };
  //===--------------------------------------------------------------------===//

  /// Returns a breadth-first iterator.
  inline breadth_first_iterator
  breadth_first_begin() {
    return breadth_first_iterator(*this, root());
  }

  inline breadth_first_iterator
  breadth_first_end() {
    return breadth_first_iterator(*this, max_node_cnt_);
  }

private:

  /// Mark the given node as a leaf node.
  inline void
  set_leaf_rec(u64 node_idx) {
    const auto recurse = is_inner_node_[node_idx];
    is_inner_node_[node_idx] = false;
    if (recurse) {
      set_leaf_rec(left_child_of(node_idx));
      set_leaf_rec(right_child_of(node_idx));
    }
  }

};
//===----------------------------------------------------------------------===//

} // namespace dtl
