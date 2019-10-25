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
    $u64 n = node_idx;
    const auto is_left_sibling = node_idx & 1ull;
    n += is_left_sibling;
    const auto is_right_sibling = (node_idx & 1ull) == 0 ? 1 : 0;
    n -= is_right_sibling;
    return n;
  }

  /// Returns the level of the given node.
  static constexpr inline u64
  level_of(u64 node_idx) {
    return log_2(node_idx + 1);
  }

  inline u64
  last_level() {
    return dtl::log_2(n_);
  }

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
    is_inner_node_.set(node_idx + offset);
    is_active_node_.set(left_child_of(node_idx) + offset);
    is_active_node_.set(right_child_of(node_idx) + offset);
  }

  /// Returns true if given node is expanded, false otherwise. A tree node is
  /// expanded iff the parent is an inner node.
  //  inline u1
  //  is_expanded(u64 node_idx) const {
  //    // TODO: why do we not allow to collapse the root?
  //    if (node_idx == root()) return true;
  //    return is_inner_node(parent_of(node_idx));
  //  }

  /// Returns true if the given node is not expanded, false otherwise.
  //  inline u1
  //  is_collapsed(u64 node_idx) const {
  //    return !is_expanded(node_idx);
  //  }

  /// Alias for 'is_expanded'.
  //  inline u1 contains(u64 node_idx) const {
  //    return is_expanded(node_idx);
  //  }


  //===--------------------------------------------------------------------===//
  struct node_t {
    $u64 idx_;
    $u64 level_;
    $u1 is_inner_;

    inline auto idx() { return idx_; }
    inline auto level() { return level_; }
    inline auto is_inner() { return is_inner_; }
    inline u1
    operator==(const node_t& other) const {
      return idx_ == other.idx_; // && level == other.level;
    }
  };
  //===--------------------------------------------------------------------===//
  struct node_simple_t {
//    const binary_tree_structure& tree_;
    $u64 idx_;

    inline auto idx() { return idx_; }
//    inline auto level() { return tree_.level_of(idx_); }
//    inline auto is_inner() { return tree_.is_inner_node(idx_); }

    inline u1
    operator==(const node_t& other) const {
      return idx_ == other.idx_; // && level == other.level;
    }
  };
  //===--------------------------------------------------------------------===//
  class breadth_first_iterator
      : public std::iterator<
            std::input_iterator_tag, // iterator_category
            $u64, // value_type
//            node_simple_t, // value_type
            u64, // difference_type
            $u64*, // pointer
//            const node_simple_t*, // pointer
            $u64 // reference
//            node_simple_t // reference
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
//      return node_simple_t { tree_, idx_ };
    }
  };
  //===--------------------------------------------------------------------===//
  class const_breadth_first_iterator
      : public std::iterator<
            std::input_iterator_tag, // iterator_category
            node_t, // value_type
            u64, // difference_type
            const node_t*, // pointer
            node_t // reference
            > {
    const binary_tree_structure& tree_;

    $u64 idx_;

    static constexpr std::size_t buf_size = 128;
    std::array<node_t, buf_size> buf_;
    std::size_t buf_read_idx_ = 0;
    std::size_t buf_end_ = 0;

  public:
    explicit inline const_breadth_first_iterator(const binary_tree_structure& tree,
        u64 start_node_idx)
        : tree_(tree), idx_(start_node_idx),
          buf_read_idx_(0), buf_end_(0) {
      assert(start_node_idx <= tree_.max_node_cnt_);
      if (start_node_idx < tree_.max_node_cnt_) {
        assert(tree_.is_active_node(start_node_idx));
        buf_[0] = node_t { start_node_idx, tree_.level_of(start_node_idx),
            tree_.is_inner_node(0) };
        buf_end_ = 1;
        ++idx_;
      }
      else {
        buf_[0] = node_t { tree_.max_node_cnt_, tree_.height_,
            false };
        buf_end_ = 1;
        ++idx_;
      }
    }

//    __forceinline__
    void __attribute__((noinline))
    next_batch() {
//      assert(level_of(idx_) > 0);
//      // Clear the buffer.
//      buf_read_idx_ = 0;
//      buf_end_ = 0;
//      ++idx_;
//      while (idx_ < tree_.max_node_cnt_ && buf_end_ < (buf_size - 64)) {
//        const auto this_level = level_of(idx_);
////        const auto parent_level = this_level - 1;
////        const auto parent_idx = tree_.parent_of(idx_);
//
//        if (this_level < 6) {
//          if (tree_.is_active_node(idx_)) {
//            u1 is_inner = tree_.is_inner_node(idx_);
//            buf_[buf_end_] = node_t { idx_, this_level, is_inner };
//            ++buf_end_;
//          }
//        }
//        else {
//          const auto remaining_bits_at_this_level =
//              tree_.first_node_idx_at_level(this_level + 1) - idx_;
//
//          const auto remaining_bits_at_parent_level =
//              tree_.first_node_idx_at_level(parent_level + 1) - parent_idx;
//
//          const auto buf_size_this_level =
//              std::min(u64(64), remaining_bits_at_this_level);
//          $u64 buf_this_level = tree_.is_inner_node_.fetch_bits(idx_ + offset,
//              idx_ + offset + buf_size_this_level);
//
//          const auto buf_size_parent_level =
//              std::min(u64(32), remaining_bits_at_parent_level);
//          $u64 buf_parent_level = tree_.is_inner_node_.fetch_bits(parent_idx + offset,
//              parent_idx + offset + buf_size_parent_level);
//
//          assert(buf_size_parent_level * 2 <= buf_size_this_level);
//
////          if (buf_parent_level != 0 || buf_this_level != 0) {
//          if (buf_parent_level != 0) {
//            for (std::size_t i = 0; i < buf_size_parent_level; ++i) {
//              u1 parent_is_inner = ((buf_parent_level >> i) & 1) != 0;
//              if (!parent_is_inner) continue;
//              u1 left_child_is_inner = ((buf_this_level >> (i * 2)) & 1) != 0;
//              u1 right_child_is_inner = ((buf_this_level >> (i * 2 + 1)) & 1) != 0;
//              buf_[buf_end_] =
//                  node_t { idx_ + (i * 2), this_level, left_child_is_inner };
//              buf_[buf_end_ + 1] =
//                  node_t { idx_ + (i * 2 + 1), this_level, right_child_is_inner };
//              buf_end_ += 2;
//            }
//          }
//          idx_ += buf_size_this_level;
//        }
//        ++idx_;
//      }
    }

    inline const_breadth_first_iterator&
    operator++() {
      ++buf_read_idx_;
      if (buf_read_idx_ >= buf_end_) {
        next_batch();
      }
      return *this;
    }

    inline const_breadth_first_iterator
    operator++(int) {
      const_breadth_first_iterator ret_val = *this;
      ++(*this);
      return ret_val;
    }

    inline bool
    operator==(const const_breadth_first_iterator& other) const {
      return (*(*this)).idx() == (*other).idx();
    }

    inline bool
    operator!=(const const_breadth_first_iterator& other) const {
      return !(*this == other);
    }

    inline reference
    operator*() const {
      if (buf_read_idx_ >= buf_end_) {
        return node_t { tree_.max_node_cnt_, tree_.height_, false };
      }
      return buf_[buf_read_idx_];
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

  /// Returns a (const) breadth-first iterator. The iterator internally uses a
  /// buffer to improve performance.
  inline const_breadth_first_iterator
  const_breadth_first_begin() const {
    return const_breadth_first_iterator(*this, root());
  }

  /// Returns a (const) breadth-first iterator.
  inline const_breadth_first_iterator
  const_breadth_first_begin(std::size_t start_node_idx) const {
    return const_breadth_first_iterator(*this, start_node_idx);
  }

  /// Returns a (const) breadth-first iterator that points one past the last node.
  inline const_breadth_first_iterator
  const_breadth_first_end() const {
    return const_breadth_first_iterator(*this, max_node_cnt_);
  }

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
