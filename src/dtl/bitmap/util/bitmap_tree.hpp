#pragma once
//===----------------------------------------------------------------------===//
#include "binary_tree_structure.hpp"
#include "plain_bitmap.hpp"
#include "rank1.hpp"
#include "rank1_logic_surf.hpp"

#include <dtl/dtl.hpp>

#include <boost/dynamic_bitset.hpp>

#include <cstring>
#include <iomanip>
#include <immintrin.h>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// Represents a bitmap as a binary tree. During construction, the bitmap tree
/// is compressed. An instance can be seen as an intermediate representation of
/// a bitmap from which a TEB is constructed.
template<i32 optimization_level_ = 3, u1 fast_path = true>
class bitmap_tree : public binary_tree_structure {
  using tree_t = dtl::binary_tree_structure;
  using bitmap_t = boost::dynamic_bitset<$u32>;

public: // TODO make protected
  /// The labels of the tree nodes.
  plain_bitmap<$u64> labels_;
protected:

  /// Represents an integer range [begin, end).
  struct range_t {
    std::size_t begin = 0;
    std::size_t end = 0;
    /// Returns true if the range is empty, false otherwise.
    inline u1 is_empty() const { return begin >= end; };
    /// Returns true if the given values is within the range.
    inline u1 contains(std::size_t i) const { return (i >= begin) && (i < end); };
    /// Returns length of the range.
    inline std::size_t length() const { return end - begin; };
  };

  //===--------------------------------------------------------------------===//
  // The following counters are required to compute the size of the resulting
  // tree-encoded bitmap. These counters are updated when the tree is modified.
  //
  // During the construction phase, the size is required to find the tree
  // instance with the lowest memory consumption. Recall that the number of
  // tree nodes does not directly correspond to the TEB size. A fully pruned
  // tree instance may result in a larger TEB than a partially pruned tree
  // instance.
  //===--------------------------------------------------------------------===//
  /// The number of inner nodes. As the tree is a full binary, the number of
  /// leaf nodes is the number of inner nodes + 1.
  std::size_t inner_node_cnt_ = 0;
  /// The number inner nodes until the first leaf occurs (in level order).
  std::size_t leading_inner_node_cnt_ = 0;
  /// The number leaf nodes after the last inner node (in level order).
  std::size_t trailing_leaf_node_cnt_ = 0;
  /// The number of 0-labels until the first 1-label occurs.
  std::size_t leading_0label_cnt_ = 0;
  /// The number of 0-labels after the last 1-label.
  std::size_t trailing_0label_cnt_ = 0;

  //===--------------------------------------------------------------------===//
  // The ranges below are used to cache the node indexes of the first and last
  // explicit node and label. This is required to incrementally update the
  // counters from above. Unlike to the counters from above, these ranges are
  // not required to estimate the size of a tree-encoded bitmap.
  //===--------------------------------------------------------------------===//
  /// The node indexes that need to be stored explicitly.
  range_t explicit_node_idxs_;
  /// The indexes of nodes with explicit labels.
  //  range_t explicit_label_idxs_; // TODO use range rather than the two variables below
  /// The first leaf node index that carries a 1-label.
  std::size_t first_node_idx_with_1label_ = 0;
  /// The last leaf node index that carries a 1-label.
  std::size_t last_node_idx_with_1label_ = 0;

  //===--------------------------------------------------------------------===//
  // The following two members cache the positions of the first and the last
  // 1-bit within the original bitmap. This allows for faster initialization
  // of the leading/trailing 0-label counters.
  //===--------------------------------------------------------------------===//
  /// The index of the first set bit in the original bitmap. The value is set
  /// to 'n' when the bits in the original bitmap are all 0.
  std::size_t first_bit_idx_ = 0;
  /// The index of the last set bit in the original bitmap. The value is set
  /// to 'n' when the bits in the original bitmap are all 0.
  std::size_t last_bit_idx_ = 0;

  /// Used to keep track of whether the counters are up to date or invalidated.
  $u1 counters_are_valid = false;

  std::size_t uncompressed_size = 0;

public:
  /// C'tor
  explicit bitmap_tree(const bitmap_t& bitmap)
      : binary_tree_structure(bitmap.size()),
        labels_(max_node_cnt_ + offset),
        inner_node_cnt_(0),
        leading_inner_node_cnt_(0),
        trailing_leaf_node_cnt_(0),
        leading_0label_cnt_(0),
        trailing_0label_cnt_(0),
        first_node_idx_with_1label_(0),
        last_node_idx_with_1label_(0),
        first_bit_idx_(0), // first and last bit idx will be initialized in init_tree()
        last_bit_idx_(0) {
    // Init the binary tree.
    init_tree(bitmap);
    uncompressed_size = estimate_encoded_size_in_bytes();

    if (!fast_path) {
      // Fused code path.
      compress();
    }
    else {
      // Classic code path.
      init_labels();
      prune_tree();

      // Run space optimizations.
      if (optimization_level_ > 1) {
        init_counters();
        run_optimize();
      }
    }
  }

protected:
  /// C'tor
  explicit bitmap_tree(const std::size_t n, u1 init = true)
      : binary_tree_structure(n, init),
        labels_(max_node_cnt_ + offset),
        inner_node_cnt_(0),
        leading_inner_node_cnt_(0),
        trailing_leaf_node_cnt_(0),
        leading_0label_cnt_(0),
        trailing_0label_cnt_(0),
        first_node_idx_with_1label_(0),
        last_node_idx_with_1label_(0),
        first_bit_idx_(0), // first and last bit idx will be initialized in init_tree()
        last_bit_idx_(0) {
  }

public:
  bitmap_tree(const bitmap_tree& other) = default;
  bitmap_tree(bitmap_tree&& other) noexcept = default;
  bitmap_tree& operator=(const bitmap_tree& other) = default;
  bitmap_tree& operator=(bitmap_tree&& other) noexcept = default;
  ~bitmap_tree() override = default;

  /// Initialize a perfect binary tree on top of the given bitmap, where
  ///  - all the inner nodes have two children
  ///  - all the leaf nodes are on the same level
  ///  - the leaf nodes are labelled with the given bitmap
  void __attribute__((noinline, hot, flatten))
  init_tree(const boost::dynamic_bitset<$u32>& bitmap) {
    // TODO: Support arbitrarily sized bitmaps.
    if (!dtl::is_power_of_two(n_)) {
      throw std::invalid_argument(
          "The length of the bitmap must be a power of two.");
    }

    u64 length = max_node_cnt_;
    u64 height = height_;
    // Copy the input bitmap to the last level of the tree.
    {
      auto i = bitmap.find_first();
      first_bit_idx_ = (i != boost::dynamic_bitset<$u32>::npos) ? i : n_;
      last_bit_idx_ = first_bit_idx_;
      while (i != boost::dynamic_bitset<$u32>::npos) {
        last_bit_idx_ = i;
        labels_.set(length / 2 + i + offset);
        i = bitmap.find_next(i);
      }
    }

    // Init the counters that are required to estimate the TEB size.
    init_counters_perfect_binary_tree();
  }

  /// Propagate the label bits along the tree (bottom-up).  The labels of an
  /// internal node is the bitwise OR of the labels of both child nodes.
  /// Pre-computing the labels of ALL tree nodes helps to avoid ad hoc
  /// computations of the labels during collapsing (or expanding) the tree
  /// structure.
  void __attribute__((noinline))
  init_labels() {
    for (auto level = last_level(); level > 0; --level) {
      const auto src_node_idx_begin = first_node_idx_at_level(level);
      const auto src_node_idx_end = first_node_idx_at_level(level + 1);
      const auto dst_node_idx_begin = first_node_idx_at_level(level - 1);
      const auto dst_node_idx_end = first_node_idx_at_level(level);
      auto src_node_idx = src_node_idx_begin;
      auto dst_node_idx = dst_node_idx_begin;
      while (src_node_idx < src_node_idx_end) {
        const auto remaining = src_node_idx_end - src_node_idx;
        assert(remaining >= 2);
        if (remaining >= 64) {
          // Process 64 nodes at a time.  For each loaded 64-bit word we write
          // a 32-bit word. For this to be efficient, we want proper alignment
          // which is why we have the +1 offset in the label and tree bitmaps.
          const auto src_label_idx = label_idx_of_node(src_node_idx);
          const auto dst_label_idx = label_idx_of_node(dst_node_idx);
          assert(src_label_idx % 64 == 0);
          assert(src_label_idx % 32 == 0);

          auto* raw_ptr = labels_.data();
          $u64* src_ptr = &raw_ptr[src_label_idx / 64];
          $u32* dst_ptr = &reinterpret_cast<$u32*>(raw_ptr)[dst_label_idx / 32];

          auto src = *src_ptr;
          src |= src >> 1;
          auto dst = _pext_u64(src, 0x5555555555555555); // FIXME non-portable code
          *dst_ptr = static_cast<$u32>(dst);

          src_node_idx += 64;
          dst_node_idx += 32;
        }
        else {
          // Process two nodes (siblings) at a time.
          u1 label_0 = label_of_node(src_node_idx);
          u1 label_1 = label_of_node(src_node_idx + 1);
          u1 parent_label = label_0 | label_1;
          const auto parent_label_idx = label_idx_of_node(dst_node_idx);
          labels_.set(parent_label_idx, parent_label);
          src_node_idx += 2;
          dst_node_idx += 1;
        }
      }
    }
  }

  /// Bottom-up pruning (loss-less).  Eliminate all sibling leaf nodes which
  /// have the same label. The algorithm terminates when all pairs of sibling
  /// leaf nodes have different labels.
  /// Note: The counters are thereby invalidated.
  void __attribute__((noinline))
  prune_tree() {
    D(init_counters();)
    for (auto level = last_level(); level > 0; --level) {
      $u64 collapse_cnt = 0;
      const auto src_node_idx_begin = first_node_idx_at_level(level);
      const auto src_node_idx_end = first_node_idx_at_level(level + 1);
      const auto dst_node_idx_begin = first_node_idx_at_level(level - 1);
      const auto dst_node_idx_end = first_node_idx_at_level(level);
      auto src_node_idx = src_node_idx_begin;
      auto dst_node_idx = dst_node_idx_begin;
      while (src_node_idx < src_node_idx_end) {
        const auto remaining = src_node_idx_end - src_node_idx;
        assert(remaining >= 2);
        if (remaining < 64) {
          u1 left_bit = labels_[src_node_idx + offset];
          u1 right_bit = labels_[src_node_idx + 1 + offset];
          u1 prune_causes_false_positives = left_bit ^ right_bit;
          u1 both_nodes_are_leaves =
              !is_inner_node(src_node_idx)
              & !is_inner_node(src_node_idx + 1);
          u1 prune = both_nodes_are_leaves & !prune_causes_false_positives;
          if (prune) {
            binary_tree_structure::set_leaf(dst_node_idx); // FIXME inefficient
            ++collapse_cnt;
          }
          src_node_idx += 2;
          dst_node_idx += 1;
        }
        else {
          // Process 64 nodes at a time.
          const auto src_idx = label_idx_of_node(src_node_idx);
          const auto dst_idx = label_idx_of_node(dst_node_idx);
          assert(src_idx % 64 == 0);
          assert(src_idx % 32 == 0);

          const auto src_word_idx = src_idx / 64;

          $u64* raw_label_ptr = labels_.data();
          u64 src_labels = raw_label_ptr[src_word_idx];
          u32 collapse_causes_false_positives = static_cast<$u32>(
              _pext_u64(src_labels ^ (src_labels >> 1), 0x5555555555555555)); // FIXME non-portable code

          $u64* raw_node_ptr = is_inner_node_.data();
          u64 src_nodes = raw_node_ptr[src_word_idx];
          u32 both_nodes_are_leaves = static_cast<$u32>(
              _pext_u64(~src_nodes & (~src_nodes >> 1), 0x5555555555555555)); // FIXME non-portable code

          u32 collapse = both_nodes_are_leaves & ~collapse_causes_false_positives;
          collapse_cnt += dtl::bits::pop_count(collapse);

          $u32* dst_nodes_ptr = &reinterpret_cast<$u32*>(raw_node_ptr)[dst_idx / 32];
          *dst_nodes_ptr = (*dst_nodes_ptr) ^ collapse;

          // Set the pruned nodes inactive.
          $u64* raw_active_node_ptr = is_active_node_.data();
          auto src_active_nodes = raw_active_node_ptr[src_word_idx];
          $u64 a = _pdep_u64(~collapse, 0x5555555555555555);
          a = a | (a << 1);
          raw_active_node_ptr[src_word_idx] = a;

          src_node_idx += 64;
          dst_node_idx += 32;
        }
      }
      D(init_counters();)

      {
        auto node_cnt_in_next_higher_level = 1ull << (level - 1);
        auto leaf_node_cnt_in_next_higher_level = collapse_cnt;
        auto inner_node_cnt_in_next_higher_level =
            node_cnt_in_next_higher_level - collapse_cnt;
        if (leaf_node_cnt_in_next_higher_level < inner_node_cnt_in_next_higher_level/2) {
          break;
        }
      }
    }
#ifndef NDEBUG
    // Validate active nodes.
    assert(is_active_node(0));
    for (std::size_t i = 1; i < max_node_cnt_; ++i) {
      if (is_inner_node(i)) {
        assert(is_active_node(i));
      }
      else {
        assert(is_active_node(i) == is_inner_node(parent_of(i)));
      }
    }
#endif
    // Invalidate the counters.
    counters_are_valid = false;
  }

  /// Determine the total number of tree nodes and which of these nodes need
  /// to be stored explicitly.
  void __attribute__((noinline))
  init_counters() {
    inner_node_cnt_ = 0;
    leading_inner_node_cnt_ = 0;
    trailing_leaf_node_cnt_ = 0;
    explicit_node_idxs_.begin = 0;
    explicit_node_idxs_.end = 0;

    leading_0label_cnt_ = 0;
    trailing_0label_cnt_ = 0;
    first_node_idx_with_1label_ = max_node_cnt_;
    last_node_idx_with_1label_ = max_node_cnt_;

    auto active_leaf_nodes = is_active_node_.and_not(is_inner_node_);

    {
      auto active_leaf_nodes_with_1_labels = active_leaf_nodes & labels_;
      first_node_idx_with_1label_ =
          active_leaf_nodes_with_1_labels.find_first() - offset;
      last_node_idx_with_1label_ =
          active_leaf_nodes_with_1_labels.find_last() - offset;

      auto active_leaf_nodes_with_0_labels = active_leaf_nodes.and_not(labels_);
      leading_0label_cnt_ =
          active_leaf_nodes_with_0_labels.count(offset, first_node_idx_with_1label_ + offset);
      trailing_0label_cnt_ =
          active_leaf_nodes_with_0_labels.count(last_node_idx_with_1label_ + offset, max_node_cnt_ + offset);
    }

    {
      auto active_inner_nodes = is_active_node_ & is_inner_node_;
      inner_node_cnt_ = active_inner_nodes.count(offset, max_node_cnt_ + offset);
      explicit_node_idxs_.end = (active_inner_nodes.find_last() - offset) + 1;
    }

    {
      explicit_node_idxs_.begin =
          active_leaf_nodes.find_first(offset, max_node_cnt_ + offset) - offset;
      trailing_leaf_node_cnt_ =
          active_leaf_nodes.count(explicit_node_idxs_.end + offset, max_node_cnt_ + offset);
    }
    leading_inner_node_cnt_ = explicit_node_idxs_.begin;
    counters_are_valid = true;
#ifndef NDEBUG
    validate_counters();
#endif
  }

  /// Determine the total number of tree nodes and which of these nodes need
  /// to be stored explicitly.
  void __attribute__((noinline))
  init_counters_naive() {
    $u1 found_leaf_node = false;
    std::size_t node_cnt = 0;
    inner_node_cnt_ = 0;
    leading_inner_node_cnt_ = 0;
    trailing_leaf_node_cnt_ = 0;
    explicit_node_idxs_.begin = 0;
    explicit_node_idxs_.end = 0;

    $u1 found_leaf_node_with_1label = false;
    leading_0label_cnt_ = 0;
    trailing_0label_cnt_ = 0;
    first_node_idx_with_1label_ = max_node_cnt_;
    last_node_idx_with_1label_ = max_node_cnt_;

    const auto it_end = breadth_first_end();
    for (auto it = breadth_first_begin(); it != it_end; ++it) {
      u64 idx = (*it);
      u64 level = level_of(idx);
      u1 is_inner = is_inner_node(idx);
      u1 is_leaf = !is_inner;

      ++node_cnt;
      inner_node_cnt_ += is_inner;

      // Count the leading inner nodes.
      if (!found_leaf_node && is_inner) {
        ++leading_inner_node_cnt_;
      }
      if (!found_leaf_node && !is_inner) {
        found_leaf_node = true;
      }

      // Count the trailing leaf nodes.
      if (!is_inner) {
        ++trailing_leaf_node_cnt_;
      }
      else {
        trailing_leaf_node_cnt_ = 0;
        explicit_node_idxs_.end = idx + 1;
      }

      if (!found_leaf_node_with_1label && is_leaf) {
        if (label_of_node(idx) == false) {
          ++leading_0label_cnt_;
        }
        else {
          found_leaf_node_with_1label = true;
          first_node_idx_with_1label_ = idx;
        }
      }

      // Count the trailing 0-labels nodes.
      if (is_leaf) {
        if (label_of_node(idx) == false) {
          ++trailing_0label_cnt_;
        }
        else {
          trailing_0label_cnt_ = 0;
          last_node_idx_with_1label_ = idx;
        }
      }
    }
    explicit_node_idxs_.begin = leading_inner_node_cnt_;
    if (explicit_node_idxs_.end < explicit_node_idxs_.begin) {
      explicit_node_idxs_.end = explicit_node_idxs_.begin;
    }

    if (!found_leaf_node_with_1label) {
      // Handle the special case where the entire bitmap is 0.
      trailing_0label_cnt_ = 0;
    }
    counters_are_valid = true;
  }

  /// Initialize the counters if necessary.
  void
  ensure_counters_are_valid() noexcept {
    if (!counters_are_valid) {
      init_counters();
    }
  }

  /// Initialize the counters. This simplified version only works with unpruned
  /// trees.
  void
  init_counters_perfect_binary_tree() {
    inner_node_cnt_ = n_ - 1;
    leading_inner_node_cnt_ = inner_node_cnt_;
    trailing_leaf_node_cnt_ = inner_node_cnt_ + 1;
    explicit_node_idxs_.begin = inner_node_cnt_;
    explicit_node_idxs_.end = inner_node_cnt_;

    first_node_idx_with_1label_ = first_bit_idx_ + n_ - 1;
    last_node_idx_with_1label_ = last_bit_idx_ + n_ - 1;
    leading_0label_cnt_ = first_bit_idx_;
    trailing_0label_cnt_ = (first_bit_idx_ != n_)
        ? n_ - last_bit_idx_ - 1
        : 0;
    counters_are_valid = true;
  }

  /// A naive implementation of the tree-based compression algorithm.
  /// Prerequisites: The tree needs to be initialized, see init_tree(), and the
  /// tree MUST NOT be pruned.
  void __attribute__((noinline))
  compress() {
    assert(explicit_node_idxs_.begin == n_ - 1);
    assert(explicit_node_idxs_.end == n_ - 1);

    if (last_level() == 0) {
      // Nothing to compress.
      return;
    }

    // TODO remove special case
    if (leading_0label_cnt_ == inner_node_cnt_ + 1
        && trailing_0label_cnt_ == inner_node_cnt_ + 1) {
      trailing_0label_cnt_ = 0;
    }

    // Identifies the tree instance with the minimum size.
    $u64 min_size_tree_idx = first_node_idx_at_level(last_level()) - 1;
    // The estimated TEB size of the smallest tree instance found so far.
    $u64 min_size = estimate_encoded_size_in_bytes();

    // The tree-based compression works bottom-up in reverse level-order.
    for (auto level = last_level() - 1; level > 0; --level) {
      // The node indexes at the current level.
      u64 node_idx_begin = first_node_idx_at_level(level);
      u64 node_idx_end = first_node_idx_at_level(level + 1);

      // Count the number of nodes that were collapsed within the current level.
      $u64 collapse_cnt = 0;

      for ($u64 idx = node_idx_end - 1; idx != node_idx_begin - 1; --idx) {
        assert(is_inner_node(idx));
        u64 left_child = left_child_of(idx);
        u64 right_child = right_child_of(idx);

        if (is_leaf_node(left_child) && is_leaf_node(right_child)) {
          u1 left_label = label_of_node(left_child);
          u1 right_label = label_of_node(right_child);
          if (left_label == right_label) {
            const auto label = left_label;
            //===----------------------------------------------------------===//
            // Prune the two child nodes, i.e. collapse the current node.
            binary_tree_structure::set_leaf(idx );
            labels_.set(idx + offset, label);
            ++collapse_cnt;
            //===----------------------------------------------------------===//
            // Update counters.
            --inner_node_cnt_;

            assert(idx < explicit_node_idxs_.begin);
            explicit_node_idxs_.begin = idx;
            leading_inner_node_cnt_ = idx;
            if (idx + 1 == explicit_node_idxs_.end) {
              // The current node was the last explicit inner node.  Due to the
              // fact that pruning happens in reverse level-order, all
              // preceding nodes are inner nodes. Thus we can simply decrement
              // end of the range.
              --explicit_node_idxs_.end;
              --trailing_leaf_node_cnt_; // -2 + 1
            }
            else if (left_child >= explicit_node_idxs_.end) {
              // The two pruned nodes were trailing implicit nodes.
              trailing_leaf_node_cnt_ -= 2;
            }

            if (label == false) {
              // 'Moving' a 0-label one level up.  The label of the current node
              // is the very first label (in level order).
              ++leading_0label_cnt_;
              if (right_child < first_node_idx_with_1label_) {
                leading_0label_cnt_ -= 2;
              }
              if (left_child > last_node_idx_with_1label_) {
                trailing_0label_cnt_ -= 2;
              }
            }
            else {
              // 'Moving' a 1-label one level up.  The label of the current node
              // is the very first label (in level order).
              leading_0label_cnt_ = 0;
              first_node_idx_with_1label_ = idx;
              if (right_child == last_node_idx_with_1label_) {
                // The last node with a 1-label has been pruned. Thus, the new
                // last node with a 1-label is in the range [idx, left_child).
                for ($u64 j = 0; j < (left_child - idx); ++j) {
                  u64 i = left_child - 1 - j;
                  if (is_leaf_node(i)) {
                    if (label_of_node(i) == true) {
                      last_node_idx_with_1label_ = i;
                      break;
                    }
                    else {
                      ++trailing_0label_cnt_;
                    }
                  }
                }
              }
            }
            //===----------------------------------------------------------===//
            // Compute the size and memorize the smallest instance.
            const auto estimated_size = estimate_encoded_size_in_bytes();
            if (estimated_size < min_size) {
              min_size_tree_idx = idx;
              min_size = estimated_size;
            }
            //===----------------------------------------------------------===//
          }
        }
      }

      if (collapse_cnt < 2) {
        // Terminate pruning. As there will be nothing to prune in the next
        // higher level. Recall, in requires at least two (sibling) nodes.
        break;
      }
      if (level < level_of(min_size_tree_idx)) {
        break;
      }
    }
#ifndef NDEBUG
    validate_counters();
#endif

    if (optimization_level_ <= 1) return;

    // Restore the tree instance with the minimum size.
    expand_until(min_size_tree_idx);

#ifndef NDEBUG
    validate_counters();
#endif
  }

  /// Expands the tree nodes which have an index less than the given node index.
  void __attribute__((noinline))
  expand_until(u64 node_idx) {
    const auto it_end = breadth_first_end();
    for (auto it = breadth_first_begin(); it != it_end; ++it) {
      u64 idx = (*it);
      if (is_inner_node(idx)) continue;
      if (right_child_of(idx) >= max_node_cnt_) break;

      if (idx >= node_idx) break;
      // Expand the leaf node to an inner node.
      // DANGER: The tree structure is modified during iterating.
      set_inner<false>(idx);
    }
  }

  // For debugging purposes.
  void __forceinline__
  validate_counters() {
#ifndef NDEBUG
    //===------------------------------------------------------------------===//
    // Validation
    //===------------------------------------------------------------------===//
    // Back up the counters.
    const std::size_t _inner_node_cnt = inner_node_cnt_;
    const std::size_t _leading_inner_node_cnt = leading_inner_node_cnt_;
    const std::size_t _trailing_leaf_node_cnt = trailing_leaf_node_cnt_;
    const std::size_t _leading_0label_cnt = leading_0label_cnt_;
    const std::size_t _trailing_0label_cnt = trailing_0label_cnt_;
    // Back up the ranges with implicit nodes and labels.
    const std::size_t _explicit_node_idxs_begin = explicit_node_idxs_.begin;
    const std::size_t _explicit_node_idxs_end = explicit_node_idxs_.end;
    const std::size_t _first_node_idx_with_1label = first_node_idx_with_1label_;
    const std::size_t _last_node_idx_with_1label = last_node_idx_with_1label_;

    init_counters_naive();
    // TODO remove special case
    if (leading_0label_cnt_ == inner_node_cnt_ + 1
        && trailing_0label_cnt_ == inner_node_cnt_ + 1) {
      trailing_0label_cnt_ = 0;
    }
    assert(_inner_node_cnt == inner_node_cnt_);
    assert(_leading_inner_node_cnt == leading_inner_node_cnt_);
    assert(_trailing_leaf_node_cnt == trailing_leaf_node_cnt_);
    assert(_explicit_node_idxs_begin == explicit_node_idxs_.begin);
    assert(_explicit_node_idxs_end == explicit_node_idxs_.end);
    assert(_leading_0label_cnt == leading_0label_cnt_);
    assert(_trailing_0label_cnt == trailing_0label_cnt_);
    assert(_first_node_idx_with_1label == first_node_idx_with_1label_);
    assert(_last_node_idx_with_1label == last_node_idx_with_1label_);
#endif
  }

  /// Estimates the size in bytes, when the bitmap tree is succinctly encoded.
  /// This function basically resembles the size_in_bytes() function of TEBs.
  std::size_t __forceinline__
  estimate_encoded_size_in_bytes() {
    u64 explicit_tree_node_cnt = inner_node_cnt_ + inner_node_cnt_ + 1
        - leading_inner_node_cnt_ - trailing_leaf_node_cnt_;

    u64 explicit_label_cnt = optimization_level_ > 2
        ? inner_node_cnt_ + 1 - leading_0label_cnt_ - trailing_0label_cnt_
        : inner_node_cnt_ + 1;
    assert(leading_0label_cnt_ <= inner_node_cnt_ + 1);
    assert(trailing_0label_cnt_ <= inner_node_cnt_ + 1);

    u64 perfect_level_cnt = dtl::log_2(leading_inner_node_cnt_ + 1) + 1;

    return estimate_encoded_size_in_bytes_v2(
        explicit_tree_node_cnt, explicit_label_cnt, perfect_level_cnt);
  }

  std::size_t __forceinline__
  estimate_encoded_size_in_bytes_v2(
      std::size_t explicit_tree_node_cnt,
      std::size_t explicit_label_cnt,
      std::size_t perfect_level_cnt) {
    constexpr u64 block_bitlength = 64; // TODO remove magic number - refer to teb storage typse
    constexpr u64 block_size = block_bitlength / 8;
    $u64 bytes = 0;

    // Bit-length of the original bitmap.
    bytes += sizeof(n_);

    // The stored length of the tree structure.
    bytes += 4;
    // The number of implicit inner nodes.
    bytes += optimization_level_ > 0 ? 4 : 0;
    // The number of implicit leaf nodes can then be computed as
    //  2n-1 - # implicit nodes - length of the tree structure bit sequence
    // The offset to the beginning of T can also be computed.
    // The height of the encoded tree (after pruning).
    bytes += 1; // actually 5 bits

    // The stored length of L.
    bytes += 4;
    // The number of implicit labels.
    bytes += optimization_level_ > 2 ? 4 : 0;
    // The offset to the beginning of L can also be computed based on the
    // size of the header, T and R.

    // Level offsets for T and L, which are required by the tree scan algorithm.
    const auto encoded_tree_height = dtl::log_2(n_) + 1; // FIXME could be lower, but its unlikely
    assert(encoded_tree_height >= perfect_level_cnt);

    if (explicit_tree_node_cnt > 1024) {
      bytes += (4 + 4) * (encoded_tree_height - perfect_level_cnt);
    }

    // Padding. We want T to be 8-byte aligned.
    bytes += 8 - (bytes % 8);

    // Tree structure
    if (explicit_tree_node_cnt > 1) {
      bytes += ((explicit_tree_node_cnt + block_bitlength - 1) / block_bitlength)
          * block_size;
    }

    // Rank helper structure
    if (explicit_tree_node_cnt > 1024) {
      bytes += dtl::rank1_logic_surf<u64>::estimate_size_in_bytes(
          explicit_tree_node_cnt);
    }

    // Labels
    bytes += ((explicit_label_cnt + block_bitlength - 1) / block_bitlength)
        * block_size;

    return bytes;
  }

  /// Returns the maximum number of tree nodes. Which is 2n-1 for perfect
  /// full binary trees.
  inline std::size_t
  max_node_cnt() const {
    return max_node_cnt_;
  }

  /// Returns the label index of the given node.
  inline auto
  label_idx_of_node(u64 node_idx) const {
    return node_idx + offset;
  }

  /// Returns the label of the given node.
  inline u1
  label_of_node(u64 node_idx) const {
    return labels_[node_idx + offset];
  }

  /// Returns the number of nodes in the tree.
  inline u32
  get_node_cnt() const noexcept {
    return get_inner_node_cnt() + get_leaf_node_cnt();
  }

  /// Returns the number of inner nodes in the tree.
  inline u32
  get_inner_node_cnt() const noexcept {
    return inner_node_cnt_;
  }

  /// Returns the number of leaf nodes in the tree.
  inline u32
  get_leaf_node_cnt() const noexcept {
    return inner_node_cnt_ + 1;
  }

  /// Returns the number of perfect levels in the tree.
  inline u32
  get_perfect_level_cnt() {
    return dtl::log_2(leading_inner_node_cnt_ + 1) + 1;
  }

  /// Returns the height of the encoded tree.
  inline u32
  get_encoded_tree_height() {
    return dtl::log_2(n_) + 1; // FIXME could be lower, but its unlikely
  }

  /// Returns the number of leading inner nodes (in level order).
  inline u32
  get_leading_inner_node_cnt() const noexcept {
    assert(explicit_node_idxs_.begin == leading_inner_node_cnt_);
    return leading_inner_node_cnt_;
  }

  /// Returns the number of trailing leaf nodes (in level order).
  inline u32
  get_trailing_leaf_node_cnt() const noexcept {
    return trailing_leaf_node_cnt_;
  }

  /// Returns the node index of the first non-implicit node. - Note that the
  /// node index refers to the index within a perfect binary tree.
  inline u32
  get_first_explicit_node_idx() const noexcept {
    return explicit_node_idxs_.begin;
  }

  /// Returns the node index of the last non-implicit node. - Note that the
  /// node index refers to the index within a perfect binary tree.
  inline u32
  get_last_explicit_node_idx() const noexcept {
    return explicit_node_idxs_.end;
  }

  /// Returns the number of leading 0-labels (in level order).
  inline u32
  get_leading_0label_cnt() const noexcept {
    return leading_0label_cnt_;
  }

  /// Returns the number of trailing leaf nodes (in level order).
  inline u32
  get_trailing_0label_cnt() const noexcept {
    return trailing_0label_cnt_;
  }

  /// Returns the index of the first leaf node that carries a 1-label.
  inline u32
  get_first_node_idx_with_1label() const noexcept {
    return first_node_idx_with_1label_;
  }

  /// Returns the index of the last leaf node that carries a 1-label.
  inline u32
  get_last_node_idx_with_1label() const noexcept {
    return last_node_idx_with_1label_;
  }

  void
  print(std::ostream& os) const noexcept {
    std::stringstream l_os;
    for (std::size_t level = 0; level <= height_; ++level) {
      std::cout << std::setw(4) << level << ": ";
      u64 node_idx_from = (1ull << level) - 1;
      u64 node_idx_to = (1ull << (level + 1)) - 1;
      for ($u64 node_idx = node_idx_from; node_idx < node_idx_to; ++node_idx) {
        const auto spaces = (level < height_)
            ? (1u << (height_ - level)) - 1
            : 0;
        if (node_idx == root()) {
          std::cout << (is_leaf_node(node_idx) ? "0" : "1");
          if (is_leaf_node(node_idx)) {
            l_os << (label_of_node(node_idx) ? "1" : "0");
          }
        }
        else {
          std::cout << (is_leaf_node(node_idx)
                  ? is_leaf_node(parent_of(node_idx)) ? " " : "0"
                  : "1");
          if (is_leaf_node(node_idx) && !is_leaf_node(parent_of(node_idx))) {
            l_os << (label_of_node(node_idx) ? "1" : "0");
          }
        }
        for (std::size_t i = 0; i < spaces; ++i) {
          std::cout << " ";
        }
      }
      std::cout << std::endl;
    }
    std::cout
        << std::setw(4) << "L" << ": "
        << l_os.str() << std::endl;
  }

private:
  /// Expands the VERY FIRST explicit node. // TODO generalize to expand arbitrary nodes
  template<u1 fast>
  void
  set_inner(u64 idx) {
    assert(is_leaf_node(idx));
    if (fast) {
      // This invalidates the active node bitmap.
      binary_tree_structure::set_inner_fast(idx);
    }
    else {
      binary_tree_structure::set_inner(idx);
    }

    const auto left_child_idx = left_child_of(idx);
    const auto right_child_idx = right_child_of(idx);

    // Update the counters.
    ++inner_node_cnt_;
    // The number of leading inner nodes increases at least by one.
    ++leading_inner_node_cnt_;
    // Check if the newly created inner node is followed by more inner nodes.
    for (auto i = idx + 1; i < max_node_cnt_ && is_inner_node(i); ++i) { // TODO maybe implement a find_next_zero()
      ++leading_inner_node_cnt_;
    }
    // Update the index of the first explicit node.
    explicit_node_idxs_.begin = leading_inner_node_cnt_;
    // Update the number of trailing leaf nodes.
    if (unlikely(explicit_node_idxs_.begin >= explicit_node_idxs_.end)) {
      // The entire tree became implicit.
      leading_inner_node_cnt_ = inner_node_cnt_;
      trailing_leaf_node_cnt_ = inner_node_cnt_ + 1;
      explicit_node_idxs_.end = explicit_node_idxs_.begin;
    }
    else {
      if (left_child_idx >= explicit_node_idxs_.end) {
        // The newly created leaf is implicit.
        ++trailing_leaf_node_cnt_;
      }
      if (right_child_idx >= explicit_node_idxs_.end) {
        // The newly created leaf is implicit.
        ++trailing_leaf_node_cnt_;
      }
    }

    // Update the number of implicit 0-labels.
    if (idx < first_node_idx_with_1label_) {
      // The current node, that has been expanded had an implicit 0-label.
      --leading_0label_cnt_;

      if (left_child_idx < first_node_idx_with_1label_) {
        // The label of the newly created leaf is implicit.
        ++leading_0label_cnt_;
      }
      if (left_child_idx > last_node_idx_with_1label_) {
        // The label of the newly created leaf is implicit.
        ++trailing_0label_cnt_;
      }
      if (left_child_idx < first_node_idx_with_1label_) {
        // The label of the newly created leaf is implicit.
        ++leading_0label_cnt_;
      }
      if (right_child_idx > last_node_idx_with_1label_) {
        // The label of the newly created leaf is implicit.
        ++trailing_0label_cnt_;
      }
    }
    else if (idx == first_node_idx_with_1label_) {
      assert(label_of_node(idx) == true);
      // The current node is the first node with a 1-label.
      // Now that this node has been expanded, the node idx with the first
      // 1-label is somewhere in (idx, left-child(idx)].
      for (std::size_t i = idx + 1; i <= left_child_idx; ++i) {
        if (is_leaf_node(i)) {
          if (label_of_node(i) == true) {
            first_node_idx_with_1label_ = i;
            break;
          }
          else {
            ++leading_0label_cnt_;
            if (i > last_node_idx_with_1label_) {
              --trailing_0label_cnt_;
            }
          }
        }
      }
      assert(first_node_idx_with_1label_ != idx);

      // The node idx with the last 1-label is now the
      // max(last_node_idx_with_1label_, right-child(idx)).
      last_node_idx_with_1label_ =
          std::max(last_node_idx_with_1label_, right_child_idx);
      if (last_node_idx_with_1label_ == right_child_idx) {
        trailing_0label_cnt_ = 0;
      }
    }
    else {
      // The current node, that has been expanded had an implicit 0-label.
      // The node idx is greater than the last node idx with a 1-label,
      // therefore, the child also have larger indexes. Thus, only the trailing
      // 0-label cnt needs to be incremented.
      ++trailing_0label_cnt_;
    }
  }

  /// Restore the initial tree state.
  void __attribute__((noinline))
  uncompress() {
    is_inner_node_.set(0, first_node_idx_at_level(last_level()) + offset);
    is_active_node_.set(0, first_node_idx_at_level(last_level() + 1) + offset);
    init_counters_perfect_binary_tree();
#ifndef NDEBUG
    validate_counters();
#endif
  }

  void __attribute__((noinline))
  run_optimize() {
    // Optimization level 2.
    if (optimization_level_ > 1) {
      // Gradual decompression.
      // Find the tree instance with the minimum size.
      auto min_idx = explicit_node_idxs_.begin;
      auto min_level = level_of(min_idx);
      auto min_size = estimate_encoded_size_in_bytes();
      const std::size_t threshold = std::max(1024ul, uncompressed_size / 10);
      auto cpy = *this; // Work with a copy.

      {
        u64 start_idx = cpy.explicit_node_idxs_.begin;
        $u64 idx = cpy.explicit_node_idxs_.begin;
        $u64 prev_idx = 0;
        $u64 idx_level = cpy.level_of(idx);
        u64 last_pruning_level = last_level();
        $u64 old_level = idx_level;
        auto prev_level_size = cpy.estimate_encoded_size_in_bytes();
        auto prev_level_min_size = min_size;
        while (idx_level < last_pruning_level
            && !(min_level + 2 < idx_level)) {
          assert(cpy.is_leaf_node(idx));

          if (idx_level > old_level) {
            {
              auto now_size = cpy.estimate_encoded_size_in_bytes();
              if (now_size > prev_level_size && old_level != min_level) {
                goto expand_done;
              }
              else {
                prev_level_size = now_size;
              }
            }
            old_level = idx_level;
            if (idx_level < cpy.last_level()) {
              auto last_explicit_node = cpy.explicit_node_idxs_.end - 1;
              auto level_of_last_explicit_node = cpy.level_of(last_explicit_node);
              if (level_of_last_explicit_node == idx_level + 1) {
                auto parent_of_last_explicit_node =
                    cpy.parent_of(level_of_last_explicit_node);
                auto x_b = cpy.first_node_idx_at_level(idx_level);
                auto x_e = parent_of_last_explicit_node;
                auto y_b = parent_of_last_explicit_node + 1;
                auto y_e = cpy.first_node_idx_at_level(idx_level + 1);
                i64 x_ones = cpy.is_inner_node_.count(x_b + offset, x_e + offset);
                i64 x_zeros = (x_e - x_b) - x_ones;
                i64 y_ones = cpy.is_inner_node_.count(y_b + offset, y_e + offset);
                i64 y_zeros = (y_e - y_b) - y_ones;
                auto r = -(x_zeros + y_zeros) + 2 * x_zeros - (x_ones + y_ones);
                if (r > 0) {
                  goto expand_done;
                }
              }
              else {
                auto b = cpy.first_node_idx_at_level(idx_level);
                auto e = cpy.first_node_idx_at_level(idx_level + 1);
                auto ones = cpy.is_inner_node_.count(
                    b + offset,
                    e + offset);
                auto zeros = (e - b) - ones;
                if (zeros/2 > ones) {
                  goto expand_done;
                }
              }
            }
          }

          cpy.template set_inner<true>(idx);
          const auto compressed_size = cpy.estimate_encoded_size_in_bytes();
          // Estimates the size of a TEB.
          if (compressed_size <= min_size) { // less than or EQUAL because we prefer more balanced trees.
            min_idx = idx;
            min_level = cpy.level_of(idx);
            min_size = compressed_size;
          }
          if (compressed_size > min_size + threshold) {
            goto expand_done;
          }
          // Note:
          // If explicit_node_idxs_.begin overtakes parent_of(first_node_idx_with_1label_)
          // L starts to increase in size. And the next time where L could
          // decrease in size is when expanding first_node_idx_with_1label_.
          // All leaf nodes in [parent_of(first_node_idx_with_1label_) + 1, first_node_idx_with_1label_)
          // carry a 0-label. Thus, the L sequence starts with a 1 and until
          // first_node_idx_with_1label_ is reached, 0's are inserted after the
          // leading 1 in L.
          // Within that range, T can increase and decrease in size.
          //  with every expand, 1-bit is added to T, however there may be
          //  following inner nodes which then become implicit.
          if (cpy.explicit_node_idxs_.begin > cpy.explicit_node_idxs_.end) {
            goto tree_implicit;
          }

          // Next
          prev_idx = idx;
          idx = cpy.explicit_node_idxs_.begin;
          idx_level = cpy.level_of(idx);
        }

        // Check if decompression stopped at the minimum sized tree.
        if (min_idx == prev_idx) {
          // Take the working copy as is.
          *this = cpy;
          // Fix active node map.
          is_active_node_.set(start_idx + offset,
              std::min(max_node_cnt_, right_child_of(idx) + 1) + offset);
          D(validate_active_nodes();)
          return;
        }

        goto expand_done;

tree_implicit:
        // Take the working copy as is.
        *this = cpy;
        // Fix active node map.
        is_active_node_.set(start_idx + offset,
            std::min(max_node_cnt_, right_child_of(idx) + 1) + offset);
        D(validate_active_nodes();)
        return;

expand_done:

        // Check whether the bitmap can be compressed.
        if (min_size > uncompressed_size) {
          // Restore the initial tree state.
          uncompress();
          return;
        }
      }

      // Restore the tree instance with the minimum size. Expand the tree down
      // to the node 'min_idx' which identifies the smallest tree instance found
      // in the previous step.
      {
        const auto start_idx = explicit_node_idxs_.begin;
        const std::size_t expand_until_idx =
            std::min(first_node_idx_at_level(height_) - 1, min_idx);
        while (explicit_node_idxs_.begin <= expand_until_idx) {
          assert(is_leaf_node(explicit_node_idxs_.begin));
          set_inner<true>(explicit_node_idxs_.begin);
        }
        // Update the active node map.
        is_active_node_.set(start_idx + offset,
            std::min(max_node_cnt_, right_child_of(expand_until_idx) + 1) + offset);
        D(validate_active_nodes();)
      }
    }
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
