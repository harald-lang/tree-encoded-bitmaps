#pragma once
//===----------------------------------------------------------------------===//
#include <iomanip>
#include <immintrin.h>

#include <boost/dynamic_bitset.hpp>
#include <dtl/dtl.hpp>
#include <dtl/bitmap/util/binary_tree_structure.hpp>

#include <dtl/bitmap/util/rank1.hpp>
#include "plain_bitmap.hpp"
#include "config.hpp"
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// Represents a bitmap as a binary tree. During construction, the bitmap tree
/// is compressed, either loss-less or lossy.
/// The template parameter controls the space optimizations.
/// An instance can be seen as an intermediate representation of a bitmap from
/// which a TEB is constructed.
template<i32 optimization_level_ = 3>
class bitmap_tree : public binary_tree_structure {

  using tree_t = dtl::binary_tree_structure;
  using bitmap_t = boost::dynamic_bitset<$u32>;

  /// The labels of the tree nodes.
  plain_bitmap<$u64> labels_;

  struct range_t {
    std::size_t begin = 0;
    std::size_t end = 0;
  };

  //===--------------------------------------------------------------------===//
  // The following counters are required to compute the size of the resulting
  // tree-encoded bitmap. These counters are updated when the tree is modified.
  //
  // During the construction phase, the size is required to find the tree
  // instance with the smallest memory consumption. Recall that the number of
  // tree nodes does not directly correspond to the TEB size. A fully pruned
  // tree instance may result in a larger TEB than a partially pruned tree
  // instance.
  //===--------------------------------------------------------------------===//
  /// The number of inner nodes.
  std::size_t inner_node_cnt_;
  /// The number of leaf nodes.
  std::size_t leaf_node_cnt_;
  /// The number inner nodes until the first leaf occurs (in level order).
  std::size_t leading_inner_node_cnt_;
  /// The number leaf nodes after the last inner node (in level order).
  std::size_t trailing_leaf_node_cnt_;
  /// The node indexes that need to be stored explicitly.
  range_t explicit_node_idxs_;
  /// The number of 0-labels until the first 1-label occurs.
  std::size_t leading_0label_cnt_;
  /// The number of 0-labels after the last 1-label.
  std::size_t trailing_0label_cnt_;
  /// The first leaf node index that carries a 1-label.
  std::size_t first_node_idx_with_1label_;
  /// The last leaf node index that carries a 1-label.
  std::size_t last_node_idx_with_1label_;
  //===--------------------------------------------------------------------===//

public:

  /// C'tor
  explicit
  bitmap_tree(const bitmap_t& bitmap, f64 fpr = 0.0)
      : binary_tree_structure(bitmap.size()),
        labels_(max_node_cnt_ + offset),
        inner_node_cnt_(0),
        leaf_node_cnt_(0),
        leading_inner_node_cnt_(0),
        trailing_leaf_node_cnt_(0),
        leading_0label_cnt_(0),
        trailing_0label_cnt_(0),
        first_node_idx_with_1label_(0),
        last_node_idx_with_1label_(0) {

    // Init the binary tree and perform bottom-up pruning.
    init_tree(bitmap);
    prune_tree();

    // Run space optimizations.
    {
      // Optimization level 2
      if (optimization_level_ > 1) {
        init_counters();
        run_optimize();
      }
    }

    // Run lossy compression algorithm.
    if (fpr > 0.0) {
      compress_lossy(bitmap, fpr); // FIXME experimental
    }
  }

  bitmap_tree(const bitmap_tree& other) = default;
  bitmap_tree(bitmap_tree&& other) noexcept = default;
  bitmap_tree& operator=(const bitmap_tree& other) = default;
  bitmap_tree& operator=(bitmap_tree&& other) noexcept = default;
  ~bitmap_tree() override = default;

  /// Initialize a perfect binary tree on top of the given bitmap, where
  ///  - all the inner nodes have two children
  ///  - all the leaf nodes are on the same level
  ///  - the leaf nodes are labelled with the given bitmap
  /// Further, the labels for ALL inner nodes are computed, to avoid ad hoc
  /// computations of the labels during collapsing (or expanding) the tree
  /// structure.
  void
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
      while (i != boost::dynamic_bitset<$u32>::npos) {
        labels_.set(length / 2 + i + offset);
        i = bitmap.find_next(i);
      }
    }

    // Propagate the label bits along the tree (bottom-up).  The labels of an
    // internal node is the bitwise OR of the labels of both child nodes.
    {
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
  }

  /// Bottom-up pruning (loss-less).  Eliminate all sibling leaf nodes which
  /// have the same label. The algorithm terminates when all pairs of sibling
  /// leaf nodes have different labels.
  void __teb_inline__
  prune_tree() {
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
        if (remaining < 64) {
          u1 left_bit = labels_[src_node_idx + offset];
          u1 right_bit = labels_[src_node_idx + 1 + offset];
          u1 prune_causes_false_positives = left_bit ^ right_bit;
          u1 both_nodes_are_leaves =
              !is_inner_node(src_node_idx)
                  & !is_inner_node(src_node_idx + 1);
          u1 prune = both_nodes_are_leaves & !prune_causes_false_positives;
          if (prune) {
            binary_tree_structure::set_leaf(dst_node_idx);
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

          $u64* raw_label_ptr = labels_.data();
          auto src_labels = raw_label_ptr[src_idx / 64];
          auto prune_causes_false_positives = static_cast<$u32>(
              _pext_u64(src_labels ^ (src_labels >> 1), 0x5555555555555555)); // FIXME non-portable code

          $u64* raw_node_ptr = is_inner_node_.data();
          auto src_nodes = raw_node_ptr[src_idx / 64];
          auto both_nodes_are_leaves = static_cast<$u32>(
              _pext_u64(~src_nodes & (~src_nodes >> 1), 0x5555555555555555)); // FIXME non-portable code

          auto prune = both_nodes_are_leaves & ~prune_causes_false_positives;

          $u32* dst_nodes_ptr = &reinterpret_cast<$u32*>(raw_node_ptr)[dst_idx / 32];
          *dst_nodes_ptr = (*dst_nodes_ptr) ^ prune;

          src_node_idx += 64;
          dst_node_idx += 32;
        }
      }
    }
  }

  /// TODO
  void __attribute__((noinline))
  init_counters() {
    inner_node_cnt_ = 0;
    leaf_node_cnt_ = 0;
    leading_inner_node_cnt_ = 0;
    trailing_leaf_node_cnt_ = 0;
    leading_0label_cnt_ = 0;
    trailing_0label_cnt_ = 0;
    first_node_idx_with_1label_ = 0;
    last_node_idx_with_1label_ = 0;
    explicit_node_idxs_.begin = 0;
    explicit_node_idxs_.end = 0;

    // Determine the total number of tree nodes and which of these nodes need
    // to be stored explicitly.
    $u1 found_leaf_node = false;
    std::size_t node_cnt = 0;
    inner_node_cnt_ = 0;
    leaf_node_cnt_ = 0;
    leading_inner_node_cnt_ = 0;
    trailing_leaf_node_cnt_ = 0;
    explicit_node_idxs_.end = 0;

    $u1 found_leaf_node_with_1label = false;
    leading_0label_cnt_ = 0;
    trailing_0label_cnt_ = 0;
    first_node_idx_with_1label_ = 0;
    last_node_idx_with_1label_ = 0;


    const auto it_end = const_breadth_first_end();
    for (auto it = const_breadth_first_begin(); it != it_end; ++it) {
      u64 idx = (*it).idx;
      u64 level = (*it).level;
      u1 is_inner = (*it).is_inner;

      ++node_cnt;
      inner_node_cnt_ += is_inner;
      leaf_node_cnt_ += !is_inner;

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
      } else {
        trailing_leaf_node_cnt_ = 0;
        explicit_node_idxs_.end = idx;
      }

      u1 is_leaf = is_leaf_node(idx);

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
  }

  void
  init_counters_old() {
    inner_node_cnt_ = 0;
    leaf_node_cnt_ = 0;
    leading_inner_node_cnt_ = 0;
    trailing_leaf_node_cnt_ = 0;
    leading_0label_cnt_ = 0;
    trailing_0label_cnt_ = 0;
    first_node_idx_with_1label_ = 0;
    last_node_idx_with_1label_ = 0;
    explicit_node_idxs_.begin = 0;
    explicit_node_idxs_.end = 0;

    // Determine the total number of tree nodes and which of these nodes need
    // to be stored explicitly.
    {
      $u1 found_leaf_node = false;
      std::size_t node_cnt = 0;
      inner_node_cnt_ = 0;
      leaf_node_cnt_ = 0;
      leading_inner_node_cnt_ = 0;
      trailing_leaf_node_cnt_ = 0;
      explicit_node_idxs_.end = 0;

      const auto it_end = breadth_first_end();
      for (auto it = breadth_first_begin(); it != it_end; ++it) {
        u64 idx = (*it).idx;
        u64 level = (*it).level;
        u1 is_inner = is_inner_node(idx);

        ++node_cnt;
        inner_node_cnt_ += is_inner;
        leaf_node_cnt_ += !is_inner;

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
          explicit_node_idxs_.end = idx;
        }
      }

      explicit_node_idxs_.begin = leading_inner_node_cnt_;
    }

    // Determine the number of labels that need to be stored explicitly.
    {
      $u1 found_leaf_node_with_1label = false;
      leading_0label_cnt_ = 0;
      trailing_0label_cnt_ = 0;
      first_node_idx_with_1label_ = 0;
      last_node_idx_with_1label_ = 0;

      const auto it_end = breadth_first_end();
      for (auto it = breadth_first_begin(); it != it_end; ++it) {
        u64 idx = (*it).idx;
        u64 level = (*it).level;
        u1 is_leaf = is_leaf_node(idx);

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

    }
  }

  /// Estimates the size in bytes, when the bitmap tree is succinctly encoded.
  /// This function basically resembles the size_in_bytes() function of TEBs.
  std::size_t
  estimate_encoded_size_in_bytes() {
    constexpr u64 block_bitlength = 64;
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

    u64 explicit_tree_node_cnt = inner_node_cnt_ + leaf_node_cnt_
        - leading_inner_node_cnt_ - trailing_leaf_node_cnt_;

    // Level offsets for T and L, which are required by the tree scan algorithm.
    const auto perfect_levels = dtl::log_2(leading_inner_node_cnt_ + 1) + 1;
    const auto encoded_tree_height = dtl::log_2(n_) + 1; // FIXME could be lower, but its unlikely
    assert(encoded_tree_height >= perfect_levels);
    if (explicit_tree_node_cnt > 1024) {
      bytes += (4 + 4) * (encoded_tree_height - perfect_levels);
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
      bytes += dtl::rank1<u64>::estimate_size_in_bytes(explicit_tree_node_cnt);
    }

    // Labels
//    assert(leading_0label_cnt_ <= leaf_node_cnt_);
//    assert(trailing_0label_cnt_ <= leaf_node_cnt_);
    u64 explicit_label_cnt = optimization_level_ > 2
        ? leaf_node_cnt_ - leading_0label_cnt_ - trailing_0label_cnt_
        : leaf_node_cnt_;
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
    return inner_node_cnt_ + leaf_node_cnt_;
  }

  /// Returns the number of leading inner nodes (in level order).
  inline u32
  get_leading_inner_node_cnt() const noexcept {
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
        }
        else {
          std::cout << (is_leaf_node(node_idx)
              ? is_leaf_node(parent_of(node_idx)) ? " " : "0"
              : "1"
          );
        }
        for (std::size_t i = 0; i < spaces; ++i) {
          std::cout << " ";
        }
      }
      std::cout << std::endl;
    }
  }

private:

  /// Expands the first explicit node.
  void __teb_inline__
  set_inner(u64 idx) {
    assert(is_leaf_node(idx));
    binary_tree_structure::set_inner(idx);
    // Update the counters.
    ++inner_node_cnt_;
    ++leaf_node_cnt_;

    // The number of leading inner nodes increases at least by one.
    ++leading_inner_node_cnt_;
    // Check if the newly created inner node is followed by more inner nodes.
    for (auto i = idx + 1; i < max_node_cnt_ && is_inner_node(i); ++i) {
      ++leading_inner_node_cnt_;
    }
    // Update the index of the first explicit node.
    explicit_node_idxs_.begin = leading_inner_node_cnt_;
    // Update the number of trailing leaf nodes.
    if (explicit_node_idxs_.begin >= explicit_node_idxs_.end) {
      // The entire tree became implicit.
      leading_inner_node_cnt_ = inner_node_cnt_;
      trailing_leaf_node_cnt_ = leaf_node_cnt_;
    }
    else {
      const auto left_child_idx = left_child_of(idx);
      if (left_child_idx >= explicit_node_idxs_.end) {
        // The newly created leaf is implicit.
        trailing_leaf_node_cnt_ += 1;
      }
      const auto right_child_idx = right_child_of(idx);
      if (right_child_idx >= explicit_node_idxs_.end) {
        // The newly created leaf is implicit.
        trailing_leaf_node_cnt_ += 1;
      }
    }

    // Update the number of implicit 0-labels.
    if (idx < first_node_idx_with_1label_) {
      // The current node, that has been expanded had an implicit 0-label.
      leading_0label_cnt_ -= 1;

      const auto left_child_idx = left_child_of(idx);
      if (left_child_idx < first_node_idx_with_1label_) {
        // The label of the newly created leaf is implicit.
        leading_0label_cnt_ += 1;
      }
      if (left_child_idx > last_node_idx_with_1label_) {
        // The label of the newly created leaf is implicit.
        trailing_0label_cnt_ += 1;
      }
      const auto right_child_idx = right_child_of(idx);
      if (left_child_idx < first_node_idx_with_1label_) {
        // The label of the newly created leaf is implicit.
        leading_0label_cnt_ += 1;
      }
      if (right_child_idx > last_node_idx_with_1label_) {
        // The label of the newly created leaf is implicit.
        trailing_0label_cnt_ += 1;
      }
    }

    if (idx == first_node_idx_with_1label_) {
      // The current node is the first node with a 1-label.
      // Now that this node has been expanded, the node idx with the first
      // 1-label is somewhere in (idx, left-child(idx)].
      for (std::size_t i = idx + 1; i <= left_child_of(idx); ++i) {
        if (is_leaf_node(i)) {
          if (label_of_node(i) == true) {
            first_node_idx_with_1label_ = i;
            break;
          }
          else {
            leading_0label_cnt_ += 1;
            if (i > last_node_idx_with_1label_) {
              trailing_0label_cnt_ -= 1;
            }
          }
        }
      }
      // The node idx with the last 1-label is now the
      // max(last_node_idx_with_1label_, right-child(idx)).
      if (right_child_of(idx) > last_node_idx_with_1label_) {
        // Update the number of trailing 0-labels.
//        for (std::size_t i = first_node_idx_with_1label_ + 1;
//             i < left_child_of(idx); ++i) {
//          if (is_leaf_node(i) && label_of_node(i) == false) {
//            trailing_0label_cnt_ -= 1;
////            std::cout << "t0lc=" << trailing_0label_cnt_ << std::endl;
//          }
//        }
        init_counters(); // FIXME performance issue
      }
      last_node_idx_with_1label_ =
          std::max(last_node_idx_with_1label_, right_child_of(idx));
    }
//    init_counters();
  }

  /// Collapse an inner node.
  void
  set_leaf(u64 idx) {
    if (label_of_node(idx) == true) {
      // Set the labels to 1 in the entire sub-tree.
      std::function<void(u64)> set_labels_to_1 = [&](u64 node_idx) { // FIXME highly inefficient
        if (node_idx >= max_node_cnt_) return;
        labels_[node_idx + offset] = true;
        set_labels_to_1(left_child_of(node_idx));
        set_labels_to_1(right_child_of(node_idx));
      };
      set_labels_to_1(idx);
    }
    binary_tree_structure::set_leaf(idx);
    init_counters(); // FIXME highly inefficient
  }
  void
  set_leaf_fast(u64 idx) {
    if (label_of_node(idx) == true) {
      // Set the labels to 1 in the entire sub-tree.
      std::function<void(u64)> set_labels_to_1 = [&](u64 node_idx) { // FIXME highly inefficient
        if (node_idx >= max_node_cnt_) return;
        labels_.set(node_idx + offset);
        set_labels_to_1(left_child_of(node_idx));
        set_labels_to_1(right_child_of(node_idx));
      };
      set_labels_to_1(idx);
    }
    binary_tree_structure::set_leaf(idx);
  }

  void __attribute__ ((noinline))
  run_optimize() {
    // Optimization level 2.
    if (optimization_level_ > 1) {
      // Estimates the size of a TEB.
      auto size = [&]() {
        return estimate_encoded_size_in_bytes();
      };

      // Gradual decompression.
      auto org_state = *this; // Memorize the fully pruned tree instance.
      // Find the tree instance with the minimum size.
      auto min_idx = root(); // Pruned all the way up to the root node.
      auto min_size = size();
      const auto it_end = breadth_first_end();
      for (auto it = breadth_first_begin(); it != it_end; ++it) {
        u64 idx = (*it).idx;
        if ((*it).is_inner) continue;
        if (right_child_of(idx) >= max_node_cnt_) break;

        // Expand the leaf node to an inner node.
        // DANGER: The tree structure is modified during iterating.
        set_inner(idx);

        const auto compressed_size = size();
        if (compressed_size <= min_size) {
          min_idx = idx;
          min_size = compressed_size;
        }
        if (min_size <= 64) break;
        if (explicit_node_idxs_.begin > explicit_node_idxs_.end) break;
        // FIXME when to stop decompression??? -------------------------------------------------------------
      }
      // Restore the tree instance with the minimum size.
      *this = org_state;
      // Expand the tree down to the node 'min_idx' which identifies the
      // smallest tree instance found in the previous step.
      for (auto it = breadth_first_begin(); it != it_end; ++it) {
        u64 idx = (*it).idx;
        if ((*it).is_inner) continue;
        if (right_child_of(idx) >= max_node_cnt_) break;

        // Expand the leaf node to an inner node.
        // DANGER: The tree structure is modified during iterating.
        set_inner(idx);

        if (idx >= min_idx) break;
        if (min_size <= 64) break;
        if (explicit_node_idxs_.begin > explicit_node_idxs_.end) break;
      }
    }
  }

  //===--------------------------------------------------------------------===//
  // Lossy compression.
  // EXPERIMENTAL CODE
  //===--------------------------------------------------------------------===//
  explicit
  bitmap_tree(const bitmap_t& bitmap, f64 fpr, f64 threshold)
      : binary_tree_structure(bitmap.size()),
        labels_(max_node_cnt_+ offset),
        inner_node_cnt_(0),
        leaf_node_cnt_(0),
        leading_inner_node_cnt_(0),
        trailing_leaf_node_cnt_(0),
        leading_0label_cnt_(0),
        trailing_0label_cnt_(0),
        first_node_idx_with_1label_(0),
        last_node_idx_with_1label_(0) {

    // Init the binary tree and perform bottom-up pruning.
    init_tree(bitmap);
    prune_tree();

    compress_lossy_pre_optimization_v2(fpr, threshold);
    init_counters();
    if (optimization_level_ > 1) {
      run_optimize();
    }
  }

  void
  compress_lossy(const bitmap_t& bitmap, f64 fpr) {
    if (fpr <= 0.0) return;

    const auto lossless_size = estimate_encoded_size_in_bytes();

    auto min_size = lossless_size;
    auto min = *this;
    auto min_threshold = 0.0;

//      $u1 monotonic = true;
//      $u1 first_lossy_pass = true;
//      std::size_t prev_lossy_size = std::numeric_limits<std::size_t>::max();
//      std::vector<std::size_t> lossy_sizes;

    for (auto threshold : {1.0/8.0, 1.0/4.0, 1.0/2.0, 1.0}) {
      bitmap_tree lossy_bt(bitmap, fpr, threshold);
      const auto size = lossy_bt.estimate_encoded_size_in_bytes();
      if (size < min_size) {
        min = lossy_bt;
        min_size = size;
        min_threshold = threshold;
      }
//        if (size <= prev_lossy_size) {
//          prev_lossy_size = size;
//        }
//        else {
//          break;
//        }
    }
    *this = min;

//      if (lossless_size > min_size) {
//
//        std::cout << ">>> saved bytes through lossy compression: "
//            << (lossless_size - min_size)
//            << " threshold = " << min_threshold
//            << std::endl;
//      }
//      else {
//        std::cout << ">>> lossless" << std::endl;
//      }
  }

  // Lossy compression.  The size of the tree structure is further reduced,
  // which causes false positive bits. The number of false positive bits is
  // limited by the given false positive rate (FPR).
  void __attribute__ ((noinline))
  compress_lossy_v2(f64 fpr) {
    // Determine maximum number of false positives.
    u64 max_fp_cnt = static_cast<u64>(n_ * fpr);
    if (max_fp_cnt == 0) return;

    // The total number of false positives introduced so far.
    $u64 total_fp_cntr = 0;

    // Pre-compute the number of associated false positives for each tree node.
    // The value fp_cntrs[i] refers to the number of false positives when
    // node i is turned into a leaf node (i.e., its sub-tree is pruned).
    std::vector<uint32_t> fp_cntrs(max_node_cnt_, 0);
    u64 length = max_node_cnt_;
    u64 height = height_;
    for ($u64 i = 0; i < length - 1; i += 2) {
      u64 left_node_idx = length - i - 2;
      u64 right_node_idx = left_node_idx + 1;

      u1 left_bit = labels_[left_node_idx + offset];
      u1 right_bit = labels_[right_node_idx + offset];

      u64 parent_node_idx = tree_t::parent_of(left_node_idx);

      const auto left_fp_cnt = fp_cntrs[left_node_idx];
      const auto right_fp_cnt = fp_cntrs[right_node_idx];
      fp_cntrs[parent_node_idx] = left_fp_cnt + right_fp_cnt;

      if (left_bit ^ right_bit) {
        // Introduce new false positives. - Note, that no new false positives
        // are introduced if the labels of the two child nodes are identical.
        u32 f = (1u << (height - tree_t::level_of(left_node_idx)));
        fp_cntrs[parent_node_idx] += f;
      }
    }

    // Pre-compute the sub-tree sizes and the number of leaf node in each
    // sub-tree.
    std::vector<uint32_t> subtree_sizes(max_node_cnt_, 0);
    std::vector<uint32_t> leaf_node_cnts(max_node_cnt_, 0);
    for ($u64 i = 0; i < length; ++i) {
      u64 node_idx = length - i - 1;
      if (is_leaf_node(node_idx)) {
        subtree_sizes[node_idx] = 1;
        leaf_node_cnts[node_idx] = 1;
      }
      else {
        const auto left_child_idx = left_child_of(node_idx);
        const auto right_child_idx = right_child_of(node_idx);
        subtree_sizes[node_idx] = subtree_sizes[left_child_idx]
            + subtree_sizes[right_child_idx] + 1;
        leaf_node_cnts[node_idx] = leaf_node_cnts[left_child_idx]
            + leaf_node_cnts[right_child_idx];
      }
    }


    {
      // Top-down pruning loop (starts at level 0).
      const auto perfect_levels = dtl::log_2(leading_inner_node_cnt_ + 1) + 1;
      std::size_t level = 0; //perfect_levels;
      auto min = *this;
      auto min_size = estimate_encoded_size_in_bytes();

      while (level <= height) {
        if (estimate_encoded_size_in_bytes() < 64) break;
        // Look for largest sub-tree to prune, so that the FPR is not exceeded.
        $u64 candidate_node_idx = 0;
        $u64 candidate_saved_byte_cnt = 0;
        auto candidate = *this;

        u64 node_idx_from = (1ull << level) - 1;
        u64 node_idx_to = (1ull << (level + 1)) - 1;

        for ($u64 node_idx = node_idx_from; node_idx < node_idx_to; ++node_idx) {
          if (is_leaf_node(node_idx)) {
            continue;
          }
          // Check whether pruning would not exceed the specified FPR.
          if (fp_cntrs[node_idx] + total_fp_cntr > max_fp_cnt) {
            continue;
          }

//          std::cout << "." << std::flush;
//          // Count the number of nodes that would be eliminated from the tree
//          // if the current node is turned into a leaf node.
////          assert(this->subtree_size(node_idx) == subtree_sizes[node_idx]);
////          const auto subtree_size = subtree_sizes[node_idx];
//          const auto subtree_size = this->subtree_size(node_idx);
//          const auto pruned_node_cnt = subtree_size - 1;
//          // Count the number of leaf nodes that would be eliminated from the
//          // tree.
////          assert(count_leaf_nodes(node_idx) == leaf_node_cnts[node_idx]);
////          const auto leaf_node_cnt = leaf_node_cnts[node_idx];
//          const auto leaf_node_cnt = count_leaf_nodes(node_idx);
//          const auto pruned_leaf_node_cnt = leaf_node_cnt - 1;
          // Estimate the amount of bits saved.
          // FIXME highly inefficient, replace by an heuristic
          auto cpy = *this;
          cpy.set_leaf(node_idx);
          cpy.run_optimize();
          const auto this_size = estimate_encoded_size_in_bytes();
          const auto cpy_size = cpy.estimate_encoded_size_in_bytes();
          const auto saved_byte_cnt = (cpy_size < this_size)
              ? this_size - cpy_size
              : std::size_t(0);

          if (saved_byte_cnt > candidate_saved_byte_cnt) {
            candidate_node_idx = node_idx;
            candidate_saved_byte_cnt = saved_byte_cnt;
            candidate = cpy;
          }
        }

        // Check if we have found a candidate to prune.
        if (candidate_saved_byte_cnt > 0) {
          *this = candidate;
//          // The actual pruning.
//          // Make the current node a leaf node. Note that we need to update
//          // the labels in the entire subtree otherwise we get side effects
//          // during gradual decompression.
//          set_leaf(candidate_node_idx);
////          run_optimize();
          total_fp_cntr += fp_cntrs[candidate_node_idx];
          std::cout << "lossy pruning. saving = " << candidate_saved_byte_cnt
              << ", remaining budget = "
              << (max_fp_cnt - total_fp_cntr)
              << "/" << max_fp_cnt
              << std::endl;
        }
        else {
          // No pruning candidate found at the current level.  Continue with the
          // next level.
          level++;
//          std::cout << "next_level" <<std::endl;
        }

      }
    }
  }

  // Lossy compression.  The size of the tree structure is further reduced,
  // which causes false positive bits. The number of false positive bits is
  // limited by the given false positive rate (FPR).
  void __attribute__ ((noinline))
  compress_lossy_pre_optimization(f64 fpr) {
    // Determine maximum number of false positives.
    u64 max_fp_cnt = static_cast<u64>(n_ * fpr);
    if (max_fp_cnt == 0) return;

    // The total number of false positives introduced so far.
    $u64 total_fp_cntr = 0;

    // Pre-compute the number of associated false positives for each tree node.
    // The value fp_cntrs[i] refers to the number of false positives when
    // node i is turned into a leaf node (i.e., its sub-tree is pruned).
    std::vector<uint32_t> fp_cntrs(max_node_cnt_, 0);
    u64 length = max_node_cnt_;
    u64 height = height_;
    for ($u64 i = 0; i < length - 1; i += 2) {
      u64 left_node_idx = length - i - 2;
      u64 right_node_idx = left_node_idx + 1;

      u1 left_bit = labels_[left_node_idx + offset];
      u1 right_bit = labels_[right_node_idx + offset];

      u64 parent_node_idx = tree_t::parent_of(left_node_idx);

      const auto left_fp_cnt = fp_cntrs[left_node_idx];
      const auto right_fp_cnt = fp_cntrs[right_node_idx];
      fp_cntrs[parent_node_idx] = left_fp_cnt + right_fp_cnt;

      if (left_bit ^ right_bit) {
        // Introduce new false positives. - Note, that no new false positives
        // are introduced if the labels of the two child nodes are identical.
        u32 f = (1u << (height - tree_t::level_of(left_node_idx)));
        fp_cntrs[parent_node_idx] += f;
      }
    }

    // Pre-compute the sub-tree sizes and the number of leaf node in each
    // sub-tree.
    std::vector<uint32_t> subtree_sizes(max_node_cnt_, 0);
    std::vector<uint32_t> leaf_node_cnts(max_node_cnt_, 0);
    for ($u64 i = 0; i < length; ++i) {
      u64 node_idx = length - i - 1;
      if (is_leaf_node(node_idx)) {
        subtree_sizes[node_idx] = 1;
        leaf_node_cnts[node_idx] = 1;
      }
      else {
        const auto left_child_idx = left_child_of(node_idx);
        const auto right_child_idx = right_child_of(node_idx);
        subtree_sizes[node_idx] = subtree_sizes[left_child_idx]
            + subtree_sizes[right_child_idx] + 1;
        leaf_node_cnts[node_idx] = leaf_node_cnts[left_child_idx]
            + leaf_node_cnts[right_child_idx];
      }
    }


    {
      // Top-down pruning loop (starts at level 0).
      std::size_t level = 0;
      while (level <= height) {
        // Look for largest sub-tree to prune, so that the FPR is not exceeded.
        $u64 candidate_node_idx = 0;
        $f64 candidate_saved_bit_cnt = 0.0;

        u64 node_idx_from = (1ull << level) - 1;
        u64 node_idx_to = (1ull << (level + 1)) - 1;

        for ($u64 node_idx = node_idx_from; node_idx < node_idx_to; ++node_idx) {
          if (is_leaf_node(node_idx)) {
            continue;
          }
          // Check whether pruning would not exceed the specified FPR.
          if (fp_cntrs[node_idx] + total_fp_cntr > max_fp_cnt) {
            continue;
          }

          // Count the number of nodes that would be eliminated from the tree
          // if the current node is turned into a leaf node.
//          assert(this->subtree_size(node_idx) == subtree_sizes[node_idx]);
//          const auto subtree_size = subtree_sizes[node_idx];
          const auto subtree_size = this->subtree_size(node_idx);
          const auto pruned_node_cnt = subtree_size - 1;
          // Count the number of leaf nodes that would be eliminated from the
          // tree.
//          assert(count_leaf_nodes(node_idx) == leaf_node_cnts[node_idx]);
//          const auto leaf_node_cnt = leaf_node_cnts[node_idx];
          const auto leaf_node_cnt = count_leaf_nodes(node_idx);
          const auto pruned_leaf_node_cnt = leaf_node_cnt - 1;
          // Compute the amount of bits saved.  One bit per node and one per
          // label.
          const auto saved_bit_cnt = pruned_node_cnt * 1.0625
              + pruned_leaf_node_cnt;

          if (saved_bit_cnt > candidate_saved_bit_cnt) {
            candidate_node_idx = node_idx;
            candidate_saved_bit_cnt = saved_bit_cnt;
          }
        }

        // Check if we have found a candidate to prune.
        if (candidate_saved_bit_cnt > 0.0) {
          // The actual pruning.
          // Set the labels to 1 in the entire sub-tree.
          visit(candidate_node_idx, [&](u64 i) {
            labels_[i + offset] = true;
          });
          // Make the current node a leaf node. Note that we need to update
          // the labels in the entire subtree otherwise we get side effects
          // during gradual decompression.
          set_leaf(candidate_node_idx);
          total_fp_cntr += fp_cntrs[candidate_node_idx];
        }
        else {
          // No pruning candidate found at the current level.  Continue with the
          // next level.
          level++;
        }

      }
    }
//
//    while (true) {
//      // Look for largest sub-tree to prune, so that the FPR is not exceeded.
//      $u64 candidate_node_idx = 0;
//      $f64 candidate_saved_bit_cnt = 0.0;
//      for (std::size_t level = 0; level < height_; ++level) {
//        u64 node_idx_from = (1ull << level) - 1;
//        u64 node_idx_to = (1ull << (level + 1)) - 1;
//        for ($u64 node_idx = node_idx_from; node_idx < node_idx_to; ++node_idx) {
//          if (is_leaf_node(node_idx)) {
//            continue;
//          }
//          // Check whether pruning would not exceed the specified FPR.
//          if (fp_cntrs[node_idx] + total_fp_cntr > max_fp_cnt) {
//            continue;
//          }
//
//          // Count the number of nodes that would be eliminated from the tree
//          // if the current node is turned into a leaf node.
//          const auto subtree_size = this->subtree_size(node_idx);
//          const auto pruned_node_cnt = subtree_size - 1;
//          // Count the number of leaf nodes that would be eliminated from the
//          // tree.
//          const auto leaf_node_cnt = count_leaf_nodes(node_idx);
//          const auto pruned_leaf_node_cnt = leaf_node_cnt - 1;
//          // Compute the amount of bits saved.  One bit per node and one per
//          // label.
//          const auto saved_bit_cnt = pruned_node_cnt * 1.0625 + pruned_leaf_node_cnt;
//
//          if (saved_bit_cnt > candidate_saved_bit_cnt) {
//            candidate_node_idx = node_idx;
//            candidate_saved_bit_cnt = saved_bit_cnt;
//          }
//        }
//        // Check, whether we found a candidate on the current tree level.
//        if (candidate_saved_bit_cnt > 0.0) break;
//      }
//
//      if (candidate_saved_bit_cnt > 0.0) {
//        // The actual pruning.
//        set_leaf(candidate_node_idx);
//        total_fp_cntr += fp_cntrs[candidate_node_idx];
//      }
//      else {
//        // Exit the pruning loop, because no candidate has been found.
//        break;
//      }
//    }
  }

  /// Lossy compression.  The size of the tree structure is further reduced,
  /// which causes false positive bits. The number of false positive bits is
  /// limited by the given false positive rate (FPR).
  void __attribute__ ((noinline))
  compress_lossy_pre_optimization_v2(f64 fpr, f64 threshold = 1.0/8.0) {
    // Determine maximum number of false positives.
    u64 max_fp_cnt = static_cast<u64>(n_ * fpr);
    if (max_fp_cnt == 0) return;

    // The total number of false positives introduced so far.
    $u64 total_fp_cntr = 0;

    std::vector<uint32_t> fp_cntrs(max_node_cnt_, 0);
    std::vector<uint32_t> intr_fp_cntrs(max_node_cnt_, 0);
    u64 length = max_node_cnt_;
    u64 height = height_;
    for ($u64 i = 0; i < length - 1; i += 2) {
      u64 left_node_idx = length - i - 2;
      u64 right_node_idx = left_node_idx + 1;

      u1 left_bit = labels_[left_node_idx + offset];
      u1 right_bit = labels_[right_node_idx + offset];

      u64 parent_node_idx = tree_t::parent_of(left_node_idx);
      binary_tree_structure::set_inner(parent_node_idx);

      const auto left_fp_cnt = fp_cntrs[left_node_idx];
      const auto right_fp_cnt = fp_cntrs[right_node_idx];
      fp_cntrs[parent_node_idx] = left_fp_cnt + right_fp_cnt;
      intr_fp_cntrs[parent_node_idx] =
          intr_fp_cntrs[left_node_idx] + intr_fp_cntrs[right_node_idx];

      u1 prune_introduces_false_positives = left_bit ^ right_bit;

      if (prune_introduces_false_positives) {
        // Introduce new false positives. - Note, that no new false positives
        // are introduced if the labels of the two child nodes are identical.
        u32 f = (1u << (height - tree_t::level_of(left_node_idx)));
        fp_cntrs[parent_node_idx] += f;
      }

      u1 both_nodes_are_leaves =
          is_leaf_node(left_node_idx) & is_leaf_node(right_node_idx);

      if (both_nodes_are_leaves && !prune_introduces_false_positives) {
        // Simple pruning step (loss less).
        binary_tree_structure::set_leaf(parent_node_idx);
      }
      else {
        // Possibly do a lossy pruning step.
        const auto span = n_ >> level_of(parent_node_idx);
        const auto fp_cnt = fp_cntrs[parent_node_idx];
        const auto intr_fp_cnt = intr_fp_cntrs[parent_node_idx];
        assert(fp_cnt >= intr_fp_cnt);

        const auto fp_weight = (fp_cnt * 1.0) / span;
        if ((fp_cnt - intr_fp_cnt + total_fp_cntr) <= max_fp_cnt
            && (fp_weight <= (threshold+0.00001))) {
//          binary_tree_structure::set_leaf(parent_node_idx);
//          labels_[parent_node_idx] = true;
          set_leaf_fast(parent_node_idx);
          const auto fp_cntr_old = total_fp_cntr;
          total_fp_cntr -= intr_fp_cntrs[parent_node_idx];
          intr_fp_cntrs[parent_node_idx] = fp_cntrs[parent_node_idx];
          total_fp_cntr += intr_fp_cntrs[parent_node_idx];

//          std::cout << "collapsing node " << parent_node_idx
//              << " with weight " << fp_weight << "."
//              << " introducing " << (total_fp_cntr - fp_cntr_old)
//              << " false positives. total = " << total_fp_cntr << " (max = "
//              << max_fp_cnt << ")"
//              << std::endl;
        }
      }
    }
  }
  //===--------------------------------------------------------------------===//

};
//===----------------------------------------------------------------------===//
} // namespace dtl