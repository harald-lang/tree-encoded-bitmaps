#pragma once

#include <iomanip>

#include <boost/dynamic_bitset.hpp>

#include <dtl/dtl.hpp>
#include <dtl/bitmap/util/binary_tree_structure.hpp>
#include <dtl/bitmap/util/rank1.hpp>

namespace dtl {

//===----------------------------------------------------------------------===//
/// Represents a bitmap as a binary tree. During construction, the bitmap tree
/// is compressed, either loss-less or lossy.
/// The template parameter controls the space optimizations.
template<i32 optimization_level_ = 2>
class bitmap_tree : public binary_tree_structure {

  using tree_t = dtl::binary_tree_structure;
  using bitmap_t = boost::dynamic_bitset<$u32>;

  /// The labels of the tree nodes.
  bitmap_t labels_;

  struct range_t {
    std::size_t begin = 0;
    std::size_t end = 0;
  };

  std::size_t inner_node_cnt_;
  std::size_t leaf_node_cnt_;
  std::size_t leading_inner_node_cnt_;
  std::size_t trailing_leaf_node_cnt_;
  range_t explicit_node_idxs_;

  std::size_t leading_0label_cnt_;
  std::size_t trailing_0label_cnt_;
  std::size_t first_node_idx_with_1label_;
  std::size_t last_node_idx_with_1label_;


public:

  /// C'tor
  explicit
  bitmap_tree(const bitmap_t& bitmap, f64 fpr = 0.0)
      : binary_tree_structure(bitmap.size()),
        labels_(max_node_cnt_),
        inner_node_cnt_(0),
        leaf_node_cnt_(0),
        leading_inner_node_cnt_(0),
        trailing_leaf_node_cnt_(0),
        leading_0label_cnt_(0),
        trailing_0label_cnt_(0),
        first_node_idx_with_1label_(0),
        last_node_idx_with_1label_(0) {

    // TODO: Support arbitrarily sized bitmaps.
    if (!dtl::is_power_of_two(n_)) {
      throw std::invalid_argument(
          "The length of the bitmap must be a power of two.");
    }

    // Initialize a perfect binary tree on top of the given bitmap.
    // ... all the inner nodes have two children
    // ... all the leaf nodes are on the same level
    // ... the leaf nodes are labelled with the given bitmap
    u64 length = max_node_cnt_;
    u64 height = height_;
    for ($u64 i = length / 2; i < length; i++) {
      labels_[i] = bitmap[i - length / 2];
    }

    // Propagate the bits along the tree (bottom-up).  The labels of an internal
    // node is the bitwise OR of the labels of both child nodes.
    for ($u64 i = 0; i < length - 1; i++) {
      u64 node_idx = length - i - 1;
      labels_[tree_t::parent_of(node_idx)] =
          labels_[tree_t::parent_of(node_idx)] | labels_[node_idx];
    }

    // Bottom-up pruning (loss-less).  Eliminate all sibling leaf nodes which
    // have the same label.
    for ($u64 i = 0; i < length - 1; i += 2) {
      u64 left_node_idx = length - i - 2;
      u64 right_node_idx = left_node_idx + 1;

      u1 left_bit = labels_[left_node_idx];
      u1 right_bit = labels_[right_node_idx];

      u64 parent_node_idx = tree_t::parent_of(left_node_idx);

      u1 prune_causes_false_positives = left_bit ^ right_bit;
      u1 both_nodes_are_leaves =
          !is_inner_node(left_node_idx)
              & !is_inner_node(right_node_idx);
      u1 prune = both_nodes_are_leaves & !prune_causes_false_positives;
      if (prune) {
        set_leaf(parent_node_idx);
      }
    }

    // Lossy compression.
    compress_lossy(fpr);

    // Determine the total number of tree nodes and which of these nodes need
    // to be stored explicitly.
    {
      $u1 found_leaf_node = false;
      std::size_t node_cnt = 0;
      explicit_node_idxs_.end = 0;

      for (auto it = breadth_first_begin(); it != breadth_first_end(); ++it) {
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
        } else {
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

      for (auto it = breadth_first_begin(); it != breadth_first_end(); ++it) {
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

    // Run space optimizations.
    {
      // Optimization level 2
      if (optimization_level_ > 1) {
        run_optimize();
      }
    }
  }

  bitmap_tree(const bitmap_tree& other) = default;
  bitmap_tree(bitmap_tree&& other) noexcept = default;
  bitmap_tree& operator=(const bitmap_tree& other) = default;
  bitmap_tree& operator=(bitmap_tree&& other) noexcept = default;
  ~bitmap_tree() override = default;

  /// Estimates the size in bytes, when the bitmap tree is succinctly encoded.
  /// This function basically resembles the size_in_bytes() function of TEBs.
  std::size_t
  estimate_encoded_size_in_bytes() {
    constexpr u64 block_bitlength = 64;
    constexpr u64 block_size = block_bitlength / 8;
    $u64 bytes = 0;

    // Bit-length of the original bitmap.
    bytes += sizeof(n_);

    // Tree structure
    u64 explicit_tree_node_cnt = inner_node_cnt_ + leaf_node_cnt_
        - leading_inner_node_cnt_ - trailing_leaf_node_cnt_;
    bytes += ((explicit_tree_node_cnt + block_bitlength - 1) / block_bitlength)
        * block_size;
    // The stored length of the tree structure.
    bytes += 4;
    // The number of implicit inner nodes.
    bytes += optimization_level_ > 0 ? 4 : 0;
    // The number of implicit leaf nodes can then be computed as
    //  2n-1 - # implicit nodes - length of the tree structure bit sequence
    // The offset to the beginning of T can also be computed.
    // The height of the encoded tree (after pruning).
    bytes += 1; // actually 5 bits

    // Rank helper structure
    bytes += dtl::rank1<u64>::estimate_size_in_bytes(explicit_tree_node_cnt);

    // Labels
    u64 explicit_label_cnt = optimization_level_ > 2
        ? leaf_node_cnt_ - leading_0label_cnt_ - trailing_0label_cnt_
        : leaf_node_cnt_;
    bytes += ((explicit_label_cnt + block_bitlength - 1) / block_bitlength)
        * block_size;
    // The stored length of L.
    bytes += 4;
    // The number of implicit labels.
    bytes += optimization_level_ > 2 ? 4 : 0;
    // The offset to the beginning of L can also be computed based on the
    // size of the header, T and R.

    // Level offsets for T and L, which are required by the tree scan algorithm.
    const auto perfect_levels = dtl::log_2(leading_inner_node_cnt_ + 1) + 1;
    const auto encoded_tree_height = dtl::log_2(n_) + 1; // FIXME could be lower, but its unlikely
    assert(encoded_tree_height >= perfect_levels);
    bytes += (4 + 4) * (encoded_tree_height - perfect_levels);

    // Padding. We want T to be 8-byte aligned.
    bytes += ((bytes + 7) / 8) * 8;

    return bytes;
  }

  /// Returns the maximum number of tree nodes. Which is 2n-1 for perfect
  /// full binary trees.
  inline std::size_t
  max_node_cnt() const {
    return max_node_cnt_;
  }

  /// Returns the label of the given node.
  inline u1
  label_of_node(u64 node_idx) const {
    return labels_[node_idx];
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

  inline u32
  get_first_node_idx_with_1label() const noexcept {
    return first_node_idx_with_1label_;
  }

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

  void __attribute__ ((noinline))
  run_optimize() {
    // Optimization level 2.
    if (optimization_level_ > 1) {
      // Estimates the size of a TEB.
      auto size = [&]() {
        return estimate_encoded_size_in_bytes();
      };

      // Gradual decompression.
      auto min = *this;
      auto min_size = size();
      for (auto it = breadth_first_begin(); it != breadth_first_end(); ++it) {
        u64 idx = (*it).idx;
        if (is_inner_node(idx)) continue;
        if (right_child_of(idx) >= max_node_cnt_) break;

        // Expand the leaf node to an inner node.
        set_inner(idx);
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
          // The current node, that has been expanded had a 0-label that was
          // implicit.
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
            for (std::size_t i = first_node_idx_with_1label_ + 1;
                 i < left_child_of(idx); ++i) {
              if (is_leaf_node(i) && label_of_node(i) == false) {
                trailing_0label_cnt_ -= 1;
              }
            }

          }
          last_node_idx_with_1label_ =
              std::max(last_node_idx_with_1label_, right_child_of(idx));
        }

        const auto compressed_size = size();
        if (compressed_size < min_size) {
          min = *this;
          min_size = compressed_size;
        }
      }
      *this = min;
    }
  }

  // Lossy compression.  The size of the tree structure is further reduced,
  // which causes false positive bits. The number of false positive bits is
  // limited by the given false positive rate (FPR).
  void
  compress_lossy(f64 fpr) {
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

      u1 left_bit = labels_[left_node_idx];
      u1 right_bit = labels_[right_node_idx];

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
          const auto subtree_size = this->subtree_size(node_idx);
          const auto pruned_node_cnt = subtree_size - 1;
          // Count the number of leaf nodes that would be eliminated from the
          // tree.
          const auto leaf_node_cnt = count_leaf_nodes(node_idx);
          const auto pruned_leaf_node_cnt = leaf_node_cnt - 1;
          // Compute the amount of bits saved.  One bit per node and one per
          // label.
          const auto saved_bit_cnt = pruned_node_cnt * 1.0625 + pruned_leaf_node_cnt;

          if (saved_bit_cnt > candidate_saved_bit_cnt) {
            candidate_node_idx = node_idx;
            candidate_saved_bit_cnt = saved_bit_cnt;
          }
        }

        // Check if we have found a candidate to prune.
        if (candidate_saved_bit_cnt > 0.0) {
          // The actual pruning.
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

};
//===----------------------------------------------------------------------===//

} // namespace dtl