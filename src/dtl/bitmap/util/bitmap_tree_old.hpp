#pragma once

#include <boost/dynamic_bitset.hpp>

#include <dtl/dtl.hpp>
#include <dtl/bitmap/util/binary_tree_structure.hpp>
#include <dtl/bitmap/util/rank1.hpp>

namespace dtl {

//===----------------------------------------------------------------------===//
/// Represents a bitmap as a binary tree. During construction, the bitmap tree
/// is compressed, either loss-less or lossy.
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

public:

  /// C'tor
  explicit
  bitmap_tree(const bitmap_t& bitmap, f64 fpr = 0.0)
      : binary_tree_structure(bitmap.size()),
        labels_(max_node_cnt_),
        inner_node_cnt_(0),
        leaf_node_cnt_(0),
        leading_inner_node_cnt_(0),
        trailing_leaf_node_cnt_(0) {

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

    // Lossy compression.  The size of the tree structure is further reduced,
    // which causes false positive bits. The number of false positive bits is
    // limited by the given false positive rate (FPR).
    if (fpr > 0.0) {
      // Determine maximum number of false positives.
      const auto tp_cnt = bitmap.count();
      const auto tn_cnt = n_ - tp_cnt;
      u64 max_fp_cnt = static_cast<u64>(tn_cnt * fpr);
      // The total number of false positives introduced so far.
      $u64 total_fp_cntr = 0;

      // Compute the number of associated false positives for each tree node.
      // The value fp_cntrs[i] refers to the number of false positives when
      // node i is turned into a leaf node (ie., the sub-tree is pruned).
      std::vector<uint32_t> fp_cntrs(max_node_cnt_, 0);
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

      // Pruning loop.
      while (true) {
        // Look for largest sub-tree to prune, so that the FPR is not exceeded.
        $u64 candidate_node_idx = 0;
        $u64 candidate_saved_bit_cnt = 0;
        for (std::size_t level = 0; level < height_; ++level) {
          u64 node_idx_from = (1ull << level) - 1;
          u64 node_idx_to = (1ull << (level + 1)) - 1;
          for ($u64 node_idx = node_idx_from; node_idx < node_idx_to;
               ++node_idx) {
            if (is_leaf_node(node_idx))
              continue;
            // Check whether pruning would not exceed the specified FPR.
            if (fp_cntrs[node_idx] + total_fp_cntr > max_fp_cnt)
              continue;

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
            const auto saved_bit_cnt = pruned_node_cnt + pruned_leaf_node_cnt;

            if (saved_bit_cnt > candidate_saved_bit_cnt) {
              candidate_node_idx = node_idx;
              candidate_saved_bit_cnt = saved_bit_cnt;
            }
          }
          // Check, whether we found a candidate on the current tree level.
          if (candidate_saved_bit_cnt > 0) break;
        }

        if (candidate_saved_bit_cnt > 0) {
          // The actual pruning.
          set_leaf(candidate_node_idx);
          total_fp_cntr += fp_cntrs[candidate_node_idx];

          // Ensure that the tree bitmap is still compressed. (bottom-up pruning)
          auto node_idx = candidate_node_idx;
          while (true) {
            if (node_idx == root()) break;
            // Go to the parent node.
            node_idx = parent_of(node_idx);
            u64 left_child_idx = left_child_of(node_idx);
            u64 right_child_idx = right_child_of(node_idx);

            u1 both_childs_are_leaves =
                !is_inner_node(left_child_idx)
                    & !is_inner_node(right_child_idx);

            u1 left_bit = labels_[left_child_idx];
            u1 right_bit = labels_[right_child_idx];
            u1 can_prune = both_childs_are_leaves && (left_bit == right_bit);

            if (can_prune) {
              set_leaf(node_idx);
            }
            else {
              break;
            }
          }
        }
        else {
          // Exit the pruning loop, because no candidate has been found.
          break;
        }
      }
    }

    // Determine the total number of tree nodes and which of these nodes need
    // to be stored explicitly.
    {
      $u1 found_leaf_node = false;
      std::size_t node_cnt = 0;
      explicit_node_idxs_.end = 0;

      for (auto it = breadth_first_begin(); it != breadth_first_end(); ++it) {
        u64 idx = *it;
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
          explicit_node_idxs_.end = *it;
        }
      }

      explicit_node_idxs_.begin = leading_inner_node_cnt_;
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
        u64 idx = *it;
        if (is_inner_node(idx)) continue;
        if (right_child_of(idx) >= max_node_cnt_) break;
        // Expand inner node.
        set_inner(idx);
        // Update the counters.
        ++inner_node_cnt_;
        ++leaf_node_cnt_;

        ++leading_inner_node_cnt_;
        for (auto i = idx + 1; i < max_node_cnt_ && is_inner_node(i); ++i) {
          ++leading_inner_node_cnt_;
        }
        explicit_node_idxs_.begin = leading_inner_node_cnt_;
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

        const auto compressed_size = size();
        if (compressed_size < min_size) {
          min = *this;
          min_size = compressed_size;
        }
      }
      *this = min;
    }
  }

  std::size_t
  estimate_encoded_size_in_bytes() {
    constexpr u64 block_bitlength = 32;
    constexpr u64 block_size = 4;
    $u64 bytes = 0;
    // Tree structure
    u64 explicit_tree_node_cnt = inner_node_cnt_ + leaf_node_cnt_
        - leading_inner_node_cnt_ - trailing_leaf_node_cnt_;
    bytes += ((explicit_tree_node_cnt + block_bitlength - 1) / block_bitlength)
        * block_size;
    // Labels
    bytes += ((leaf_node_cnt_ + block_bitlength - 1) / block_bitlength)
        * block_size;
    // Rank support
    bytes += dtl::rank1::estimate_size_in_bytes(explicit_tree_node_cnt);
    // Bit-length of the original bitmap.
    bytes += sizeof(n_);
    // The number of implicit inner nodes.
    bytes += 4;
    // FIXME: + number of tree bits + number of label bits
    return bytes;
  }

//  std::size_t
//  determine_explicit_node_cnt() {
//    $u1 found_leaf_node = false;
//    std::size_t inner_node_cnt = 0;
//    std::size_t leaf_node_cnt = 0;
//    std::size_t leading_inner_node_cnt = 0;
//    std::size_t trailing_leaf_node_cnt = 0;
//
//    for (auto it = breadth_first_begin(); it != breadth_first_end(); ++it) {
//      u64 idx = *it;
//      u1 is_inner = is_inner_node(idx);
//      inner_node_cnt += is_inner;
//      leaf_node_cnt += !is_inner;
//
//      // Count the leading inner nodes.
//      if (!found_leaf_node && is_inner) {
//        ++leading_inner_node_cnt;
//      }
//      if (!found_leaf_node && !is_inner) {
//        found_leaf_node = true;
//      }
//
//      // Count the trailing leaf nodes.
//      if (!is_inner) {
//        ++trailing_leaf_node_cnt;
//      } else {
//        trailing_leaf_node_cnt = 0;
//      }
//    }
//    std::cout
//        << "[bt]  total: " << (inner_node_cnt + leaf_node_cnt)
//        << ", leading: " << leading_inner_node_cnt
//        << ", explicit: " << (inner_node_cnt + leaf_node_cnt
//          - leading_inner_node_cnt - trailing_leaf_node_cnt)
//        << ", trailing: " << trailing_leaf_node_cnt
//        << ", labels: " << leaf_node_cnt
//        << std::endl;
//    return inner_node_cnt + leaf_node_cnt
//        - leading_inner_node_cnt - trailing_leaf_node_cnt;
//  }

//  struct counters_t {
//    std::size_t inner_node_cnt = 0;
//    std::size_t leaf_node_cnt = 0;
//    std::size_t leading_inner_node_cnt = 0;
//    std::size_t trailing_leaf_node_cnt = 0;
//    range_t explicit_node_idx_ {0,0};
//  };
//
//  counters_t
//  _count() {
//    counters_t c;
//    $u1 found_leaf_node = false;
//
//    for (auto it = breadth_first_begin(); it != breadth_first_end(); ++it) {
//      u64 idx = *it;
//      u1 is_inner = is_inner_node(idx);
//      c.inner_node_cnt += is_inner;
//      c.leaf_node_cnt += !is_inner;
//
//      // Count the leading inner nodes.
//      if (!found_leaf_node && is_inner) {
//        ++c.leading_inner_node_cnt;
//      }
//      if (!found_leaf_node && !is_inner) {
//        found_leaf_node = true;
//      }
//
//      // Count the trailing leaf nodes.
//      if (!is_inner) {
//        ++c.trailing_leaf_node_cnt;
//      } else {
//        c.trailing_leaf_node_cnt = 0;
//        c.explicit_node_idx_.end = idx;
//      }
//    }
//    c.explicit_node_idx_.begin = c.leading_inner_node_cnt;
//    return c;
//  }

  /// Returns the maximum number of tree nodes. Which is 2n-1 for binary trees.
  inline std::size_t
  max_node_cnt() const {
    return max_node_cnt_;
  }

  inline u1
  label_of_node(u64 node_idx) const {
    return labels_[node_idx];
  }

  inline u32
  get_node_cnt() const noexcept {
    return inner_node_cnt_ + leaf_node_cnt_;
  }

  inline u32
  get_leading_inner_node_cnt() const noexcept {
    return leading_inner_node_cnt_;
  }

  inline u32
  get_trailing_leaf_node_cnt() const noexcept {
    return trailing_leaf_node_cnt_;
  }

  inline u32
  get_first_explicit_node_idx() const noexcept {
    return explicit_node_idxs_.begin;
  }

  inline u32
  get_last_explicit_node_idx() const noexcept {
    return explicit_node_idxs_.end;
  }

};
//===----------------------------------------------------------------------===//

} // namespace dtl