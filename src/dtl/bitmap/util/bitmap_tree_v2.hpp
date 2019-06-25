#pragma once

#include <iomanip>

#include <boost/dynamic_bitset.hpp>

#include <dtl/dtl.hpp>
#include <dtl/bitmap/util/binary_tree_structure.hpp>
#include <dtl/bitmap/util/rank1.hpp>
#include <dtl/color.hpp>

namespace dtl {
// FIXME implementation still incomplete
//===----------------------------------------------------------------------===//
/// Represents a bitmap as a binary tree. During construction, the bitmap tree
/// is compressed, either loss-less or lossy.
/// The template parameter controls the space optimizations.
template<i32 optimization_level_ = 3>
class bitmap_tree : public binary_tree_structure {

  using tree_t = dtl::binary_tree_structure;
  using bitmap_t = boost::dynamic_bitset<$u32>;

  /// The labels of the tree nodes.
  bitmap_t labels_;

  struct range_t {
    std::size_t begin = 0;
    std::size_t end = 0;
    u1 empty() { return end <= begin; }
  };

  std::size_t inner_node_cnt_;
  std::size_t leaf_node_cnt_;
  std::size_t implicit_leading_inner_node_cnt_;
  std::size_t implicit_trailing_leaf_node_cnt_;
//  std::size_t first_explicit_node_idx_;
//  std::size_t last_explicit_node_idx_;
  range_t explicit_node_idxs_;

  std::size_t leading_0label_cnt_;
  std::size_t trailing_0label_cnt_;
//  std::size_t first_node_idx_with_1label_;
//  std::size_t last_node_idx_with_1label_;
  range_t explicit_leaf_idxs_;

  std::size_t pop_count;

public:

  /// C'tor
  explicit
  bitmap_tree(const bitmap_t& bitmap, f64 fpr = 0.0)
      : binary_tree_structure(bitmap.size()),
        labels_(max_node_cnt_),
        inner_node_cnt_(0),
        leaf_node_cnt_(0),
        implicit_leading_inner_node_cnt_(0),
        implicit_trailing_leaf_node_cnt_(0),
        explicit_node_idxs_(),
        leading_0label_cnt_(0),
        trailing_0label_cnt_(0),
        explicit_leaf_idxs_(),
        pop_count(0) {

    // Init the binary tree.
    init_tree(bitmap);
    reduce_tree();

//
//    // Run space optimizations.
////    {
////      // Optimization level 2
////      if (optimization_level_ > 1) {
////        init_counters();
////        run_optimize();
////      }
////    }
//    const auto lossless_size = estimate_encoded_size_in_bytes();
//
//    if (fpr > 0.0) {
//      auto min_size = lossless_size;
//      auto min = *this;
//      auto min_threshold = 0.0;
//
////      $u1 monotonic = true;
////      $u1 first_lossy_pass = true;
////      std::size_t prev_lossy_size = std::numeric_limits<std::size_t>::max();
////      std::vector<std::size_t> lossy_sizes;
//
//      for (auto threshold : {1.0/8.0, 1.0/4.0, 1.0/2.0, 1.0}) {
//        bitmap_tree lossy_bt(bitmap, fpr, threshold);
//        const auto size = lossy_bt.estimate_encoded_size_in_bytes();
//        if (size < min_size) {
//          min = lossy_bt;
//          min_size = size;
//          min_threshold = threshold;
//        }
////        if (size <= prev_lossy_size) {
////          prev_lossy_size = size;
////        }
////        else {
////          break;
////        }
//      }
//      *this = min;
//
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
//    }
  }

  bitmap_tree(const bitmap_tree& other) = default;
  bitmap_tree(bitmap_tree&& other) noexcept = default;
  bitmap_tree& operator=(const bitmap_tree& other) = default;
  bitmap_tree& operator=(bitmap_tree&& other) noexcept = default;
  ~bitmap_tree() override = default;

private:

  explicit
  bitmap_tree(const bitmap_t& bitmap, f64 fpr, f64 threshold)
      : binary_tree_structure(bitmap.size()),
        labels_(max_node_cnt_),
        inner_node_cnt_(0),
        leaf_node_cnt_(0),
        implicit_leading_inner_node_cnt_(0),
        implicit_trailing_leaf_node_cnt_(0),
        explicit_node_idxs_(),
        leading_0label_cnt_(0),
        trailing_0label_cnt_(0),
        explicit_leaf_idxs_() {

    // Init the binary tree.
    init_tree(bitmap);
    reduce_tree();

//    compress_lossy_pre_optimization_v2(fpr, threshold);
//    init_counters();
//    if (optimization_level_ > 1) {
//      run_optimize();
//    }
  }

  /// Initialize a perfect binary tree on top of the given bitmap.
  ///  - all the inner nodes have two children
  ///  - all the leaf nodes are on the same level
  ///  - the leaf nodes are labelled with the given bitmap
  void
  init_tree(const bitmap_t& bitmap) {
    if (!dtl::is_power_of_two(n_)) {
      throw std::invalid_argument(
          "The length of the bitmap must be a power of two.");
    }

    u64 length = max_node_cnt_;
    u64 height = height_;
    for ($u64 i = length / 2; i < length; i++) {
      labels_[i] = bitmap[i - length / 2];
    }

    // Propagate the bits along the tree (bottom-up).  The labels of an internal
    // node is the bitwise OR of the labels of both child nodes.
    for ($u64 i = 0; i < length - 1; i++) {
      u64 node_idx = length - i - 1;
      labels_[parent_of(node_idx)] =
          labels_[parent_of(node_idx)] | labels_[node_idx];
    }

    // Initialize the counters.
    inner_node_cnt_ = n_ - 1;
    leaf_node_cnt_ = n_;
    pop_count = 0;

    implicit_leading_inner_node_cnt_ = inner_node_cnt_;
    implicit_trailing_leaf_node_cnt_ = leaf_node_cnt_;
    // Defaults to the first leaf node.
    explicit_node_idxs_.begin = inner_node_cnt_;
    explicit_node_idxs_.end = explicit_node_idxs_.begin;

    leading_0label_cnt_ = 0;
    trailing_0label_cnt_ = 0;

    explicit_leaf_idxs_.begin = inner_node_cnt_;
    explicit_leaf_idxs_.end = explicit_leaf_idxs_.begin;
    {
      $u1 found_1label = false;
      for (auto idx = first_node_idx_at_level(last_level());
           idx < first_node_idx_at_level(last_level() + 1); ++idx) {

        u1 label = label_of_node(idx);
        pop_count += label;

        if (!found_1label) {
          explicit_leaf_idxs_.begin = idx;
          explicit_leaf_idxs_.end = idx;
          if (label) {
            // First 1-label.
            found_1label = true;
          }
          else {
            // Leading 0-label(s).
            ++leading_0label_cnt_;
          }
        }

        found_1label |= label;

        if (found_1label) {
          if (label) {
            trailing_0label_cnt_ = 0;
            explicit_leaf_idxs_.end = idx + 1;
          }
          else {
            ++trailing_0label_cnt_;
          }
        }
      }
    }
    counters_are_valid();
  }

  /// Handle the special case where the entire bitmap is zero.
  void
  set_false() {
    for (std::size_t i = 1; i < first_node_idx_at_level(last_level() + 1); ++i) {
      set_leaf(i);
    }
    set_leaf(root());
    inner_node_cnt_ = 0;
    leaf_node_cnt_ = 1;
    implicit_leading_inner_node_cnt_ = 0;
    implicit_trailing_leaf_node_cnt_ = 1;
    explicit_node_idxs_.begin = explicit_node_idxs_.end = 1;

    labels_[root()] = false;
    leading_0label_cnt_ = 1;
    trailing_0label_cnt_ = 0;
    explicit_leaf_idxs_.begin = explicit_leaf_idxs_.end = 1;
  }

  /// Bottom-up pruning (lossless).  Eliminate all sibling leaf nodes which
  /// have the same label.
  void
  reduce_tree() {

    if (pop_count == 0) {
      // Handle the special case where the entire bitmap is zero.
      set_false();
      std::cout << "-----------------------------" << std::endl;
      std::cout << (*this) << std::endl;
      return;
    }

    explicit_node_idxs_.begin = first_node_idx_at_level(last_level());
    explicit_node_idxs_.end = explicit_node_idxs_.begin;

    //===------------------------------------------------------------------===//
    // The minimum encoded bitmap size observed so far.
    auto min_size = estimate_encoded_size_in_bytes();
    // As the tree is reduced starting from the highest inner node index in
    // descending order, we can identify each reduced tree instance by the
    // node idx that has been reduced last.
    // Thus, to keep track of the smallest TEB that has been observed so far,
    // we only need to memorize the node index of the last collapsed node.
    $u64 min_collapsed_node_idx = max_node_cnt_ - n_ + 1;
    //===------------------------------------------------------------------===//

    std::cout << (*this) << std::endl;

    // Iterate over all inner nodes in descending order.
    u64 length = max_node_cnt_ - n_;
    for ($u64 i = 0; i < length - 1; ++i) {
      // The current node, which is currently the highest order inner node.
      u64 node_idx = length - i - 1;
      if (node_idx == root()) {
        // Do not collapse the root node.
        break;
      }

      u64 left_child_idx = left_child_of(node_idx);
      u64 right_child_idx = left_child_idx + 1;

      u1 left_bit = label_of_node(left_child_idx);
      u1 right_bit = label_of_node(right_child_idx);

      u1 prune_causes_false_positives = left_bit ^ right_bit;
      u1 both_child_nodes_are_leaves =
          is_leaf_node(left_child_idx) & is_leaf_node(right_child_idx);
      u1 prune = both_child_nodes_are_leaves & !prune_causes_false_positives;
      if (prune) {
        // Prune the two child nodes and turn the current node a leaf.
        set_leaf(node_idx);

        //===--------------------------------------------------------------===//
        // Adjust the counters (which are required to compute the TEB size).
        inner_node_cnt_ -= 1;
        leaf_node_cnt_ -= 1;

        // --- Implicit / explicit nodes ---
        // All nodes with a lower index than the current are inner nodes.
        implicit_leading_inner_node_cnt_ = node_idx;

        if (explicit_node_idxs_.begin == node_idx + 1
            && explicit_node_idxs_.end == node_idx + 1) {
          explicit_node_idxs_.begin = node_idx;
          explicit_node_idxs_.end = node_idx;
        }
        else {
          explicit_node_idxs_.begin = node_idx;
        }
        implicit_trailing_leaf_node_cnt_ = inner_node_cnt_ + leaf_node_cnt_
            - implicit_leading_inner_node_cnt_
            - (explicit_node_idxs_.end - explicit_node_idxs_.begin);
//        if (explicit_node_idxs_.end == node_idx) {
//          explicit_node_idxs_.end = node_idx;
//          implicit_trailing_leaf_node_cnt_ -= 2;
//        }

        // --- Implicit / explicit labels ---
        if (!right_bit) {
          // We are pruning two leaves with 0-labels.

          // Adjust the trailing 0-label count if necessary.

          if (explicit_leaf_idxs_.begin == right_child_idx + 1
              && explicit_leaf_idxs_.end == right_child_idx + 1) {
            // Pruning leading 0-labels.
            explicit_leaf_idxs_.begin = left_child_idx;
            explicit_leaf_idxs_.end = left_child_idx;
            leading_0label_cnt_ -= 1;
          }
          else {
            if (left_child_idx >= explicit_leaf_idxs_.end) {
              // Pruning trailing 0-labels.
              trailing_0label_cnt_ -= 2;
            }
            else {
              // Pruning explicit 0 labels.
            }
          }

          // Adjust the leading 0-label count. The begin of the explicit leaves
          // range does not change either.
          leading_0label_cnt_ += 1;
        }
        else {
          // We are pruning two leaves with 1-labels.

          if (explicit_leaf_idxs_.end == right_child_idx + 1) {
            // Pruning the the last explicit leaves.
            // The new last explicit leaf is now somewhere within the
            // range [node_idx, left_child_idx). // TODO quadratic runtime?
            for (std::size_t i = left_child_idx - 1; i >= node_idx; --i) {
              if (is_leaf_node(i) && label_of_node(i)) break;
              explicit_leaf_idxs_.end = i;
              ++trailing_0label_cnt_;
            }

          }

          // The current node becomes the first explicit leaf, as it has a
          // 1-label and all prior nodes are inner nodes.
          explicit_leaf_idxs_.begin = node_idx;
          leading_0label_cnt_ = 0;
        }
        //===--------------------------------------------------------------===//

        // TODO compute the size and compare to the minimum
        std::cout << "-----------------------------" << std::endl;
        std::cout << (*this) << std::endl;
        counters_are_valid();
      }
    }
  }

  // Lossy compression.  The size of the tree structure is further reduced,
  // which causes false positive bits. The number of false positive bits is
  // limited by the given false positive rate (FPR).
  void __attribute__ ((noinline))
  compress_lossy(f64 fpr, f64 threshold = 1.0/8.0) {
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

      u1 left_bit = labels_[left_node_idx];
      u1 right_bit = labels_[right_node_idx];

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
          set_leaf(parent_node_idx);
          const auto fp_cntr_old = total_fp_cntr;
          total_fp_cntr -= intr_fp_cntrs[parent_node_idx];
          intr_fp_cntrs[parent_node_idx] = fp_cntrs[parent_node_idx];
          total_fp_cntr += intr_fp_cntrs[parent_node_idx];
        }
      }
    }
  }

public:

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
        - implicit_leading_inner_node_cnt_ - implicit_trailing_leaf_node_cnt_;

    // Level offsets for T and L, which are required by the tree scan algorithm.
    const auto perfect_levels = dtl::log_2(implicit_leading_inner_node_cnt_ + 1) + 1;
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
    return implicit_leading_inner_node_cnt_;
  }

  /// Returns the number of trailing leaf nodes (in level order).
  inline u32
  get_trailing_leaf_node_cnt() const noexcept {
    return implicit_trailing_leaf_node_cnt_;
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
    return explicit_leaf_idxs_.begin;
  }

  inline u32
  get_last_node_idx_with_1label() const noexcept {
    return explicit_leaf_idxs_.end;
  }

  void
  print(std::ostream& os) const noexcept {
    std::stringstream t_out;
    std::stringstream l_out;
    t_out << "T:" << std::endl;
    l_out << "L:" << std::endl;

    dtl::color_modifier col_default(dtl::color::gray);
    dtl::color_modifier col_explicit(dtl::color::yellow);
    std::string implicit_inner_node_char = ".";
    std::string explicit_inner_node_char =
        col_explicit.str() + "." + col_default.str();
    std::string implicit_leaf_node_char = "*";
    std::string explicit_leaf_node_char =
        col_explicit.str() + "*" + col_default.str();

    std::string implicit_0label_char = "0";
    std::string explicit_0label_char =
        col_explicit.str() + "0" + col_default.str();
    std::string explicit_1label_char =
        col_explicit.str() + "1" + col_default.str();

    for (std::size_t level = 0; level <= height_; ++level) {
      t_out << std::setw(4) << level << ": ";
      l_out << std::setw(4) << level << ": ";
      u64 node_idx_from = (1ull << level) - 1;
      u64 node_idx_to = (1ull << (level + 1)) - 1;
      for ($u64 node_idx = node_idx_from; node_idx < node_idx_to; ++node_idx) {

        const auto node_is_explicit = node_idx >= explicit_node_idxs_.begin
            && node_idx < explicit_node_idxs_.end;
        const auto inner_char = node_is_explicit
            ? explicit_inner_node_char : implicit_inner_node_char;
        const auto leaf_char = node_is_explicit
            ? explicit_leaf_node_char : implicit_leaf_node_char;

        const auto label_is_explicit = node_idx >= explicit_leaf_idxs_.begin
            && node_idx < explicit_leaf_idxs_.end;
        const auto zero_label = label_is_explicit
            ? explicit_0label_char : implicit_0label_char;
        const auto one_label = explicit_1label_char;


        const auto spaces = (level < height_)
            ? (1u << (height_ - level)) - 1
            : 0;
        if (node_idx == root()) {
          t_out << (is_leaf_node(node_idx) ? leaf_char : inner_char);
          l_out << (is_leaf_node(node_idx)
              ? label_of_node(node_idx) ? one_label : zero_label
              : ".");
        }
        else {

          t_out << (is_leaf_node(node_idx)
              ? is_leaf_node(parent_of(node_idx)) ? " " : leaf_char
              : inner_char
          );
          l_out << (is_leaf_node(node_idx) & is_inner_node(parent_of(node_idx))
              ? label_of_node(node_idx) ? one_label : zero_label
              : is_inner_node(node_idx) ? "." : " "
          );
        }
        for (std::size_t i = 0; i < spaces; ++i) {
          t_out << " ";
          l_out << " ";
        }
      }
      t_out << std::endl;
      l_out << std::endl;
    }
    t_out << "cntrs: "
        << "#inner=" << inner_node_cnt_
        << ", #leaves=" << leaf_node_cnt_
        << ", leading=" << implicit_leading_inner_node_cnt_
        << ", trailing=" << implicit_trailing_leaf_node_cnt_
        << std::endl;
    l_out << "cntrs: "
        << "leading=" << leading_0label_cnt_
        << ", trailing=" << trailing_0label_cnt_
        << std::endl;
    os << t_out.str() << l_out.str();
  }


  // FOR DEBUGGING PURPOSES
  u1
  counters_are_valid() {
    const auto inner_node_cnt_old = inner_node_cnt_;
    const auto leaf_node_cnt_old = leaf_node_cnt_;
    const auto implicit_leading_inner_node_cnt_old = implicit_leading_inner_node_cnt_;
    const auto implicit_trailing_leaf_node_cnt_old = implicit_trailing_leaf_node_cnt_;
    const auto explicit_node_idxs_old = explicit_node_idxs_;
    const auto leading_0label_cnt_old = leading_0label_cnt_;
    const auto trailing_0label_cnt_old = trailing_0label_cnt_;
    const auto explicit_leaf_idxs_old = explicit_leaf_idxs_;
    
    inner_node_cnt_ = 0;
    leaf_node_cnt_ = 0;
    implicit_leading_inner_node_cnt_ = 0;
    implicit_trailing_leaf_node_cnt_ = 0;
    explicit_node_idxs_.begin = 0;
    explicit_node_idxs_.end = 0;
    
    leading_0label_cnt_ = 0;
    trailing_0label_cnt_ = 0;
    explicit_leaf_idxs_.begin = 0;
    explicit_leaf_idxs_.end = 0;

    // Determine the total number of tree nodes and which of these nodes need
    // to be stored explicitly.
    $u1 found_leaf_node = false;
    std::size_t node_cnt = 0;

    $u1 found_leaf_node_with_1label = false;

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
        ++implicit_leading_inner_node_cnt_;
      }
      if (!found_leaf_node && !is_inner) {
        found_leaf_node = true;
      }

      // Count the trailing leaf nodes.
      if (!is_inner) {
        ++implicit_trailing_leaf_node_cnt_;
      } else {
        implicit_trailing_leaf_node_cnt_ = 0;
        explicit_node_idxs_.end = idx;
      }

      u1 is_leaf = is_leaf_node(idx);

      if (is_leaf) {
        u1 label = label_of_node(idx);
        if (!found_leaf_node_with_1label) {
          if (!label) {
            ++leading_0label_cnt_;
          }
          else {
            found_leaf_node_with_1label = true;
            explicit_leaf_idxs_.begin = idx;
            explicit_leaf_idxs_.end = idx + 1;
          }
        }
        else {
          if (!label) {
            ++trailing_0label_cnt_;
          }
          else {
            trailing_0label_cnt_ = 0;
            explicit_leaf_idxs_.end = idx + 1;
          }
        }
      }
    }
    explicit_node_idxs_.begin = implicit_leading_inner_node_cnt_;
    if (explicit_node_idxs_.end < explicit_node_idxs_.begin) {
      explicit_node_idxs_.end = explicit_node_idxs_.begin;
    }

    assert(inner_node_cnt_old == inner_node_cnt_);
    assert(leaf_node_cnt_old == leaf_node_cnt_);
    assert(implicit_leading_inner_node_cnt_old == implicit_leading_inner_node_cnt_);
    assert(implicit_trailing_leaf_node_cnt_old == implicit_trailing_leaf_node_cnt_);
    assert(explicit_node_idxs_old.begin == explicit_node_idxs_.begin);
    assert(explicit_node_idxs_old.end == explicit_node_idxs_.end);
    assert(leading_0label_cnt_old == leading_0label_cnt_);
    assert(trailing_0label_cnt_old == trailing_0label_cnt_);
    assert(explicit_leaf_idxs_old.begin == explicit_leaf_idxs_.begin);
    assert(explicit_leaf_idxs_old.end == explicit_leaf_idxs_.end);
  }

};
//===----------------------------------------------------------------------===//

} // namespace dtl