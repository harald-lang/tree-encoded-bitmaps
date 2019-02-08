#pragma once

#include <bitset>
#include <list>
#include <queue>
#include <stack>
#include <vector>

#include "boost/dynamic_bitset.hpp"

#include <dtl/bits.hpp>
#include <dtl/bitmap/util/rank1.hpp>
#include <dtl/bitmap/util/rank1_surf.hpp>
#include <dtl/bitmap/util/rank1_surf_cached.hpp>
#include <dtl/bitmap/util/bitmap_tree.hpp>
#include <dtl/bitmap/util/binary_tree_structure.hpp>
#include <dtl/static_stack.hpp>
#include <dtl/static_stack2.hpp>
#include <dtl/bitmap/util/rank1_naive.hpp>
#include <dtl/bitmap/util/rank1_interceptor.hpp>
#include <dtl/bitmap/util/rank1_super_fast.hpp>

namespace dtl {


#if !defined(__teb_inline__)
#if defined(NDEBUG)
// Release build.
#define __teb_inline__ inline __attribute__((always_inline))
//#define __teb_inline__ __attribute__((noinline))
#else
#define __teb_inline__
#endif
#endif

//===----------------------------------------------------------------------===//
/// Encodes a bitmap of length n as a binary tree.
///
/// The implementation supports the optimization levels 0, 1, and 2.
///   0 = implements the core idea of TEBs
///   1 = uses the implicit inner node optimization
///   2 = additionally applies the gradual decompression optimization.
/// Please note, that the different optimization levels can be chosen due to
/// academic reasons only. I.e., to analyse and compare the different
/// optimizations levels.
///
/// Terminology:
/// - node index: the nodes identifier, which is its position within the
///       bit sequence of the encoded tree.
/// - path: also identifies a tree node, but it encodes the path to the node,
///       starting at the root
/// - perfect levels: the number of tree levels in which the tree is perfect.
/// - top nodes: refer to the nodes within the last perfect level.
///
template<i32 optimization_level_ = 2>
class teb {

public:

  // The fundamental storage type. The size of a TEB is a multiple of
  // sizeof(_block_type).
  using _block_type = $u64;
  using word_type = _block_type;

  using position_t = $u32;
  using bitmap_t = boost::dynamic_bitset<word_type>;

  /// The number of bits in the bitmap.
  $u64 n_;

  /// The tree encoded bitmap.
  bitmap_t structure_;
  bitmap_t labels_;

  /// Support data structure for rank1 operations on the tree structure.
  static constexpr u1 inclusive = true;
  static constexpr u1 non_inclusive = false;
//  using rank_support = dtl::rank1_naive<word_type>;
//  using rank_support = dtl::rank1_surf_cached<word_type>;
//  using rank_support = dtl::rank1<word_type>;
//  using rank_support = dtl::rank1_interceptor<dtl::rank1_surf_cached<word_type, non_inclusive>>;
//  using rank_support = dtl::rank1_interceptor<dtl::rank1_surf_cached<word_type, inclusive>>;
//  using rank_support = dtl::rank1_interceptor<dtl::rank1_surf<word_type, inclusive>>;
//  using rank_support = dtl::rank1_surf_cached<word_type, inclusive>;
//  using rank_support = dtl::rank1_surf<word_type, inclusive>;
  using rank_support = dtl::rank1_super_fast<word_type, inclusive>;
  rank_support rank_;

  /// The number of implicit inner nodes in the tree structure.
  $u32 implicit_inner_node_cnt_;

  /// For testing purposes only.
  $u32 implicit_leaf_node_cnt_;

public:

  /// Tree-encode the given bitmap with an optional false positive rate.
  explicit
  teb(const boost::dynamic_bitset<$u32>& bitmap, f64 fpr = 0.0)
      : n_(bitmap.size()) {

    // Construct a binary tree that represents the given bitmap.
    // Space-optimizations are performed in the (non-encoded) bitmap tree.
    dtl::bitmap_tree<optimization_level_> bitmap_tree(bitmap, fpr);

    // Encode the tree into level-order.
    implicit_inner_node_cnt_ = 0;
    implicit_leaf_node_cnt_ = 0;
    if (optimization_level_ > 0) {
      // Implicit tree nodes are only considered with -o1 or higher.
      implicit_inner_node_cnt_ = bitmap_tree.get_leading_inner_node_cnt();
      implicit_leaf_node_cnt_ = bitmap_tree.get_trailing_leaf_node_cnt();
    }
    std::size_t node_cntr = 0;
    for (auto it = bitmap_tree.breadth_first_begin();
         it != bitmap_tree.breadth_first_end();
         ++it) {
      ++node_cntr;

      // Omit the implicit nodes.
      if (node_cntr <= implicit_inner_node_cnt_) {
        continue;
      }
      u64 idx = *it;

      // Emit a 1-bit if the current node is an inner node, a 0-bit otherwise.
      u1 is_inner = bitmap_tree.is_inner_node(idx);
      if (is_inner
          || optimization_level_ == 0
          || (optimization_level_ > 0
              && idx <= bitmap_tree.get_last_explicit_node_idx())) {
        structure_.push_back(is_inner);
      }
      if (!is_inner) {
        // Add the label of the leaf node.
        labels_.push_back(bitmap_tree.label_of_node(idx));
      }
    }

    // Hack: ensure structure is not null
    if (structure_.size() == 0) {
      structure_.push_back(false);
    }

    // Init rank1 support data structure.
    rank_.init(structure_);
  }

  teb(const teb& other) = default;
  teb(teb&& other) = default;
  teb& operator=(const teb& other) = default;
  teb& operator=(teb&& other) = default;
  ~teb() = default;

  inline void
  run_optimize() {
    // Optimization level 2.
    if (optimization_level_ > 1) {
      auto min = *this;
      auto min_size = size_in_byte();
      // Gradual decompression.
      $u64 global_cntr = 0;
      $u64 local_cntr = 0;
      while (decompress()) {
        global_cntr++;
        local_cntr++;
        const auto compressed_size = size_in_byte();
        if (compressed_size < min_size) {
          min = *this;
          min_size = compressed_size;
          local_cntr = 0;
        }
      }
      *this = min;
    }
  }

  /// For testing purposes only.
  u1
  decompress() {
    if (structure_.size() == 0) {
      return false;
    }
    assert(structure_[0] == false);
    assert(structure_[structure_.size() - 1] == true);

    // The leaf node, that is turned into an inner node.
    u64 node_idx = implicit_inner_node_cnt_;
    assert(is_leaf_node(node_idx));
    u1 node_label = get_label(node_idx);
    // The position of the nodes' label.
    u64 node_label_idx = 0;

    // The new child nodes. Note, we do not need the rank function here,
    // because all prior nodes are inner nodes.
    u64 abs_left_child_idx = 2 * node_idx + 1;
    u64 abs_right_child_idx = abs_left_child_idx + 1;
    u64 left_child_idx = (2 * node_idx + 1) - implicit_inner_node_cnt_;
    u64 right_child_idx = left_child_idx + 1;

    // The new structure of the gradually decompressed tree.
    bitmap_t s(structure_.size() + 2);
    // The first bit is the bit that corresponds to the affected node, which
    // is turned into an inner node.
    s[0] = true;
    // Copy the structure, but leave space for the two new leaf nodes.
    for (std::size_t i = 1; i < std::min(left_child_idx, s.size()); ++i) {
      assert(i < s.size());
      s[i] = i < structure_.size() ? structure_[i] : false;
    }
    if (left_child_idx < s.size()) {
      s[left_child_idx] = false;
    }
    if (right_child_idx < s.size()) {
      s[right_child_idx] = false;
    }
    for (std::size_t i = right_child_idx + 1; i < s.size(); i++) {
      assert(i < s.size());
      s[i] = (i - 2) < structure_.size() ? structure_[i - 2] : false;
    }

    // Count the number of trailing 1-bits.
    std::size_t trailing_1bits = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
      if (!s[i]) {
        break;
      }
      trailing_1bits++;
    }
    // Remove the trailing 1-bits from the tree structure.
    s >>= trailing_1bits;
    s.resize(s.size() - trailing_1bits);

    // Update the number of implicit inner nodes.
    implicit_inner_node_cnt_ += static_cast<u32>(trailing_1bits);

    // Remove the most significant 0-bits from the tree structure.
    while (!s.empty() && !s[s.size() - 1]) {
      s.pop_back();
      implicit_leaf_node_cnt_++;
    }
    s.resize(s.size()); // shrink to fit
    std::swap(structure_, s);

    // Init rank1 support data structure.
    rank_.init(structure_);

    // The new labels of the gradually decompressed tree.
    bitmap_t l(labels_.size() + 1);

    // The positions of the new label.
    u64 left_child_label_idx = get_label_idx(abs_left_child_idx);
    u64 right_child_label_idx = get_label_idx(abs_right_child_idx);

    for (std::size_t i = 0; i < left_child_label_idx; ++i) {
      assert(i < l.size());
      assert((i + 1) < labels_.size());
      l[i] = labels_[i + 1];
    }
    l[left_child_label_idx] = node_label;
    l[right_child_label_idx] = node_label;
    for (std::size_t i = right_child_label_idx + 1; i < l.size(); ++i) {
      assert(i < l.size());
      assert((i - 1) < labels_.size());
      l[i] = labels_[i - 1];
    }
    l.resize(l.size()); // shrink to fit
    std::swap(labels_, l);

    return true;
  }

  /// Decodes the level-order encoding to a bitmap.
  boost::dynamic_bitset<$u32> __teb_inline__
  to_bitset() const noexcept {
    boost::dynamic_bitset<$u32> ret(n_); // the resulting bitmap

    // special case if the tree is only a root node
    if(implicit_inner_node_cnt_ == 0 && structure_.size() <= 1){
      //std::cout << "to_bitset() only root case" << std::endl;
      if(labels_[0]) {
        ret.flip();
      }
      return ret;
    }

    // normal cases: construct the tree using the lo_construction and a level counter
    $u64 level = 0; // current height
    $u64 write_pointer = 0;
    u64 tree_height = dtl::log_2(n_);

    // a doubly linked list<node_idx, level> to manage the decoding for all nodes
    std::list<std::pair<$u64, $u64>> nodes;

    std::function<void(std::pair<u64, u64>)> decode_tree = [&](std::pair<u64, u64> n) {
      u64 idx = n.first;
      u64 level = n.second;
      if (is_inner_node(idx)) {
        u64 l_child = left_child(idx);
        u64 r_child = right_child(idx);
        u64 child_level = level+1;

        nodes.push_front(std::make_pair(r_child, child_level)); // first push the right child to the front
        nodes.push_front(std::make_pair(l_child, child_level)); // then push the left child to the front
        // resulting list: {left child, right child, ...}

      } else {

        // Write the label to the bitset.
        u1 label = get_label(idx);
        u64 to_write = 1 << (tree_height - level) ;// number of tuples represented by the label

        for(auto i = 0; i < to_write; ++i) { // TODO: optimize
          ret[i + write_pointer] = label;
        }
        write_pointer += to_write;
      }
    };

    nodes.push_front(std::make_pair(0,0)); // push the root node to the list

    while(!nodes.empty()){
      /*
      for(auto e : nodes){
        std::cout << e.first << "|" << e.second << "  ,  ";
      }*/
      //std::cout << std::endl;
      // the function is always called for the first node in the list
      //remove the node for which this function is executed
      std::pair<$u64,$u64> node = nodes.front();
      nodes.pop_front();
      decode_tree(node);
    }

    return ret;
  }

  /// Return the size in bytes.
  std::size_t __teb_inline__
  size_in_byte() const noexcept {
    constexpr u64 block_bitlength = sizeof(_block_type) * 8;
    constexpr u64 block_size = sizeof(_block_type);
    $u64 bytes = 0;
    // Tree structure
    bytes += ((structure_.size() + block_bitlength - 1) / block_bitlength) * block_size;
    // Labels
    bytes += ((labels_.size() + block_bitlength - 1) / block_bitlength) * block_size;
    // Rank support
    bytes += rank_.size_in_bytes();
    // Bit-length of the original bitmap.
    bytes += sizeof(n_);
//    if (optimization_level_ > 0) {
      // The number of implicit inner nodes.
      bytes += sizeof(implicit_inner_node_cnt_);
//    }
    // FIXME: + number of tree bits + number of label bits
    return bytes;
  }

  u1 __teb_inline__
  operator!=(teb& other) const noexcept {
    return !(*this == other);
  }

  u1 __teb_inline__
  operator==(teb& other) const noexcept {
    return implicit_inner_node_cnt_ == other.implicit_inner_node_cnt_
        && structure_ == other.structure_
        && labels_ == other.labels_;
  }

  void
  print(std::ostream& os) const noexcept {
    os << "implicit nodes (internal/external) = "
       << implicit_inner_node_cnt_
       << "/"
       << implicit_leaf_node_cnt_
       << ", perfect levels = "
       << determine_perfect_tree_levels(implicit_inner_node_cnt_)
       << ", tree bits = " << structure_.size()
       << ", label bits = " << labels_.size()
       << ", opt level = " << optimization_level_
       << ", n = " << n_
       << ", rank size = " << rank_.size_in_bytes()
       << ", size = " << size_in_byte()
       << "\n | ";

    if (implicit_inner_node_cnt_ > 0) {
//      for ($i64 i = 0; i < implicit_inner_node_cnt_; i++) {
//        os << "1";
//      }
      os << "'";
    }
    for ($i64 i = 0; i < structure_.size(); i++) {
      os << (structure_[i] ? "1" : "0");
    }
//    os << "'";
//    for ($i64 i = implicit_inner_node_cnt_ + structure_.size();
//         i < n_; i++) {
//      os << "0";
//    }
    os << "\n | ";
    for ($i64 i = 0; i < labels_.size(); i++) {
      os << (labels_[i] ? "1" : "0");
    }
  }

  /// Return the name of the implementation.
  static std::string
  name() noexcept {
    return "teb";
  }

//  /// Returns the value of the bit at the position pos.
//  u1 __teb_inline__
//  test(const std::size_t pos) const noexcept {
//    auto n_log2 = dtl::log_2(n_);
//    $u64 node_idx = 0;
//    if (is_leaf_node(node_idx)) { // FIXME: eliminate special case!!!
//      return get_label(node_idx);
//    }
//    for ($u64 i = n_log2 - 1; i < n_log2; i--) {
//      u1 bit = dtl::bits::bit_test(pos, i);
//      auto r = rank(node_idx + 1);
//      node_idx = 2 * r - 1 + bit; // right child if bit is set, left child otherwise
//      if (is_leaf_node(node_idx)) {
//        u64 label_idx = node_idx - rank(node_idx); // FIXME: do not call rank() twice!!!
//        auto label = labels_[label_idx];
//        return label;
//      }
//    }
//    std::cout << "BÄM" << std::endl;
//    std::exit(42);
//  }

  /// Returns true if all bits are set, false otherwise.
  u1 __teb_inline__
  all() noexcept {
    // FIXME: this works only if the tree mask is in a compressed state
    return is_leaf_node(0) // root is the only node (a leaf)
        && get_label(0) == true; // and the label is 1
  }

  /// Returns true if all bits are zero, false otherwise.
  u1 __teb_inline__
  none() noexcept {
    // FIXME: this works only if the tree mask is in a compressed state
    return is_leaf_node(0) // root is the only node (a leaf)
        && get_label(0) == false; // and the label is 0
  }

  /// Returns the length of the original bitmap.
  std::size_t
  size() const noexcept {
    return n_;
  }

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

  /// Computes the height of the tree based on n. - Note that a tree, consisting
  /// of a single (root) node has a height of 0.
  static inline u64
  determine_tree_height(u64 n) noexcept {
    return dtl::log_2(n);
  }

  /// Determines the path to the common ancestor of the two nodes specified
  /// by there paths.
  static inline u64
  determine_common_ancestor_path(u64 src_path, u64 dst_path) {
    //TODO should use positions instead of paths (at least for the second argument)
    assert(src_path != dst_path);
    const auto a = src_path << (dtl::bits::lz_count(src_path) + 1);
    const auto b = dst_path << (dtl::bits::lz_count(dst_path) + 1);
    assert(a < b);
    const auto src_path_len = sizeof(src_path) * 8
        - (dtl::bits::lz_count(src_path) + 1);
    const auto common_prefix_length = dtl::bits::lz_count(a ^ b);
    const auto common_ancestor_path =
        src_path >> (src_path_len - common_prefix_length);
    return common_ancestor_path;
  }

  static inline u64
  determine_level_of(u64 path) {
    const auto lz_cnt_path = dtl::bits::lz_count(path);
    const auto level = sizeof(u64) * 8 - 1 - lz_cnt_path;
    return level;
  }
  //===--------------------------------------------------------------------===//

  //===--------------------------------------------------------------------===//
  /// 1-fill iterator, with skip support.
  class iter {

    /// The fundamental type to encode paths within the tree.
    using path_t = $u64;

    struct stack_entry {
      $u64 node_idx;
      path_t path;
      $u32 rank;
      $u32 level;
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
    const teb& teb_;
    /// The height of the tree structure. // TODO refers to max height
    u64 tree_height_;
    /// The number of perfect levels in the tree structure.
    u64 perfect_levels_;
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
    //===------------------------------------------------------------------===//

  public:

    /// Constructs an iterator for the given TEB instance. After construction,
    /// the iterator points to the first 1-fill.
    explicit __teb_inline__
    iter(const teb& teb) noexcept :
        teb_(teb),
        tree_height_(determine_tree_height(teb.n_)),
        perfect_levels_(
            determine_perfect_tree_levels(teb_.implicit_inner_node_cnt_)),
        top_node_idx_begin_((1ull << (perfect_levels_ - 1)) - 1),
        top_node_idx_end_((1ull << perfect_levels_) - 1),
        top_node_idx_current_(top_node_idx_begin_),
        node_idx_(top_node_idx_current_),
        path_(path_t(1) << (perfect_levels_ - 1))
    {
      // Initialize the stack.
      stack_entry entry; // TODO emplace back
      entry.node_idx = top_node_idx_begin_;
      entry.path = path_t(1) << (perfect_levels_ - 1);
      entry.rank = teb_.rank_inclusive(top_node_idx_begin_);
      // TODO populate level
//      entry.level = 0;
      stack_.push(entry);
      next();
    }

    __teb_inline__
    iter(iter&&) noexcept = default;

    /// Forwards the iterator to the next 1-fill (if any).
    /// Use the functions pos() and length() to get the 1-fill the iterator
    /// is currently pointing to.
    void __teb_inline__
    next() noexcept __attribute__ ((flatten, hot)) {
      D(std::cout << "next()" << std::endl;)
      while (top_node_idx_current_ < top_node_idx_end_) {
        D(std::cout << "next top node" << std::endl;)
outer_loop_begin:
        while (!stack_.empty()) {
          stack_entry node_info = stack_.top();
          D(std::cout << "pop: " << node_info << std::endl;)
          D(std::cout << "rank: " << teb_.rank_inclusive(node_info.node_idx) << std::endl;)
          assert(node_info.rank == teb_.rank_inclusive(node_info.node_idx));
          stack_.pop();

          $u1 label;
          while (teb_.is_inner_node(node_info.node_idx)) {
loop_begin:
            label = false;
            // Determine left and right child. - Both exist, because the tree
            // is full binary.
            // Note: If the current node is an inner node, the rank is always
            //       required.
            u64 right_child_idx =
                (node_info.node_idx < teb_.implicit_inner_node_cnt_)
                  ? 2 * node_info.node_idx + 2
                  : 2 * node_info.rank;
            assert(right_child_idx == teb_.right_child(node_info.node_idx));
            u64 left_child_idx = right_child_idx - 1;
            const auto right_child_is_inner = 0 + teb_.is_inner_node(right_child_idx);
            const auto left_child_is_inner = 0 + teb_.is_inner_node(left_child_idx);
            const auto children_are_inner =
                (left_child_is_inner << 1) | right_child_is_inner; // TODO fetch both bits in one go

            // Eagerly fetch the labels.
            u64 left_child_rank = teb_.rank_inclusive(left_child_idx);
            u64 left_child_label_idx = left_child_idx - left_child_rank
                + left_child_is_inner; // prevent underflow
//            u1 left_child_label = teb_.labels_[left_child_label_idx];

            u64 right_child_label_idx = left_child_label_idx + 1
                - left_child_is_inner; // adjust index if necessary
//            u1 right_child_label = teb_.labels_[right_child_label_idx];

            // Switch over the different cases.
            switch (children_are_inner) {
              case 0b00: {
                //===------------------------------------------------------===//
                // Both childs are leaf nodes.
                //===------------------------------------------------------===//
                // Compute the rank for one child, and derive the rank of the
                // other one.
                // Rank is required to compute the label index.
                u1 left_child_label = teb_.labels_[left_child_label_idx];
                u1 right_child_label = teb_.labels_[right_child_label_idx];
                // Fetch the labels.
                //
                // Note for UN-optimized TEBs the following holds:
                //   One child has a 1-label the other has
                //   a 0-label, which is guaranteed by the bottom-up pruning.
                //   Thus, right_label == !left_label.
                //
                // However, gradual decompression (optimization level >= 2) may
                // expand nodes and replicate labels. Thus, the sub-tree is no
                // longer guaranteed to be 'compressed', and therefore, both
                // labels need to be inspected.
                if (optimization_level_ < 2) {
                  // Go to the node which has the 1-label.
                  node_info.node_idx =
                      left_child_label ? left_child_idx : right_child_idx;
                  node_info.path = (node_info.path << 1) | !left_child_label;
                  node_info.rank = left_child_rank;
                  goto produce_output;
                }
                else {
                  // The price we pay for better compression ratios.
//                  u1 right_child_label = teb_.labels_[right_child_label_idx];
                  u64 both_labels = 2 * left_child_label + right_child_label;
                  switch (both_labels) {
                    case 0b00: {
                      goto outer_loop_begin;
                    }
                    case 0b01: {
                      node_info.node_idx = right_child_idx;
                      node_info.path = (node_info.path << 1) | 1;
                      node_info.rank = left_child_rank + 1; // TODO double check
                      goto produce_output;
                    }
                    case 0b10: {
                      node_info.node_idx = left_child_idx;
                      node_info.path <<= 1;
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
                u1 left_child_label = teb_.labels_[left_child_label_idx];

                // Derive rank of the right child.  The following works, because
                // the left child is leaf (0-bit) and the right is inner (1-bit).
                u64 right_child_rank = left_child_rank + 1; // TODO paper

                if (left_child_label) {
                  // Produce the output for the left child iff it has a 1-label,
                  // otherwise it can be ignored.
                  node_info.node_idx = left_child_idx;
                  label = true; // TODO remove
                  node_info.path <<= 1;
                  node_info.rank = left_child_rank;
                  // Push the right child on the stack.
                  stack_entry right_child_info;
                  right_child_info.node_idx = right_child_idx;
                  right_child_info.path = node_info.path | 1;
                  right_child_info.rank = right_child_rank;
                  stack_.push(right_child_info);
                  // TODO also memorize the rank of the node
                  goto produce_output;
                }
                // Else, go to right child, ignoring the left child.
                node_info.node_idx = right_child_idx;
                node_info.path = (node_info.path << 1) | 1;
                node_info.rank = right_child_rank;
                goto loop_begin;
              }
              case 0b10: {
                //===------------------------------------------------------===//
                // Left child is an inner, right child is a leaf node.
                //===------------------------------------------------------===//

                // Determine whether the right child produces an output
                // (label = 1).
                u1 right_child_label = teb_.labels_[right_child_label_idx];
                if (right_child_label) { // FIXME DEP
                  // Push the right child on the stack iff it has a 1-label,
                  // otherwise the right child is ignored.
                  stack_entry right_child_info;
                  right_child_info.node_idx = right_child_idx;
                  right_child_info.path = (node_info.path << 1) | 1;
                  // Rank of the right child is equal to the rank of the left
                  // child.
                  u64 right_child_rank = left_child_rank; // TODO paper
                  right_child_info.rank = right_child_rank;
                  stack_.push(right_child_info);
                  // TODO also memorize the rank of the node
                }
                // Go to left child.
                node_info.node_idx = left_child_idx;
                node_info.path <<= 1;
                node_info.rank = left_child_rank;
                goto loop_begin;
              }
              case 0b11: {
                //===------------------------------------------------------===//
                // Both childs are an inner nodes.
                //===------------------------------------------------------===//

                u64 right_child_rank = left_child_rank + 1; // TODO paper (maybe)
                // Go to left child.
                node_info.node_idx = left_child_idx;
                node_info.path <<= 1;
                node_info.rank = left_child_rank;
                // Push the right child on the stack.
                stack_entry right_child_info;
                right_child_info.node_idx = right_child_idx;
                right_child_info.path = node_info.path | 1;
                right_child_info.rank = right_child_rank;
                stack_.push(right_child_info);
                goto loop_begin;
              }
              default:
                __builtin_unreachable();
            }
          }
          // Reached a leaf node.
          label = teb_.get_label(node_info.node_idx);
          D(std::cout << "reached leaf " << node_info.node_idx  << std::endl;)
//          D(std::cout << "label = " << label  << std::endl;)
          if (label) {
produce_output:
            // Produce output (a 1-fill).
            const auto lz_cnt_path = dtl::bits::lz_count(node_info.path);
            const auto level = sizeof(path_t) * 8 - 1 - lz_cnt_path;
            // Toggle sentinel bit (= highest bit set) and add offset.
            pos_ = (node_info.path ^ (1ull << level)) << (tree_height_ - level);
            // The length of the 1-fill.
            length_ = teb_.n_ >> level;
            node_idx_ = node_info.node_idx;
            path_ = node_info.path;
            D(std::cout << "return [" << pos_ << ","
                        << (pos_ + length_) << ")"
                        << std::endl;)
            return;
          }
        }

        // TODO bypass the stack here
        ++top_node_idx_current_;
        stack_entry next_node;
        next_node.node_idx = top_node_idx_current_;
        next_node.path = path_t(top_node_idx_current_ - top_node_idx_begin_);
        // Set the sentinel bit.
        next_node.path |= path_t(1) << (perfect_levels_ - 1);
        next_node.rank = teb_.rank_inclusive(top_node_idx_current_);
        stack_.push(next_node);
      }
      pos_ = teb_.n_;
      length_ = 0;
    }

    /// Navigate to the desired position, starting from the trees' root node.
    void __teb_inline__
    nav_from_root_to(const std::size_t to_pos) noexcept {
      //===----------------------------------------------------------------===//
      // (Re-)initialize the iterator state.
      stack_.clear();

      // Determine the top-node path.
      $u64 level = perfect_levels_ - 1;
      const auto foo = to_pos >> (tree_height_ - level);
      path_ = (path_t(1) << level) | foo;

      // Determine the top-node idx.
      top_node_idx_current_ = top_node_idx_begin_ + foo;
      const auto bar = determine_level_of(path_);
      //===----------------------------------------------------------------===//

      node_idx_ = top_node_idx_current_;
      nav_downwards(to_pos);
      return;
//      // Walk down the tree to the desired position.
//      // TODO possibly inline the nav_downwards function to save a few instructions.
//      std::size_t i = tree_height_ - perfect_levels_;
//      while (true) {
//        // First check, if this is already a leaf node.
//        if (teb_.is_leaf_node(node_idx_)) {
//          // Reached the desired position.
//          if (teb_.get_label(node_idx_)) {
//            // Found the corresponding leaf node.
//            // Toggle sentinel bit (= highest bit set) and add offset.
//            pos_ = (path_ ^ (1ull << level)) << (tree_height_ - level);
//            // The length of the 1-fill.
//            length_ = teb_.n_ >> level;
//            // Adjust the current position and fill-length.
//            length_ -= to_pos - pos_;
//            pos_ = to_pos;
//            return;
//          }
//          else {
//            // Search forward to the next 1-fill.
//            next();
//            return;
//          }
//        }
//
//        // Navigate downwards the tree.
//        // 0 -> go to left child, 1 -> go to right child
//        u1 direction_bit = dtl::bits::bit_test(to_pos, i);
//        i--;
//        const auto right_child = teb_.right_child(node_idx_);
//        const auto left_child = right_child - 1;
//        level++;
//        if (!direction_bit) {
//          // Go to left child.
//          stack_.push((right_child << 32) | ((path_ << 1) | 1));
//          path_ <<= 1;
//          node_idx_ = left_child;
//        }
//        else {
//          // Go to right child.
//          path_ = (path_ << 1) | 1;
//          node_idx_ = right_child;
//        }
//      }
    }

    /// Navigate downwards the tree to the desired position, starting from the
    /// current node. The destination position must be part of the current sub-
    /// tree, otherwise the behavior is undefined.
    void __teb_inline__
    nav_downwards(const std::size_t to_pos) noexcept {
      $u64 level = determine_level_of(path_);
      std::size_t i = tree_height_ - level - 1;
      while (true) {
        // First check, if this is already a leaf node.
        if (teb_.is_leaf_node(node_idx_)) {
          // Reached the desired position.
          if (teb_.get_label(node_idx_)) {
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
        const auto right_child_idx = teb_.right_child(node_idx_);
        const auto left_child_idx = right_child_idx - 1;
        level++;
        if (!direction_bit) {
          // Push the right child only if necessary.
          if (teb_.is_inner_node(right_child_idx)
              || teb_.get_label(right_child_idx)) {
            stack_entry right_child_info;
            right_child_info.node_idx = right_child_idx;
            right_child_info.path = (path_ << 1) | 1;
            right_child_info.rank = teb_.rank_inclusive(right_child_idx); // FIXME re-use previously computed rank (get_label())
            stack_.push(right_child_info);
          }
          // Go to left child.
          path_ <<= 1;
          node_idx_ = left_child_idx;
        }
        else {
          // Go to right child.
          path_ = (path_ << 1) | 1;
          node_idx_ = right_child_idx;
        }
      }
    }

    /// Forwards the iterator to the given position. The function first
    /// determines the common ancestor node and then decides whether it is
    /// cheaper to navigate starting from the current node or from the root
    /// node.
    void __teb_inline__
    nav_to(const std::size_t to_pos) {
      assert(to_pos >= pos_ + length_);

      const path_t to_path = to_pos | path_t(1) << tree_height_;
      const path_t from_path = path_;

      const auto common_ancestor_path =
          determine_common_ancestor_path(from_path, to_path);
      const auto right_child_of_common_ancestor_path =
          (common_ancestor_path << 1) | 1ull;

      // The common ancestor must be in the perfect tree part. - The reason is,
      // that the perfect levels are skipped during downward traversal and
      // therefore no nodes from these levels will ever be pushed on the stack.
      // TODO change the condition to something like this:
      // TODO level_of(ca) - (perfect_levels_ - 1) < level_of(current node) - level_of(ca)
      if (determine_level_of(common_ancestor_path) <= perfect_levels_) {
        nav_from_root_to(to_pos);
        return;
      }

      // Walk up the tree to the common ancestor.
      stack_entry node;
      while (path_ != right_child_of_common_ancestor_path) {
        if (stack_.empty()) {
          // Note: It is no longer guaranteed, that the common ancestor is on
          //       the stack since we push only inner nodes and leaf nodes with
          //       1-labels on the stack.
          nav_from_root_to(to_pos);
          return;
        }
        node = stack_.top();
        path_ = node.path;
        stack_.pop();
      }
      node_idx_ = node.node_idx;
      nav_downwards(to_pos);
    }

    /// Fast-forwards the iterator to the given position.
    void __teb_inline__
    skip_to(const std::size_t to_pos) noexcept {
      nav_to(to_pos);
    }

    /// Returns true if the iterator reached the end, false otherwise.
    u1 __teb_inline__
    end() const noexcept {
      return pos_ == teb_.n_;
    }

    /// Returns the starting position of the current 1-fill.
    u64 __teb_inline__
    pos() const noexcept {
      return pos_;
    }

    /// Returns the length of the current 1-fill.
    u64 __teb_inline__
    length() const noexcept {
      return length_;
    }

    /// Returns the path of the current tree node.
    u64 __teb_inline__
    path() const noexcept {
      return path_;
    }

  };

  iter __teb_inline__
  it() const noexcept {
    return std::move(iter(*this));
  }
  //===--------------------------------------------------------------------===//

  /// Returns the name of the instance including the most important parameters
  /// in JSON.
  std::string
  info() const noexcept {
    auto determine_compressed_tree_height = [&]() {
      auto i = it();
      $u64 height = 0;
      while (!i.end()) {
        const auto h = determine_level_of(i.path());
        height = std::max(height, h);
        i.next();
      }
      return height;
    };
    return "{\"name\":\"" + name() + "\""
        + ",\"n\":" + std::to_string(n_)
        + ",\"size\":" + std::to_string(size_in_byte())
        + ",\"tree_bits\":" + std::to_string(structure_.size())
        + ",\"label_bits\":" + std::to_string(labels_.size())
        + ",\"implicit_inner_nodes\":"
          + std::to_string(implicit_inner_node_cnt_)
        + ",\"implicit_leaf_nodes\":"
          + std::to_string(implicit_leaf_node_cnt_)
        + ",\"tree_height\":"
          + std::to_string(determine_compressed_tree_height())
        + ",\"perfect_levels\":"
          + std::to_string(determine_perfect_tree_levels(implicit_inner_node_cnt_))
        + ",\"opt_level\":" + std::to_string(optimization_level_)
        + ",\"rank\":" + rank_.info()
        + "}";
  }

private:

  u1 __teb_inline__
  is_inner_node(u64 node_idx) const {
    const std::size_t implicit_1bit_cnt = implicit_inner_node_cnt_;
    if (node_idx < implicit_1bit_cnt) {
      // Implicit inner node.
      return true;
    }
    if ((node_idx - implicit_1bit_cnt) >= structure_.size()) {
      // Implicit leaf node.
      return false;
    }
    return structure_[node_idx - implicit_1bit_cnt];
  }

  u1 __teb_inline__
  is_leaf_node(u64 node_idx) const {
    return !is_inner_node(node_idx);
  }

  u64 __teb_inline__
  left_child(u64 node_idx) const {
    const std::size_t implicit_1bit_cnt = implicit_inner_node_cnt_;
    $u64 left_child_idx;
    if (node_idx < implicit_1bit_cnt) {
      left_child_idx = 2 * node_idx + 1;
    }
    else {
//      left_child_idx = 2 * (rank(node_idx) + 1) - 1;
      left_child_idx = 2 * rank_inclusive(node_idx) - 1;
    }
    return left_child_idx;
//    return right_child(node_idx) - 1;
  }

  u64 __teb_inline__
  right_child(u64 node_idx) const {
    const std::size_t implicit_1bit_cnt = implicit_inner_node_cnt_;
    $u64 right_child_idx;
    if (node_idx < implicit_1bit_cnt) {
      right_child_idx = 2 * node_idx + 2;
    }
    else {
      right_child_idx = 2 * rank_inclusive(node_idx);
    }
    return right_child_idx;
  }

  std::size_t __teb_inline__
  get_label_idx(u64 node_idx) const {
//    assert(is_leaf_node(node_idx));
    return node_idx - rank_inclusive(node_idx);
//    return node_idx - rank_inclusive(node_idx);
//    return node_idx - rank(node_idx);
//    return node_idx - rank(node_idx + 1);
  }

  u1 __teb_inline__
  get_label(u64 node_idx) const {
    u64 label_idx = get_label_idx(node_idx);
    return labels_[label_idx];
  }

//  u64 __teb_inline__
//  rank(u64 node_idx) const {
//    const std::size_t implicit_1bit_cnt = implicit_inner_node_cnt_;
//    if (node_idx < implicit_1bit_cnt) {
//      return node_idx;
//    }
//    const auto i = std::min(node_idx - implicit_1bit_cnt, structure_.size());
//    const auto r = rank_(i, structure_.m_bits.data());
//    const auto ret_val = implicit_1bit_cnt + r;
//    return ret_val;
//  }

//  u64 __teb_inline__
//  rank_inclusive(u64 node_idx) const {
//    const std::size_t implicit_1bit_cnt = implicit_inner_node_cnt_;
//    if (node_idx < implicit_1bit_cnt) {
//      return node_idx;
//    }
//    const auto i = std::min(node_idx - implicit_1bit_cnt, structure_.size());
//    const auto r = rank_(i, structure_.m_bits.data())
//        + (i < structure_.size() ? structure_[i] : 0); // inclusive rank // TODO optimize
//    const auto ret_val = implicit_1bit_cnt + r;
//    return ret_val;
//  }

  u64 __teb_inline__
  rank_inclusive(u64 node_idx) const {
    const std::size_t implicit_1bit_cnt = implicit_inner_node_cnt_;
    if (node_idx < implicit_1bit_cnt) {
      return node_idx + 1;
    }
    const auto i = std::min(node_idx - implicit_1bit_cnt, structure_.size());
    const auto r = rank_(i, structure_.m_bits.data());
    const auto ret_val = implicit_1bit_cnt + r;
    return ret_val;
  }

};
//===----------------------------------------------------------------------===//

}; // namespace dtl
