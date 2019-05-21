#pragma once

// TODO remove
#define __teb_inline__ __attribute__((noinline))

#if !defined(__teb_inline__)
#if defined(NDEBUG)
// Release build.
#define __teb_inline__ inline __attribute__((always_inline))
//#define __teb_inline__ __attribute__((noinline))
#else
#define __teb_inline__
#endif
#endif

#include <bitset>
#include <list>
#include <queue>
#include <stack>
#include <vector>

#include "boost/dynamic_bitset.hpp"

#include <dtl/dtl.hpp>
#include <dtl/bits.hpp>
#include <dtl/bitmap/util/rank1.hpp>
#include <dtl/bitmap/util/rank1_surf.hpp>
#include <dtl/bitmap/util/rank1_surf_cached.hpp>
#include <dtl/bitmap/util/bit_buffer.hpp>
#ifdef __AVX2__
#include <dtl/bitmap/util/bit_buffer_avx2.hpp>
#endif // __AVX2__
#ifdef __AVX512BW__
#include <dtl/bitmap/util/bit_buffer_avx512.hpp>
#endif // __AVX512BW__
#include <dtl/bitmap/util/bitmap_tree.hpp>
#include <dtl/bitmap/util/binary_tree_structure.hpp>
#include <dtl/static_stack.hpp>
#include <dtl/static_stack2.hpp>
#include <dtl/bitmap/util/rank1_naive.hpp>
#include <dtl/bitmap/util/rank1_interceptor.hpp>
#include <dtl/bitmap/util/rank1_super_fast.hpp>
#include <dtl/bitmap/util/bitmap_view.hpp>
#include <dtl/math.hpp>

#include "teb_scan_util.hpp"

namespace dtl {
//===----------------------------------------------------------------------===//
/// Encodes a bitmap of length n as a binary tree.
///
/// The implementation supports the optimization levels 0, 1, and 2.
///   0 = implements the core idea of TEBs
///   1 = uses the implicit inner node optimization
///   2 = additionally applies the gradual decompression optimization.
///   3 = uses the implicit labels optimization.
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
static constexpr i32 teb_default_opt_level = 3;
template<i32 optimization_level_ = teb_default_opt_level>
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

  bitmap_view<word_type> T_;
  $u64 structure_bit_cnt_;
  bitmap_view<word_type> L_;
  $u64 label_bit_cnt_;

  /// Support data structure for rank1 operations on the tree structure.
  static constexpr u1 inclusive = true;
  using rank_support = dtl::rank1_surf<word_type, inclusive>;
  rank_support rank_;

  /// The number of implicit inner nodes in the tree structure.
  $u32 implicit_inner_node_cnt_;

  /// For testing purposes only.
  $u32 implicit_leaf_node_cnt_;

  /// The number of implicit leading 0-labels.
  $u32 implicit_leading_label_cnt_ = 0;

  /// For testing purposes only.
  $u32 implicit_trailing_label_cnt_ = 0;

  std::array<std::size_t, 32> level_offsets_structure_;
  std::array<std::size_t, 32> level_offsets_labels_;

  $u32 encoded_tree_height_ = 0;

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
    if (optimization_level_ > 2) {
      // Implicit labels are only considered with -o3 or higher.
      implicit_leading_label_cnt_ = bitmap_tree.get_leading_0label_cnt();
      implicit_trailing_label_cnt_ = bitmap_tree.get_trailing_0label_cnt();
    }
    std::size_t node_cntr = 0;
    std::size_t leaf_node_cntr = 0;
    std::size_t current_level = ~0ull;
    for (auto it = bitmap_tree.breadth_first_begin();
         it != bitmap_tree.breadth_first_end();
         ++it) {
      ++node_cntr;

      u64 idx = (*it).idx;
      u64 level = (*it).level;

      if (current_level != level) {
        level_offsets_structure_[level] = node_cntr - 1;
        level_offsets_labels_[level] = leaf_node_cntr;
        current_level = level;
      }

      // Omit the implicit nodes.
      if (node_cntr <= implicit_inner_node_cnt_) {
        continue;
      }

      u1 is_inner = bitmap_tree.is_inner_node(idx);
      leaf_node_cntr += !is_inner;

      // Emit a 1-bit if the current node is an inner node, a 0-bit otherwise.
      if (is_inner
          || optimization_level_ == 0
          || (optimization_level_ > 0
              && idx <= bitmap_tree.get_last_explicit_node_idx())) {
        structure_.push_back(is_inner);
      }
      if (!is_inner) {
        if (optimization_level_ > 2) {
          // Add the label of the leaf node (if necessary).
          if (idx >= bitmap_tree.get_first_node_idx_with_1label()
              && idx <= bitmap_tree.get_last_node_idx_with_1label()) {
            labels_.push_back(bitmap_tree.label_of_node(idx));
          }
        }
        else {
          // Add the label of the leaf node.
          labels_.push_back(bitmap_tree.label_of_node(idx));
        }
      }
    }

    // TODO do we really need to store this?
    encoded_tree_height_ = current_level + 1;

    // Hack: ensure structure is not null
    if (structure_.size() == 0) {
      structure_.push_back(false);
    }

    // Hack: add an additional word to the tree structure, which helps to
    //       eliminate a branch in the scan iterator.
//    for (std::size_t i = 0; i < (sizeof(word_type) * 8); ++i) {
//      structure_.push_back(false);
//    }

    // Init rank1 support data structure.
    rank_.init(structure_);

    T_.init(data_view<const word_type> {
        structure_.m_bits.data(),
        structure_.m_bits.data() + structure_.m_bits.size(),
    });
    structure_bit_cnt_ = structure_.size();
    L_.init(data_view<const word_type> {
        labels_.m_bits.data(),
        labels_.m_bits.data() + labels_.m_bits.size(),
    });
    label_bit_cnt_ = labels_.size();
  }

  teb(const teb& other) = default;
  teb(teb&& other) noexcept = default;
  teb& operator=(const teb& other) = default;
  teb& operator=(teb&& other) = default;
  ~teb() = default;

  /// Return the size in bytes.
  std::size_t __teb_inline__
  size_in_byte() const noexcept {
    constexpr u64 block_bitlength = sizeof(_block_type) * 8;
    constexpr u64 block_size = sizeof(_block_type);
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
    const auto perfect_levels =
        determine_perfect_tree_levels(implicit_inner_node_cnt_);
    const auto encoded_tree_height = dtl::log_2(n_) + 1; // FIXME could be lower, but its unlikely
    assert(encoded_tree_height >= perfect_levels);
    bytes += (4 + 4) * (encoded_tree_height - perfect_levels);

    // Padding. We want T to be 8-byte aligned.
    bytes += 8 - (bytes % 8);

    // Tree structure
    bytes += ((structure_.size() + block_bitlength - 1) / block_bitlength)

        * block_size;
    // Rank helper structure
    bytes += rank_.size_in_bytes();

    // Labels
    bytes += ((labels_.size() + block_bitlength - 1) / block_bitlength)
        * block_size;

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
       << ", implicit labels (leading/trailing) = "
       << implicit_leading_label_cnt_
       << "/"
       << implicit_trailing_label_cnt_
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
      os << "'";
    }
    for ($i64 i = 0; i < structure_.size(); i++) {
      os << (structure_[i] ? "1" : "0");
    }
    os << "\n | ";
    if (implicit_leading_label_cnt_ > 0) {
      os << "'";
    }
    for ($i64 i = 0; i < labels_.size(); i++) {
      os << (labels_[i] ? "1" : "0");
    }
  }

  /// Return the name of the implementation.
  static std::string
  name() noexcept {
    if (optimization_level_ == teb_default_opt_level) {
      return "teb";
    }
    else {
      return "teb_o" + std::to_string(optimization_level_);
    }
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
  determine_common_ancestor_path(u64 src_path, u64 dst_path) noexcept {
    //TODO should use positions instead of paths (at least for the second argument)
    assert(src_path != dst_path);
    const auto t0 = dtl::bits::lz_count(src_path) + 1;
    const auto a = src_path << t0;
    const auto b = dst_path << (dtl::bits::lz_count(dst_path) + 1);
    assert(a < b);
    const auto src_path_len = sizeof(src_path) * 8 - t0;
    const auto common_prefix_length = dtl::bits::lz_count(a ^ b);
    const auto common_ancestor_path =
        src_path >> (src_path_len - common_prefix_length);
    return common_ancestor_path;
  }

  static inline void
  determine_common_ancestor_path(u64 src_path, u64 dst_path,
      $u64& out_common_ancestor_path, $u64& out_common_ancestor_level) noexcept {
    //TODO should use positions instead of paths (at least for the second argument)
    assert(src_path != dst_path);
    const auto t0 = dtl::bits::lz_count(src_path) + 1;
    const auto a = src_path << t0;
    const auto b = dst_path << (dtl::bits::lz_count(dst_path) + 1);
    assert(a < b);
    const auto src_path_len = sizeof(src_path) * 8 - t0;
    const auto common_prefix_length = dtl::bits::lz_count(a ^ b);
    out_common_ancestor_path =
        src_path >> (src_path_len - common_prefix_length);
    out_common_ancestor_level = common_prefix_length;
  }

  static inline void
  determine_common_ancestor_path2(u64 src_path, u64 dst_pos, u64 tree_height,
      $u64& out_common_ancestor_path, $u64& out_common_ancestor_level) noexcept {
    const auto t0 = dtl::bits::lz_count(src_path) + 1;
    const auto a = src_path << t0;
    const auto b = dst_pos << (sizeof(dst_pos) * 8 - tree_height);
    assert(a < b);
    const auto src_path_len = sizeof(src_path) * 8 - t0;
    const auto common_prefix_length = dtl::bits::lz_count(a ^ b);
    out_common_ancestor_path =
        src_path >> (src_path_len - common_prefix_length);
    out_common_ancestor_level = common_prefix_length;
  }

  static inline u64
  determine_level_of(u64 path) noexcept {
    const auto lz_cnt_path = dtl::bits::lz_count(path);
    const auto level = sizeof(u64) * 8 - 1 - lz_cnt_path;
    return level;
  }
  //===--------------------------------------------------------------------===//

  //===--------------------------------------------------------------------===//
  /// 1-fill iterator, with skip support.
  class iter;

  iter __teb_inline__
  it() const noexcept {
    return std::move(iter(*this));
  }
  //===--------------------------------------------------------------------===//

  //===--------------------------------------------------------------------===//
  /// Scan iterator, with INEFFICIENT skip support.
  class scan_iter;

  scan_iter __teb_inline__
  scan_it() const noexcept {
    return std::move(scan_iter(*this));
  }
  //===--------------------------------------------------------------------===//

  using scan_iter_type = scan_iter;
  using skip_iter_type = iter;

  /// Returns the name of the instance including the most important parameters
  /// in JSON.
  std::string
  info() const noexcept {
//    auto determine_compressed_tree_height = [&]() {
//      auto i = it();
//      $u64 height = 0;
//      while (!i.end()) {
//        const auto h = determine_level_of(i.path());
//        height = std::max(height, h);
//        i.next();
//      }
//      return height;
//    };
    return "{\"name\":\"" + name() + "\""
        + ",\"n\":" + std::to_string(n_)
        + ",\"size\":" + std::to_string(size_in_byte())
        + ",\"tree_bits\":" + std::to_string(structure_.size())
        + ",\"label_bits\":" + std::to_string(labels_.size())
        + ",\"implicit_inner_nodes\":"
          + std::to_string(implicit_inner_node_cnt_)
        + ",\"implicit_leaf_nodes\":"
          + std::to_string(implicit_leaf_node_cnt_)
//        + ",\"tree_height\":"
//          + std::to_string(determine_compressed_tree_height())
        + ",\"perfect_levels\":"
          + std::to_string(determine_perfect_tree_levels(implicit_inner_node_cnt_))
        + ",\"opt_level\":" + std::to_string(optimization_level_)
        + ",\"rank\":" + rank_.info()
        + ",\"leading_zero_labels\":" + std::to_string(implicit_leading_label_cnt_)
        + ",\"trailing_zero_labels\":" + std::to_string(implicit_trailing_label_cnt_)
        + "}";
  }

private:

  u1 __teb_inline__
  is_inner_node(u64 node_idx) const noexcept {
    const std::size_t implicit_1bit_cnt = implicit_inner_node_cnt_;
    const std::size_t implicit_leaf_begin =
        structure_bit_cnt_ + implicit_inner_node_cnt_;
    if (node_idx < implicit_1bit_cnt) {
      // Implicit inner node.
      return true;
    }
    if (node_idx >= implicit_leaf_begin) {
      // Implicit leaf node.
      return false;
    }
    return T_[node_idx - implicit_1bit_cnt];
  }

  u1 __teb_inline__
  is_inner_node_branchfree(u64 node_idx) const noexcept { // SLOW!!!
    const std::size_t implicit_1bit_cnt = implicit_inner_node_cnt_;
    auto node_offset = node_idx - implicit_1bit_cnt;
    const std::size_t implicit_leaf_begin =
        structure_bit_cnt_ + implicit_inner_node_cnt_;
    u1 is_implicit_inner = node_idx < implicit_1bit_cnt;
    u1 is_implicit_leaf = node_idx >= implicit_leaf_begin;

    u1 is_implicit = is_implicit_inner | is_implicit_leaf;
    u1 implicit_ret_val = is_implicit_inner;
    node_offset -= (is_implicit) * node_offset;
    u1 read_val = T_[node_offset] & !is_implicit;

    u1 ret_val = is_implicit
        ? implicit_ret_val
        : read_val;
    return ret_val;
  }

  u1 __teb_inline__
  is_leaf_node(u64 node_idx) const noexcept {
    return !is_inner_node(node_idx);
  }

  u64 __teb_inline__
  left_child(u64 node_idx) const noexcept {
    return 2 * rank_inclusive(node_idx) - 1;
  }

  u64 __teb_inline__
  right_child(u64 node_idx) const noexcept {
    return 2 * rank_inclusive(node_idx);
  }

  std::size_t __teb_inline__
  get_label_idx(u64 node_idx) const noexcept {
    return node_idx - rank_inclusive(node_idx);
  }

  u1 __teb_inline__
  get_label_by_idx(u64 label_idx) const noexcept {
    if (optimization_level_ > 2) {
      const auto implicit_leading_label_cnt = implicit_leading_label_cnt_;
      const std::size_t implicit_trailing_0labels_begin =
          label_bit_cnt_ + implicit_leading_label_cnt;
      if (label_idx < implicit_leading_label_cnt) {
        // An implicit leading 0-label.
        return false;
      }
      if (label_idx >= implicit_trailing_0labels_begin) {
        // An implicit trailing 0-label.
        return false;
      }
      return L_[label_idx - implicit_leading_label_cnt];
    }
    else {
      return labels_[label_idx];
    }
  }

  u1 __teb_inline__
  get_label(u64 node_idx) const noexcept {
    u64 label_idx = get_label_idx(node_idx);
    return get_label_by_idx(label_idx);
  }

  u64 __teb_inline__
  rank_inclusive(u64 node_idx) const noexcept {
    const std::size_t implicit_1bit_cnt = implicit_inner_node_cnt_;
    if (node_idx < implicit_1bit_cnt) {
      return node_idx + 1;
    }
    assert(structure_bit_cnt_ > 0);
    const auto i = std::min(node_idx - implicit_1bit_cnt, structure_bit_cnt_ - 1);
    const auto r = rank_(i, T_.data_.begin());
    const auto ret_val = implicit_1bit_cnt + r;
    return ret_val;
  }

};
//===----------------------------------------------------------------------===//
template<i32 optimization_level_>
class teb<optimization_level_>::iter {

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
      partition_shift_(tree_height_ - (perfect_levels_ - 1)),
      top_node_idx_begin_((1ull << (perfect_levels_ - 1)) - 1),
      top_node_idx_end_((1ull << perfect_levels_) - 1),
      top_node_idx_current_(top_node_idx_begin_),
      node_idx_(top_node_idx_current_),
      path_(path_t(1) << (perfect_levels_ - 1))
  {
    // Initialize the stack.
    stack_entry& entry = stack_.push();
    entry.node_idx = top_node_idx_begin_;
    entry.path = path_t(1) << (perfect_levels_ - 1);
    entry.rank = teb_.rank_inclusive(top_node_idx_begin_);
    entry.level = perfect_levels_ - 1;
    next();
  }

  __teb_inline__
  iter(iter&&) noexcept = default;

  /// Forwards the iterator to the next 1-fill (if any).
  /// Use the functions pos() and length() to get the 1-fill the iterator
  /// is currently pointing to.
  void __teb_inline__
  next() noexcept __attribute__ ((flatten, hot)) {
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
          const auto right_child_is_inner = 0 + teb_.is_inner_node(right_child_idx);
          const auto left_child_is_inner = 0 + teb_.is_inner_node(left_child_idx);
          const auto children_are_inner =
              (left_child_is_inner << 1) | right_child_is_inner; // TODO fetch both bits in one go

          // Compute the rank for one child, and derive the rank of the
          // other one.
          // Rank is required to compute the label index.
          u64 left_child_rank = teb_.rank_inclusive(left_child_idx);
          u64 left_child_label_idx = left_child_idx - left_child_rank
              + left_child_is_inner; // prevent underflow
          u64 right_child_label_idx = left_child_label_idx + 1
              - left_child_is_inner; // adjust index if necessary

//          // TODO Eagerly fetch the labels.
//          u1 left_child_label = teb_.L_[left_child_label_idx];
//          u1 right_child_label = teb_.L_[right_child_label_idx];



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
                node_info.node_idx =
                    left_child_label ? left_child_idx : right_child_idx;
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
//      next_node.level = determine_level_of(next_node.path);
      next_node.level = perfect_levels_ - 1;
      next_node.rank = teb_.rank_inclusive(top_node_idx_current_);
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
    //===----------------------------------------------------------------===//

    node_idx_ = top_node_idx_current_;
    nav_downwards(to_pos);
    return;
  }

  /// Navigate downwards the tree to the desired position, starting from the
  /// current node. The destination position must be part of the current sub-
  /// tree, otherwise the behavior is undefined.
  void __teb_inline__
  nav_downwards(const std::size_t to_pos) noexcept {
    $u64 level = determine_level_of(path_);
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
    determine_common_ancestor_path2(from_path, to_pos, tree_height_,
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
    const auto right_child_of_common_ancestor_path =
        (common_ancestor_path << 1) | 1ull;
    const auto right_child_of_common_ancestor_level = common_ancestor_level + 1;

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
//      stack_.rewind(i);
//      next();
//      return;
//    }
//
//    if (found) {
//      stack_.rewind(i - 1);
//      node_idx_ = stack_[i - 1].node_idx;
//      path_ = path;
//      nav_downwards(to_pos);
//      return;
//    }
//
//    nav_from_root_to(to_pos);
//    return;
// --- slower ---

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
    if (to_pos >= teb_.n_) {
      pos_ = teb_.n_;
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
  u1 __forceinline__ //__teb_inline__
  end() const noexcept {
    return pos_ == teb_.n_;
  }

  /// Returns the starting position of the current 1-fill.
  u64 __forceinline__ //__teb_inline__
  pos() const noexcept {
    return pos_;
  }

  /// Returns the length of the current 1-fill.
  u64 __forceinline__ //__teb_inline__
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
    return determine_level_of(path_);
  }

  /// Returns the number of perfect tree levels.
  u64 __teb_inline__
  perfect_levels() const noexcept {
    return perfect_levels_;
  }

};
//===----------------------------------------------------------------------===//
/// 1-fill iterator, with skip support.
template<i32 optimization_level_>
class teb<optimization_level_>::scan_iter {

  /// The fundamental type to encode paths within the tree.
  using path_t = $u64;

  static constexpr u64 DEFAULT_BATCH_SIZE = 128;

  struct range_t {
    $u64 pos;
    $u64 length;
  };

  struct scanner_state_t {
    $u64 node_idx_;
    $u64 label_idx_;

    void __teb_inline__
    init(u64 node_idx, u64 label_offset) {
      node_idx_ = node_idx;
      label_idx_ = label_offset;
    }

    void
    print(std::ostream& os) const noexcept {
      os << "scan_state["
         << "node_idx=" << node_idx_
         << ",label_idx=" << label_idx_ << "]";
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

  /// One scanner (or stream) per level.
  std::array<scanner_state_t, 32> scanner_states_;

  /// A batch of results.
  std::vector<range_t> results_;
  $u64 result_cnt_;
  u64 batch_size_;
  $u64 result_read_pos_;

  $u64 scan_path_;
  $u64 scan_path_level_;

  $u32 alpha_;

  std::size_t first_1label_idx_ = 0;
  //===--------------------------------------------------------------------===//


  //===--------------------------------------------------------------------===//
  // Helper functions.
  //===--------------------------------------------------------------------===//
  auto __teb_inline__
  determine_level(u32 alpha) {
    // LSB is always set. Thus, lz_count(alpha) is never undefined.
    return sizeof(alpha) * 8 - dtl::bits::lz_count(alpha);
  };

  auto __teb_inline__
  is_left_child(path_t path) {
    return (path & path_t(1)) == 0;
  };
  //===--------------------------------------------------------------------===//

public:

  /// Constructs an iterator for the given TEB instance. After construction,
  /// the iterator points to the first 1-fill.
  explicit __teb_inline__
  scan_iter(const teb& teb, u64 batch_size = DEFAULT_BATCH_SIZE) noexcept
      : teb_(teb),
        tree_height_(determine_tree_height(teb.n_)),
        perfect_levels_(
            determine_perfect_tree_levels(teb_.implicit_inner_node_cnt_)),
        top_node_idx_begin_((1ull << (perfect_levels_ - 1)) - 1),
        top_node_idx_end_((1ull << perfect_levels_) - 1),
        top_node_idx_current_(top_node_idx_begin_),
        results_(batch_size),
        result_cnt_(0),
        result_read_pos_(0),
        batch_size_(batch_size),
        alpha_(0) {
    results_.reserve(batch_size);
    if (teb_.encoded_tree_height_ == 1) {
      if (teb_.labels_[0]) {
        results_[0].pos = 0;
        results_[0].length = teb_.n_;
        results_[1].pos = teb_.n_;
        results_[1].length = 0;
        result_cnt_ = 2;
      }
      else {
        results_[0].pos = teb_.n_;
        results_[0].length = 0;
        result_cnt_ = 1;
      }
      return;
    }
    // Initialize the scanners.
    for (std::size_t level = 0; level < teb_.encoded_tree_height_; ++level) {
      const auto node_offset = teb_.level_offsets_structure_[level];
      const auto label_offset = teb_.level_offsets_labels_[level];
      scanner_states_[level].init(node_offset, label_offset);
    }
    // Initialize path and level variable.
    scan_path_ = path_t(0);
    $u32 level = 0;
    for (; level < teb_.encoded_tree_height_; ++level) {
      u64 node_idx = teb_.level_offsets_structure_[level];
      const auto node_bit = teb_.is_inner_node(node_idx);
      if (!node_bit) break;
    }
    scan_path_level_ = level;
    // Iterator is now positioned at the first leaf node.
#ifdef __AVX512BW__
    next_batch_avx512();
#else
#ifdef __AVX2__
    next_batch_avx2();
#else
    next_batch();
#endif
#endif
  }

  __teb_inline__
  scan_iter(scan_iter&&) noexcept = default;
  scan_iter(const scan_iter& other) = default;
  scan_iter& operator=(const scan_iter& other) = default;
  scan_iter& operator=(scan_iter&& other) = default;
  ~scan_iter() = default;

  /// Forwards the iterator to the next 1-fill (if any).
  /// Use the functions pos() and length() to get the 1-fill the iterator
  /// is currently pointing to.
  void __teb_inline__
  next() noexcept __attribute__ ((flatten, hot)) {
    assert(!end());
    ++result_read_pos_;
    if (result_read_pos_ == result_cnt_) {
      result_read_pos_ = 0;
      result_cnt_ = 0;
#ifdef __AVX512BW__
      next_batch_avx512();
#else
#ifdef __AVX2__
      next_batch_avx2();
#else
      next_batch();
#endif
#endif
    }
  }

  void
  next_batch() noexcept __attribute__((flatten, hot, noinline)) {
    const auto h = tree_height_;
    const auto n = teb_.size();
    // Note: The path variable does NOT contain a sentinel bit. Instead, we
    //       keep track of the paths' level using a separate variable.
    register auto path = scan_path_;
    register auto path_level = scan_path_level_;
    auto pos = path << (h - path_level);
    auto length = n >> path_level;
    auto result_cnt = result_cnt_;
    auto alpha = 0u;

    const data_view<const word_type> T {
        teb_.structure_.m_bits.data(),
        teb_.structure_.m_bits.data() + teb_.structure_.m_bits.size(),
    };
    const data_view<const word_type> L {
        teb_.labels_.m_bits.data(),
        teb_.labels_.m_bits.data() + teb_.labels_.m_bits.size(),
    };

    // Bit buffers for the tree structure.
    dtl::bit_buffer<8> t_bb0;
    dtl::bit_buffer<8> t_bb1;
    dtl::bit_buffer<8> t_bb2;
    dtl::bit_buffer<8> t_bb3;

    // Bit buffers for the labels.
    dtl::bit_buffer<8> l_bb0;
    dtl::bit_buffer<8> l_bb1;
    dtl::bit_buffer<8> l_bb2;
    dtl::bit_buffer<8> l_bb3;

    // Perfect levels. (all bits set to 1)
    for (std::size_t level = 0; level < perfect_levels_ - 1; ++level) {
      u64 buffer_idx = level / 8;
      u64 slot_idx = level % 8;
      switch (buffer_idx) {
        case 0: t_bb0.set(slot_idx, ~0ul); l_bb0.set(slot_idx, ~0ul); break;
        case 1: t_bb1.set(slot_idx, ~0ul); l_bb1.set(slot_idx, ~0ul); break;
        case 2: t_bb2.set(slot_idx, ~0ul); l_bb2.set(slot_idx, ~0ul); break;
        case 3: t_bb3.set(slot_idx, ~0ul); l_bb3.set(slot_idx, ~0ul); break;
      }
    }
    // --

    // Advance the tree structure bit buffers.
    auto advance_tree_scanners = [&](std::size_t level_begin,
                                     std::size_t level_end) {
      static constexpr u64 p = (1ul << 8) - 1;
      u64 a = (~0ul) << level_begin;
      u64 b = (~0ul) >> (64 - level_end);
      u64 m = a & b;
      switch (level_begin) {
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
          t_bb0.increment((m >> (8 * 0)) & p);
        case 8:
        case 9:
        case 10:
        case 11:
        case 12:
        case 13:
        case 14:
        case 15:
          t_bb1.increment((m >> (8 * 1)) & p);
        case 16:
        case 17:
        case 18:
        case 19:
        case 20:
        case 21:
        case 22:
        case 23:
          t_bb2.increment((m >> (8 * 2)) & p);
        default:
          t_bb3.increment((m >> (8 * 3)) & p);
      }
      // Update the alpha vector.
      u64 r0 = t_bb0.read();
      u64 r1 = t_bb1.read();
      u64 r2 = t_bb2.read();
      u64 r3 = t_bb3.read();
      alpha = (r3 << (8 * 3)) | (r2 << (8 * 2)) | (r1 << (8 * 1)) | r0;
    };

    // Advance the label bit buffers.
    auto advance_label_scanner = [&](u32 level) {
      u64 buffer_idx = level / 8;
      u64 slot_idx = level % 8;
      u64 slot_mask = 1ul << slot_idx;
      $u64 label = 0;
      switch (buffer_idx) {
        case 0: l_bb0.increment(slot_mask); break;
        case 1: l_bb1.increment(slot_mask); break;
        case 2: l_bb2.increment(slot_mask); break;
        case 3: l_bb3.increment(slot_mask); break;
      }
    };

    while (result_cnt < batch_size_ - 8) {

      // Initialize the bit buffers.
      t_bb0.reset_read_mask();
      t_bb1.reset_read_mask();
      t_bb2.reset_read_mask();
      t_bb3.reset_read_mask();
      l_bb0.reset_read_mask();
      l_bb1.reset_read_mask();
      l_bb2.reset_read_mask();
      l_bb3.reset_read_mask();
      for (std::size_t level = perfect_levels_ - 1;
           level < teb_.encoded_tree_height_; ++level) {
        u64 buffer_idx = level / 8;
        u64 slot_idx = level % 8;
        i64 delta = - static_cast<i64>(teb_.implicit_inner_node_cnt_);
        switch (buffer_idx) {
          case 0: t_bb0.set(slot_idx, fetch_n_bits(T, static_cast<i64>(scanner_states_[level].node_idx_) + delta, 8)); l_bb0.set(slot_idx, fetch_n_bits(L, static_cast<i64>(scanner_states_[level].label_idx_), 8)); break;
          case 1: t_bb1.set(slot_idx, fetch_n_bits(T, static_cast<i64>(scanner_states_[level].node_idx_) + delta, 8)); l_bb1.set(slot_idx, fetch_n_bits(L, static_cast<i64>(scanner_states_[level].label_idx_), 8)); break;
          case 2: t_bb2.set(slot_idx, fetch_n_bits(T, static_cast<i64>(scanner_states_[level].node_idx_) + delta, 8)); l_bb2.set(slot_idx, fetch_n_bits(L, static_cast<i64>(scanner_states_[level].label_idx_), 8)); break;
          case 3: t_bb3.set(slot_idx, fetch_n_bits(T, static_cast<i64>(scanner_states_[level].node_idx_) + delta, 8)); l_bb3.set(slot_idx, fetch_n_bits(L, static_cast<i64>(scanner_states_[level].label_idx_), 8)); break;
        }
      }

      // Initialize alpha. (which reads the first bits in the bit buffers)
      {
        u64 r0 = t_bb0.read();
        u64 r1 = t_bb1.read();
        u64 r2 = t_bb2.read();
        u64 r3 = t_bb3.read();
        alpha = (r3 << (8 * 3)) | (r2 << (8 * 2)) | (r1 << 8) | r0;
      }

      // Do multiple iteration to consume the bit buffers.
      for ($i32 b = 0; b < 7; ++b) {
        u64 length = n >> path_level;

        assert(pos < n);
        assert(pos + length <= n);
        assert(path_level >= 1);
        assert(length > 0);

        // Read the label of the current leaf node.
        u64 buffer_idx = path_level / 8;
        u64 slot_idx = path_level % 8;
        u64 label_mask = 1ul << slot_idx;
        $u64 label = 0;
        switch (buffer_idx) {
          case 0: label = l_bb0.read(label_mask) != 0; break;
          case 1: label = l_bb1.read(label_mask) != 0; break;
          case 2: label = l_bb2.read(label_mask) != 0; break;
          case 3: label = l_bb3.read(label_mask) != 0; break;
        }
        // Produce output (a 1-fill).
        results_[result_cnt].pos = pos;
        results_[result_cnt].length = length;
        result_cnt += label;


        // Increment the current position.
        pos += length;
        if (pos == n) goto done; // TODO eliminate branch

        // Walk upwards until a left child is found. (might be the current one)
        const auto advance_end = path_level + 1;
        const auto up_steps = dtl::bits::tz_count(~path);
        path >>= up_steps;
        path_level -= up_steps;
        // Go to right sibling.
        path = path | 1;

        // Advance the scanners and update the alpha vector.
        const auto advance_begin = path_level;
        advance_tree_scanners(advance_begin, advance_end);
        advance_label_scanner(advance_end - 1);
        // Walk downwards to the left-most leaf in that sub-tree.
        const auto down_steps = dtl::bits::tz_count(~(alpha >> path_level));
        path_level += down_steps;
        path <<= down_steps;
      }

      // Update the scanner states.
      auto get_tree_scanner_read_pos = [&](u32 level) {
        u64 buffer_idx = level / 8;
        u64 slot_idx = level % 8;
        $i64 ret_val = 0;
        switch (buffer_idx) {
          case 0: ret_val = t_bb0.get_read_pos(slot_idx); break;
          case 1: ret_val = t_bb1.get_read_pos(slot_idx); break;
          case 2: ret_val = t_bb2.get_read_pos(slot_idx); break;
          case 3: ret_val = t_bb3.get_read_pos(slot_idx); break;
        }
        return ret_val;
      };
      auto get_label_scanner_read_pos = [&](u32 level) {
        u64 buffer_idx = level / 8;
        u64 slot_idx = level % 8;
        $i64 ret_val = 0;
        switch (buffer_idx) {
          case 0: ret_val = l_bb0.get_read_pos(slot_idx); break;
          case 1: ret_val = l_bb1.get_read_pos(slot_idx); break;
          case 2: ret_val = l_bb2.get_read_pos(slot_idx); break;
          case 3: ret_val = l_bb3.get_read_pos(slot_idx); break;
        }
        return ret_val;
      };
      for ($u32 level = perfect_levels_ - 1;
           level < teb_.encoded_tree_height_; ++level) {
        scanner_states_[level].node_idx_ += get_tree_scanner_read_pos(level);
        scanner_states_[level].label_idx_ += get_label_scanner_read_pos(level);
      }
    }

done:
    if (result_cnt < batch_size_
        && pos == n) {
      results_[result_cnt].pos = n;
      results_[result_cnt].length = 0;
      ++result_cnt;
      length = 0;
    }
    scan_path_ = path;
    scan_path_level_ = path_level;
    result_cnt_ = result_cnt;
    alpha_ = alpha;
  }

#ifdef __AVX2__
  void
  next_batch_avx2() noexcept __attribute__((flatten, hot, noinline)) {
    const auto h = tree_height_;
    const auto n = teb_.size();
    // Note: The path variable does NOT contain a sentinel bit. Instead, we
    //       keep track of the paths' level using a separate variable.
    register auto path = scan_path_;
    register auto path_level = scan_path_level_;
    auto pos = path << (h - path_level);
    auto length = n >> path_level;
    auto result_cnt = result_cnt_;
    auto alpha = 0u;
    auto alpha_la = 0u;
    auto beta = 0u;

    dtl::r256 tmp_t { .u8 = {
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0 }};
    dtl::r256 tmp_l { .u8 = {
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0 }};

    const data_view<const word_type> T {
        teb_.structure_.m_bits.data(),
        teb_.structure_.m_bits.data() + teb_.structure_.m_bits.size(),
    };
    const data_view<const word_type> L {
        teb_.labels_.m_bits.data(),
        teb_.labels_.m_bits.data() + teb_.labels_.m_bits.size(),
    };

    // Bit buffers for the tree structure.
    dtl::bit_buffer_avx2<> tree_bit_buffer;

    // Bit buffers for the labels.
    dtl::bit_buffer_avx2<> label_bit_buffer;

    // Perfect levels. (all bits set to 1)
    for (std::size_t level = 0; level < perfect_levels_ - 1; ++level) {
      const auto slot_idx = level;
      tmp_t.u8[slot_idx] = ~u8(0);
      tmp_l.u8[slot_idx] = ~u8(0);
    }
    tree_bit_buffer.set_raw(tmp_t.i);
    label_bit_buffer.set_raw(tmp_l.i);

    assert(batch_size_ > 8);
    while (result_cnt < batch_size_ - 8
        && pos < n) {
      // Initialize the bit buffers.
      tree_bit_buffer.reset_read_mask();
      label_bit_buffer.reset_read_mask();
      for (std::size_t level = perfect_levels_ - 1;
           level < teb_.encoded_tree_height_; ++level) {
        const auto slot_idx = level;
        i64 delta = - static_cast<i64>(teb_.implicit_inner_node_cnt_);
        tmp_t.u8[slot_idx] = fetch_n_bits(T, static_cast<i64>(scanner_states_[level].node_idx_) + delta, 8);
        tmp_l.u8[slot_idx] = fetch_n_bits(L, static_cast<i64>(scanner_states_[level].label_idx_), 8);
      }
      tree_bit_buffer.set_raw(tmp_t.i);
      label_bit_buffer.set_raw(tmp_l.i);

      // Initialize alpha. (which reads the first bits in the bit buffers)
      alpha = tree_bit_buffer.read();
      alpha_la = tree_bit_buffer.read_ahead();

      // Initialize beta.
      beta = label_bit_buffer.read();


      // Do multiple iteration to consume the bit buffers.
      $u1 is_not_done = true;
      for ($i32 b = 0; b < 6 && is_not_done; ++b) {
        u64 length = n >> path_level;

        assert(pos < n);
        assert(pos + length <= n);
        assert(path_level >= 1);
        assert(length > 0);

        // Produce output (a 1-fill).
        results_[result_cnt].pos = pos;
        results_[result_cnt].length = length;
        // Read the label of the current leaf node.
        u64 label = (beta >> path_level) & 1;
        result_cnt += label;

        // Increment the current position.
        pos += length;
        is_not_done = pos != n;

        // Advance the label scanner.
        label_bit_buffer.increment(u32(1) << path_level);
        beta = label_bit_buffer.read();

        // Walk upwards until a left child is found. (might be the current one)
        const auto advance_end = path_level + 1;
        u32 advance_mask_hi = (~u32(0)) >> (32 - advance_end);
        const auto up_steps = dtl::bits::tz_count(~path);
        path >>= up_steps;
        path_level -= up_steps;
        // Go to right sibling.
        path = path | 1;

        // Advance the scanners and update the alpha vector.
        const auto advance_begin = path_level;
        u32 advance_mask_lo = (~u32(0)) << advance_begin;
        u32 advance_mask = advance_mask_lo & advance_mask_hi;

        const auto alpha_next =
            (alpha_la & advance_mask) | (alpha & (~advance_mask));
        tree_bit_buffer.increment(advance_mask);

        // Walk downwards to the left-most leaf in that sub-tree.
        const auto down_steps =
            dtl::bits::tz_count(~(alpha_next >> path_level));
        path_level += down_steps;
        path <<= down_steps;

        // Update the alpha vector.
        alpha = alpha_next;
        alpha_la = tree_bit_buffer.read_ahead();
      }

      // Update the scanner states.
      dtl::r256 tmp_tm {.i = tree_bit_buffer.get_read_mask()};
      dtl::r256 tmp_lm {.i = label_bit_buffer.get_read_mask()};
      for ($u32 level = perfect_levels_ - 1;
           level < teb_.encoded_tree_height_; ++level) {
        scanner_states_[level].node_idx_ += dtl::bits::tz_count(tmp_tm.u8[level]);
        scanner_states_[level].label_idx_ += dtl::bits::tz_count(tmp_lm.u8[level]);
      }
    }

done:
    if (result_cnt < batch_size_
        && pos == n) {
      results_[result_cnt].pos = n;
      results_[result_cnt].length = 0;
      ++result_cnt;
    }
    scan_path_ = path;
    scan_path_level_ = path_level;
    result_cnt_ = result_cnt;
    alpha_ = alpha;
  }
#endif // __AVX2__

#ifdef __AVX512BW__
  void //__teb_inline__
  next_batch_avx512() noexcept __attribute__ ((flatten, hot, noinline)) {
    const auto h = tree_height_;
    const auto eth = teb_.encoded_tree_height_;
    const auto u = perfect_levels_ - 1;
    const auto n = teb_.size();
    // Note: The path variable does NOT contain a sentinel bit. Instead, we
    //       keep track of the paths' level using a separate variable.
    register auto path = scan_path_;
    register auto path_level = scan_path_level_;
    auto pos = path << (h - path_level);
    auto length = n >> path_level;
    auto result_cnt = result_cnt_;
    auto alpha = 0u;
    auto alpha_la = 0u;
    auto beta = 0u;

    dtl::r512 tmp_t { .u16 = {
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0 }};
    dtl::r512 tmp_l { .u16 = {
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0,
        0,0,0,0 }};

    const data_view<const word_type> T {
        teb_.structure_.m_bits.data(),
        teb_.structure_.m_bits.data() + teb_.structure_.m_bits.size(),
    };
    const data_view<const word_type> L {
        teb_.labels_.m_bits.data(),
        teb_.labels_.m_bits.data() + teb_.labels_.m_bits.size(),
    };

    // Bit buffers for the tree structure.
    dtl::bit_buffer_avx512<> tree_bit_buffer;

    // Bit buffers for the labels.
    dtl::bit_buffer_avx512<> label_bit_buffer;

    // Perfect levels. (all bits set to 1)
    for (std::size_t level = 0; level < u; ++level) {
      const auto slot_idx = level;
      tmp_t.u16[slot_idx] = ~u16(0);
      tmp_l.u16[slot_idx] = ~u16(0);
    }
    tree_bit_buffer.set_raw(tmp_t.i);
    label_bit_buffer.set_raw(tmp_l.i);

    const auto batch_size = DEFAULT_BATCH_SIZE;
    while (result_cnt < batch_size - 16 && pos < n) {

      // Initialize the bit buffers.
      tree_bit_buffer.reset_read_mask();
      label_bit_buffer.reset_read_mask();
      for (std::size_t level = u; level < eth; ++level) {
        const auto slot_idx = level;
        i64 delta = - static_cast<i64>(teb_.implicit_inner_node_cnt_);
        tmp_t.u16[slot_idx] = fetch_n_bits(T, static_cast<i64>(scanner_states_[level].node_idx_) + delta, 16);
        tmp_l.u16[slot_idx] = fetch_n_bits(L, static_cast<i64>(scanner_states_[level].label_idx_), 16);
      }
      tree_bit_buffer.set_raw(tmp_t.i);
      label_bit_buffer.set_raw(tmp_l.i);

      // Initialize alpha. (which reads the first bits in the bit buffers)
      alpha = tree_bit_buffer.read();
      alpha_la = tree_bit_buffer.read_ahead();

      // Initialize beta.
      beta = label_bit_buffer.read();

      // Do multiple iteration to consume the bit buffers.
      $u1 is_not_done = true;
      for ($i32 b = 0; b < 15 && is_not_done; ++b) {
        u64 length = n >> path_level;

        assert(pos < n);
        assert(pos + length <= n);
        assert(path_level >= 1);
        assert(length > 0);

        // Produce output (a 1-fill).
        results_[result_cnt].pos = pos;
        results_[result_cnt].length = length;
        // Read the label of the current leaf node.
        u64 label = (beta >> path_level) & 1;
        result_cnt += label;

        // Increment the current position.
        pos += length;
        is_not_done = pos != n;

        // Advance the label scanner.
        label_bit_buffer.increment(__mmask32(1) << path_level);
        beta = label_bit_buffer.read();

        // Walk upwards until a left child is found. (might be the current one)
        const auto advance_end = path_level + 1;
        const __mmask32 advance_mask_hi = (~__mmask32(0)) >> (32 - advance_end);
        const auto up_steps = dtl::bits::tz_count(~path);
        path >>= up_steps;
        path_level -= up_steps;
        // Go to right sibling.
        path = path | 1;

        // Advance the scanners and update the alpha vector.
        const auto advance_begin = path_level;
        const __mmask32 advance_mask_lo = (~__mmask32(0)) << advance_begin;
        const __mmask32 advance_mask = advance_mask_lo & advance_mask_hi;

        const auto alpha_next =
            (alpha_la & advance_mask) | (alpha & (~advance_mask));
        tree_bit_buffer.increment(advance_mask);

        // Walk downwards to the left-most leaf in that sub-tree.
        const auto down_steps =
            dtl::bits::tz_count(~(alpha_next >> path_level));
        path_level += down_steps;
        path <<= down_steps;

        // Update the alpha vector.
        alpha = alpha_next;
        alpha_la = tree_bit_buffer.read_ahead();
      }

      // Update the scanner states.
      dtl::r512 tmp_tm {.i = tree_bit_buffer.get_read_mask()};
      dtl::r512 tmp_lm {.i = label_bit_buffer.get_read_mask()};
      for ($u32 level = u; level < eth; ++level) {
        scanner_states_[level].node_idx_ += dtl::bits::tz_count(tmp_tm.u16[level]);
        scanner_states_[level].label_idx_ += dtl::bits::tz_count(tmp_lm.u16[level]);
      }
    }

done:
    if (result_cnt < batch_size && pos == n) {
      results_[result_cnt].pos = n;
      results_[result_cnt].length = 0;
      ++result_cnt;
    }
    scan_path_ = path;
    scan_path_level_ = path_level;
    result_cnt_ = result_cnt;
    alpha_ = alpha;
  }
#endif // __AVX512BW__

  /// Forwards the iterator to the given position. - Note: This functions
  /// simply calls next() until the desired position has been reached.
  void __teb_inline__
  skip_to(const std::size_t to_pos) noexcept {
    if (to_pos >= teb_.n_) {
      results_[0].pos = teb_.n_;
      results_[0].length = 0;
      result_cnt_ = 1;
      result_read_pos_ = 0;
      return;
    }

    while (!end() && pos() + length() <= to_pos) {
      next();
    }
    // Adjust the current position and fill-length.
    if (!end() && pos() < to_pos) {
      results_[result_read_pos_].length -=
          to_pos - results_[result_read_pos_].pos;
      results_[result_read_pos_].pos = to_pos;
    }
  }

  /// Returns true if the iterator reached the end, false otherwise.
  u1 __forceinline__ //__teb_inline__
  end() const noexcept {
    return results_[result_read_pos_].pos == teb_.n_;
  }

  /// Returns the starting position of the current 1-fill.
  u64 __forceinline__ //__teb_inline__
  pos() const noexcept {
    assert(result_read_pos_ < result_cnt_);
    return results_[result_read_pos_].pos;
  }

  /// Returns the length of the current 1-fill.
  u64 __forceinline__ //__teb_inline__
  length() const noexcept {
    assert(result_read_pos_ < result_cnt_);
    return results_[result_read_pos_].length;
  }

};
//===----------------------------------------------------------------------===//
}; // namespace dtl
