#pragma once

#include <bitset>
#include <list>
#include <queue>
#include <stack>
#include <vector>

#include "boost/dynamic_bitset.hpp"

#include <dtl/bits.hpp>
#include <dtl/bitmap/util/rank1.hpp>
#include <dtl/bitmap/util/bitmap_tree.hpp>
#include <dtl/bitmap/util/binary_tree_structure.hpp>
#include <dtl/static_stack.hpp>

namespace dtl {

//===----------------------------------------------------------------------===//
/// Encodes a bitmap of length N as a binary tree.
///
/// The implementation supports the optimization levels 0, 1, and 2.
///   0 = implements the core idea of TEBs
///   1 = uses the implicit inner node optimization
///   2 = additionally applies the gradual decompression optimization.
/// Please note, that the different optimization levels can be chosen due to
/// academic reasons only. I.e., to analyse and compare the different
/// optimizations levels.
template<i32 optimization_level_ = 2>
class teb {

public:

  // The fundamental storage type. The size of a TEB is a multiple of
  // sizeof(_block_type).
  using _block_type = $u32;

  using position_t = $u32;
  using bitmap_t = boost::dynamic_bitset<$u32>;

  /// The number of bits in the bitmap.
  $u64 n_;

  /// The tree encoded bitmap.
  bitmap_t structure_;
  bitmap_t labels_;

  /// Support data structure for rank1 operations on the tree structure.
  dtl::rank1 rank_;

  /// The number of implicit inner nodes in the tree structure.
  $u32 implicit_inner_node_cnt_;

  /// For testing purposes only.
  $u32 implicit_leaf_node_cnt_;

public:

  /// Tree-encode the given bitmap with an optional false positive rate.
  explicit
  teb(const bitmap_t& bitmap, f64 fpr = 0.0)
      : n_(bitmap.size()) {

    // Construct a binary tree that represents the given bitmap.
    dtl::bitmap_tree bitmap_tree(bitmap, fpr);

    // Encode the tree into level-order.
    for (auto it = bitmap_tree.breadth_first_begin();
         it != bitmap_tree.breadth_first_end();
         ++it) {
      u64 idx = *it;
      // Emit a 1-bit if the current node is an inner node, a 0-bit otherwise.
      u1 is_inner = bitmap_tree.is_inner_node(idx);
      structure_.push_back(is_inner);
      if (!is_inner) {
        // Add the label of the leaf node.
        labels_.push_back(bitmap_tree.label_of_node(idx));
      }
    }

    // Optimization level 1.
    implicit_inner_node_cnt_ = 0;
    implicit_leaf_node_cnt_ = 0;
    if (optimization_level_ > 0) {
      // Count the number of consecutive least significant 1-bits.
      std::size_t trailing_1bits = 0;
      for (std::size_t i = 0; i < structure_.size(); ++i) {
        if (!structure_[i]) {
          break;
        }
        trailing_1bits++;
      }
      implicit_inner_node_cnt_ = static_cast<u32>(trailing_1bits);

      // Remove the trailing 1-bits from the tree structure.
      structure_ >>= implicit_inner_node_cnt_;
      structure_.resize(structure_.size() - implicit_inner_node_cnt_);

      // Remove the most significant 0-bits from the tree structure.
      while (structure_.size() > 0 && !structure_[structure_.size() - 1]) {
        structure_.pop_back();
        implicit_leaf_node_cnt_++;
      }
    }

    // Init rank1 support data structure.
    rank_.init(structure_);

    // Optimization level 2
    if (optimization_level_ > 1) {
      run_optimize();
    }
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
  boost::dynamic_bitset<$u32> __forceinline__
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
  std::size_t __forceinline__
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

  u1 __forceinline__
  operator!=(teb& other) const noexcept {
    return !(*this == other);
  }

  u1 __forceinline__
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

  /// Returns the value of the bit at the position pos.
  u1 __forceinline__
  test(const std::size_t pos) const noexcept {
    auto n_log2 = dtl::log_2(n_);
    $u64 node_idx = 0;
    if (is_leaf_node(node_idx)) { // FIXME: eliminate special case!!!
      return get_label(node_idx);
    }
    for ($u64 i = n_log2 - 1; i < n_log2; i--) {
      u1 bit = dtl::bits::bit_test(pos, i);
      auto r = rank(node_idx + 1);
      node_idx = 2 * r - 1 + bit; // right child if bit is set, left child otherwise
      if (is_leaf_node(node_idx)) {
        u64 label_idx = node_idx - rank(node_idx); // FIXME: do not call rank() twice!!!
        auto label = labels_[label_idx];
        return label;
      }
    }
    std::cout << "BÄM" << std::endl;
    std::exit(42);
  }

  /// Returns true if all bits are set, false otherwise.
  u1 __forceinline__
  all() noexcept {
    // FIXME: this works only if the tree mask is in a compressed state
    return is_leaf_node(0) // root is the only node (a leaf)
        && get_label(0) == true; // and the label is 1
  }

  /// Returns true if all bits are zero, false otherwise.
  u1 __forceinline__
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

  /// Computes the number of perfect levels in the tree structure based on the
  /// number of implicit inner nodes. The number of perfect levels is at least
  /// one, because there is at least one node in the tree structure.
  static inline u64
  determine_perfect_tree_levels(u64 implicit_inner_node_cnt) noexcept {
    return implicit_inner_node_cnt == 0
           ? 1
           : dtl::log_2(implicit_inner_node_cnt + 1) + 1;
  }

  /// Computes the height of the tree based on n. - Note that a tree, consisting
  /// of a single (root) node has a height of 0.
  static inline u64
  determine_tree_height(u64 n) noexcept {
    return dtl::log_2(n);
  }

  //===--------------------------------------------------------------------===//
  /// 1-fill iterator, with skip support.
  class iter {

    /// The fundamental type to encode paths within the tree.
    using path_t = $u64;

    /// Reference to the TEB instance.
    const teb& teb_;
    /// The height of the tree structure.
    u64 tree_height_;
    /// The number of perfect levels in the tree structure.
    u64 perfect_levels_;
    /// First node index in the last perfect level.
    u64 top_node_idx_end_;
    /// Last node index + 1 in the last perfect level.
    u64 top_node_idx_begin_;
    /// The current top node.
    $u64 top_node_idx_current_;
    /// The stack contains the tree nodes that need to be visited.
    static_stack<$u64, 32> stack_;
    /// Points to the beginning of the current 1-fill.
    $u64 pos_;
    /// The length of the current 1-fill.
    $u64 length_;
    //===------------------------------------------------------------------===//

  public:

    /// Constructs an iterator for the given TEB instance. After construction,
    /// the iterator points to the first 1-fill.
    explicit
    iter(const teb& teb) noexcept :
        teb_(teb),
        tree_height_(determine_tree_height(teb.n_)),
        perfect_levels_(
            determine_perfect_tree_levels(teb_.implicit_inner_node_cnt_)),
        top_node_idx_end_((1ull << perfect_levels_) - 1),
        top_node_idx_begin_(top_node_idx_end_ - (1ll << (perfect_levels_ - 1))),
        top_node_idx_current_(top_node_idx_begin_) {
      u64 root_node_idx = 0;
      if (teb.is_leaf_node(root_node_idx)) {
        // Handle the special case, where the root node is a leaf node.
        u1 label = teb.get_label(root_node_idx);
        if (label) {
          pos_ = 0;
          length_ = teb.n_;
        }
        else {
          pos_ = teb.n_;
          length_ = 0;
        }
      }
      else {
        // Initialize the stack.
        path_t path = path_t(1) << (perfect_levels_ - 1);
        stack_.push((top_node_idx_begin_ << 32) | path);
        next();
      }
    }

    iter(iter&&) noexcept = default;

    /// Forwards the iterator to the next 1-fill (if any).
    /// Use the functions pos() and length() to get the 1-fill the iterator
    /// is currently pointing to.
    void __forceinline__
    next() noexcept {
      while (top_node_idx_current_ < top_node_idx_end_) {
        while (!stack_.empty()) {
          u64 pair = stack_.top();
          $u64 node_idx = pair >> 32;
          $u64 path = pair & ((u64(1) << 32) - 1);
          stack_.pop();
          while (teb_.is_inner_node(node_idx)) {
            // Push right child on the stack and go to left child.
            node_idx = teb_.left_child(node_idx);
            u64 right_child_idx = node_idx + 1;
            path = path << 1;
            path_t right_child_path = path | 1;
            stack_.push((right_child_idx << 32) | right_child_path);
          }
          // Reached a leaf node.
          u1 label = teb_.get_label(node_idx);
          if (label) {
            // Produce output (a 1-fill).
            const auto lz_cnt_path = dtl::bits::lz_count(path);
            const auto level = sizeof(path_t) * 8 - 1 - lz_cnt_path;
            // Toggle sentinel bit (= highest bit set) and add offset.
            pos_ = (path ^ (1ull << level)) << (tree_height_ - level);
            // The length of the 1-fill.
            length_ = teb_.n_ >> level;
            return;
          }
        }

        ++top_node_idx_current_;
        auto path = path_t(top_node_idx_current_ - top_node_idx_begin_);
        // Set the sentinel bit.
        path |= path_t(1) << (perfect_levels_ - 1);
        stack_.push((top_node_idx_current_ << 32) | path);
      }
      pos_ = teb_.n_;
      length_ = 0;
    }

    /// Navigate to the desired position, starting from the trees' root node.
    void __forceinline__
    nav_to(const std::size_t to_pos) noexcept {
      //===----------------------------------------------------------------===//
      // (Re-)initialize the iterator state.
      stack_.clear();

      // Determine the top-node path.
      $u64 level = perfect_levels_ - 1;
      const auto foo = to_pos >> (tree_height_ - level);
      path_t path = (path_t(1) << level) | foo;

      // Determine the top-node idx.
      top_node_idx_current_ = top_node_idx_begin_ + foo;
      //===----------------------------------------------------------------===//

      $u64 node_idx = top_node_idx_current_;
      // Walk down the tree to the desired position.
      std::size_t i = tree_height_ - perfect_levels_;
      while (true) {
        // First check, if this is already a leaf node.
        if (teb_.is_leaf_node(node_idx)) {
          // Reached the desired position.
          if (teb_.get_label(node_idx)) {
            // done
            const auto lz_cnt_path = dtl::bits::lz_count(path);
            // Toggle sentinel bit (= highest bit set) and add offset.
            pos_ = (path ^ (1ull << level)) << (tree_height_ - level);
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
        // 0 -> goto left child, 1 -> goto right child
        u1 direction_bit = dtl::bits::bit_test(to_pos, i);
        i--;
        const auto r = teb_.rank(node_idx + 1);
        const auto right_child = 2 * r;
        const auto left_child = right_child - 1;
        level++;
        if (!direction_bit) {
          // Go to left child.
          stack_.push((right_child << 32) | ((path << 1) | 1));
          path <<= 1;
          node_idx = left_child;
        }
        else {
          // Go to right child.
          path = (path << 1) | 1;
          node_idx = right_child;
        }
      }
    }

    /// Fast-forwards the iterator to the given position.
    void __forceinline__
    skip_to(const std::size_t to_pos) noexcept {
      nav_to(to_pos);
    }

    /// Returns true if the iterator reached the end, false otherwise.
    u1 __forceinline__
    end() const noexcept {
      return pos_ == teb_.n_;
    }

    /// Returns the starting position of the current 1-fill.
    u64 __forceinline__
    pos() const noexcept {
      return pos_;
    }

    /// Returns the length of the current 1-fill.
    u64 __forceinline__
    length() const noexcept {
      return length_;
    }

  };

  iter __forceinline__
  it() const noexcept {
    return std::move(iter(*this));
  }

private:

  u1 __forceinline__
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

  u1 __forceinline__
  is_leaf_node(u64 node_idx) const {
    return !is_inner_node(node_idx);
  }

  /// Important: rank() calculates the rank of the prefix -> we need idx + 1
  u64 __forceinline__
  left_child(u64 node_idx) const {
    const std::size_t implicit_1bit_cnt = implicit_inner_node_cnt_;
    if (node_idx < implicit_1bit_cnt) {
      return 2 * node_idx + 1;
    }
    return 2 * rank(node_idx + 1) - 1;
  }

  u64 __forceinline__
  right_child(u64 node_idx) const {
    const std::size_t implicit_1bit_cnt = implicit_inner_node_cnt_;
    if (node_idx < implicit_1bit_cnt) {
      return 2 * node_idx + 2;
    }
    return 2 * rank(node_idx + 1);
  }

  std::size_t __forceinline__
  get_label_idx(u64 node_idx) const {
    return node_idx - rank(node_idx);
  }

  u1 __forceinline__
  get_label(u64 node_idx) const {
    u64 label_idx = get_label_idx(node_idx);
    return labels_[label_idx];
  }

  u64 __forceinline__
  rank_scan(u64 node_idx) const {
    const std::size_t implicit_1bit_cnt = implicit_inner_node_cnt_;
    if (node_idx < implicit_1bit_cnt) {
      return node_idx;
    }
    const auto i = std::min(node_idx - implicit_1bit_cnt, structure_.size());
    auto ret_val = implicit_1bit_cnt;
    for (std::size_t j = 0; j < i; ++j) {
      ret_val += structure_[j];
    }
    return ret_val;
  }

  u64 __forceinline__
  rank(u64 node_idx) const {
    const std::size_t implicit_1bit_cnt = implicit_inner_node_cnt_;
    if (node_idx < implicit_1bit_cnt) {
      return node_idx;
    }
    const auto i = std::min(node_idx - implicit_1bit_cnt, structure_.size());
    const auto ret_val = implicit_1bit_cnt + rank_(i);
    return ret_val;
  }

};
//===----------------------------------------------------------------------===//

}; // namespace dtl
