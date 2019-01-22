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

#define VERBOSE_OUT true

namespace dtl {

//===----------------------------------------------------------------------===//
/// Encodes a bitmap of length N as a binary tree.
/// The tree structure is encoded in level-order.
class teb {
//public:

  using bitmap_t = boost::dynamic_bitset<$u32>;

  /// The number of bits in the bitmap.
  u64 n_;

  /// The tree encoded bitmap.
  bitmap_t structure_;
  bitmap_t labels_;

  /// Support data structure for rank1 operations on the tree structure.
  dtl::rank1 rank_;

  /// The number of perfect levels in the tree structure.
  uint8_t perfect_levels_;

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

    // Count the number of consecutive least significant 1-bits.
    std::size_t trailing_1bits = 0;
    for (std::size_t i = 0; i < structure_.size(); ++i) {
      if (!structure_[i]) {
        break;
      }
      trailing_1bits++;
    }

    // Determine the number of "perfect levels" in the binary tree.
    perfect_levels_ = 0;
    if (dtl::is_power_of_two(trailing_1bits + 1)) {
      perfect_levels_ = static_cast<uint8_t>(dtl::log_2(trailing_1bits + 1));
    }
    else {
      perfect_levels_ = static_cast<uint8_t>(dtl::log_2(
          dtl::prev_power_of_two(trailing_1bits + 1)));
    }

    // Remove the 1-bits of the perfect levels from the tree structure.
    const std::size_t shift_amount = (1ull << perfect_levels_) - 1;
    structure_ >>= shift_amount;
    structure_.resize(structure_.size() - shift_amount);

    // Remove the most significant 0-bits from the tree structure.
    while (structure_.size() > 0 && !structure_[structure_.size() - 1]) {
      structure_.pop_back();
    }

    // Init rank1 support data structure.
    rank_.init(structure_);
  }

//  explicit
//  teb(u64 n, const bitmap_t& structure, const bitmap_t& labels)
//      : n_(n), structure_(structure), labels_(labels) {
//    rank_.init(structure_);
//  }

  /// Decodes the level-order encoding to a bitmap.
  boost::dynamic_bitset<$u32> __forceinline__
  to_bitset() const {
    boost::dynamic_bitset<$u32> ret(n_); // the resulting bitmap

    // special case if the tree is only a root node
    if(perfect_levels_ == 0 && structure_.size() <= 1){
      //std::cout << "to_bitset() only root case" << std::endl;
      if(labels_[0]){
        ret.flip();
      }
      return ret;
    }

    // normal cases: construct the tree using the lo_construction and a level counter
    $u64 level = 0; //current height
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
  size_in_byte() const {
    u64 lo_struct_size = (structure_.size() + 7) / 8;
    u64 lo_labels_size = (labels_.size() + 7) / 8;
    u64 rank_supp_bytes = rank_.size_in_bytes();
    return lo_struct_size + lo_labels_size + rank_supp_bytes;
  }

  /// Return the size in bytes.
  std::size_t __forceinline__
  serialized_size_in_byte() const {
    u64 lo_struct_size = (structure_.size() + 7) / 8;
    u64 lo_labels_size = (labels_.size() + 7) / 8;
    return 4 + lo_struct_size + lo_labels_size;
  }

  u1 __forceinline__
  operator!=(teb& other) const {
    return !(*this == other);
  }

  u1 __forceinline__
  operator==(teb& other) const {
    return perfect_levels_ == other.perfect_levels_
        && structure_ == other.structure_
        && labels_ == other.labels_;
  }

  void
  print(std::ostream& os) const {
    os << "perfect levels=" << static_cast<int32_t>(perfect_levels_) << " | ";
    const std::size_t implicit_1bit_cnt = (1ull << perfect_levels_) - 1;

    if (perfect_levels_ > 0) {
      for ($i64 i = 0; i < implicit_1bit_cnt; i++) {
        os << "1";
      }
      os << "'";
    }
    for ($i64 i = 0; i < structure_.size(); i++) {
      os << (structure_[i] ? "1" : "0");
    }
    os << " | ";
    for ($i64 i = 0; i < labels_.size(); i++) {
      os << (labels_[i] ? "1" : "0");
    }
  }

  /// Return the name of the implementation.
  static std::string
  name() {
    return "teb";
  }

  /// Returns the value of the bit at the position pos.
  u1 __forceinline__
  test(const std::size_t pos) const {
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
  all() {
    // FIXME: this works only if the tree mask is in a compressed state
    return is_leaf_node(0) // root is the only node (a leaf)
        && get_label(0) == true; // and the label is 1
  }

  /// Returns true if all bits are zero, false otherwise.
  u1 __forceinline__
  none() {
    // FIXME: this works only if the tree mask is in a compressed state
    return is_leaf_node(0) // root is the only node (a leaf)
        && get_label(0) == false; // and the label is 0
  }

  std::size_t
  size() const {
    return n_;
  }

  //===----------------------------------------------------------------------===//
  /// 1-fill iterator, with skip support.
  class iter {

    using path_t = uint64_t;
    static constexpr path_t path_msb = path_t(1) << (sizeof(path_t) * 8 - 1);

    const teb& tm_;
    u64 tree_height = dtl::log_2(tm_.n_);


    //===----------------------------------------------------------------------===//
    // Iterator state
    //===----------------------------------------------------------------------===//

    static_stack<$u64, 32> stack_;

    /// encodes the path to the current node (the highest set bit is a sentinel bit)
    path_t path_ = 1;
    /// the level of the current tree node
    $u64 level_ = 0; // FIXME somewhat redundant with path and length
    /// points to the beginning of a 1-fill
    $u64 pos_;
    /// the length of the current 1-fill
    $u64 length_;
    //===----------------------------------------------------------------------===//

  public:

    void __forceinline__
    next() {
      while (!stack_.empty()) {
        u64 pair = stack_.top();
        u64 node_idx = pair >> 32;
        u64 path = pair & ((u64(1) << 32) - 1);
        stack_.pop();

        //std::cout << "structure size: " << tm_.structure_.size() << " label_size: " << tm_.labels_.size() << std::endl;
        //std::cout << "node index: " << node_idx << std::endl;
        if(node_idx >= tm_.structure_.size()){
          break;
        }

        if (!tm_.is_leaf_node(node_idx)) {
          // goto left child
          const auto r = tm_.rank(node_idx + 1);
          const auto right_child = 2 * r;
          const auto left_child = right_child - 1;
          const auto left_child_path = path << 1;
          const auto right_child_path = left_child_path | 1;
          stack_.push((right_child << 32) | right_child_path);
          stack_.push((left_child << 32) | left_child_path);
        }
        else {
          u1 label = tm_.get_label(node_idx);
          if (label) {
            // produce output (a 1-fill)
            const auto lz_cnt_path = dtl::bits::lz_count(path);
            level_ = sizeof(path_t) * 8 - 1 - lz_cnt_path;
            pos_ = (path ^ (path_msb >> lz_cnt_path)) << (tree_height - level_); // toggle sentinel bit (= highest bit set) and add offset
            length_ = tm_.n_ >> level_; // the length of the 1-fill
            path_ = path;
            return;
          }
        }
      }
      pos_ = tm_.n_;
      length_ = 0;
    }

    explicit
    iter(const teb& tm) : tm_(tm) {
      const auto n_log2 = dtl::log_2(tm_.n_);
      u64 root_node_idx = 0;
      if (tm.is_leaf_node(root_node_idx)) {
        u1 label = tm.get_label(root_node_idx);
        if (label) {
          pos_ = 0;
          length_ = tm.n_;
          level_ = 0;
        }
        else {
          pos_ = tm.n_;
          length_ = 0;
        }
        return;
      }
//      stack_.push(std::make_pair(root_node_idx, path_t(1)));
      stack_.push((root_node_idx << 32) | 1);

      next();
    }

    iter(iter&&) = default;

    static path_t __forceinline__
    toggle_msb(path_t i) {
      return i ^ (path_t(1) << (sizeof i * 8 - dtl::bits::lz_count(i) - 1));
    }

    void __forceinline__
    nav_to(const std::size_t to_pos) {
      if (to_pos >= tm_.n_) {
        pos_ = tm_.n_;
        length_ = 0;
        return;
      }
      level_ = 0;
      stack_.clear();
      $u64 node_idx = 0;
      path_ = 1;
//      std::cout << "to_pos=" << std::bitset<64>(to_pos) << std::endl;
      // walk down the tree to the desired position
      std::size_t i = tree_height - 1;
      while (true) {
//        std::cout << "path=" << std::bitset<64>(path_) << std::endl;
        // first check, if this is already a leaf node
        if (tm_.is_leaf_node(node_idx)) {
          // reached the desired position
          if (tm_.get_label(node_idx)) {
            // done
            const auto lz_cnt_path = dtl::bits::lz_count(path_);
            pos_ = (path_ ^ (path_msb >> lz_cnt_path)) << (tree_height - level_); // toggle sentinel bit (= highest bit set) and add offset
            length_ = tm_.n_ >> level_; // the length of the 1-fill
            // adjust the current position and fill-length
            length_ -= to_pos - pos_;
            pos_ = to_pos;
            return;
          }
          else {
            // search forward to the next 1-fill
            next();
            return;
          }
        }

        // navigate downwards the tree
        u1 bit = dtl::bits::bit_test(to_pos, i); // 0 -> goto left child, 1 -> goto right child
        i--;
        const auto r = tm_.rank(node_idx + 1);
        const auto right_child = 2 * r;
        const auto left_child = right_child - 1;
        level_++;
        if (!bit) {
          // goto left child
          stack_.push((right_child << 32) | ((path_ << 1) | 1));
          path_ <<= 1;
          node_idx = left_child;
        }
        else {
          // goto right child
          path_ = (path_ << 1) | 1;
          node_idx = right_child;
        }
      }
    }

    void __forceinline__
    skip_to(const std::size_t to_pos) {
      nav_to(to_pos);
    }

    // FIXME: buggy
    void __forceinline__
    skip_to_OFF(const std::size_t to_pos) {
      assert(to_pos >= pos_ + length_);
      if (to_pos >= tm_.n_) {
        pos_ = tm_.n_;
        return;
      }
      if (to_pos == pos_ + length_) {
        next();
        return;
      }

      // determine the common ancestor
      const auto shift_amount = ((sizeof(path_t) * 8) - tree_height);
//      std::cout << "path=" << std::bitset<64>(path_) << std::endl;
      const auto a = toggle_msb(path_) << shift_amount;
//      std::cout << "   a=" << std::bitset<64>(a) << std::endl;
      const auto b = to_pos << shift_amount;
//      std::cout << "   b=" << std::bitset<64>(b) << std::endl;
      const auto a_xor_b = a ^ b;
      const auto common_prefix_len = a_xor_b == 0 ? 0 : dtl::bits::lz_count(a_xor_b);

      // walk up the tree to the common ancestor
//      stack_.pop(); // requires the TM to be compressed?
      const auto level_of_common_ancestor = common_prefix_len;
      while (true) {
        u64 pair = stack_.top();
        u64 node_idx = pair >> 32;
        u64 path = pair & ((u64(1) << 32) - 1);
//        $u64 node_idx = stack_.top().first;
//        $u64 path = stack_.top().second;
        const auto lz_cnt_path = dtl::bits::lz_count(path);
        const auto level = sizeof(path_t) * 8 - 1 - lz_cnt_path;
        if (level_of_common_ancestor + 1 == level) {
          level_ = level;
          break;
        }
        stack_.pop();
        if (stack_.empty()) {
          // end
          pos_ = tm_.n_;
          break;
        }
      }

      // common ancestor
      u64 pair = stack_.top();
      $u64 node_idx = pair >> 32;
      $u64 path = pair & ((u64(1) << 32) - 1);
//      $u64 node_idx = stack_.top().first;
//      $u64 path = stack_.top().second;
      stack_.pop();

      // walk down the tree to the desired position
      std::size_t i = tree_height - level_ - 1;
      while (true) {

        // first check, if this is already a leaf node
        if (tm_.is_leaf_node(node_idx)) {
          // reached the desired position
          if (tm_.get_label(node_idx)) {
            // done
            const auto lz_cnt_path = dtl::bits::lz_count(path);
            pos_ = (path ^ (path_msb >> lz_cnt_path)) << (tree_height - level_); // toggle sentinel bit (= highest bit set) and add offset
            length_ = tm_.n_ >> level_; // the length of the 1-fill
            // adjust the current position and fill-length
            length_ -= to_pos - pos_;
            pos_ = to_pos;
            return;
          }
          else {
            // search forward to the next 1-fill
            next();
            return;
          }
        }

        // navigate downwards the tree
        u1 bit = dtl::bits::bit_test(to_pos, i--); // 0 -> goto left child, 1 -> goto right child
        const auto r = tm_.rank(node_idx + 1);
        const auto right_child = 2 * r;
        const auto left_child = right_child - 1;
        level_++;
        if (!bit) {
          // goto left child
//          stack_.push(std::make_pair(right_child, (path << 1) | 1));
          stack_.push((right_child << 32) | ((path << 1) | 1));
          path <<= 1;
          node_idx = left_child;
        }
        else {
          // goto right child
          path = (path << 1) | 1;
          node_idx = right_child;
        }

      }
    }

    u1 __forceinline__
    end() const noexcept {
      return pos_ == tm_.n_;
    }

    u64 __forceinline__
    pos() const noexcept {
      return pos_;
    }

    u64 __forceinline__
    length() const noexcept {
      return (length_ > tm_.n_ ? 0 : length_);
    }

  };

  iter __forceinline__
  it() const {
    return std::move(iter(*this));
  }

private:

  u1 __forceinline__
  is_inner_node(u64 node_idx) const {
    const std::size_t implicit_1bit_cnt = (1ull << perfect_levels_) - 1;
    if (node_idx < implicit_1bit_cnt) return true;
    if ((node_idx - implicit_1bit_cnt) >= structure_.size()) {
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
    const std::size_t implicit_1bit_cnt = (1ull << perfect_levels_) - 1;
    if (node_idx < implicit_1bit_cnt) {
      return 2 * node_idx + 1;
    }
    return 2 * rank(node_idx + 1) - 1;
  }

  u64 __forceinline__
  right_child(u64 node_idx) const {
    const std::size_t implicit_1bit_cnt = (1ull << perfect_levels_) - 1;
    if (node_idx < implicit_1bit_cnt) {
      return 2 * node_idx + 2;
    }
    return 2 * rank(node_idx + 1);
  }

  u1 __forceinline__
  get_label(u64 node_idx) const {
    u64 label_idx = node_idx - rank(node_idx);
    return labels_[label_idx];
  }

  u64 __forceinline__
  rank(u64 node_idx) const {
    const std::size_t implicit_1bit_cnt = (1ull << perfect_levels_) - 1;
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

