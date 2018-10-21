#pragma once

#include <bitset>
#include <list>
#include <queue>
#include <vector>

#include <dtl_storage/tree_mask.hpp>
#include <stack>

#include "boost/dynamic_bitset.hpp"

#include "rank.hpp"
#include "dynamic_tree.hpp"
#include "static_stack.hpp"

#define VERBOSE_OUT true

namespace dtl {

/// Encodes a bitmap of length N as a binary tree.
/// The tree structure is encoded in level-order.
class dynamic_tree_mask_lo {
public:

  using bitmap_t = boost::dynamic_bitset<$u32>;


//  //===----------------------------------------------------------------------===//
//  /// Helper structure to navigate within the tree structure.
//  class traversal {
////  private:
//  public:
//
//    using path_t = uint64_t;
//    static constexpr path_t path_msb = path_t(1) << (sizeof(path_t) * 8 - 1);
//
//    const bitmap_t& structure_;
//    const bitmap_t& labels_;
//    const rank1& rank_;
//
//    std::size_t s_pos_;
//    std::size_t l_pos_;
//
//    path_t path_; // encodes the path to the current node (the highest set bit is a sentinel bit)
//    $u32 level_ = 0;
//
//  public:
//
//    /// C'tor
//    __forceinline__
//    traversal(const bitmap_t& structure,
//              const bitmap_t& labels,
//              const rank1& rank)
//        : structure_(structure), labels_(labels), rank_(rank), s_pos_(0), l_pos_(0), path_(1) {}
//
//    __forceinline__ explicit
//    traversal(const dynamic_tree_mask_lo& tm) : traversal(tm.structure_, tm.labels_, tm.rank_) {}
//
//    __forceinline__
//    ~traversal() = default;
//
//    __forceinline__
//    traversal(const traversal& other) = default;
//
//    __forceinline__
//    traversal(traversal&& other) noexcept = delete;
//
//    __forceinline__
//    traversal&
//    operator=(const traversal& other) = delete;
//
//    __forceinline__
//    traversal&
//    operator=(traversal&& other) noexcept = delete;
//
//    /// Return 'true' if the current node is an inner node, 'false' otherwise.
//    __forceinline__ u1
//    is_inner_node() const { return structure_[s_pos_]; }
//
//    /// Return 'true' if the current node is leaf node, 'false' otherwise.
//    __forceinline__ u1
//    is_leaf_node() const { return !is_inner_node(); }
//
//    /// Return 'true' if the current node a left child, 'false' otherwise.
//    __forceinline__ u1
//    is_left_child() const { return (path_ & 1ull) == 0; }
//
//    /// Return 'true' if the current node a right child, 'false' otherwise.
//    __forceinline__ u1
//    is_right_child() const { return !is_left_child(); }
//
//    /// Return the label of the current leaf node.
//    /// The result is undefined, if the current node is an inner node.
//    __forceinline__ u1
//    get_label() const { assert(is_leaf_node()); return labels_[l_pos_]; }
//
//    /// Navigate to the next node (in level-order).
//    __forceinline__ u1
//    next() {
//      if (end()) return false;
//      assert(l_pos_ < labels_.size());
//      if (structure_[s_pos_] /* is inner node */) {
//        // go to left child
//        path_ <<= 1;
//        level_++;
//      }
//      else {
//        // is leaf node
//        path_++;
//        l_pos_++;
//        const auto tzc = dtl::bits::tz_count(path_);
//        path_ >>= tzc;
//        level_ -= tzc;
//      }
//      s_pos_++;
//      return true;
//    }
//
//    __forceinline__ u1
//    end() const {
//      return s_pos_ == structure_.size();
//    }
//
//    /// Return the level of the current node.
//    __forceinline__ auto
//    get_level() const {
//      return level_;
//    }
//
//    /// Navigate to the left child.
//    /// The current node must be an inner node.
//    __forceinline__ void
//    goto_left_child() {
//      assert(is_inner_node());
//      s_pos_++;
//      level_++;
//    }
//
//    /// Navigate to the right child. (Naive implementation)
//    /// The current node must be an inner node.
//    __forceinline__ void
//    goto_right_child() {
//      if (end()) return;
//      assert(is_inner_node());
//      // navigate to the left hand side child
//      next();
//      // determine the level of the left child
//      const auto level = get_level();
//      // traverse left sub tree until the level is reached again
//      next();
//      while (get_level() != level) {
//        next();
//      }
//    }
//
//    /// Compute and return the level-order node index of the current node.
//    __forceinline__ std::size_t
//    get_node_idx() {
//      return tree_mask_po::get_node_idx(path_);
//    }
//  };
//  //===----------------------------------------------------------------------===//


  // the number of bits
  u64 N;

  // level-order encoding
  bitmap_t structure_;
  bitmap_t labels_;

  dtl::rank1 rank_;

  /// C'tor
  explicit
  dynamic_tree_mask_lo(const boost::dynamic_bitset<$u32>& bitmask)
      : N(bitmask.size()) {

    if (!dtl::is_power_of_two(N)) {
      throw std::invalid_argument("The length of the bitmask must be a power of two.");
    }

    using tree_t = dtl::dynamic_binary_tree_structure;

    tree_t tree_structure(N);
    u64 length = tree_structure.max_node_cnt;
    u64 height = tree_structure.height;
    bitmap_t labels(length, false);

    // initialize a complete binary tree
    // ... all the inner nodes have two children
    // ... the leaf nodes are labelled with the given bitmask
    for ($u64 i = length / 2; i < length; i++) {
      labels[i] = bitmask[i - length / 2];
    }
    // propagate the mask bits along the tree (bottom-up)
    for ($u64 i = 0; i < length - 1; i++) {
      u64 node_idx = length - i - 1;
      labels[tree_t::parent_of(node_idx)] = labels[tree_t::parent_of(node_idx)] | labels[node_idx];
    }

    // bottom-up pruning (loss-less)
    for ($u64 i = 0; i < length - 1; i += 2) {
      u64 left_node_idx = length - i - 2;
      u64 right_node_idx = left_node_idx + 1;

      u1 left_bit = labels[left_node_idx];
      u1 right_bit = labels[right_node_idx];

      u64 parent_node_idx = tree_t::parent_of(left_node_idx);

      u1 prune_causes_false_positives = left_bit ^ right_bit;
      u1 both_nodes_are_leaves = !tree_structure.is_inner_node(left_node_idx) & !tree_structure.is_inner_node(right_node_idx);
      u1 prune = both_nodes_are_leaves & !prune_causes_false_positives;
      if (prune) {
        tree_structure.set_leaf(parent_node_idx);
      }
    }

    // Encode the tree into level-order

    std::queue<$u64> fifo; // explored nodes set

    // Allocate enough space for the structure and labels
    structure_.resize(tree_structure.max_node_cnt);
    labels_.resize(labels.size());
    $u64 struct_cnt = 0;
    $u64 label_cnt = 0;

    std::function<void(u64)> add_node = [&](u64 idx){
      u1 is_inner = tree_structure.is_inner_node(idx);
      if(is_inner){
        // add the children of the current node to the structure
        u64 l_child = tree_structure.left_child_of(idx);
        u64 r_child = tree_structure.right_child_of(idx);

        u1 l_child_is_inner = tree_structure.is_inner_node(l_child);
        u1 r_child_is_inner = tree_structure.is_inner_node(r_child);

        structure_[struct_cnt++] = l_child_is_inner;
        structure_[struct_cnt++] = r_child_is_inner;

        // push them to the fifo queue to traverse them later
        fifo.push(l_child);
        fifo.push(r_child);

        // done for this node

      } else { // leaf node
        // add the label of the leaf node
        labels_[label_cnt++] = labels[idx];

        // no children, nothing to add to the fifo queue, done
      }
    };

    // Special case for the root node: if the tree is only the root, structure is 0
    {
      u1 root_is_inner = tree_structure.is_inner_node(0);

      // add the root to the tree structure
      structure_[struct_cnt++] = root_is_inner;

      if(root_is_inner) {
        // add the root to the fifo to add the rest of the tree
        fifo.push(0);
      }
      else {
        labels_[label_cnt++] = labels[0];
        // tree is only the root, we are done
      }
    }

    while(!fifo.empty()){
      u64 next_node = fifo.front();
      add_node(next_node);
      fifo.pop();
    }

    structure_.resize(struct_cnt);
    labels_.resize(label_cnt);

//    rank_support.set_vector(&structure_);
    rank_.init(structure_);
  }

  explicit
  dynamic_tree_mask_lo(u64 N, const bitmap_t& structure, const bitmap_t& labels)
      : N(N) {

    structure_ = structure;
    labels_ = labels;
    rank_.init(structure_);
  }

  __forceinline__ u1
  is_inner_node(u64 node_idx) const {
    return structure_[node_idx];
  }

  __forceinline__ u1
  is_leaf_node(u64 node_idx) const {
    return !structure_[node_idx];
  }

  __forceinline__ u1
  is_root_node(u64 node_idx) const {
    return (node_idx == 0);
  }

  __forceinline__ u1
  is_left_child(u64 node_idx) const {
    return (node_idx != 0) && (node_idx & 1);
  }

  __forceinline__ u1
  is_right_child(u64 node_idx) const {
    return (node_idx != 0) && !(node_idx & 1);
  }

  // naive rank for tests
  __forceinline__ u64
  rank_off(u64 node_idx) const {
    $u64 rank = 0;
    for(auto i = 0; i < node_idx; i++) {
      rank += structure_[i];
    }
    return rank;
  }

  // naive rank for tests
  __forceinline__ u64
  rank(u64 node_idx) const {
    return rank_(node_idx);
  }

  /// Important: rank() calculates the rank of the prefix -> we need idx + 1
  __forceinline__ u64
  left_child(u64 node_idx) const {
    return 2 * rank(node_idx + 1) - 1;
  }

  __forceinline__ u64
  right_child(u64 node_idx) const {
    return 2 * rank(node_idx + 1);
  }

  __forceinline__ u1
  get_label(u64 node_idx) const {
    u64 label_idx = node_idx - rank(node_idx);
    return this->labels_[label_idx];
  }

  /// Decodes the level-order encoding to a bitmap.
  __forceinline__ boost::dynamic_bitset<$u32>
  to_bitset() const {
    boost::dynamic_bitset<$u32> ret(N, false); // the resulting bitmap

    // special case if the tree is only a root node
    if(structure_.size() == 1){
      //std::cout << "to_bitset() only root case" << std::endl;
      if(labels_[0]){
        return boost::dynamic_bitset<$u32>(N, true);
      } else {
        return ret;
      }
    }

    // normal cases: construct the tree using the lo_construction and a level counter
    $u64 level = 0; //current height
    $u64 write_pointer = 0;
    u64 tree_height = dtl::log_2(N);

    // a doubly linked list<node_idx, level> to manage the decoding for all nodes
    std::list<std::pair<$u64, $u64>> nodes;

    std::function<void(std::pair<u64, u64>)> decode_tree = [&](std::pair<u64, u64> n){
      u64 idx = n.first;
      u64 level = n.second;

      if(is_inner_node(idx)){
        //std::cout << "Inner node idx:" << idx << std::endl;

        u64 l_child = this->left_child(idx);
        u64 r_child = this->right_child(idx);
        u64 child_level = level+1;

        //std::cout << "Left child: " << l_child << " Right_child: " << r_child << std::endl;

        nodes.push_front(std::make_pair(r_child, child_level)); // first push the right child to the front
        nodes.push_front(std::make_pair(l_child, child_level)); // then push the left child to the front
        // resulting list: {left child, right child, ...}

      } else {
        //std::cout << "Leaf node" << std::endl;
        // write the label to the bitset
//        u1 label = labels_[idx - rank_support.rank(idx+1)];
        u1 label = labels_[idx - rank(idx + 1)];
        u64 to_write = 1 << (tree_height - level) ;// number of tuples represented by the label

        for(auto i = 0; i < to_write; ++i){
          ret[i + write_pointer] = label;
        }

        //std::cout << "Written: " << ret << " Level: " << level << " Tree_Height: " << tree_height << std::endl;

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
  __forceinline__ std::size_t
  size_in_byte() const {
    u64 lo_struct_size = (structure_.size() + 7) / 8;
    u64 lo_labels_size = (labels_.size() + 7) / 8;
    u64 rank_supp_bytes = rank_.size_in_bytes();
    return lo_struct_size + lo_labels_size + rank_supp_bytes;
  }

  u1 operator!=(dynamic_tree_mask_lo& other) const {
    return (this->structure_ != other.structure_ || this->labels_ != other.labels_);
  }

  u1 operator==(dynamic_tree_mask_lo& other) const {
    return (this->structure_ == other.structure_ && this->labels_ == other.labels_);
  }

  /// Bitwise XOR without compression of the resulting tree
  dynamic_tree_mask_lo
  operator^(const dynamic_tree_mask_lo& other) const{

    struct node {
        $u64 node_idx;
        $u64 node_pos; // position in the more full of both trees
        $u1 xor_bit;

        node() : node_idx(0), node_pos(0), xor_bit(false) {};

        node(u64 idx, u64 pos , u1 bit) : node_idx(idx), node_pos(pos), xor_bit(bit) {};

        void
        operator=(node& other) {
          this->node_idx = other.node_idx;
          this->node_pos = other.node_pos;
          this->xor_bit = other.xor_bit;
        }
    };

    bitmap_t structure;
    bitmap_t labels;
    std::vector<std::pair<$u64, $u64>> level_offset; // to keep track of the different level begins

    // TODO optimize: determine and skip common prefix

    // a XOR b : a = this, b = other
    std::queue<node> fifo_a; // nodes of the next level of a
    std::queue<node> fifo_b; // nodes of the next level of b

    $u64 node_pos = 1; // needed to keep track of the corresponding nodes of tree a & b

    // push the root-nodes into the queues
    fifo_a.push(node(0, node_pos, this->is_inner_node(0) ? false : this->get_label(0)));
    fifo_b.push(node(0, node_pos, other.is_inner_node(0) ? false : other.get_label(0)));
    //frontier_level_nodes = fifo_a.size(), fifo_b.size();
    node_pos++;

    while(fifo_a.size() > 0 || fifo_b.size() > 0){

      //std::cout << std::endl << "Size queue_a: " << fifo_a.size() << " queue_b: " << fifo_b.size() << std::endl;

      node curr_a;
      node curr_b;

      if(fifo_a.size() > 0){
        curr_a = fifo_a.front();
      }

      if(fifo_b.size() > 0){
        curr_b = fifo_b.front();
      }

      if(curr_a.node_pos == curr_b.node_pos && curr_a.node_pos != 0){

        fifo_a.pop();
        fifo_b.pop();

        // both nodes are no fillers
        u8 c = static_cast<u8>(this->is_inner_node(curr_a.node_idx) |
                               (static_cast<u8>(other.is_inner_node(curr_b.node_idx)) << 1));
        $u1 bit;

        switch (c) {
          case 0b00: // a and b: leaf
            //std::cout << "Case: 0b00" << std::endl;

            bit = this->get_label(curr_a.node_idx) ^ other.get_label(curr_b.node_idx);

            structure.push_back(false); // insert the leaf
            labels.push_back(bit);
            break;

          case 0b01: // a: inner, b: leaf
          {
            //std::cout << "Case: 0b01" << std::endl;

            // add the inner node from a
            structure.push_back(true);

            u1 bit_b = other.get_label(curr_b.node_idx);

            // add the children of a to the queue
            node left_child(this->left_child(curr_a.node_idx), node_pos++, bit_b);
            node right_child(this->right_child(curr_a.node_idx), node_pos++, bit_b);

            fifo_a.push(left_child);
            fifo_a.push(right_child);

            break;
          }
          case 0b10: // a: leaf, b: inner
          {
            //std::cout << "Case: 0b10" << std::endl;

            // add the inner node from b
            structure.push_back(true);

            u1 bit_a = this->get_label(curr_a.node_idx);

            // add the children of b to the queue
            node left_child(other.left_child(curr_b.node_idx), node_pos++, bit_a);
            node right_child(other.right_child(curr_b.node_idx), node_pos++, bit_a);

            fifo_b.push(left_child);
            fifo_b.push(right_child);

            break;
          }
          case 0b11: // a and b: inner
          {
            //std::cout << "Case: 0b11" << std::endl;

            // add the inner node
            structure.push_back(true);

            // add the left child of a and b
            node left_child_a(this->left_child(curr_a.node_idx), node_pos, false);
            node left_child_b(other.left_child(curr_b.node_idx), node_pos, false);
            node_pos++;

            node right_child_a(this->right_child(curr_a.node_idx), node_pos, false);
            node right_child_b(other.right_child(curr_b.node_idx), node_pos, false);
            node_pos++;

            // add the children of a
            fifo_a.push(left_child_a);
            fifo_a.push(right_child_a);

            // add the children of b
            fifo_b.push(left_child_b);
            fifo_b.push(right_child_b);

            break;
          }
        }

      } else { // if we are in a subtree of one of the both nodes

        //std::cout << "Node a: " << curr_a.node_idx << " | " << curr_a.node_pos << " | " << curr_a.xor_bit << std::endl;
        //std::cout << "Node b: " << curr_b.node_idx << " | " << curr_b.node_pos << " | " << curr_b.xor_bit << std::endl;

        if((curr_a.node_pos < curr_b.node_pos && curr_a.node_pos != 0) || curr_b.node_pos == 0){ // add part of the subtree of a

          //std::cout << "Case: subtree_a" << std::endl;

          fifo_a.pop();

          if(this->is_inner_node(curr_a.node_idx)){ // current node is an inner node -> add children to the queue

            //std::cout << "Is inner node" << std::endl;
            // add the inner node to the structure
            structure.push_back(true);

            // add the children of a to the queue
            node left_child(this->left_child(curr_a.node_idx), node_pos++, curr_a.xor_bit);
            node right_child(this->right_child(curr_a.node_idx), node_pos++, curr_a.xor_bit);

            fifo_a.push(left_child);
            fifo_a.push(right_child);

          } else { // node is a leaf, add it to the structure

            //std::cout << "Is leaf node" << std::endl;
            // add the leaf to the structure
            structure.push_back(false);
            labels.push_back(this->get_label(curr_a.node_idx) ^ curr_a.xor_bit);
          }

        } else { // add part of the subtree of b

          //std::cout << "Case: subtree_b" << std::endl;

          fifo_b.pop();

          if(other.is_inner_node(curr_b.node_idx)){ // current node is an inner node -> add children to the queue

            //std::cout << "Is inner node" << std::endl;

            // add the inner node to the structure
            structure.push_back(true);

            // add the children of b to the queue
            node left_child(other.left_child(curr_b.node_idx), node_pos++, curr_b.xor_bit);
            node right_child(other.right_child(curr_b.node_idx), node_pos++, curr_b.xor_bit);

            fifo_b.push(left_child);
            fifo_b.push(right_child);

          } else { // node is a leaf, add it to the structure

            //std::cout << "Is leaf node" << std::endl;

            // add the leaf to the structure
            structure.push_back(false);
            labels.push_back(other.get_label(curr_b.node_idx) ^ curr_b.xor_bit);
          }
        }

        if(curr_a.node_pos == 0 && fifo_a.size() == 1){
          fifo_a.pop();
        }

        if(curr_b.node_pos == 0 && fifo_b.size() == 1 ){
          fifo_b.pop();
        }
      }

      //dynamic_tree_mask_lo<N> tmp(structure, labels);
      //std::cout << tmp << std::endl;
    }

    dynamic_tree_mask_lo ret(N, structure, labels);
    return ret;
  }

  /// Bitwise XOR (range encoding)
  dynamic_tree_mask_lo
  xor_re(const dynamic_tree_mask_lo& other) const {
    return *this ^ other; // TODO: find a way to exploit range encoding properties
  }

  void
  print(std::ostream& os) const {

    for ($i64 i = 0; i < structure_.size(); i++) {
      os << (structure_[i] ? "1" : "0");
    }
    os << " | ";
    for ($i64 i = 0; i < labels_.size(); i++) {
      os << (labels_[i] ? "1" : "0");
    }
  }



  /// Bitwise AND without compression of the resulting tree
  dynamic_tree_mask_lo
  operator&(const dynamic_tree_mask_lo& other) const {

    // the output
    bitmap_t out_structure;
    bitmap_t out_labels;

    struct node {
      $u64 node_idx;
      $u64 node_pos; // position in the more full of both trees

      node() : node_idx(0), node_pos(0) {};

      node(u64 idx, u64 pos) : node_idx(idx), node_pos(pos) {};

      void operator=(node& other){
        this->node_idx = other.node_idx;
        this->node_pos = other.node_pos;
      }
    };


    std::queue<node> fifo_a; // nodes of the next level of a
    std::queue<node> fifo_b; // nodes of the next level of b

    $u64 node_pos = 1;

    // push the root-nodes into the queues
    fifo_a.push(node(0, node_pos));
    fifo_b.push(node(0, node_pos));
    node_pos++;

    while (!fifo_a.empty() || !fifo_b.empty()) {

      node curr_a;
      node curr_b;

      if (!fifo_a.empty()) {
        curr_a = fifo_a.front();
      }

      if (!fifo_b.empty()){
        curr_b = fifo_b.front();
      }

      if(curr_a.node_pos == curr_b.node_pos && curr_a.node_pos != 0){

        fifo_a.pop();
        fifo_b.pop();

        // both nodes are no fillers
        u8 c = static_cast<u8>(this->is_inner_node(curr_a.node_idx) |
            (static_cast<u8>(other.is_inner_node(curr_b.node_idx)) << 1));

        $u1 bit;

        switch (c) {
          case 0b00: // a and b: leaf
            bit = this->get_label(curr_a.node_idx) & other.get_label(curr_b.node_idx);

            out_structure.push_back(false); // insert the leaf
            out_labels.push_back(bit);
            break;

          case 0b01: // a: inner, b: leaf
          {
            //std::cout << "Case: 0b01" << std::endl;

            if (!other.get_label(curr_b.node_idx)) { // a & 0 -> leaf 0
              out_structure.push_back(false);
              out_labels.push_back(false);
            } else {
              // add the inner node from a
              out_structure.push_back(true);

              // add the children of a to the queue
              node left_child(this->left_child(curr_a.node_idx), node_pos++);
              node right_child(this->right_child(curr_a.node_idx), node_pos++);

              fifo_a.push(left_child);
              fifo_a.push(right_child);
            }

            break;
          }
          case 0b10: // a: leaf, b: inner
          {
            //std::cout << "Case: 0b10" << std::endl;

            if (!this->get_label(curr_a.node_idx)) { // 0 & b -> leaf 0
              out_structure.push_back(false);
              out_labels.push_back(false);
            } else {
              // add the inner node from a
              out_structure.push_back(true);

              // add the children of a to the queue
              node left_child(other.left_child(curr_b.node_idx), node_pos++);
              node right_child(other.right_child(curr_b.node_idx), node_pos++);

              fifo_b.push(left_child);
              fifo_b.push(right_child);
            }

            break;
          }
          case 0b11: // a and b: inner
          {
            //std::cout << "Case: 0b11" << std::endl;

            // add the inner node
            out_structure.push_back(true);

            // add the left child of a and b
            node left_child_a(this->left_child(curr_a.node_idx), node_pos);
            node left_child_b(other.left_child(curr_b.node_idx), node_pos);
            node_pos++;

            node right_child_a(this->right_child(curr_a.node_idx), node_pos);
            node right_child_b(other.right_child(curr_b.node_idx), node_pos);
            node_pos++;

            // add the children of a
            fifo_a.push(left_child_a);
            fifo_a.push(right_child_a);

            // add the children of b
            fifo_b.push(left_child_b);
            fifo_b.push(right_child_b);

            break;
          }



        }
      }
      else { // if we are in a subtree of one of the both nodes

        //std::cout << "Node a: " << curr_a.node_idx << " | " << curr_a.node_pos << std::endl;
        //std::cout << "Node b: " << curr_b.node_idx << " | " << curr_b.node_pos << std::endl;

        if ((curr_a.node_pos < curr_b.node_pos && curr_a.node_pos != 0) || curr_b.node_pos == 0) { // add part of the subtree of a

          //std::cout << "Case: subtree_a" << std::endl;

          fifo_a.pop();

          if (this->is_inner_node(curr_a.node_idx)) { // current node is an inner node -> add children to the queue

            //std::cout << "Is inner node" << std::endl;
            // add the inner node to the structure
            out_structure.push_back(true);

            // add the children of a to the queue
            node left_child(this->left_child(curr_a.node_idx), node_pos++);
            node right_child(this->right_child(curr_a.node_idx), node_pos++);

            fifo_a.push(left_child);
            fifo_a.push(right_child);

          } else { // node is a leaf, add it to the structure

            //std::cout << "Is leaf node" << std::endl;
            // add the leaf to the structure
            out_structure.push_back(false);
            out_labels.push_back(this->get_label(curr_a.node_idx));
          }

        }
        else { // add part of the subtree of b

          //std::cout << "Case: subtree_b" << std::endl;

          fifo_b.pop();

          if (other.is_inner_node(curr_b.node_idx)) { // current node is an inner node -> add children to the queue

            //std::cout << "Is inner node" << std::endl;

            // add the inner node to the structure
            out_structure.push_back(true);

            // add the children of b to the queue
            node left_child(other.left_child(curr_b.node_idx), node_pos++);
            node right_child(other.right_child(curr_b.node_idx), node_pos++);

            fifo_b.push(left_child);
            fifo_b.push(right_child);

          }
          else { // node is a leaf, add it to the structure

            //std::cout << "Is leaf node" << std::endl;

            // add the leaf to the structure
            out_structure.push_back(false);
            out_labels.push_back(other.get_label(curr_b.node_idx));
          }
        }

        if(curr_a.node_pos == 0 && fifo_a.size() == 1){
          fifo_a.pop();
        }

        if(curr_b.node_pos == 0 && fifo_b.size() == 1 ){
          fifo_b.pop();
        }
      }

      //tree_mask_lo tmp(structure, labels);
      //std::cout << tmp << std::endl;
    }

    dynamic_tree_mask_lo ret(N, out_structure, out_labels);
    return ret;
  }

  /// Bitwise AND (range encoding)
  dynamic_tree_mask_lo
  and_re(const dynamic_tree_mask_lo& other) const {
    return *this & other; // TODO: find a way to exploit range encoding properties
  }

  /// Computes (a XOR b) & this
  /// Note: this, a and b must be different instances. Otherwise, the behavior is undefined.
  dynamic_tree_mask_lo&
  fused_xor_and(const dynamic_tree_mask_lo& a, const dynamic_tree_mask_lo& b) const {

    // get (a XOR b)
    dynamic_tree_mask_lo tree_mask_xor = a ^ b;

    // get this & (a XOR b)
    dynamic_tree_mask_lo tree_mask_and = *this & tree_mask_xor;

    return tree_mask_and; // FIXME: reference to local variable ‘tree_mask_and’ returned
  }

  /// Return the name of the implementation.
  static std::string
  name() {
    return "dynamic_tree_mask_lo";
  }

  /// Returns the value of the bit at the position pos.
  u1
  test(const std::size_t pos) const {
    auto n_log2 = dtl::log_2(N);
    $u64 node_idx = 0;
    if (is_leaf_node(node_idx)) {
      return get_label(node_idx);
    }
    for ($u64 i = n_log2 - 1; i < n_log2; i--) {
      u1 bit = dtl::bits::bit_test(pos, i);
      auto r = rank(node_idx + 1);
      node_idx = 2 * r - 1 + bit; // right child if bit is set, left child otherwise
      if (is_leaf_node(node_idx)) {
        u64 label_idx = node_idx - rank(node_idx);
        auto label = labels_[label_idx];
        return label;
      }
    }
    std::cout << "BÄM" << std::endl;
    std::exit(42);
  }

  u1
  all() {
    return structure_[0] == false // root is the only node (a leaf)
        && labels_[0] == true; // and the label is 1
  }

  u1
  none() {
    return structure_[0] == false // root is the only node (a leaf)
        && labels_[0] == false; // and the label is 0
  }


  //===----------------------------------------------------------------------===//
  /// 1-fill iterator, with skip support.
  class iter {

    using path_t = uint64_t;
    static constexpr path_t path_msb = path_t(1) << (sizeof(path_t) * 8 - 1);

    const dynamic_tree_mask_lo& tm_;
    u64 tree_height = dtl::log_2(tm_.N);


    //===----------------------------------------------------------------------===//
    // Iterator state
    //===----------------------------------------------------------------------===//

//    std::stack<std::pair<$u64, path_t>> stack_;
//    std::stack<$u64> stack_;
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

    void
    next() {
      while (!stack_.empty()) {
        u64 pair = stack_.top();
        u64 node_idx = pair >> 32;
        u64 path = pair & ((u64(1) << 32) - 1);
        stack_.pop();
        if (likely(!tm_.is_leaf_node(node_idx))) {
          // goto left child
          const auto r = tm_.rank(node_idx + 1);
          const auto right_child = 2 * r;
          const auto left_child = right_child - 1;
          const auto left_child_path = path << 1;
          const auto right_child_path = left_child_path | 1;
          stack_.push((right_child << 32) | right_child_path);
          stack_.push((left_child << 32) | left_child_path);
//          stack_.push(std::make_pair(right_child, right_child_path));
//          stack_.push(std::make_pair(left_child, left_child_path));
        }
        else {
          u1 label = tm_.get_label(node_idx);
          if (label) {
            // produce output (a 1-fill)
            const auto lz_cnt_path = dtl::bits::lz_count(path);
            level_ = sizeof(path_t) * 8 - 1 - lz_cnt_path;
            pos_ = (path ^ (path_msb >> lz_cnt_path)) << (tree_height - level_); // toggle sentinel bit (= highest bit set) and add offset
            length_ = tm_.N >> level_; // the length of the 1-fill
            path_ = path;
            return;
          }
        }
      }
      pos_ = tm_.N;
    }

    explicit
    iter(const dynamic_tree_mask_lo& tm) : tm_(tm) {
      const auto n_log2 = dtl::log_2(tm_.N);
      u64 root_node_idx = 0;
      if (tm.is_leaf_node(root_node_idx)) {
        u1 label = tm.get_label(root_node_idx);
        if (label) {
          pos_ = 0;
          length_ = tm.N;
          level_ = 0;
        }
        else {
          pos_ = tm.N;
          length_ = 0;
        }
        return;
      }
//      stack_.push(std::make_pair(root_node_idx, path_t(1)));
      stack_.push((root_node_idx << 32) | 1);

      next();
    }

    iter(iter&&) = default;

    void
    skip_to(const std::size_t to_pos) {
      assert(to_pos >= pos_ + length_);

      // determine the common ancestor
      const auto shift_amount = ((sizeof(path_t) * 8) - level_);
      const auto a = path_ << shift_amount;
      const auto b = to_pos << shift_amount;
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
          pos_ = tm_.N;
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
            length_ = tm_.N >> level_; // the length of the 1-fill
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

    u1
    end() const noexcept {
      return pos_ == tm_.N;
    }

    u64
    pos() const noexcept {
      return pos_;
    }

    u64
    length() const noexcept {
      return length_;
    }

  };

  iter
  it() const {
    return iter(*this);
  }

  struct range {
            // [first;last)
            $u64 first = 0;
            $u64 last = 0;

            range() = default;

            range(u64 start, u64 length) : first(start), last(start + length) {}

            range(u64 start_1, u64 length_1, u64 start_2, u64 length_2){
              u64 tmp_first = std::max(start_1, start_2);
              u64 tmp_last = std::min(start_1 + length_1, start_2 + length_2);
              if (tmp_first >= tmp_last) {
                first = 0;
                last = 0;
              } else {
                first = tmp_first;
                last = tmp_last;
              }
            }

            __forceinline__ void
            set(u64 start, u64 length) noexcept {
              first = start;
              last = first + length;
            }

            __forceinline__ void
            set_first_last(u64 start_, u64 last_) noexcept {
              first = start_;
              last = last_;
            }

            void
            print() const noexcept {
              std::cout << "[" << first << ";" << last << ")" << std::endl;
            }

            __forceinline__ u1
            is_empty() const noexcept {
              return first == last;
            }

            __forceinline__ range
            range_and(u64 start_1, u64 length_1, u64 start_2, u64 length_2) const noexcept {
              u64 tmp_first = std::max(start_1, start_2);
              u64 tmp_last = std::min(start_1 + length_1, start_2 + length_2);
              range r;

              if (tmp_first >= tmp_last)
                return r;

              r.set_first_last(tmp_first, tmp_last);
              return r;
            }

            __forceinline__ range
            operator&(const range& other) const noexcept {
              range r;

              if(is_empty())
                return r;

              if(other.is_empty())
                return r;

              r.set_first_last(std::max(first, other.first),std::min(last, other.last));

              if(r.first >= r.last) // empty range
                r.set(0,0);

              return r;
            }
        };

  // simple iterator without skip_to, only using next
  // always returns the first next possible range -> i.e.: [4;9), [9;13)...
  // (is skip_to sufficient? test it...)
  class iter_and_simple {

      using tm_lo = dynamic_tree_mask_lo;

      iter it_1;
      iter it_2;
      bool end_it_1 = false;
      bool end_it_2 = false;

      $u1 side = false; // false -> consume from it_1, true -> consume from it_2

      range r;
      range r_1;
      range r_2;

  public:

      explicit
      iter_and_simple(const tm_lo& tm_1, const tm_lo& tm_2) : it_1(tm_1), it_2(tm_2) {
        next();
      };

      void next() {

        while(!end()){

          if(!side | r_1.is_empty()) {
            end_it_1 = it_1.end();

            if(!end_it_1) {
              r_1.set(it_1.pos(), it_1.length());
              it_1.next();
            } else
              r_1.set(it_1.pos(),0);
          }

          if(side | r_2.is_empty()) {
            end_it_2 = it_2.end();
            if(!end_it_2) {
              r_2.set(it_2.pos(), it_2.length());
              it_2.next();
            } else
              r_2.set(it_2.pos(), 0);
          }

          r = r_1 & r_2;

          if(r_1.last < r_2.last)
            side = false;
          else
            side = true;


          // we have an overlapping interval
          if(!r.is_empty())
            return; // end search

        }

        r.set(it_1.pos(),0);
      }

      u1
      end() const noexcept {
        return end_it_1 | end_it_2;
      }

      const range
      matches() const noexcept {
        return r;
      }

  };

  // simple iterator with skip_to
  // always returns the first next possible range -> i.e.: [4;9), [9;13)...
  // TODO does not work, as skip_to is not working properly
  class iter_and_with_skip_to {

      using tm_lo = dynamic_tree_mask_lo;

      iter it_1;
      iter it_2;
      bool end_it_1 = false;
      bool end_it_2 = false;

      $u1 side = false; // false -> consume from it_1, true -> consume from it_2

      range r;
      range r_1;
      range r_2;

  public:

      explicit
      iter_and_with_skip_to(const tm_lo& tm_1, const tm_lo& tm_2) : it_1(tm_1), it_2(tm_2) {
        next();
      };

      void next() {

        while(!end()){

          if(!side | r_1.is_empty()) {
            end_it_1 = it_1.end();

            if(!end_it_1) {
              r_1.set(it_1.pos(), it_1.length());
            }else
              r_1.set(it_1.pos(),0);
          }

          if(side | r_2.is_empty()) {
            end_it_2 = it_2.end();

            if(!end_it_2) {
              r_2.set(it_2.pos(), it_2.length());
            } else
              r_2.set(it_2.pos(), 0);
          }

          r = r_1 & r_2;

          #if VERBOSE_OUT
            std::cout << "     r_1: "; r_1.print();
            std::cout << "     r_2: "; r_2.print();
            std::cout << "     r  : "; r.print();
            std::cout << std::endl;
          #endif

          if(!r.is_empty()) {
            std::cout << "Normal case" << std::endl;
            if (r_1.last < r_2.last) {
              it_1.next();
              side = false;
            } else {
              it_2.next();
              side = true;
            }
            std::cout << "Current state: it_1.end: " << (it_1.end() ? "True" : "False") <<
            " it_2.end: " << (it_2.end() ? "True" : "False") << std::endl;

            return;
          } else {

            if(it_1.end() | it_2.end()){
              return;
            }

            std::cout << "Check skipping" << std::endl;


            if(r_1.last <= r_2.first){
              std::cout << "it_1 is smaller" << std::endl;
              if(!it_2.end()){
                std::cout << "it_1 skip to: " << r_2.first << std::endl;
                it_1.skip_to(r_2.first);
                std::cout << "New pos: " << it_1.pos() << std::endl;
              }
              side = false;
            } else {

              if(r_2.last <= r_1.first){
                std::cout << "it_2 is smaller" << std::endl;

                if(!it_1.end()){
                  std::cout << "it_2 skip to: " << r_1.first << std::endl;
                  it_2.skip_to(r_1.first);
                  std::cout << "New pos: " << it_2.pos() << std::endl;
                }
                side = true;
              }
            }
          }
        }

        r.set(it_1.pos(),0);
      }

      u1
      end() const noexcept {
        return end_it_1 | end_it_2;
      }

      const range
      matches() const noexcept {
        return r;
      }

  };


  // iterator without skip_to, only using next
  // returns the largest possible range -> i.e.: [4;9), [9;13) = [4;13)
  class iter_and {

            using tm_lo = dynamic_tree_mask_lo;

            iter it_1;
            iter it_2;
            bool end_it_1 = false;
            bool end_it_2 = false;

            $u1 side = false; // false -> consume from it_1, true -> consume from it_2

            range r;
            range r_1;
            range r_2;

        public:

            explicit
            iter_and(const tm_lo& tm_1, const tm_lo& tm_2) : it_1(tm_1), it_2(tm_2) {
              next();
            };

            void next() {

              while(!end()){

                if(!side | r_1.is_empty()) {
                  end_it_1 = it_1.end();

                  if(!end_it_1) {
                      r_1.set(it_1.pos(), it_1.length());
                      it_1.next();

                      while(!it_1.end() && it_1.pos() == r_1.last){
                        r_1.set_first_last(r_1.first, it_1.pos() + it_1.length());
                        it_1.next();
                      }
                  } else
                    r_1.set(it_1.pos(),0);
                }

                if(side | r_2.is_empty()) {
                  end_it_2 = it_2.end();

                  if(!end_it_2) {
                    r_2.set(it_2.pos(), it_2.length());
                    it_2.next();

                    while(!it_2.end() && it_2.pos() == r_2.last){
                      r_2.set_first_last(r_2.first, it_2.pos() + it_2.length());
                      it_2.next();
                    }
                  } else
                    r_2.set(it_2.pos(), 0);
                }

                r = r_1 & r_2;

                if(r_1.last < r_2.last)
                  side = false;
                else
                  side = true;


                // we have an overlapping interval
                if(!r.is_empty())
                  return; // end search

              }

              r.set(it_1.pos(),0);
            }

            u1
            end() const noexcept {
              return end_it_1 | end_it_2;
            }

            const range
            matches() const noexcept {
              return r;
            }

        };

//  struct iter_alex {
//
//      const dynamic_tree_mask_lo& tm;
//      u64 tree_height = dtl::log_2(tm.N);
//
//      //===----------------------------------------------------------------------===//
//      // Iterator state
//      //===----------------------------------------------------------------------===//
//
//      /// <level, node position in struct>
//      using node = std::pair<$u64, u64>;
//
//      std::list<node> state;
//
//      /// points to the beginning of a 1-fill
//      $u64 pos = 0; //
//
//      /// the length of the current 1-fill
//      $u64 length = tm.N;
//
//
//      void next() {
//
//        while (!state.empty()) {
//          u64 level = state.front().first;
//          u64 node_idx = state.front().second;
//          state.pop_front();
//
//          if (!tm.is_leaf_node(node_idx)) {
//
//            const auto rank = tm.rank(node_idx);
//            const auto l_child = 2 * rank - 1;
//            const auto r_child = 2 * rank;
//
//            state.push_front(std::make_pair(level + 1, r_child));
//            state.push_front(std::make_pair(level + 1, l_child));
//          } else {
//
//            u1 label = tm.get_label(node_idx);
//
//            if (label) { // true
//
//              //TODO
//              length = 1 << (tree_height - level); // // todo nachrechnen ob das stimmt
//
//              return;
//            } else {
//              // increment the pos and proceed
//              pos +=  1 << (tree_height - level); // todo nachrechnen ob das stimmt
//            }
//          }
//
//        }
//
//        pos = tm.N;
//      }
//
//      explicit
//      iter(const dynamic_tree_mask_lo& tm){
//      }
//
//
//      void skip_to() {
//
//      }
//
//      u1 end() {
//        return pos == tm.N;
//      }
//
//  };

};

//TODO:
// 1. XOR
// 2. Fused XOR
// 3. BP

}; // namespace dtl