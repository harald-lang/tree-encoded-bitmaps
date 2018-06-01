#pragma once

#include <queue>

#include <boost/dynamic_bitset.hpp>

#include <dtl/tree_mask.hpp>

#include <sdsl/int_vector.hpp>
// Boost
#include "boost/dynamic_bitset.hpp"
#include <queue>
#include <list>

namespace dtl{
/// Encodes a bitmap of length N as a binary tree.
/// The tree structure is encoded in level-order.

template<std::size_t _N>
class tree_mask_lo{
public:
  static constexpr auto N = _N;

  // level-order encoding
  sdsl::bit_vector lo_struc;
  sdsl::bit_vector lo_label;

  sdsl::rank_support_v5<> rank_support;

  /// C'tor
  explicit
  tree_mask_lo(const std::bitset<N>& bitmask){

    using tree_t = dtl::binary_tree_structure<N>;

    static constexpr u64 length = tree_t::max_node_cnt;
    static constexpr u64 height = tree_t::height;

    tree_t tree_structure;
    std::bitset<length> labels;

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
    lo_struc.resize(tree_structure.max_node_cnt);
    lo_label.resize(labels.size());
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

        lo_struc[struct_cnt++] = l_child_is_inner;
        lo_struc[struct_cnt++] = r_child_is_inner;

        // push them to the fifo queue to traverse them later
        fifo.push(l_child);
        fifo.push(r_child);

        // done for this node

      } else { // leaf node
        // add the label of the leaf node
        lo_label[label_cnt++] = labels[idx];

        // no children, nothing to add to the fifo queue, done
      }
    };

      // Special case for the root node: if the tree is only the root, structure is 0
      {
        bool root_is_inner = tree_structure.is_inner_node(0);

      // add the root to the tree structure
      lo_struc[struct_cnt++] = root_is_inner;

      if(root_is_inner){
        // add the root to the fifo to add the rest of the tree
        fifo.push(0);
      } else {
        lo_label[label_cnt++] = labels[0];
        // tree is only the root, we are done
      }
    }

    while(!fifo.empty()){
      u64 next_node = fifo.front();
      add_node(next_node);
      fifo.pop();
    }

    lo_struc.resize(struct_cnt);
    lo_label.resize(label_cnt);

    rank_support.set_vector(&lo_struc);
  }

  __forceinline__ u1
  is_inner_node(u64 node_idx){
    return lo_struc[node_idx] == 1 ? true : false;
  }

  __forceinline__ u1
  is_leaf_node(u64 node_idx){
    return lo_struc[node_idx] == 0 ? true : false;
  }

  // naive rank for tests
  u64 rank(u64 node_idx){
    $u64 rank = 0;

    for(auto i = 0; i <= node_idx; i++)
      if(lo_struc[i])
        rank++;

    return rank;
  }

  /// Important: rank_supportv5.rank() calculates the rank of the prefix -> we need idx + 1
  __forceinline__ u64
  left_child(u64 node_idx){
    return 2* rank_support.rank(node_idx+1) -1;
  }

  __forceinline__ u64
  right_child(u64 node_idx){
    return 2* rank_support.rank(node_idx+1);
  }


  /// decodes the level-order encoding to a bitmap
  __forceinline__ std::bitset<N>
  to_bitset(){

    std::bitset<N> ret; // the resulting bitset

    // special case if the tree is only a root node
    if(lo_struc.size() == 1){
      if(lo_label[0]){
        ret.set();
      } else {
        ret.reset();
      }

      return ret;
    }

    // normal cases: construct the tree using the lo_construction and a level counter
    $u64 level = 0; //current height
    $u64 write_pointer = 0;
    constexpr u64 tree_height = dtl::ct::log_2<N>::value;

    // a doubly linked list<node_idx, level> to manage the decoding for all nodes
    std::list<std::pair<$u64, $u64>> nodes;

    std::function<void(std::pair<u64, u64>)> decode_tree = [&](std::pair<u64, u64> n){
      u64 idx = n.first;
      u64 level = n.second;

      if(is_inner_node(idx)){

        u64 l_child = left_child(idx);
        u64 r_child = right_child(idx);

        nodes.push_front(std::make_pair(r_child, level+1)); // first push the right child to the front
        nodes.push_front(std::make_pair(l_child, level+1)); // then push the left child to the front
        // resulting list: {left child, right child, ...}

      } else {
        // write the label to the bitset
        u1 label = lo_label[idx - rank_support.rank(idx+1)];
        u64 to_write = 1 << (tree_height - level) ;// number of tuples represented by the label

        for(auto i = 0; i < to_write; ++i){
          ret[i + write_pointer] = label;
        }

        write_pointer += to_write;
      }
    };

    nodes.push_front(std::make_pair(0,0)); // push the root node to the list

    while(!nodes.empty()){
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
  size_in_byte() {
    u64 lo_struct_size = lo_struc.bit_size();
    u64 lo_labels_size = lo_label.bit_size();

    // the required space of an int_vector with n bits: 64*ceil(n/64+1) bit = 8*ceil(n/64+1) byte
    u64 lo_struct_bytes = 8 * std::ceil((lo_struct_size*1.0 / 64) +1);
    u64 lo_labels_bytes = 8 * std::ceil((lo_labels_size*1.0 / 64) +1);

    // the additional required space for the rank_support_v5 is: 0.0625*n bit
    u64 rank_supp_bytes = ((rank_support.size() * 0.0625) + 7) / 8;

    // std::cout << "Struct Bits: " << lo_struct_size << " Size: " << lo_struct_bytes << std::endl;
    // std::cout << "Labels Bits: " << lo_labels_size << " Size: " << lo_labels_bytes << std::endl;
    // std::cout << "R_Supp Size: " << rank_supp_bytes << std::endl;

    return lo_struct_bytes + lo_labels_bytes + rank_supp_bytes;
  }

  /// Bitwise XOR
  tree_mask_lo
  operator^(const tree_mask_lo& other) const{
    // TODO
  }

  /// Bitwise XOR (Works only in combination with range encoding (RE),
  /// i.e., the following must hold: this[i] == true => other[i] == true)
  tree_mask_lo
  xor_re(const tree_mask_lo& other) const {
    // TODO
  }

  /// Computes (a XOR b) & this
  /// Note: this, a and b must be different instances. Otherwise, the behavior is undefined.
  tree_mask_po&
  fused_xor_and(const tree_mask_po& a, const tree_mask_po& b) {
    // TODO
  }

};

//TODO:
// 1. XOR
// 2. XOR RE
// 3. Fused XOR
// 4. BP

}; // namespace dtl