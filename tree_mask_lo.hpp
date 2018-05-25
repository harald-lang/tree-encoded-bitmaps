#pragma once

#include <queue>

#include <boost/dynamic_bitset>
#include <boost/dynamic_bitset.hpp>

#include <dtl/tree_mask.hpp>

#include <sdsl/int_vector.hpp>

namespace dtl{
  /// Encodes a bitmap of length N as a binary tree.
  /// The tree structure is encoded in level-order.

  template<std::size_t N>
  class tree_mask_lo{

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

      // TODO combine pruned structure with encoding -> Analyse the labels after pruning and add them
      // level-order encoding

      boost::dynamic_bitset<> structure;
      boost::dynamic_bitset<> labels;
      std::queue fifo; // explored nodes set (fifo)

      std::function<void(u64)> add_node = [&](u64 idx){
        bool is_inner = tree.is_inner_node(idx);

        if(is_inner){
          // add the children of the current node to the structure
          u64 left_child = tree.left_child_of(idx);
          u64 right_child = tree.right_child_of(idx);

          bool left_child_is_inner = tree.is_inner_node(left_child);
          bool right_child_is_inner = tree.is_inner_node(right_child);

          structure.push_back(left_child_is_inner);
          structure.push_back(right_child_is_inner);

          // push them to the fifo queue to traverse them later
          fifo.push(left_child);
          fifo.push(right_child);

          // done for this node

        } else { // leaf node
          // add the label of the leaf node
          // TODO add labels: labels.push_back(labels[idx]);

          // add the two external nodes for the leaf
          structure.push_back(0);
          structure.push_back(0);

          // no children, nothing to add to the fifo queue, done
        }
      };

      // Special case for the root node: if the tree is only the root, structure is 0
      {
        bool root_is_inner = tree.is_inner_node(0);

        // add the root to the tree structure
        structure.push_back(root_is_inner);

        if(root_is_inner){
          // add the root to the fifo to add the rest of the tree
          fifo.push(0);
        } else {
          // tree is only the root, we are done
        }
      }

      while(!fifo.empty()){
        u64 next_node = fifo.front();
        add_node(next_node);
        fifo.pop();
      }
    }

  };

}; // namespace dtl