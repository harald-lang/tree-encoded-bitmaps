#pragma once

#include <bitset>

#include <dtl/dtl.hpp>
#include <dtl/math.hpp>

namespace dtl {

  class dynamic_binary_tree_structure {

  public:

    /// The number of leaf nodes.
    u64 N;

    /// The max number of nodes.
    u64 max_node_cnt;

    /// The tree height.
    u64 height;

    /// Indicates whether a node is an internal or a leaf node.
    std::vector<$u1> _is_inner_node;

  public:

    explicit
    dynamic_binary_tree_structure(u64 N)
        :N(N), max_node_cnt(2 * N - 1), height(dtl::log_2(N)) {

      _is_inner_node.resize(max_node_cnt, false);

      // initialize a complete binary tree
      // ... all the inner nodes have two children
      for ($u64 i = 0; i < max_node_cnt / 2; i++) {
        _is_inner_node[i] = true;
      }
    }

    static inline u64
    root() {
      return 0;
    }

    static inline u64
    parent_of(u64 node_idx) {
      return (node_idx - 1) / 2;
    }

    static inline u64
    left_child_of(u64 node_idx) {
      return 2 * node_idx + 1;
    }

    static inline u64
    right_child_of(u64 node_idx) {
      return 2 * node_idx + 2;
    }

    static inline u64
    level_of(u64 node_idx) {
      return log_2(node_idx + 1);
    }

    inline u1
    is_inner_node(u64 node_idx) const {
      return _is_inner_node[node_idx];
    }

    inline u1
    is_leaf_node(u64 node_idx) const {
      return ! _is_inner_node[node_idx];
    }

    inline void
    set_leaf(u64 node_idx) {
      _is_inner_node[node_idx] = false;
    }

    inline void
    set_inner(u64 node_idx) {
      _is_inner_node[node_idx] = true;
    }

  };


  using dynamic_full_binary_tree = dynamic_binary_tree_structure;


} // namespace dtl
