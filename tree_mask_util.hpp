#pragma once

#include <bitset>

#include <dtl/dtl.hpp>
#include <dtl/tree.hpp>

namespace dtl {

//===----------------------------------------------------------------------===//
/// Check whether the given tree mask is compressed.
template<typename tree_mask_t>
u1
is_compressed(const tree_mask_t& tm) {
  constexpr auto N = tree_mask_t::N;
  using traversal_t = typename tree_mask_t::traversal;

  dtl::binary_tree_structure<N> bt;
  std::bitset<N> bt_labels;

  // decode the tree mask into an implicit representation
  traversal_t tmt(tm);
  do {
    const auto node_idx = tmt.get_node_idx();
    if (tmt.is_inner_node()) {
      bt.set_inner(node_idx);
    }
    else {
      bt.set_leaf(node_idx);
      bt_labels[node_idx] = tmt.get_label();
    }
    tmt.next();
  } while (!tmt.end());

  std::vector<$u64> stack;
  stack.push_back(0);
  while (!stack.empty()) {
    const auto current_node = stack.back();
    stack.pop_back();
    if (bt.is_inner_node(current_node)) {
      const auto left_child = bt.left_child_of(current_node);
      const auto right_child = bt.right_child_of(current_node);
      if (bt.is_leaf_node(left_child) && bt.is_leaf_node(right_child)) {
        // both child must have different labels
        if (bt_labels[left_child] == bt_labels[right_child]) {
          return false; // the given tree mask is not valid.
        }
      }
      if (bt.is_inner_node(right_child)) {
        stack.push_back(right_child);
      }
      if (bt.is_inner_node(left_child)) {
        stack.push_back(left_child);
      }
      continue;
    }
  }
  return true;
}
//===----------------------------------------------------------------------===//


template<typename tree_mask_t>
u1
is_valid(const tree_mask_t& tm) {

}


} // namespace dtl