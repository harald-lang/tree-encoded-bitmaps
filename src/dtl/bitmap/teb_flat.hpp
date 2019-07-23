#pragma once

#include "teb.hpp"

namespace dtl {
//===----------------------------------------------------------------------===//
class teb_flat {

  using word_type = dtl::teb<3>::word_type;
  static constexpr auto word_size = sizeof(word_type);
  static constexpr auto word_bitlength = word_size * 8;
  
  using rank_type = dtl::teb<3>::rank_support;
  using size_type = rank_type::size_type;
  
  using bitmap_fn = dtl::bitmap_fun<word_type>;

  const word_type* const ptr_;

  size_type n_;
  size_type implicit_inner_node_cnt_;
  size_type implicit_leading_label_cnt_;
  size_type perfect_levels_;
  size_type encoded_tree_height_;

  const word_type* structure_ptr_;
  size_type structure_bit_cnt_;

  const word_type* label_ptr_;
  size_type label_bit_cnt_;

  const size_type* rank_lut_ptr_;
  rank_type rank_;

public:

  explicit
  teb_flat(const word_type* const ptr) : ptr_(ptr) {
    auto hdr = reinterpret_cast<u32* const>(ptr);
    // Read the header.
    n_ = hdr[0];
    implicit_inner_node_cnt_ = hdr[1];
    implicit_leading_label_cnt_ = hdr[2];
    structure_bit_cnt_ = hdr[3];
    label_bit_cnt_ = hdr[4];
    perfect_levels_ = hdr[5] & 0b11111;
    encoded_tree_height_ = (hdr[5] >> 5) & 0b11111;

    // Initialize pointers to the tree structure, the labels, and the rank
    // helper structure.
    const auto hdr_bytes = 6 * 4;
    const auto hdr_word_cnt = (hdr_bytes
        + (sizeof(word_type) - 1)) / word_size;
    const auto t_begin = &ptr[hdr_word_cnt];
    const auto l_begin = t_begin + (structure_bit_cnt_ + word_bitlength - 1)
        / word_bitlength;
    const auto r_begin = l_begin + (label_bit_cnt_ + word_bitlength - 1)
        / word_bitlength;

    structure_ptr_ = structure_bit_cnt_ > 0 ? t_begin : nullptr;
    label_ptr_ = label_bit_cnt_ > 0 ? l_begin : nullptr;

    if (structure_bit_cnt_ > 1024) {
      rank_lut_ptr_ = reinterpret_cast<const size_type* const>(r_begin);
    }
    else {
      // Initialize rank helper structure.
      rank_.init(structure_ptr_, label_ptr_);
      rank_lut_ptr_ = rank_.lut.data();
    }
  }

  /// Returns the value of the bit at the given position.
  u1 __teb_inline__
  test(const std::size_t pos) const noexcept {
    const auto tree_height = dtl::log_2(n_);

    size_type level = perfect_levels_ - 1;
    const auto foo = pos >> (tree_height - level);

    // Determine the top-node idx.
    const auto top_node_idx_begin = (1ull << (perfect_levels_ - 1)) - 1;
    size_type node_idx = top_node_idx_begin + foo;

    auto i = tree_height - 1 - level;
    while (!is_leaf_node(node_idx)) {
      u1 direction_bit = dtl::bits::bit_test(pos, i);
      node_idx = 2 * rank_inclusive(node_idx) - 1; // left child
      node_idx += direction_bit; // right child if bit is set, left child otherwise
      --i;
    }
    return get_label(node_idx);
  }

private:

  size_type __teb_inline__
  rank_inclusive(size_type node_idx) const noexcept {
    const auto implicit_1bit_cnt = implicit_inner_node_cnt_;
    if (node_idx < implicit_1bit_cnt) {
      return node_idx + 1;
    }
    assert(structure_bit_cnt_ > 0);
    const auto i = std::min(node_idx - implicit_1bit_cnt, structure_bit_cnt_ - 1);
    const auto r = rank_type::get(rank_lut_ptr_, i, structure_ptr_);
    const auto ret_val = implicit_1bit_cnt + r;
    return ret_val;
  }

  u1 __teb_inline__
  is_inner_node(size_type node_idx) const noexcept {
    const auto implicit_1bit_cnt = implicit_inner_node_cnt_;
    const auto implicit_leaf_begin =
        structure_bit_cnt_ + implicit_inner_node_cnt_;
    if (node_idx < implicit_1bit_cnt) {
      // Implicit inner node.
      return true;
    }
    if (node_idx >= implicit_leaf_begin) {
      // Implicit leaf node.
      return false;
    }
    return bitmap_fn::test(structure_ptr_, node_idx - implicit_1bit_cnt);
  }

  u1 __teb_inline__
  is_leaf_node(size_type node_idx) const noexcept {
    return !is_inner_node(node_idx);
  }

  size_type __teb_inline__
  get_label_idx(size_type node_idx) const noexcept {
    return node_idx - rank_inclusive(node_idx);
  }

  u1 __teb_inline__
  get_label_by_idx(size_type label_idx) const noexcept {
    const auto implicit_leading_label_cnt = implicit_leading_label_cnt_;
    const auto implicit_trailing_labels_begin =
        label_bit_cnt_ + implicit_leading_label_cnt;
    if (label_idx < implicit_leading_label_cnt) {
      // An implicit leading 0-label.
      return false;
    }
    if (label_idx >= implicit_trailing_labels_begin) {
      // An implicit trailing 0-label.
      return false;
    }
    return bitmap_fn::test(label_ptr_, label_idx - implicit_leading_label_cnt);
  }

  u1 __teb_inline__
  get_label(size_type node_idx) const noexcept {
    const auto label_idx = get_label_idx(node_idx);
    return get_label_by_idx(label_idx);
  }

};
//===----------------------------------------------------------------------===//
} // namespace dtl