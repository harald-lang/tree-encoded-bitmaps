#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/dtl.hpp>
#include <dtl/bitmap/util/bitmap_tree.hpp>
#include "boost/dynamic_bitset.hpp"
#include "util/bitmap_view.hpp"
#include "teb_flat.hpp"
#include "teb_types.hpp"
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
class teb_builder {

  using size_type = teb_size_type;
  using word_type = teb_word_type;
  static constexpr auto word_size = sizeof(word_type);
  static constexpr auto word_bitlength = word_size * 8;

  /// The intermediate representation of the bitmap.
  dtl::bitmap_tree<> bitmap_tree_;

public:

  /// C'tor
  explicit
  teb_builder(const boost::dynamic_bitset<$u32>& bitmap, f64 fpr = 0.0);

  /// Returns the serialized size in number of words.
  std::size_t
  serialized_size_in_words();

  /// Serializes the TEB to the given destination address.
  void
  serialize(word_type* ptr);

private:

  /// Determine the length of the encoded tree structure.
  inline std::size_t
  explicit_node_cnt() {
    return bitmap_tree_.get_node_cnt()
        - bitmap_tree_.get_leading_inner_node_cnt()
        - bitmap_tree_.get_trailing_leaf_node_cnt();
  }

  /// Determine the number of labels.
  inline std::size_t
  explicit_label_cnt() {
    return bitmap_tree_.get_leaf_node_cnt()
        - bitmap_tree_.get_leading_0label_cnt()
        - bitmap_tree_.get_trailing_0label_cnt();
  }

  /// Returns the length of the tree in number of words.
  inline std::size_t
  tree_word_cnt() {
    return (explicit_node_cnt() + word_bitlength - 1) / word_bitlength;
  }

  /// Returns the length of the rank lookup table in number of words.
  inline std::size_t
  rank_word_cnt() {
    const auto tree_bits = explicit_node_cnt();
    return (tree_bits > 0)
        ? (teb_rank_type::estimate_size_in_bytes(tree_bits) + word_size - 1) / word_size
        : 0;
  }

  /// Returns the length of the label bitmap in number of words.
  inline std::size_t
  label_word_cnt() {
    return (explicit_label_cnt() + word_bitlength - 1) / word_bitlength;
  }

  /// Returns the length of the additional meta data in number of words.
  inline std::size_t
  metadata_word_cnt() {
    const auto entry_cnt = bitmap_tree_.get_encoded_tree_height()
        - bitmap_tree_.get_perfect_level_cnt();
    const auto size_in_bytes = 2 * sizeof(size_type) * entry_cnt;
    const auto word_cnt = (size_in_bytes + word_size - 1) / word_size;
    return (explicit_node_cnt() > 1024) ? word_cnt : 0; // TODO remove magic number
  }

};
//===----------------------------------------------------------------------===//
teb_builder::teb_builder(const boost::dynamic_bitset<$u32>& bitmap, f64 fpr)
    : bitmap_tree_(bitmap, fpr) {
}
//===----------------------------------------------------------------------===//
std::size_t
teb_builder::serialized_size_in_words() {
  std::size_t word_cnt = 0;
  word_cnt += sizeof(teb_header) / word_size;
  word_cnt += tree_word_cnt();
  word_cnt += rank_word_cnt();
  word_cnt += label_word_cnt();
  word_cnt += metadata_word_cnt();
  return word_cnt;
}
//===----------------------------------------------------------------------===//
void
teb_builder::serialize(word_type* dst) {
  // Prepare the header.
  teb_header hdr;
  hdr.n = static_cast<u32>(bitmap_tree_.n_);
  hdr.tree_bit_cnt = static_cast<u32>(explicit_node_cnt());
  hdr.implicit_inner_node_cnt = static_cast<u32>(bitmap_tree_.get_leading_inner_node_cnt());
  hdr.label_bit_cnt = static_cast<u32>(explicit_label_cnt());
  hdr.implicit_leading_label_cnt = static_cast<u32>(bitmap_tree_.get_leading_0label_cnt());
  hdr.perfect_level_cnt = static_cast<u8>(bitmap_tree_.get_perfect_level_cnt());
  hdr.encoded_tree_height = static_cast<u8>(bitmap_tree_.get_encoded_tree_height());
  hdr.has_level_offsets = hdr.tree_bit_cnt > 1024 ? u8(1) : u8(0); // TODO remove magic number

  word_type* ptr = dst;

  // Write the header.
  auto* p = reinterpret_cast<teb_header*>(ptr);
  *p = hdr;

  // Write the tree and the labels.
  auto* tree_ptr = teb_flat::get_tree_ptr(ptr);
  auto tree_word_cnt = teb_flat::get_tree_word_cnt(ptr);
  dtl::bitmap_view<word_type> tree_data(tree_ptr, tree_ptr + tree_word_cnt);
  std::size_t tree_data_write_pos = 0;

  auto* label_ptr = teb_flat::get_label_ptr(ptr);
  auto label_word_cnt = teb_flat::get_label_word_cnt(ptr);
  dtl::bitmap_view<word_type> label_data(label_ptr, label_ptr + label_word_cnt);
  std::size_t label_data_write_pos = 0;

  std::array<std::size_t, 32> level_offsets_tree;
  std::array<std::size_t, 32> level_offsets_labels;

  // Encode the tree into level-order.
  std::size_t node_cntr = 0;
  std::size_t leaf_node_cntr = 0;
  std::size_t current_level = ~0ull;
  for (auto it = bitmap_tree_.breadth_first_begin();
       it != bitmap_tree_.breadth_first_end(); ++it) {

    u64 idx = (*it).idx;
    u64 level = (*it).level;
    u1 is_inner = (*it).is_inner;

    ++node_cntr;
    leaf_node_cntr += !is_inner;

    // Memorize the level offsets.
    if (current_level != level) {
      level_offsets_tree[level] = node_cntr - 1;
      level_offsets_labels[level] = leaf_node_cntr;
      current_level = level;
    }

    // Skip the implicit leading inner nodes.
    if (node_cntr <= hdr.implicit_inner_node_cnt) {
      continue;
    }

    // Emit a 1-bit if the current node is an inner node, a 0-bit otherwise.
    if (idx < bitmap_tree_.get_last_explicit_node_idx()) {
      assert(tree_data_write_pos < hdr.tree_bit_cnt);
      tree_data.set(tree_data_write_pos, is_inner);
      ++tree_data_write_pos;
    }
    // Add the label of the leaf node (if necessary).
    if (!is_inner) {
      if (idx >= bitmap_tree_.get_first_node_idx_with_1label()
          && idx <= bitmap_tree_.get_last_node_idx_with_1label()) {
        assert(label_data_write_pos < hdr.label_bit_cnt);
        label_data.set(label_data_write_pos, bitmap_tree_.label_of_node(idx));
        ++label_data_write_pos;
      }
    }
  }
  // Fill the remaining bits with 0s.
  for (std::size_t i = tree_data_write_pos; i < hdr.tree_bit_cnt; ++i) {
    tree_data.set(i, false);
  }
  for (std::size_t i = label_data_write_pos; i < hdr.label_bit_cnt; ++i) {
    label_data.set(i, false);
  }

  // Write the rank LuT.
  auto* rank_ptr = teb_flat::get_rank_ptr(ptr);
  if (rank_ptr != nullptr) {
    teb_rank_type::init_inplace(tree_ptr, tree_ptr + tree_word_cnt, rank_ptr);
  }

  // Write the additional meta data.
  if (hdr.has_level_offsets) {
    // Level offset LuT.
    const auto entry_cnt = hdr.encoded_tree_height - hdr.perfect_level_cnt;
    auto* metadata_ptr = teb_flat::get_metadata_ptr(ptr);
    auto* ofs_tree = reinterpret_cast<size_type*>(metadata_ptr);
    auto* ofs_labels = ofs_tree + entry_cnt;
    for (std::size_t i = 0; i < entry_cnt; ++i) {
      ofs_tree[i] =
          static_cast<size_type>(level_offsets_tree[i + hdr.perfect_level_cnt]);
      ofs_labels[i] =
          static_cast<size_type>(level_offsets_labels[i + hdr.perfect_level_cnt]);
    }
  }

}
//===----------------------------------------------------------------------===//
}; // namespace dtl
