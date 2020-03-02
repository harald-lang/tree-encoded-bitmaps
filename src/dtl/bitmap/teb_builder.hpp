#pragma once
//===----------------------------------------------------------------------===//
#include "teb_flat.hpp"
#include "teb_types.hpp"
#include "util/bitmap_tree.hpp"
#include "util/bitmap_view.hpp"

#include <dtl/dtl.hpp>

#include <boost/dynamic_bitset.hpp>
#include <boost/concept_check.hpp>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// Converts a plain bitmap into a serialized TEB.
class teb_builder {
  using size_type = teb_size_type;
  using word_type = teb_word_type;
  static constexpr auto word_size = sizeof(word_type);
  static constexpr auto word_bitlength = word_size * 8;

  /// The intermediate representation of the bitmap.
  dtl::bitmap_tree<> bitmap_tree_;

public:
  /// C'tor
  explicit teb_builder(const boost::dynamic_bitset<$u32>& bitmap)
      : bitmap_tree_(bitmap) {
    bitmap_tree_.ensure_counters_are_valid();
  }

  explicit teb_builder(const bitmap_tree<>&& bitmap_tree)
      : bitmap_tree_(std::move(bitmap_tree)) {
    bitmap_tree_.ensure_counters_are_valid();
  }

  /// Returns the serialized size in number of words.
  inline std::size_t
  serialized_size_in_words() {
    std::size_t word_cnt = 0;
    word_cnt += sizeof(teb_header) / word_size;
    word_cnt += tree_word_cnt();
    word_cnt += rank_word_cnt();
    word_cnt += label_word_cnt();
    word_cnt += metadata_word_cnt();
    return word_cnt;
  }

  /// Serializes the TEB to the given destination address.
  void
  serialize(word_type* dst);

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
inline void
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

  auto* label_ptr = teb_flat::get_label_ptr(ptr);
  auto label_word_cnt = teb_flat::get_label_word_cnt(ptr);
  dtl::bitmap_view<word_type> label_data(label_ptr, label_ptr + label_word_cnt);

  std::array<std::size_t, 32> level_offsets_tree; // NOLINT
  std::array<std::size_t, 32> level_offsets_labels; // NOLINT

  const std::size_t tree_bit_cnt = tree_word_cnt * word_bitlength;
  const std::size_t label_bit_cnt = label_word_cnt * word_bitlength;

#ifndef NDEBUG
  // Slow-path. Left there for validation. // TODO remove
  // Encode the tree into level-order.
  std::size_t tree_data_write_pos = 0;
  std::size_t label_data_write_pos = 0;

  std::size_t node_cntr = 0;
  std::size_t leaf_node_cntr = 0;
  std::size_t current_level = ~0ull;

  // Skip over the perfect levels.
  node_cntr = bitmap_tree_.get_leading_inner_node_cnt();

  const auto it_end = bitmap_tree_.breadth_first_end();
  for (auto it = bitmap_tree_.breadth_first_begin(node_cntr);
      it != it_end; ++it) {
    u64 idx = (*it);
    u64 level = bitmap_tree_.level_of(idx);

    ++node_cntr;

    // Memorize the level offsets.
    if (current_level != level) {
      level_offsets_tree[level] = node_cntr - 1;
      level_offsets_labels[level] = leaf_node_cntr;
      current_level = level;
    }

    assert(node_cntr > hdr.implicit_inner_node_cnt);

    u1 is_inner = bitmap_tree_.is_inner_node(idx);
    u1 is_leaf = !is_inner;
    leaf_node_cntr += is_leaf;

    // Emit a 1-bit if the current node is an inner node, a 0-bit otherwise.
    if (idx < bitmap_tree_.get_last_explicit_node_idx()) {
      assert(tree_data_write_pos < hdr.tree_bit_cnt);
      tree_data.set(tree_data_write_pos, is_inner);
      ++tree_data_write_pos;
    }
    // Add the label of the leaf node (if necessary).
    if (is_leaf) {
      if (idx >= bitmap_tree_.get_first_node_idx_with_1label()
          && idx <= bitmap_tree_.get_last_node_idx_with_1label()) {
        assert(label_data_write_pos < hdr.label_bit_cnt);
        label_data.set(label_data_write_pos, bitmap_tree_.label_of_node(idx));
        ++label_data_write_pos;
      }
    }
  }
  // Fill the remaining bits with 0s.
  tree_data.clear(tree_data_write_pos, hdr.tree_bit_cnt);
  label_data.clear(label_data_write_pos, hdr.label_bit_cnt);
#endif

  {
    auto* is_active_ptr = bitmap_tree_.is_active_node_.data();
    auto* is_inner_ptr = bitmap_tree_.is_inner_node_.data();
    auto* label_ptr = bitmap_tree_.labels_.data();
    const auto word_cnt = bitmap_tree_.is_active_node_.data_end()
        - bitmap_tree_.is_active_node_.data_begin();

    auto succinct_tree_writer = tree_data.writer(0, hdr.tree_bit_cnt);
    auto succinct_labels_writer = label_data.writer(0, hdr.label_bit_cnt);

    const std::size_t c = bitmap_tree_.get_leading_inner_node_cnt() + bitmap_tree_.offset;
    const std::size_t z = bitmap_tree_.get_leading_0label_cnt();

    auto append_to_tree = [&](word_type bits_to_append, std::size_t bit_cnt) {
      succinct_tree_writer.write(bits_to_append, bit_cnt);
    };

    // Keep track of the total number of observed labels.
    std::size_t label_cntr = 0;
    auto append_to_labels = [&](word_type bits_to_append, const std::size_t bit_cnt) {
      if (likely(label_cntr >= z)) {
        succinct_labels_writer.write(bits_to_append, bit_cnt);
      }
      else if (label_cntr + bit_cnt > z) {
        auto o = (label_cntr + bit_cnt) - z;
        assert(o < word_bitlength);
        auto a = bits_to_append >> (bit_cnt - o);
        auto c = bit_cnt - (bit_cnt - o);
        succinct_labels_writer.write(a, c);
      }
      // Update the counter.
      label_cntr += bit_cnt;
    };

    // Iterate over the (implicit) tree word-wise in level order. As we skip
    // over the perfect levels, the first relevant node may not be word aligned,
    // thus we handle the first read differently.
    {
      const auto w = c / word_bitlength;
      const auto o = c % word_bitlength;
      auto is_active_word = is_active_ptr[w] >> o;
      auto is_inner_word = is_inner_ptr[w] >> o;
      auto label_word = label_ptr[w] >> o;

      // Extract the relevant tree bits from the current word.
      u64 tree_bits_extract_mask = is_active_word;
      u64 tree_bits_to_append = _pext_u64(is_inner_word, tree_bits_extract_mask);

      // Extract the relevant label bits.
      u64 label_bits_extract_mask = ~is_inner_word & is_active_word;
      u64 label_bits_to_append = _pext_u64(label_word, label_bits_extract_mask);

      // Determine the number of extracted bits.
      u64 extracted_tree_bit_cnt = dtl::bits::pop_count(tree_bits_extract_mask);
      u64 extracted_label_bit_cnt = dtl::bits::pop_count(label_bits_extract_mask);

      append_to_tree(tree_bits_to_append, extracted_tree_bit_cnt);
      append_to_labels(label_bits_to_append, extracted_label_bit_cnt);
    }
    // Read the rest of the tree (word aligned).
    for (std::size_t w = (c / word_bitlength) + 1; w < word_cnt; ++w) {
      // The current part of the tree.
      auto is_active_word = is_active_ptr[w];
      auto is_inner_word = is_inner_ptr[w];
      auto label_word = label_ptr[w];

      // Extract the relevant tree bits from the current word.
      u64 tree_bits_extract_mask = is_active_word;
      u64 tree_bits_to_append = _pext_u64(is_inner_word, tree_bits_extract_mask);

      // Extract the relevant label bits.
      u64 label_bits_extract_mask = ~is_inner_word & is_active_word;
      u64 label_bits_to_append = _pext_u64(label_word, label_bits_extract_mask);

      // Determine the number of extracted bits.
      u64 extracted_tree_bit_cnt = dtl::bits::pop_count(tree_bits_extract_mask);
      u64 extracted_label_bit_cnt = dtl::bits::pop_count(label_bits_extract_mask);

      append_to_tree(tree_bits_to_append, extracted_tree_bit_cnt);
      append_to_labels(label_bits_to_append, extracted_label_bit_cnt);
    }
    succinct_tree_writer.flush();
    succinct_labels_writer.flush();

    // Compute the level offsets.
    if (hdr.has_level_offsets) {
      // Used to compute the offsets for the labels.
      const auto active_leaf_nodes =
          bitmap_tree_.is_active_node_.and_not(bitmap_tree_.is_inner_node_);
      // The number of entries in the table.
      const auto entry_cnt = hdr.encoded_tree_height - hdr.perfect_level_cnt;

      std::size_t offset_tree = 0;
      std::size_t offset_labels = 0;
      {
        // Init offsets.
        assert(hdr.perfect_level_cnt > 0);
        const std::size_t level = hdr.perfect_level_cnt - 1;
        offset_tree = bitmap_tree_.is_active_node_.count(
            bitmap_tree_.offset,
            bitmap_tree_.first_node_idx_at_level(level) + bitmap_tree_.offset);
        offset_labels = active_leaf_nodes.count(
            bitmap_tree_.offset,
            bitmap_tree_.first_node_idx_at_level(level) + bitmap_tree_.offset);
      }

      for (std::size_t i = 0; i < entry_cnt; ++i) {
        const std::size_t level = (i + hdr.perfect_level_cnt);
        offset_tree += bitmap_tree_.is_active_node_.count(
          bitmap_tree_.first_node_idx_at_level(level - 1) + bitmap_tree_.offset,
          bitmap_tree_.first_node_idx_at_level(level) + bitmap_tree_.offset);
        offset_labels += active_leaf_nodes.count(
            bitmap_tree_.first_node_idx_at_level(level - 1) + bitmap_tree_.offset,
            bitmap_tree_.first_node_idx_at_level(level) + bitmap_tree_.offset);
        level_offsets_tree[level] = offset_tree;
        level_offsets_labels[level] = offset_labels;
      }
    }
  }

  // Write the rank LuT.
  auto* rank_ptr = teb_flat::get_rank_ptr(ptr);
  if (rank_ptr != nullptr) {
    teb_rank_logic_type::init_inplace(
        tree_ptr, tree_ptr + tree_word_cnt, rank_ptr);
  }

  // Write the additional meta data.
  if (hdr.has_level_offsets) {
    // Level offset LuT.
    const auto entry_cnt = hdr.encoded_tree_height - hdr.perfect_level_cnt;
    auto* metadata_ptr = teb_flat::get_metadata_ptr(ptr);
    auto* ofs_tree = reinterpret_cast<size_type*>(metadata_ptr);
    auto* ofs_labels = ofs_tree + entry_cnt;
    for (std::size_t i = 0; i < entry_cnt; ++i) {
      ofs_tree[i] = static_cast<size_type>(level_offsets_tree[i + hdr.perfect_level_cnt]);
      ofs_labels[i] = static_cast<size_type>(level_offsets_labels[i + hdr.perfect_level_cnt]);
    }
  }
}
//===----------------------------------------------------------------------===//
}; // namespace dtl
