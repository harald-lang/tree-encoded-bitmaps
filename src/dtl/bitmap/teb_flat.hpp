#pragma once
//===----------------------------------------------------------------------===//
#include <string>

#include "teb.hpp"
#include "teb_types.hpp"
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
class teb_iter;
class teb_scan_iter;
class teb_wrapper;
//===----------------------------------------------------------------------===//
/// The TEB logic that is used to access a serialized TEB.
class teb_flat {

  friend class teb_iter;
  friend class teb_scan_iter;
  friend class teb_wrapper;

  using word_type = teb_word_type;
  static constexpr auto word_size = sizeof(word_type);
  static constexpr auto word_bitlength = word_size * 8;
  
  using rank_type = teb_rank_type;
  using size_type = teb_size_type;
  using bitmap_fn = dtl::bitmap_fun<word_type>;

  const word_type* const ptr_;
  const teb_header* const hdr_;

  const size_type n_;
  const size_type tree_height_;
  const size_type implicit_inner_node_cnt_;
  const size_type implicit_leading_label_cnt_;
  const size_type perfect_level_cnt_;
  const size_type encoded_tree_height_;

  const word_type* tree_ptr_;
  const size_type tree_bit_cnt_;
  const dtl::bitmap_view<const word_type> T_;

  const word_type* label_ptr_;
  const size_type label_bit_cnt_;
  const dtl::bitmap_view<const word_type> L_;

  const size_type* rank_lut_ptr_;
  rank_type rank_;

  const size_type* level_offsets_tree_lut_ptr_;
  const size_type* level_offsets_labels_lut_ptr_;

public:

  /// Returns the pointer to the header.
  static constexpr const teb_header* const
  get_header_ptr(const word_type* const ptr) {
    return reinterpret_cast<const teb_header* const>(ptr);
  }

  /// Returns the length of the header in number of words.
  static constexpr std::size_t
  get_header_word_cnt(const word_type* const ptr) {
    return sizeof(teb_header) / word_size;
  }

  /// Returns the pointer to the tree structure, or NULL if no tree exists.
  static constexpr const word_type* const
  get_tree_ptr(const word_type* const ptr) {
    const auto* hdr = get_header_ptr(ptr);
    return (hdr->tree_bit_cnt != 0)
        ? ptr + get_header_word_cnt(ptr)
        : nullptr;
  }
  static constexpr word_type*
  get_tree_ptr(word_type* ptr) {
    auto* hdr = get_header_ptr(ptr);
    return (hdr->tree_bit_cnt != 0)
        ? ptr + get_header_word_cnt(ptr)
        : nullptr;
  }

  /// Returns the length of the tree in number of words.
  static constexpr std::size_t
  get_tree_word_cnt(const word_type* const ptr) {
    const auto* hdr = get_header_ptr(ptr);
    return (hdr->tree_bit_cnt + word_bitlength - 1) / word_bitlength;
  }

  /// Returns the pointer to the rank structure, or NULL if it does not exist.
  static constexpr const size_type* const
  get_rank_ptr(const word_type* const ptr) {
    const auto* hdr = get_header_ptr(ptr);
    const auto hdr_word_cnt = get_header_word_cnt(ptr);
    return (hdr->tree_bit_cnt > 1024)
        ? reinterpret_cast<const size_type* const>(ptr + get_header_word_cnt(ptr) + get_tree_word_cnt(ptr))
        : nullptr;
  }
  static constexpr size_type*
  get_rank_ptr(word_type* ptr) {
    auto* hdr = get_header_ptr(ptr);
    auto hdr_word_cnt = get_header_word_cnt(ptr);
    return (hdr->tree_bit_cnt > 1024)
        ? reinterpret_cast<size_type*>(ptr + get_header_word_cnt(ptr) + get_tree_word_cnt(ptr))
        : nullptr;
  }

  /// Returns the length of the rank helper structure in number of words.
  static constexpr std::size_t
  get_rank_word_cnt(const word_type* const ptr) {
    const auto* hdr = get_header_ptr(ptr);
    const auto rank_size_bytes = rank_type::estimate_size_in_bytes(hdr->tree_bit_cnt);
    return (hdr->tree_bit_cnt > 1024)
        ? (rank_size_bytes + word_size - 1) / word_size
        : 0;
  }

  /// Returns the pointer to the labels, or NULL if it does not exist.
  static constexpr const word_type* const
  get_label_ptr(const word_type* const ptr) {
    const auto* hdr = get_header_ptr(ptr);
    const auto rank_word_cnt = get_rank_word_cnt(ptr);
    return (hdr->label_bit_cnt != 0)
        ? ptr + get_header_word_cnt(ptr) + get_tree_word_cnt(ptr) + get_rank_word_cnt(ptr)
        : nullptr;
  }
  static constexpr word_type*
  get_label_ptr(word_type* ptr) {
    auto* hdr = get_header_ptr(ptr);
    auto rank_word_cnt = get_rank_word_cnt(ptr);
    return (hdr->label_bit_cnt != 0)
        ? ptr + get_header_word_cnt(ptr) + get_tree_word_cnt(ptr) + get_rank_word_cnt(ptr)
        : nullptr;
  }

  /// Returns the length of the labels in number of words.
  static constexpr std::size_t
  get_label_word_cnt(const word_type* const ptr) {
    const auto* hdr = get_header_ptr(ptr);
    return (hdr->label_bit_cnt + word_bitlength - 1) / word_bitlength;
  }

  /// Returns the pointer to the additional meta data, or NULL if it does not exist.
  static constexpr const word_type* const
  get_metadata_ptr(const word_type* const ptr) {
    const auto* hdr = get_header_ptr(ptr);
    return (hdr->has_level_offsets != 0)
        ? ptr + get_header_word_cnt(ptr) + get_tree_word_cnt(ptr)
            + get_rank_word_cnt(ptr) + get_label_word_cnt(ptr)
        : nullptr;
  }
  static constexpr word_type*
  get_metadata_ptr(word_type* ptr) {
    auto* hdr = get_header_ptr(ptr);
    return (hdr->has_level_offsets != 0)
        ? ptr + get_header_word_cnt(ptr) + get_tree_word_cnt(ptr)
            + get_rank_word_cnt(ptr) + get_label_word_cnt(ptr)
        : nullptr;
  }

  /// Returns the length of the additional meta data in number of words.
  static constexpr std::size_t
  get_metadata_word_cnt(const word_type* const ptr) {
    const auto* hdr = get_header_ptr(ptr);
    const auto entry_cnt = hdr->encoded_tree_height - hdr->perfect_level_cnt;
    const auto size_in_bytes = 2 * sizeof(size_type) * entry_cnt;
    const auto word_cnt = (size_in_bytes + word_size - 1) / word_size;
    return hdr->has_level_offsets ? word_cnt : 0;
  }

  /// Construct a TEB instance that provides all the logic to work with a
  /// serialized TEB.
  explicit
  teb_flat(const word_type* const ptr)
      : ptr_(ptr),
        hdr_(get_header_ptr(ptr)),
        n_(get_header_ptr(ptr)->n),
        tree_height_(dtl::log_2(get_header_ptr(ptr)->n)),
        implicit_inner_node_cnt_(get_header_ptr(ptr)->implicit_inner_node_cnt),
        implicit_leading_label_cnt_(get_header_ptr(ptr)->implicit_leading_label_cnt),
        perfect_level_cnt_(get_header_ptr(ptr)->perfect_level_cnt),
        encoded_tree_height_(get_header_ptr(ptr)->encoded_tree_height),
        tree_ptr_(get_tree_ptr(ptr) != nullptr ? get_tree_ptr(ptr) : &teb_null_word),
        T_((get_tree_ptr(ptr) != nullptr)
                ? get_tree_ptr(ptr)
                : &teb_null_word,
           (get_tree_ptr(ptr) != nullptr)
                ? get_tree_ptr(ptr) + teb_flat::get_tree_word_cnt(ptr)
                : &teb_null_word + 1),
        tree_bit_cnt_((get_tree_ptr(ptr) != nullptr)
                ? get_header_ptr(ptr)->tree_bit_cnt
                : 1),
        rank_lut_ptr_(get_rank_ptr(ptr)),
        label_ptr_(get_label_ptr(ptr)),
        L_(get_label_ptr(ptr), get_label_ptr(ptr) + teb_flat::get_label_word_cnt(ptr)),
        label_bit_cnt_(get_header_ptr(ptr)->label_bit_cnt),
        level_offsets_tree_lut_ptr_(
            reinterpret_cast<const size_type*>(get_metadata_ptr(ptr))),
        level_offsets_labels_lut_ptr_(
            reinterpret_cast<const size_type*>(get_metadata_ptr(ptr))
                + get_header_ptr(ptr)->encoded_tree_height
                - get_header_ptr(ptr)->perfect_level_cnt) {
        // TODO maybe it is not necessary to copy the header data

    // Initialize rank helper structure.
    if (rank_lut_ptr_ == nullptr) {
      // The TEB instance does not contain a rank LuT. This typically happens
      // when the tree structure is very small. In that case the LuT is
      // initialized on the fly.
      rank_.init(tree_ptr_, tree_ptr_ + get_tree_word_cnt(ptr));
      rank_lut_ptr_ = rank_.lut.data();
    }
  }

  teb_flat(const teb_flat& other) = default;
  teb_flat(teb_flat&& other) noexcept = default;
  teb_flat& operator=(const teb_flat& other) = default;
  teb_flat& operator=(teb_flat&& other) = default;
  ~teb_flat() = default;

  /// Return the name of the implementation.
  static std::string
  name() noexcept {
    return "teb_flat";
  }

  /// Returns the length of the original bitmap.
  std::size_t
  size() const noexcept {
    return hdr_->n;
  }

  /// Returns the value of the bit at the given position.
  u1 __teb_inline__
  test(const std::size_t pos) const noexcept {
    size_type level = perfect_level_cnt_ - 1;
    const auto foo = pos >> (tree_height_ - level);

    // Determine the top-node idx.
    const auto top_node_idx_begin = (1ull << (perfect_level_cnt_ - 1)) - 1;
    size_type node_idx = top_node_idx_begin + foo;

    auto i = tree_height_ - 1 - level;
    while (!is_leaf_node(node_idx)) {
      u1 direction_bit = dtl::bits::bit_test(pos, i);
      node_idx = 2 * rank_inclusive(node_idx) - 1; // left child
      node_idx += direction_bit; // right child if bit is set, left child otherwise
      --i;
    }
    return get_label(node_idx);
  }

  /// Return the size in bytes.
  std::size_t __teb_inline__
  size_in_byte() const noexcept {
    std::size_t word_cnt = 0;
    word_cnt += get_header_word_cnt(ptr_);
    word_cnt += get_tree_word_cnt(ptr_);
    word_cnt += get_rank_word_cnt(ptr_);
    word_cnt += get_label_word_cnt(ptr_);
    word_cnt += get_metadata_word_cnt(ptr_);
    return word_cnt * word_size;
  }

  /// For debugging purposes.
  void __forceinline__
  print(std::ostream& os) const noexcept {
    os << "implicit inner nodes = "
        << implicit_inner_node_cnt_
        << ", implicit leading labels = "
        << implicit_leading_label_cnt_
        << ", perfect levels = "
        << perfect_level_cnt_
        << ", tree bits = " << tree_bit_cnt_
        << ", label bits = " << label_bit_cnt_
        << ", n = " << n_
        << ", encoded tree height = " << encoded_tree_height_
        << ", rank size = " << (get_rank_word_cnt(ptr_) * word_size)
        << ", size = " << size_in_byte()
        << "\n | ";

    if (implicit_inner_node_cnt_ > 0) {
      os << "'";
    }
    for ($i64 i = 0; i < tree_bit_cnt_; i++) {
      os << (T_[i] ? "1" : "0");
    }
    os << "\n | ";
    if (implicit_leading_label_cnt_ > 0) {
      os << "'";
    }
    for ($i64 i = 0; i < label_bit_cnt_; i++) {
      os << (L_[i] ? "1" : "0");
    }
    os << "\n";
    if (has_level_offsets()) {
      os << "Level offsets (tree):   ";
      for (size_type l = 0; l < encoded_tree_height_; ++l) {
        os << l << "=" << get_level_offset_tree(l) << " ";
      }
      os << "\n";
      os << "Level offsets (labels): ";
      for (size_type l = 0; l < encoded_tree_height_; ++l) {
        os << l << "=" << (l < perfect_level_cnt_ ? 0 : get_level_offset_labels(l)) << " ";
      }
    }
  }

private:

  /// Computes the (inclusive) rank of the given tree node.
  size_type __teb_inline__
  rank_inclusive(size_type node_idx) const noexcept {
    const auto implicit_1bit_cnt = implicit_inner_node_cnt_;
    if (node_idx < implicit_1bit_cnt) {
      return node_idx + 1;
    }
    assert(tree_bit_cnt_ > 0);
    const auto i = std::min(node_idx - implicit_1bit_cnt, tree_bit_cnt_ - 1);
    const auto r = rank_type::get(rank_lut_ptr_, i, tree_ptr_);
    const auto ret_val = implicit_1bit_cnt + r;
    return ret_val;
  }

  /// Returns true if the given node is an inner node, false otherwise.
  u1 __teb_inline__
  is_inner_node(size_type node_idx) const noexcept {
    const auto implicit_1bit_cnt = implicit_inner_node_cnt_;
    const auto implicit_leaf_begin =
        tree_bit_cnt_ + implicit_inner_node_cnt_;
    if (node_idx < implicit_1bit_cnt) {
      // Implicit inner node.
      return true;
    }
    if (node_idx >= implicit_leaf_begin) {
      // Implicit leaf node.
      return false;
    }
    return bitmap_fn::test(tree_ptr_, node_idx - implicit_1bit_cnt);
  }

  /// Returns true if the given node is a leaf node, false otherwise.
  u1 __teb_inline__
  is_leaf_node(size_type node_idx) const noexcept {
    return !is_inner_node(node_idx);
  }

  /// Returns the label index for the given node.
  size_type __teb_inline__
  get_label_idx(size_type node_idx) const noexcept {
    return node_idx - rank_inclusive(node_idx);
  }

  /// Returns the label at index.
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

  /// Returns the label of the given node.
  u1 __teb_inline__
  get_label(size_type node_idx) const noexcept {
    const auto label_idx = get_label_idx(node_idx);
    return get_label_by_idx(label_idx);
  }

  /// Returns true if this instance contains the offsets of the individual
  /// tree levels. This is required by the tree scan iterator. Typically, the
  /// offsets aren't present when the tree structure is very small.
  u1 __teb_inline__
  has_level_offsets() const noexcept {
    return hdr_->has_level_offsets;
  }

  /// Returns the position in T where the first tree node of level is stored.
  size_type __teb_inline__
  get_level_offset_tree(size_type level) const noexcept {
    assert(has_level_offsets());
    assert(level < encoded_tree_height_);
    if (level < perfect_level_cnt_) {
      return dtl::binary_tree_structure::first_node_idx_at_level(level);
    }
    else {
      return level_offsets_tree_lut_ptr_[level - perfect_level_cnt_];
    }
  }

  /// Returns the position in L where the first label of level is stored.
  size_type __teb_inline__
  get_level_offset_labels(size_type level) const noexcept {
    assert(has_level_offsets());
    assert(level < encoded_tree_height_);
    if (level < perfect_level_cnt_) {
      return 0;
    }
    else {
      return level_offsets_labels_lut_ptr_[level - perfect_level_cnt_];
    }
  }

};
//===----------------------------------------------------------------------===//
} // namespace dtl
