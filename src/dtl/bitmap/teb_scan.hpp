#pragma once

#include <bitset>
#include <list>
#include <queue>
#include <stack>
#include <vector>

#include "boost/dynamic_bitset.hpp"

#include <dtl/dtl.hpp>
#include <dtl/bits.hpp>
#include <dtl/bitmap/util/rank1.hpp>
#include <dtl/bitmap/util/rank1_surf.hpp>
#include <dtl/bitmap/util/rank1_surf_cached.hpp>
#include <dtl/bitmap/util/bit_buffer.hpp>
#include <dtl/bitmap/util/bitmap_tree.hpp>
#include <dtl/bitmap/util/binary_tree_structure.hpp>
#include <dtl/static_stack.hpp>
#include <dtl/static_stack2.hpp>
#include <dtl/bitmap/util/rank1_naive.hpp>
#include <dtl/bitmap/util/rank1_interceptor.hpp>
#include <dtl/bitmap/util/rank1_super_fast.hpp>
#include <dtl/bitmap/util/bitmap_view.hpp>
#include <dtl/math.hpp>

#include "teb_scan_util.hpp"

namespace dtl {

// TODO remove
//#define __teb_inline__ __attribute__((noinline))

#if !defined(__teb_inline__)
#if defined(NDEBUG)
// Release build.
#define __teb_inline__ inline __attribute__((always_inline))
//#define __teb_inline__ __attribute__((noinline))
#else
#define __teb_inline__
#endif
#endif

//===----------------------------------------------------------------------===//
/// Encodes a bitmap of length n as a binary tree.
///
/// The implementation supports the optimization levels 0, 1, and 2.
///   0 = implements the core idea of TEBs
///   1 = uses the implicit inner node optimization
///   2 = additionally applies the gradual decompression optimization.
/// Please note, that the different optimization levels can be chosen due to
/// academic reasons only. I.e., to analyse and compare the different
/// optimizations levels.
///
/// Terminology:
/// - node index: the nodes identifier, which is its position within the
///       bit sequence of the encoded tree.
/// - path: also identifies a tree node, but it encodes the path to the node,
///       starting at the root
/// - perfect levels: the number of tree levels in which the tree is perfect.
/// - top nodes: refer to the nodes within the last perfect level.
///
template<i32 optimization_level_ = 2>
class teb_scan {

public:

  // The fundamental storage type. The size of a TEB is a multiple of
  // sizeof(_block_type).
  using _block_type = $u64;
  using word_type = _block_type;

  using position_t = $u32;
  using bitmap_t = boost::dynamic_bitset<word_type>;

  /// The number of bits in the bitmap.
  $u64 n_;

  /// The tree encoded bitmap.
  bitmap_t structure_;
  bitmap_t labels_;

  bitmap_view<word_type> T_;
  $u64 structure_bit_cnt_;
  bitmap_view<word_type> L_;

  /// Support data structure for rank1 operations on the tree structure.
  static constexpr u1 inclusive = true;
  static constexpr u1 non_inclusive = false;
//  using rank_support = dtl::rank1_naive<word_type>;
//  using rank_support = dtl::rank1_surf_cached<word_type>;
//  using rank_support = dtl::rank1<word_type>;
//  using rank_support = dtl::rank1_interceptor<dtl::rank1_surf_cached<word_type, non_inclusive>>;
//  using rank_support = dtl::rank1_interceptor<dtl::rank1_surf_cached<word_type, inclusive>>;
//  using rank_support = dtl::rank1_interceptor<dtl::rank1_surf<word_type, inclusive>>;
//  using rank_support = dtl::rank1_surf_cached<word_type, inclusive>;
//  using rank_support = dtl::rank1_super_fast<word_type, inclusive>;
//  using rank_support = dtl::rank1_naive<word_type>;
  using rank_support = dtl::rank1_surf<word_type, inclusive>;
  rank_support rank_;

  /// The number of implicit inner nodes in the tree structure.
  $u32 implicit_inner_node_cnt_;

  /// For testing purposes only.
  $u32 implicit_leaf_node_cnt_;

  std::array<std::size_t, 32> level_offsets_structure_;
  std::array<std::size_t, 32> level_offsets_labels_;

  $u32 encoded_tree_height_ = 0;

public:

  /// Tree-encode the given bitmap with an optional false positive rate.
  explicit
  teb_scan(const boost::dynamic_bitset<$u32>& bitmap, f64 fpr = 0.0)
      : n_(bitmap.size()) {

    // Construct a binary tree that represents the given bitmap.
    // Space-optimizations are performed in the (non-encoded) bitmap tree.
    dtl::bitmap_tree<optimization_level_> bitmap_tree(bitmap, fpr);

    // Encode the tree into level-order.
    implicit_inner_node_cnt_ = 0;
    implicit_leaf_node_cnt_ = 0;
    if (optimization_level_ > 0) {
      // Implicit tree nodes are only considered with -o1 or higher.
      implicit_inner_node_cnt_ = bitmap_tree.get_leading_inner_node_cnt();
      implicit_leaf_node_cnt_ = bitmap_tree.get_trailing_leaf_node_cnt();
    }
    std::size_t node_cntr = 0;
    std::size_t leaf_node_cntr = 0;
    std::size_t current_level = ~0ull;
    for (auto it = bitmap_tree.breadth_first_begin();
         it != bitmap_tree.breadth_first_end();
         ++it) {
      ++node_cntr;

      u64 idx = (*it).idx;
      u64 level = (*it).level;

      if (current_level != level) {
        level_offsets_structure_[level] = node_cntr - 1;
        level_offsets_labels_[level] = leaf_node_cntr;
        current_level = level;
      }

      // Omit the implicit nodes.
      if (node_cntr <= implicit_inner_node_cnt_) {
        continue;
      }

      u1 is_inner = bitmap_tree.is_inner_node(idx);
      leaf_node_cntr += !is_inner;

      // Emit a 1-bit if the current node is an inner node, a 0-bit otherwise.
      if (is_inner
          || optimization_level_ == 0
          || (optimization_level_ > 0
              && idx <= bitmap_tree.get_last_explicit_node_idx())) {
        structure_.push_back(is_inner);
      }
      if (!is_inner) {
        // Add the label of the leaf node.
        labels_.push_back(bitmap_tree.label_of_node(idx));
      }
    }

    // TODO do we really need to store this?
    encoded_tree_height_ = current_level + 1;

    // Hack: ensure structure is not null
    if (structure_.size() == 0) {
      structure_.push_back(false);
    }

    // Hack: add an additional word to the tree structure, which helps to
    //       eliminate a branch in the scan iterator.
//    for (std::size_t i = 0; i < (sizeof(word_type) * 8); ++i) {
//      structure_.push_back(false);
//    }

    // Init rank1 support data structure.
    rank_.init(structure_);

    T_.init(data_view<const word_type> {
        structure_.m_bits.data(),
        structure_.m_bits.data() + structure_.m_bits.size(),
    });
    structure_bit_cnt_ = structure_.size();
    L_.init(data_view<const word_type> {
        labels_.m_bits.data(),
        labels_.m_bits.data() + labels_.m_bits.size(),
    });
  }

  teb_scan(const teb_scan& other) = default;
  teb_scan(teb_scan&& other) noexcept = default;
  teb_scan& operator=(const teb_scan& other) = default;
  teb_scan& operator=(teb_scan&& other) = default;
  ~teb_scan() = default;

  /// Decodes the level-order encoding to a bitmap.
  boost::dynamic_bitset<$u32> __teb_inline__
  to_bitset() const noexcept {
    boost::dynamic_bitset<$u32> ret(n_); // the resulting bitmap

    // special case if the tree is only a root node
    if(implicit_inner_node_cnt_ == 0 && structure_.size() <= 1){
      //std::cout << "to_bitset() only root case" << std::endl;
      if(labels_[0]) {
        ret.flip();
      }
      return ret;
    }

    // normal cases: construct the tree using the lo_construction and a level counter
    $u64 level = 0; // current height
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
  std::size_t __teb_inline__ // TODO include level offsets
  size_in_byte() const noexcept {
    constexpr u64 block_bitlength = sizeof(_block_type) * 8;
    constexpr u64 block_size = sizeof(_block_type);
    $u64 bytes = 0;
    // Tree structure
    bytes += ((structure_.size() + block_bitlength - 1) / block_bitlength) * block_size;
    // Labels
    bytes += ((labels_.size() + block_bitlength - 1) / block_bitlength) * block_size;
    // Rank support
    bytes += rank_.size_in_bytes();
    // Bit-length of the original bitmap.
    bytes += sizeof(n_);
//    if (optimization_level_ > 0) {
      // The number of implicit inner nodes.
      bytes += sizeof(implicit_inner_node_cnt_);
//    }
    // FIXME: + number of tree bits + number of label bits
    return bytes;
  }

  u1 __teb_inline__
  operator!=(teb_scan& other) const noexcept {
    return !(*this == other);
  }

  u1 __teb_inline__
  operator==(teb_scan& other) const noexcept {
    return implicit_inner_node_cnt_ == other.implicit_inner_node_cnt_
        && structure_ == other.structure_
        && labels_ == other.labels_;
  }

  void
  print(std::ostream& os) const noexcept {
    os << "implicit nodes (internal/external) = "
       << implicit_inner_node_cnt_
       << "/"
       << implicit_leaf_node_cnt_
       << ", perfect levels = "
       << determine_perfect_tree_levels(implicit_inner_node_cnt_)
       << ", tree bits = " << structure_.size()
       << ", label bits = " << labels_.size()
       << ", opt level = " << optimization_level_
       << ", n = " << n_
       << ", rank size = " << rank_.size_in_bytes()
       << ", size = " << size_in_byte()
       << "\n | ";

    if (implicit_inner_node_cnt_ > 0) {
//      for ($i64 i = 0; i < implicit_inner_node_cnt_; i++) {
//        os << "1";
//      }
      os << "'";
    }
    for ($i64 i = 0; i < structure_.size(); i++) {
      os << (structure_[i] ? "1" : "0");
    }
//    os << "'";
//    for ($i64 i = implicit_inner_node_cnt_ + structure_.size();
//         i < n_; i++) {
//      os << "0";
//    }
    os << "\n | ";
    for ($i64 i = 0; i < labels_.size(); i++) {
      os << (labels_[i] ? "1" : "0");
    }
  }

  /// Return the name of the implementation.
  static std::string
  name() noexcept {
    return "teb_scan";
  }

  /// Returns true if all bits are set, false otherwise.
  u1 __teb_inline__
  all() noexcept {
    // FIXME: this works only if the tree mask is in a compressed state
    return is_leaf_node(0) // root is the only node (a leaf)
        && get_label(0) == true; // and the label is 1
  }

  /// Returns true if all bits are zero, false otherwise.
  u1 __teb_inline__
  none() noexcept {
    // FIXME: this works only if the tree mask is in a compressed state
    return is_leaf_node(0) // root is the only node (a leaf)
        && get_label(0) == false; // and the label is 0
  }

  /// Returns the length of the original bitmap.
  std::size_t
  size() const noexcept {
    return n_;
  }

  //===--------------------------------------------------------------------===//
  // Helper functions.
  //===--------------------------------------------------------------------===//
  /// Computes the number of perfect levels in the tree structure based on the
  /// number of implicit inner nodes. The number of perfect levels is at least
  /// one, because there is at least one node in the tree structure.
  static inline u64
  determine_perfect_tree_levels(u64 implicit_inner_node_cnt) noexcept {
    return dtl::log_2(implicit_inner_node_cnt + 1) + 1;
  }

  /// Computes the height of the tree based on n. - Note that a tree, consisting
  /// of a single (root) node has a height of 0.
  static inline u64
  determine_tree_height(u64 n) noexcept {
    return dtl::log_2(n);
  }

  /// Determines the path to the common ancestor of the two nodes specified
  /// by there paths.
  static inline u64
  determine_common_ancestor_path(u64 src_path, u64 dst_path) {
    //TODO should use positions instead of paths (at least for the second argument)
    assert(src_path != dst_path);
    const auto a = src_path << (dtl::bits::lz_count(src_path) + 1);
    const auto b = dst_path << (dtl::bits::lz_count(dst_path) + 1);
    assert(a < b);
    const auto src_path_len = sizeof(src_path) * 8
        - (dtl::bits::lz_count(src_path) + 1);
    const auto common_prefix_length = dtl::bits::lz_count(a ^ b);
    const auto common_ancestor_path =
        src_path >> (src_path_len - common_prefix_length);
    return common_ancestor_path;
  }

  static inline u64
  determine_level_of(u64 path) {
    const auto lz_cnt_path = dtl::bits::lz_count(path);
    const auto level = sizeof(u64) * 8 - 1 - lz_cnt_path;
    return level;
  }
  //===--------------------------------------------------------------------===//

  //===--------------------------------------------------------------------===//
  /// 1-fill iterator, with skip support.
  class iter {

    /// The fundamental type to encode paths within the tree.
    using path_t = $u64;

    static constexpr u64 DEFAULT_BATCH_SIZE = 128;

    struct range_t {
      $u64 pos;
      $u64 length;
    };

    struct scanner_state_t {
      $u64 node_idx_;
      $u64 label_idx_;

      void __teb_inline__
      init(const iter& /*iter*/, u64 node_idx, u64 label_offset) {
        node_idx_ = node_idx;
        label_idx_ = label_offset;
      }

        void
      print(std::ostream& os) const noexcept {
        os << "scan_state["
           << "node_idx=" << node_idx_
           << ",label_idx=" << label_idx_ << "]";
      }
    };

    /// Reference to the TEB instance.
    const teb_scan& teb_;
    /// The height of the tree structure. // TODO refers to max height
    u64 tree_height_;
    /// The number of perfect levels in the tree structure.
    u64 perfect_levels_;
    /// First node index in the last perfect level.
    u64 top_node_idx_begin_;
    /// Last node index + 1 in the last perfect level.
    u64 top_node_idx_end_;
    /// The current top node.
    $u64 top_node_idx_current_;

    /// One scanner (or stream) per level.
    std::array<scanner_state_t, 32> scanner_states_;

    /// A batch of results.
    std::vector<range_t> results_;
    $u64 result_cnt_;
    u64 batch_size_;
    $u64 result_read_pos_;

    $u64 scan_pos_;
    $u64 scan_length_;
    $u64 scan_path_;
    $u64 scan_path_level_;

    $u32 alpha_;
    //===------------------------------------------------------------------===//


    //===------------------------------------------------------------------===//
    // Helper functions.
    //===------------------------------------------------------------------===//
    auto __teb_inline__
    determine_level(u32 alpha) {
      // LSB is always set. Thus, lz_count(alpha) is never undefined.
      return sizeof(alpha) * 8 - dtl::bits::lz_count(alpha);
    };

    auto __teb_inline__
    is_left_child(path_t path) {
      return (path & path_t(1)) == 0;
    };
    //===------------------------------------------------------------------===//

  public:

    /// Constructs an iterator for the given TEB instance. After construction,
    /// the iterator points to the first 1-fill.
    explicit __teb_inline__
    iter(const teb_scan& teb_scan, u64 batch_size = DEFAULT_BATCH_SIZE) noexcept
        : teb_(teb_scan),
          tree_height_(determine_tree_height(teb_scan.n_)),
          perfect_levels_(
             determine_perfect_tree_levels(teb_.implicit_inner_node_cnt_)),
          top_node_idx_begin_((1ull << (perfect_levels_ - 1)) - 1),
          top_node_idx_end_((1ull << perfect_levels_) - 1),
          top_node_idx_current_(top_node_idx_begin_),
          results_(batch_size),
          result_cnt_(0),
          result_read_pos_(0),
          batch_size_(std::min(batch_size, teb_.size())),
          alpha_(0) {
      results_.reserve(batch_size);
      if (teb_.encoded_tree_height_ == 1) {
        scan_pos_ = teb_.n_;
        scan_length_ = 0;
        result_cnt_ = 1;
        if (teb_.labels_[0]) {
          results_[0].pos = 0;
          results_[0].length = teb_.n_;
        }
        else {
          results_[0].pos = teb_.n_;
          results_[0].length = 0;
        }
        return;
      }
      // Initialize the scanners.
      for (std::size_t level = 0; level < teb_.encoded_tree_height_; ++level) {
        scanner_states_[level].init(*this,
                              teb_.level_offsets_structure_[level],
                              teb_.level_offsets_labels_[level]);
      }
      // Initialize alpha vector.
      alpha_ = 0;
      for (std::size_t level = 0; level < teb_.encoded_tree_height_; ++level) {
        u64 node_idx = teb_.level_offsets_structure_[level];
        u32 node_bit = teb_.is_inner_node(node_idx) ? 1u : 0u;
        alpha_ |= node_bit << level;
      }
      // Populate pos and length.
      u32 level =  dtl::bits::tz_count(~alpha_);
      scan_pos_ = 0;
      scan_length_ = teb_.n_ >> level;
      scan_path_ = path_t(0);
      scan_path_level_ = level;
      // Iterator is now positioned at the first leaf node.
      const auto first_leaf_node_idx = teb_.level_offsets_structure_[level];
      // Determine the label.
      u1 first_leaf_label = teb_.get_label(first_leaf_node_idx);
      if (first_leaf_label) {
        results_[0].pos = 0;
        results_[0].length = scan_length_;
        ++result_cnt_;
      }
      next_batch();
    }

    __teb_inline__
    iter(iter&&) noexcept = default;
    iter(const iter& other) = default;
    iter& operator=(const iter& other) = default;
    iter& operator=(iter&& other) = default;
    ~iter() = default;

    /// Forwards the iterator to the next 1-fill (if any).
    /// Use the functions pos() and length() to get the 1-fill the iterator
    /// is currently pointing to.
    void __teb_inline__
    next() noexcept __attribute__ ((flatten, hot)) {
      assert(!end());
      ++result_read_pos_;
      if (result_read_pos_ == result_cnt_) {
        result_read_pos_ = 0;
        result_cnt_ = 0;
        next_batch();
      }
    }

    void //__teb_inline__
    next_batch() noexcept __attribute__ ((flatten, hot, noinline)) {
      const auto h = tree_height_;
      const auto n = teb_.size();
      auto pos = scan_pos_;
      auto length = scan_length_;
      if (pos + length >= n) {
        results_[result_cnt_].pos = n;
        results_[result_cnt_].length = 0;
        ++result_cnt_;
        return;
      }
      // Note: The path variable does NOT contain a sentinel bit. Instead, we
      //       keep track of the paths' level using a separate variable.
      register auto path = scan_path_;
      register auto path_level = scan_path_level_;
      auto result_cnt = result_cnt_;
      auto alpha = alpha_;

      // Initialize the scanners.
      const data_view<const word_type> T {
          teb_.structure_.m_bits.data(),
          teb_.structure_.m_bits.data() + teb_.structure_.m_bits.size(),
      };
      const data_view<const word_type> L {
          teb_.labels_.m_bits.data(),
          teb_.labels_.m_bits.data() + teb_.labels_.m_bits.size(),
      };

      // Bit buffers for the tree structure.
      dtl::bit_buffer<8> t_bb0;
      dtl::bit_buffer<8> t_bb1;
      dtl::bit_buffer<8> t_bb2;
      dtl::bit_buffer<8> t_bb3;

      // Bit buffers for the labels.
      dtl::bit_buffer<8> l_bb0;
      dtl::bit_buffer<8> l_bb1;
      dtl::bit_buffer<8> l_bb2;
      dtl::bit_buffer<8> l_bb3;

      // Perfect levels. (all bits set to 1)
      for (std::size_t level = 0; level < perfect_levels_ - 1; ++level) {
        u64 buffer_idx = level / 8;
        u64 slot_idx = level % 8;
        switch (buffer_idx) {
          case 0: t_bb0.set(slot_idx, ~0ul); l_bb0.set(slot_idx, ~0ul); break;
          case 1: t_bb1.set(slot_idx, ~0ul); l_bb1.set(slot_idx, ~0ul); break;
          case 2: t_bb2.set(slot_idx, ~0ul); l_bb2.set(slot_idx, ~0ul); break;
          case 3: t_bb3.set(slot_idx, ~0ul); l_bb3.set(slot_idx, ~0ul); break;
        }
      }
      // --

      // Advance the tree structure bit buffers.
      auto advance_tree_scanners = [&](std::size_t level_begin,
                                       std::size_t level_end) {
        static constexpr u64 p = (1ul << 8) - 1;
        u64 a = (~0ul) << level_begin;
        u64 b = (~0ul) >> (64 - level_end);
        u64 m = a & b;
        D(std::cout << std::bitset<32>(m) << std::endl;)
        u64 m0 = (m >> (8 * 0)) & p; t_bb0.increment(m0); u64 r0 = t_bb0.read();
        u64 m1 = (m >> (8 * 1)) & p; t_bb1.increment(m1); u64 r1 = t_bb1.read();
        u64 m2 = (m >> (8 * 2)) & p; t_bb2.increment(m2); u64 r2 = t_bb2.read();
        u64 m3 = (m >> (8 * 3)) & p; t_bb3.increment(m3); u64 r3 = t_bb3.read();
        // Update the alpha vector.
        alpha = (r3 << (8 * 3)) | (r2 << (8 * 2)) | (r1 << (8 * 1)) | r0;
        D(std::cout << "alpha: " << std::bitset<32>(alpha) << std::endl;)
      };

      // Advance the label bit buffers.
      auto advance_label_scanner = [&](u32 level) {
        u64 buffer_idx = level / 8;
        u64 slot_idx = level % 8;
        u64 slot_mask = 1ul << slot_idx;
        $u64 label = 0;
        std::cout << "label++" << std::endl;
        switch (buffer_idx) {
          case 0: l_bb0.increment(slot_mask); break;
          case 1: l_bb1.increment(slot_mask); break;
          case 2: l_bb2.increment(slot_mask); break;
          case 3: l_bb3.increment(slot_mask); break;
        }
      };

      while (result_cnt < batch_size_ - 8 && pos < n) {

        // Initialize the bit buffers.
        std::cout << "--fetch--" << std::endl;
        for (std::size_t level = perfect_levels_ - 1;
             level < teb_.encoded_tree_height_; ++level) {
          u64 buffer_idx = level / 8;
          u64 slot_idx = level % 8;
          i64 delta = - static_cast<i64>(teb_.implicit_inner_node_cnt_);
          switch (buffer_idx) {
            case 0: t_bb0.set(slot_idx, fetch_n_bits(T, static_cast<i64>(scanner_states_[level].node_idx_) + delta, 8)); l_bb0.set(slot_idx, fetch_n_bits(L, static_cast<i64>(scanner_states_[level].label_idx_), 8)); break;
            case 1: t_bb1.set(slot_idx, fetch_n_bits(T, static_cast<i64>(scanner_states_[level].node_idx_) + delta, 8)); l_bb1.set(slot_idx, fetch_n_bits(L, static_cast<i64>(scanner_states_[level].label_idx_), 8)); break;
            case 2: t_bb2.set(slot_idx, fetch_n_bits(T, static_cast<i64>(scanner_states_[level].node_idx_) + delta, 8)); l_bb2.set(slot_idx, fetch_n_bits(L, static_cast<i64>(scanner_states_[level].label_idx_), 8)); break;
            case 3: t_bb3.set(slot_idx, fetch_n_bits(T, static_cast<i64>(scanner_states_[level].node_idx_) + delta, 8)); l_bb3.set(slot_idx, fetch_n_bits(L, static_cast<i64>(scanner_states_[level].label_idx_), 8)); break;
          }
        }

        // Do eight iterations. (which is the number of buffered bits per level)
        for ($i32 b = 0; b < 8; ++b) {
          assert(pos <= n);
          assert(pos + length <= n);
          assert(path_level >= 1);
          assert(length > 0);

          // Increment the current position.
          pos += length;
          if (pos == n) goto done; // TODO eliminate branch
          assert(pos < n);

          const auto advance_end = path_level + 1;
          // Walk upwards until a left child is found. (might be the current one)
          const auto up_steps = dtl::bits::tz_count(~path);
          path >>= up_steps;
          path_level -= up_steps;
          // Go to right sibling.
          path = path | 1;

          // Advance the scanners and update the alpha vector.
          const auto advance_begin = path_level;
          advance_tree_scanners(advance_begin, advance_end);
          advance_label_scanner(advance_end - 1);
          // Walk downwards to the left-most leaf in that sub-tree.
          const auto down_steps = dtl::bits::tz_count(~(alpha >> path_level));
          path_level += down_steps;
          path <<= down_steps;

          // Toggle sentinel bit (= highest bit set) and add offset.
          pos = path << (h - path_level); // FIXME remove. pos already set
          // The length of the 1-fill.
          length = n >> path_level;

          // Read the label.
          u64 buffer_idx = path_level / 8;
          u64 slot_idx = path_level % 8;
          u64 label_mask = 1ul << slot_idx;
          $u64 label = 0;
          switch (buffer_idx) {
            case 0: label = l_bb0.read(label_mask) != 0; break;
            case 1: label = l_bb1.read(label_mask) != 0; break;
            case 2: label = l_bb2.read(label_mask) != 0; break;
            case 3: label = l_bb3.read(label_mask) != 0; break;
          }
          // Produce output (a 1-fill).
          std::cout << "pos=" << pos
                    << ", length=" << length
                    << ", label=" << label
                    << std::endl;
          results_[result_cnt].pos = pos;
          results_[result_cnt].length = length;
          result_cnt += label;
        }

        // Update the scanner states.
        auto get_tree_scanner_read_pos = [&](u32 level) {
          u64 buffer_idx = level / 8;
          u64 slot_idx = level % 8;
          $i64 ret_val = 0;
          switch (buffer_idx) {
            case 0: ret_val = t_bb0.get_read_pos(slot_idx); break;
            case 1: ret_val = t_bb1.get_read_pos(slot_idx); break;
            case 2: ret_val = t_bb2.get_read_pos(slot_idx); break;
            case 3: ret_val = t_bb3.get_read_pos(slot_idx); break;
          }
          return ret_val;
        };
        auto get_label_scanner_read_pos = [&](u32 level) {
          u64 buffer_idx = level / 8;
          u64 slot_idx = level % 8;
          $i64 ret_val = 0;
          switch (buffer_idx) {
            case 0: ret_val = l_bb0.get_read_pos(slot_idx); break;
            case 1: ret_val = l_bb1.get_read_pos(slot_idx); break;
            case 2: ret_val = l_bb2.get_read_pos(slot_idx); break;
            case 3: ret_val = l_bb3.get_read_pos(slot_idx); break;
          }
          return ret_val;
        };
        for ($u32 level = perfect_levels_ - 1;
             level < teb_.encoded_tree_height_; ++level) {
          scanner_states_[level].node_idx_ += get_tree_scanner_read_pos(level);
          scanner_states_[level].label_idx_ += get_label_scanner_read_pos(level);
        }
      }

done:
      if (result_cnt < batch_size_
          && pos == n) {
        results_[result_cnt].pos = n;
        results_[result_cnt].length = 0;
        ++result_cnt;
        length = 0;
      }
      scan_pos_ = pos;
      scan_length_ = length;
      scan_path_ = path;
      scan_path_level_ = path_level;
      result_cnt_ = result_cnt;
      alpha_ = alpha;
    }

    /// Forwards the iterator to the given position. - Note: This functions
    /// simply calls next() until the desired position has been reached.
    void __teb_inline__
    skip_to(const std::size_t to_pos) noexcept {
      if (to_pos >= teb_.n_) {
        results_[0].pos = teb_.n_;
        results_[0].length = 0;
        result_cnt_ = 1;
        result_read_pos_ = 0;
        return;
      }

      while (!end() && pos() + length() <= to_pos) {
        next();
      }
      // Adjust the current position and fill-length.
      if (!end() && pos() < to_pos) {
        results_[result_read_pos_].length -=
            to_pos - results_[result_read_pos_].pos;
        results_[result_read_pos_].pos = to_pos;
      }
    }

    /// Returns true if the iterator reached the end, false otherwise.
    u1 __teb_inline__
    end() const noexcept {
      return results_[result_read_pos_].pos == teb_.n_;
    }

    /// Returns the starting position of the current 1-fill.
    u64 __teb_inline__
    pos() const noexcept {
      assert(result_read_pos_ < result_cnt_);
      return results_[result_read_pos_].pos;
    }

    /// Returns the length of the current 1-fill.
    u64 __teb_inline__
    length() const noexcept {
      assert(result_read_pos_ < result_cnt_);
      return results_[result_read_pos_].length;
    }

//    /// Returns the path of the current tree node.
//    u64 __teb_inline__
//    path() const noexcept {
//      return path_;
//    }

  };

  iter __teb_inline__
  it() const noexcept {
    return std::move(iter(*this));
  }
  //===--------------------------------------------------------------------===//

  /// Returns the name of the instance including the most important parameters
  /// in JSON.
  std::string
  info() const noexcept {
//    auto determine_compressed_tree_height = [&]() {
//      auto i = it();
//      $u64 height = 0;
//      while (!i.end()) {
//        const auto h = determine_level_of(i.path());
//        height = std::max(height, h);
//        i.next();
//      }
//      return height;
//    };
    return "{\"name\":\"" + name() + "\""
        + ",\"n\":" + std::to_string(n_)
        + ",\"size\":" + std::to_string(size_in_byte())
        + ",\"tree_bits\":" + std::to_string(structure_.size())
        + ",\"label_bits\":" + std::to_string(labels_.size())
        + ",\"implicit_inner_nodes\":"
          + std::to_string(implicit_inner_node_cnt_)
        + ",\"implicit_leaf_nodes\":"
          + std::to_string(implicit_leaf_node_cnt_)
//        + ",\"tree_height\":"
//          + std::to_string(determine_compressed_tree_height())
        + ",\"perfect_levels\":"
          + std::to_string(determine_perfect_tree_levels(implicit_inner_node_cnt_))
        + ",\"opt_level\":" + std::to_string(optimization_level_)
        + ",\"rank\":" + rank_.info()
        + "}";
  }

private:

  u1 __teb_inline__
  is_inner_node(u64 node_idx) const {
    const std::size_t implicit_1bit_cnt = implicit_inner_node_cnt_;
    const std::size_t implicit_leaf_begin =
        structure_bit_cnt_ + implicit_inner_node_cnt_;
    if (node_idx < implicit_1bit_cnt) {
      // Implicit inner node.
      return true;
    }
    if (node_idx >= implicit_leaf_begin) {
      // Implicit leaf node.
      return false;
    }
    return T_[node_idx - implicit_1bit_cnt];
  }

  u1 __teb_inline__
  is_leaf_node(u64 node_idx) const {
    return !is_inner_node(node_idx);
  }

  u64 __teb_inline__
  left_child(u64 node_idx) const {
    return 2 * rank_inclusive(node_idx) - 1;
  }

  u64 __teb_inline__
  right_child(u64 node_idx) const {
    return 2 * rank_inclusive(node_idx);
  }

  std::size_t __teb_inline__
  get_label_idx(u64 node_idx) const {
    return node_idx - rank_inclusive(node_idx);
  }

  u1 __teb_inline__
  get_label(u64 node_idx) const {
    u64 label_idx = get_label_idx(node_idx);
    return labels_[label_idx];
  }

  u64 __teb_inline__
  rank_inclusive(u64 node_idx) const {
    const std::size_t implicit_1bit_cnt = implicit_inner_node_cnt_;
    if (node_idx < implicit_1bit_cnt) {
      return node_idx + 1;
    }
    assert(structure_bit_cnt_ > 0);
    const auto i = std::min(node_idx - implicit_1bit_cnt, structure_bit_cnt_ - 1); // FIXME DAMN this inclusive rank >:/
    const auto r = rank_(i, T_.data_.begin());
    const auto ret_val = implicit_1bit_cnt + r;
    return ret_val;
  }

};
//===----------------------------------------------------------------------===//
}; // namespace dtl
