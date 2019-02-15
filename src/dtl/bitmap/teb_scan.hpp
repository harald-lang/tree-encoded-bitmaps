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
#include <dtl/bitmap/util/bitmap_tree.hpp>
#include <dtl/bitmap/util/binary_tree_structure.hpp>
#include <dtl/static_stack.hpp>
#include <dtl/static_stack2.hpp>
#include <dtl/bitmap/util/rank1_naive.hpp>
#include <dtl/bitmap/util/rank1_interceptor.hpp>
#include <dtl/bitmap/util/rank1_super_fast.hpp>
#include <dtl/bitmap/util/bitmap_view.hpp>
#include <dtl/math.hpp>

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

    struct scan_t {
//      using cache_word_type = $u64;
      using cache_word_type = word_type;
      static constexpr u64 cache_word_bitlength = (sizeof(cache_word_type) * 8);

      $u64 idx;
      $u64 label_idx;
      $u64 bitmap_cache_idx;

      cache_word_type structure_cache;
      $i64 structure_block_idx;
      cache_word_type label_cache;
      $i64 label_block_idx;

      void
      print(std::ostream& os) const noexcept {
        os << "scan[node_idx=" << idx
           << ",label_idx=" << label_idx
           << ",bitmap_cache_idx=" << bitmap_cache_idx
           << ",structure_block_idx=" << structure_block_idx
           << ",structure_cache=";
        for (std::size_t i = bitmap_cache_idx; i < cache_word_bitlength; ++i) {
          std::cout << (dtl::bits::bit_test(structure_cache, i) ? "1" : "0");
        }
        std::cout << "]";
      }

      void __teb_inline__
      init(const iter& iter, u64 node_idx, u64 label_offset) {
        // Set the node index.
        idx = node_idx;
        label_idx = label_offset;

        // Fetch the corresponding word from the tree structure.
        const std::size_t implicit_1bit_cnt =
            iter.teb_.implicit_inner_node_cnt_;
        const std::size_t implicit_leaf_begin =
            iter.teb_.structure_bit_cnt_ + implicit_1bit_cnt;

        i64 phys_node_idx = i64(node_idx) - implicit_1bit_cnt;
        if (phys_node_idx >= 0) {
          structure_block_idx = phys_node_idx / i64(cache_word_bitlength);
          bitmap_cache_idx = phys_node_idx % cache_word_bitlength;
        }
        else {
          structure_block_idx = (phys_node_idx / i64(cache_word_bitlength)) - 1;
          bitmap_cache_idx = phys_node_idx % cache_word_bitlength;
          if (phys_node_idx % cache_word_bitlength == 0) {
            ++structure_block_idx;
          }
        }


        if (structure_block_idx < 0) {
          // Implicit inner node.
          structure_cache = ~word_type(0);
        }
        else if (structure_block_idx >= iter.teb_.structure_.m_bits.size()) {
          // Implicit leaf node.
          structure_cache = word_type(0);
        }
        else {
          structure_cache = iter.teb_.structure_.m_bits[structure_block_idx];
        }
        assert(bitmap_cache_idx < cache_word_bitlength);
      }

      void __teb_inline__
      advance_fetch(const iter& iter) {
        // The cached bitmap fragment has been fully consumed.
        ++structure_block_idx;
        D(std::cout << "fetch " << structure_block_idx << std::flush;)
        if (structure_block_idx < 0) {
          // Implicit inner nodes.
          D(std::cout << " implicit inner" << std::endl;)
          structure_cache = ~word_type(0);
        }
        else if (structure_block_idx >= iter.teb_.structure_.m_bits.size()) {
          // Implicit leaf nodes.
          D(std::cout << " implicit leaf" << std::endl;)
          structure_cache = word_type(0);
        }
        else {
          structure_cache = iter.teb_.structure_.m_bits[structure_block_idx];
        }
        bitmap_cache_idx = 0;
      }

      void __teb_inline__
      advance(const iter& iter, u1 current_node_is_leaf) {
//        label_idx += !is_inner_node(iter);
        ++idx;
        label_idx += current_node_is_leaf;
        ++bitmap_cache_idx;
        if (unlikely(bitmap_cache_idx == cache_word_bitlength)) {
          advance_fetch(iter);
        }
      }

      /// Returns true if the scanner currently points to an inner node,
      /// false otherwise.
      u1 __teb_inline__
      is_inner_node(const iter& iter) {
        assert(bitmap_cache_idx < (sizeof(word_type) * 8));
        assert(dtl::bits::bit_test(structure_cache, bitmap_cache_idx) == iter.teb_.is_inner_node(idx));
        return dtl::bits::bit_test(structure_cache, bitmap_cache_idx);
//        return iter.teb_.is_inner_node(idx);
      }

      /// Returns the label of the current node. This function must only be
      /// called, when the current node is a leaf node.
      u1 __teb_inline__
      get_label(const iter& iter) {
        assert(!is_inner_node(iter));
        return iter.teb_.labels_[label_idx];
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
    /// One scanner / stream per level.
    std::array<scan_t, 32> scanners_;

    /// A batch of results.
    std::vector<range_t> results_;
    $u64 result_cnt_;
    u64 batch_size_;
    $u64 result_read_pos_;

    $u64 scan_pos_;
    $u64 scan_length_;
    $u64 scan_path_;

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
        scanners_[level].init(*this,
                              teb_.level_offsets_structure_[level],
                              teb_.level_offsets_labels_[level]);
      }
      // Initialize alpha vector.
      alpha_ = 0;
      for (std::size_t level = 0; level < teb_.encoded_tree_height_; ++level) {
        u32 node_bit = scanners_[level].is_inner_node(*this) ? 1u : 0u;
        alpha_ |= node_bit << level;
      }

      // Populate pos and length.
      u32 level =  dtl::bits::tz_count(~alpha_);
      scan_pos_ = 0;
      scan_length_ = teb_.n_ >> level;
      scan_path_ = path_t(1) << level;
      // Iterator is positioned at the first leaf node.
      const auto first_leaf_node_idx = scanners_[level].idx;
      // Determine the label.
      u1 first_leaf_label = scanners_[level].get_label(*this);
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
//      D(std::cout << "next_batch()" << std::endl;)
      D(constexpr auto O = optimization_level_;)
      const auto h = tree_height_;
      const auto n = teb_.size();
      auto pos = scan_pos_;
      auto length = scan_length_;
      register auto path = scan_path_;
      register auto path_level = determine_level_of(path);
      auto result_cnt = result_cnt_;
      auto alpha = alpha_;

      auto advance_scanner = [&](std::size_t i) {
        u1 was_inner_node = dtl::bits::bit_test(alpha, i);
        scanners_[i].advance(*this, !was_inner_node);
        u1 is_inner_node = scanners_[i].is_inner_node(*this);
        u1 x = was_inner_node ^ is_inner_node;
        alpha ^= u32(x) << i;
        D(std::cout <<"["<<O<<"]"<< "scanner " << i << ": " << scanners_[i] << std::endl;)
      };

      while (result_cnt < batch_size_ && pos < n) {
        assert(pos <= n);
        assert(pos + length <= n);
        assert(path >= 1);
        assert(length > 0);

        // Increment the current position.
        pos += length;
        if (pos == n) break; // TODO eliminate branch
        assert(pos < n);

        const auto advance_end = path_level + 1;

        D(std::cout <<"["<<O<<"]"<< "before up: " << std::bitset<32>(path) << std::endl;)
        // Walk upwards until a left child is found. (might be the current one)
        const auto up_steps = dtl::bits::tz_count(~path);
        path >>= up_steps;
        path_level -= up_steps;
        // Go to right sibling.
        D(std::cout <<"["<<O<<"]"<< "after up:  " << std::bitset<32>(path) << std::endl;)
        path = path | 1;
        D(std::cout <<"["<<O<<"]"<< "right sib: " << std::bitset<32>(path) << std::endl;)


        // Advance the scanners and update the alpha vector.
        D(std::cout <<"["<<O<<"]"<< "alpha bef: " << std::bitset<32>(alpha) << std::endl;)

        const auto advance_begin = path_level;

        for (std::size_t i = 0; i < advance_begin; ++i) {
          D(std::cout <<"["<<O<<"]"<< "scanner " << i << ": " << scanners_[i] << std::endl;)
        }
        D(std::cout <<"["<<O<<"]"<< "advancing: [" << advance_begin << ", " << advance_end << ")" << std::endl;)
        for (std::size_t i = advance_begin; i < advance_end; ++i) {
          advance_scanner(i);
        }
        for (std::size_t i = advance_end; i < teb_.encoded_tree_height_; ++i) {
          D(std::cout <<"["<<O<<"]"<< "scanner " << i << ": " << scanners_[i] << std::endl;)
        }

        D(std::cout <<"["<<O<<"]"<< "alpha aft: " << std::bitset<32>(alpha) << std::endl;)

        // Walk downwards to the left-most leaf in that sub-tree.
        const auto down_steps = dtl::bits::tz_count(~(alpha >> path_level));
        path_level += down_steps;
        path <<= down_steps;
        D(std::cout <<"["<<O<<"]"<< "after down:" << std::bitset<32>(path) << std::endl;)
        D(std::cout <<"["<<O<<"]"<< "level:     " << path_level << std::endl;)


        // Toggle sentinel bit (= highest bit set) and add offset.
        pos = (path ^ (1ull << path_level)) << (h - path_level);
        // The length of the 1-fill.
        length = n >> path_level;
        D(std::cout <<"["<<O<<"]"<< "pos=" << pos << ", len=" << length << ", level=" << path_level << std::endl;)

        u64 label = scanners_[path_level].get_label(*this);

        if (unlikely(label)) {
          // Produce output (a 1-fill).
          results_[result_cnt].pos = pos;
          results_[result_cnt].length = length;
          ++result_cnt;
        }
      }

      if (result_cnt < batch_size_
          && pos == n) {
        results_[result_cnt].pos = n;
        results_[result_cnt].length = 0;
        ++result_cnt;
      }
      scan_pos_ = pos;
      scan_length_ = length;
      scan_path_ = path;
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
