#pragma once
//===----------------------------------------------------------------------===//
#include "teb.hpp"
#include "teb_flat.hpp"
#include "teb_iter.hpp"
#include "teb_scan_util.hpp"
#include "util/bit_buffer.hpp"
#include "util/bit_buffer_avx2.hpp"
#include "util/bit_buffer_avx512.hpp"

#include <dtl/dtl.hpp>
#include <dtl/iterator.hpp>

#include <array>
#include <cassert>
#include <iostream>
#include <memory>
#include <ostream>
#include <string>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// 1-fill iterator, WITHOUT efficient skip support. The scan iterator falls
/// back to a regular iterator if the level offsets are not present.
class teb_scan_iter {
  static constexpr u32 optimization_level_ = 3; // TODO remove

  using word_type = teb_word_type;

  /// The fundamental type to encode paths within the tree.
  using path_t = $u64;

  static constexpr u64 DEFAULT_BATCH_SIZE = 128 + 1;

  struct range_t { // TODO move to teb_types.hpp
    $u64 pos;
    $u64 length;
  };

  struct scanner_state_t {
    $u64 node_idx_;
    $u64 label_idx_;

    void __teb_inline__
    init(u64 node_idx, u64 label_offset) {
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
  const teb_flat& teb_;
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

  $u64 scan_path_;
  $u64 scan_path_level_;

  $u32 alpha_;

  std::size_t first_1label_idx_ = 0;

  u1 fallback_to_default_iter;
  std::unique_ptr<teb_iter> default_iter;
  //===--------------------------------------------------------------------===//

  //===--------------------------------------------------------------------===//
  // Helper functions.
  //===--------------------------------------------------------------------===//
  auto __teb_inline__
  determine_level(u32 alpha) {
    // LSB is always set. Thus, lz_count(alpha) is never undefined.
    return sizeof(alpha) * 8 - dtl::bits::lz_count(alpha);
  };

  auto __teb_inline__
  is_left_child(path_t path) {
    return (path & path_t(1)) == 0;
  };
  //===--------------------------------------------------------------------===//

public:
  /// Constructs an iterator for the given TEB instance. After construction,
  /// the iterator points to the first 1-fill.
  explicit __teb_inline__
  teb_scan_iter(const teb_flat& teb, u64 batch_size = DEFAULT_BATCH_SIZE) noexcept
      : teb_(teb),
        tree_height_(dtl::teb<>::determine_tree_height(teb.n_)),
        perfect_levels_(teb.perfect_level_cnt_),
        top_node_idx_begin_((1ull << (teb.perfect_level_cnt_ - 1)) - 1),
        top_node_idx_end_((1ull << teb.perfect_level_cnt_) - 1),
        top_node_idx_current_((1ull << (teb.perfect_level_cnt_ - 1)) - 1),
        results_(batch_size),
        result_cnt_(0),
        result_read_pos_(0),
        batch_size_(batch_size),
        alpha_(0),
        fallback_to_default_iter(!teb.has_level_offsets()),
        default_iter(nullptr) {
    if (fallback_to_default_iter) {
      default_iter = std::make_unique<teb_iter>(teb);
      results_[0].pos = default_iter->pos();
      results_[0].length = default_iter->length();
      result_cnt_ = 1;
      return;
    }

    // Hack for very sparse/unclustered bitmaps.
    first_1label_idx_ = teb_.L_.find_first(); // TODO always 0?

    results_.reserve(batch_size);
    if (teb_.encoded_tree_height_ == 1) {
      if (teb_.L_[0]) {
        results_[0].pos = 0;
        results_[0].length = teb_.n_;
        results_[1].pos = teb_.n_;
        results_[1].length = 0;
        result_cnt_ = 2;
      }
      else {
        results_[0].pos = teb_.n_;
        results_[0].length = 0;
        result_cnt_ = 1;
      }
      return;
    }
    // Initialize the scanners.
    for (std::size_t level = 0; level < teb_.encoded_tree_height_; ++level) {
      const auto node_offset = teb_.get_level_offset_tree(level);
      const auto label_offset = teb_.get_level_offset_labels(level);
      scanner_states_[level].init(node_offset, label_offset);
    }
    // Initialize path and level variable.
    scan_path_ = path_t(0);
    $u32 level = 0;
    for (; level < teb_.encoded_tree_height_; ++level) {
      u64 node_idx = teb_.get_level_offset_tree(level);
      const auto node_bit = teb_.is_inner_node(node_idx);
      if (!node_bit) break;
    }
    scan_path_level_ = level;
    // Iterator is now positioned at the first leaf node.
    get_next_batch();
  }

  teb_scan_iter(teb_scan_iter&&) noexcept = default;
  teb_scan_iter(const teb_scan_iter& other) = default;
  teb_scan_iter& operator=(const teb_scan_iter& other) = default;
  teb_scan_iter& operator=(teb_scan_iter&& other) = default;
  ~teb_scan_iter() = default;

  /// Forwards the iterator to the next 1-fill (if any).
  /// Use the functions pos() and length() to get the 1-fill the iterator
  /// is currently pointing to.
  void __teb_inline__
  next() noexcept __attribute__((flatten, hot)) {
    assert(!end());
    ++result_read_pos_;
    if (result_read_pos_ == result_cnt_) {
      result_read_pos_ = 0;
      result_cnt_ = 0;
      // Get the next batch of results.
      get_next_batch();
      //compact_results(); // disabled, because it lowers performance
    }
  }

  /// Merge contiguous 1-runs.
  void __forceinline__
  compact_results() {
    if (result_cnt_ == 0) return;
    u1 is_last_batch = results_[result_cnt_ - 1].length == 0;
    if (is_last_batch) --result_cnt_;
    if (result_cnt_ > 0) {
      std::size_t write_pos = 0;
      for (std::size_t i = 1; i < result_cnt_; ++i) {
        const auto e = results_[write_pos].pos + results_[write_pos].length;
        if (results_[i].pos == e) {
          results_[write_pos].length += results_[i].length;
        }
        else {
          ++write_pos;
          results_[write_pos] = results_[i];
        }
      }
      result_cnt_ = write_pos + 1;
    }
    if (is_last_batch) {
      results_[result_cnt_].pos = teb_.size();
      results_[result_cnt_].length = 0;
      ++result_cnt_;
    }
  }

  $u1 first_ = true;
  void
  get_next_batch() {
    if (fallback_to_default_iter) {
      // Fall back to the default iterator, as the given TEB instance does not
      // provide the necessary meta data to perform a tree scan.
      next_batch_from_default_iter();
      return;
    }
    const auto tree_levels = tree_height_ + 1;
    // Check, if a special implementation can be used to produce the next batch.
    if (perfect_levels_ == tree_levels) {
      // The bitmap is not compressed.
      next_batch_uncompressed();
    }
    else if (perfect_levels_ == tree_levels - 1) {
      next_batch_2levels();
    }
    else if (perfect_levels_ == tree_levels - 2) {
      next_batch_3levels();
    }
    else {
      // clang-format off
      // Use the common tree-scan algorithm.  The actual implementation is
      // chosen at compile time based on the target architecture.
      #ifdef __AVX512BW__
      next_batch_avx512(); // A SIMD implementation for AVX-512 architectures.
      #else
//      #ifdef __AVX2__
//      next_batch_avx2(); // A SIMD implementation for AVX2 architectures.
//      #else
//      next_batch_swar(); // A SIMD-within-a-register implementation.
//      if (first_) {first_ = false; std::cout << "STP" << std::endl;}
      next_batch_scalar_stack_perfect();
//      #endif
      #endif
    }

//    if (first_) {first_ = false; std::cout << "S" << std::endl;}
//    next_batch_scalar();

//    if (first_) {first_ = false; std::cout << "ST" << std::endl;}
//    next_batch_scalar_stack();

//    if (first_) {first_ = false; std::cout << "STP" << std::endl;}
//  -->  next_batch_scalar_stack_perfect();

//    if (first_) {first_ = false; std::cout << "SWAR" << std::endl;}
//    next_batch_swar();

//    #ifdef __AVX512BW__
//    next_batch_avx512(); // A SIMD implementation for AVX-512 architectures.
//    #else
//    #ifdef __AVX2__
//    next_batch_avx2(); // A SIMD implementation for AVX2 architectures.
//    #else
//    next_batch_swar(); // A SIMD-within-a-register implementation.
//    #endif
//    #endif
    // clang-format on
  }

  void next_batch_from_default_iter() noexcept __attribute__((flatten, hot, noinline)) {
    auto* it = default_iter.get();
    assert(!it->end());
    const auto n = teb_.size();
    std::size_t pos = 0;
    while (result_cnt_ < batch_size_ - 1) {
      it->next();
      if (it->end()) {
        results_[result_cnt_].pos = n;
        results_[result_cnt_].length = 0;
        ++result_cnt_;
        break;
      }
      else {
        pos = it->pos();
        results_[result_cnt_].pos = pos;
        results_[result_cnt_].length = it->length();
        ++result_cnt_;
      }
    }
  }

  /// Handle the special case where the labels are equal to the uncompressed
  /// bitmap.
  void
  next_batch_uncompressed() noexcept __attribute__((flatten, hot, noinline)) {
    const auto h = tree_height_ + 1;
    const auto n = teb_.size();
    // Note: The path variable does NOT contain a sentinel bit. Instead, we
    //       keep track of the paths' level using a separate variable.
    auto pos = scan_path_; // path << (h - path_level);
    const auto length = 1; // constant in that case
    auto result_cnt = result_cnt_;

    const data_view<const word_type> L {
      // TODO use the bitmap view instead
      teb_.L_.data_.begin(),
      teb_.L_.data_.end(),
    };

    i64 delta_l = 0 - static_cast<i64>(teb_.implicit_leading_label_cnt_);
    while (result_cnt < batch_size_ - 64) {
      // Initialize the bit buffer.
      u64 bb = fetch_n_bits(L, static_cast<i64>(pos) + delta_l, 64, false, false);

      // Do multiple iterations to consume the bit buffer.
      for ($i32 b = 0; b < 63; ++b) { // TODO 64?
        const auto label = (bb >> b) & 1;
        // Produce output (a 1-fill).
        results_[result_cnt].pos = pos;
        results_[result_cnt].length = length;
        result_cnt += label;

        ++pos;
        if (pos == n) goto done;
      }
    }

  done:
    if (result_cnt < batch_size_
        && pos == n) {
      results_[result_cnt].pos = n;
      results_[result_cnt].length = 0;
      ++result_cnt;
    }
    scan_path_ = pos;
    result_cnt_ = result_cnt;
  }

  void
  next_batch_2levels() noexcept __attribute__((flatten, hot, noinline)) {
    const auto h = tree_height_ + 1;
    const auto n = teb_.size();
    const auto u = perfect_levels_ - 1;
    // Note: The path variable does NOT contain a sentinel bit. Instead, we
    //       keep track of the paths' level using a separate variable.
    register auto path = scan_path_;
    register auto path_level = u;
    auto pos = scan_path_;
    auto length = n >> path_level;
    auto result_cnt = result_cnt_;

    const data_view<const word_type> T {
      // TODO use the bitmap view instead
      teb_.T_.data_.begin(),
      teb_.T_.data_.end(),
    };
    const data_view<const word_type> L {
      // TODO use the bitmap view instead
      teb_.L_.data_.begin(),
      teb_.L_.data_.end(),
    };

    const std::size_t level = h - 1;
    while (result_cnt < batch_size_ - 64) {
      // Initialize the bit buffers. - Bit buffers try to reduce the amount of
      // load instructions from T and L.

      // The last tree level is h - 1, but on that level, all tree nodes are
      // leaves (0-bit). Thus we only need to read the tree structure bits
      // at level h - 2.
      i64 delta_t = 0 - static_cast<i64>(teb_.implicit_inner_node_cnt_);
      u64 bb_t = fetch_n_bits(T, static_cast<i64>(scanner_states_[h - 2].node_idx_) + delta_t, 64);

      // The above does not apply for the labels. Here we need to read the
      // label bits at level h - 2 and h - 1.
      i64 delta_l = 0 - static_cast<i64>(teb_.implicit_leading_label_cnt_);
      u64 bb_l0 = fetch_n_bits(L, static_cast<i64>(scanner_states_[h - 2].label_idx_) + delta_l, 64, false, false);
      u64 bb_l1 = fetch_n_bits(L, static_cast<i64>(scanner_states_[h - 1].label_idx_) + delta_l, 64, false, false);

      // Do multiple iterations to consume the bit buffer.
      $i32 t_idx = 0;
      $i32 l0_idx = 0;
      $i32 l1_idx = 0;
      for (; t_idx < 31; ++t_idx) { // TODO 32?
        u1 t = ((bb_t >> t_idx) & 1) == 1;

        if (t) {
          // Observed an inner node at level h - 2.
          // Consume two labels from level h - 1.
          const auto label0 = (bb_l1 >> l1_idx) & 1;
          const auto label1 = (bb_l1 >> (l1_idx + 1)) & 1;
          l1_idx += 2;

          // Produce output (a 1-fill).
          results_[result_cnt].pos = pos + (!label0);
          const auto len = 0 + label0 + label1;
          results_[result_cnt].length = len;
          result_cnt += (len > 0);
        }
        else {
          // Observed a leaf node at level h - 2.
          // Consume one label from level h - 2.
          const auto label = (bb_l0 >> l0_idx) & 1;
          l0_idx += 1;
          // Produce output (a 1-fill).
          results_[result_cnt].pos = pos;
          results_[result_cnt].length = 2;
          result_cnt += label;
        }

        pos += 2;
        if (pos == n)
          goto done;
        assert(pos < n);
      }

      // Update the scanner states.
      scanner_states_[h - 2].node_idx_ += t_idx;
      scanner_states_[h - 2].label_idx_ += l0_idx;
      scanner_states_[h - 1].label_idx_ += l1_idx;
    }

  done:
    if (result_cnt < batch_size_
        && pos == n) {
      results_[result_cnt].pos = n;
      results_[result_cnt].length = 0;
      ++result_cnt;
    }
    scan_path_ = pos;
    result_cnt_ = result_cnt;
  }

  void
  next_batch_3levels() noexcept __attribute__((flatten, hot, noinline)) {
    const auto h = tree_height_ + 1;
    const auto n = teb_.size();
    const auto u = perfect_levels_ - 1;
    // Note: The path variable does NOT contain a sentinel bit. Instead, we
    //       keep track of the paths' level using a separate variable.
    register auto path = scan_path_;
    register auto path_level = u;
    auto pos = scan_path_;
    auto length = n >> path_level;
    auto result_cnt = result_cnt_;

    const data_view<const word_type> T {
      // TODO use the bitmap view instead
      teb_.T_.data_.begin(),
      teb_.T_.data_.end(),
    };
    const data_view<const word_type> L {
      // TODO use the bitmap view instead
      teb_.L_.data_.begin(),
      teb_.L_.data_.end(),
    };

    const std::size_t level = h - 1;
    while (result_cnt < batch_size_ - 1 - 64) {
      // Initialize the bit buffers. - Bit buffers try to reduce the amount of
      // load instructions from T and L.

      i64 delta = 0 - static_cast<i64>(teb_.implicit_inner_node_cnt_);
      u64 bb_t0 = fetch_n_bits(T, static_cast<i64>(scanner_states_[h - 3].node_idx_) + delta, 16);
      u64 bb_t1 = fetch_n_bits(T, static_cast<i64>(scanner_states_[h - 2].node_idx_) + delta, 32);

      i64 delta_l = 0 - static_cast<i64>(teb_.implicit_leading_label_cnt_);
      u64 bb_l0 = fetch_n_bits(L, static_cast<i64>(scanner_states_[h - 3].label_idx_) + delta_l, 16, false, false);
      u64 bb_l1 = fetch_n_bits(L, static_cast<i64>(scanner_states_[h - 2].label_idx_) + delta_l, 32, false, false);
      u64 bb_l2 = fetch_n_bits(L, static_cast<i64>(scanner_states_[h - 1].label_idx_) + delta_l, 64, false, false);

      // Do multiple iterations to consume the bit buffer.
      $i32 t0_idx = 0;
      $i32 t1_idx = 0;
      $i32 l0_idx = 0;
      $i32 l1_idx = 0;
      $i32 l2_idx = 0;
      for (; t0_idx < 15; ++t0_idx) { // TODO 32?
        u1 t0 = ((bb_t0 >> t0_idx) & 1) == 1;

        if (t0) {
          // Observed an inner node ...
          // Go to the two child nodes.
          for (std::size_t j = 0; j < 2; ++j) {
            u1 t1 = ((bb_t1 >> t1_idx) & 1) == 1;
            t1_idx += 1;
            if (t1) {
              // Observed an inner node ...
              // Consume two labels from level ....
              const auto label0 = (bb_l2 >> l2_idx) & 1;
              const auto label1 = (bb_l2 >> (l2_idx + 1)) & 1;
              l2_idx += 2;
              // Produce output (a 1-fill).
              results_[result_cnt].pos = pos + (!label0);
              const auto len = 0 + label0 + label1;
              results_[result_cnt].length = len;
              result_cnt += len;
            }
            else {
              // Observed a leaf node at level ....
              // Consume one label from level ....
              const auto label = (bb_l1 >> l1_idx) & 1;
              l1_idx += 1;
              // Produce output (a 1-fill).
              results_[result_cnt].pos = pos;
              results_[result_cnt].length = 2;
              result_cnt += label;
            }
            pos += 2;
          }
        }
        else {
          // Observed a leaf node at level ....
          // Consume one label from level ....
          const auto label = (bb_l0 >> l0_idx) & 1;
          l0_idx += 1;
          // Produce output (a 1-fill).
          results_[result_cnt].pos = pos;
          results_[result_cnt].length = 4;
          result_cnt += label;
          pos += 4;
        }

        if (pos == n)
          goto done;
        assert(pos < n);
      }

      // Update the scanner states.
      scanner_states_[h - 3].node_idx_ += t0_idx;
      scanner_states_[h - 2].node_idx_ += t1_idx;
      scanner_states_[h - 3].label_idx_ += l0_idx;
      scanner_states_[h - 2].label_idx_ += l1_idx;
      scanner_states_[h - 1].label_idx_ += l2_idx;
    }

  done:
    if (result_cnt < batch_size_
        && pos == n) {
      results_[result_cnt].pos = n;
      results_[result_cnt].length = 0;
      ++result_cnt;
    }
    scan_path_ = pos;
    result_cnt_ = result_cnt;
  }

  void
  next_batch_swar() noexcept __attribute__((flatten, hot, noinline)) {
    const auto h = tree_height_;
    const auto n = teb_.size();
    // Note: The path variable does NOT contain a sentinel bit. Instead, we
    //       keep track of the paths' level using a separate variable.
    register auto path = scan_path_;
    register auto path_level = scan_path_level_;
    auto pos = path << (h - path_level);
    auto length = n >> path_level;
    auto result_cnt = result_cnt_;
    auto alpha = 0u;

    const data_view<const word_type> T {
      // TODO use the bitmap view instead
      teb_.T_.data_.begin(),
      teb_.T_.data_.end(),
    };
    const data_view<const word_type> L {
      // TODO use the bitmap view instead
      teb_.L_.data_.begin(),
      teb_.L_.data_.end(),
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
        case 0:
          t_bb0.set(slot_idx, ~0ul);
          l_bb0.set(slot_idx, ~0ul);
          break;
        case 1:
          t_bb1.set(slot_idx, ~0ul);
          l_bb1.set(slot_idx, ~0ul);
          break;
        case 2:
          t_bb2.set(slot_idx, ~0ul);
          l_bb2.set(slot_idx, ~0ul);
          break;
        case 3:
          t_bb3.set(slot_idx, ~0ul);
          l_bb3.set(slot_idx, ~0ul);
          break;
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
      switch (level_begin) {
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
          t_bb0.increment((m >> (8 * 0)) & p);
        case 8:
        case 9:
        case 10:
        case 11:
        case 12:
        case 13:
        case 14:
        case 15:
          t_bb1.increment((m >> (8 * 1)) & p);
        case 16:
        case 17:
        case 18:
        case 19:
        case 20:
        case 21:
        case 22:
        case 23:
          t_bb2.increment((m >> (8 * 2)) & p);
        default:
          t_bb3.increment((m >> (8 * 3)) & p);
      }
      // Update the alpha vector.
      u64 r0 = t_bb0.read();
      u64 r1 = t_bb1.read();
      u64 r2 = t_bb2.read();
      u64 r3 = t_bb3.read();
      alpha = (r3 << (8 * 3)) | (r2 << (8 * 2)) | (r1 << (8 * 1)) | r0;
    };

    // Advance the label bit buffers.
    auto advance_label_scanner = [&](u32 level) {
      u64 buffer_idx = level / 8;
      u64 slot_idx = level % 8;
      u64 slot_mask = 1ul << slot_idx;
      $u64 label = 0;
      switch (buffer_idx) {
        case 0: l_bb0.increment(slot_mask); break;
        case 1: l_bb1.increment(slot_mask); break;
        case 2: l_bb2.increment(slot_mask); break;
        case 3: l_bb3.increment(slot_mask); break;
      }
    };

    while (result_cnt < batch_size_ - 8) {
      // Initialize the bit buffers.
      t_bb0.reset_read_mask();
      t_bb1.reset_read_mask();
      t_bb2.reset_read_mask();
      t_bb3.reset_read_mask();
      l_bb0.reset_read_mask();
      l_bb1.reset_read_mask();
      l_bb2.reset_read_mask();
      l_bb3.reset_read_mask();
      for (std::size_t level = perfect_levels_ - 1;
           level < teb_.encoded_tree_height_;
           ++level) {
        u64 buffer_idx = level / 8;
        u64 slot_idx = level % 8;
        i64 delta = -static_cast<i64>(teb_.implicit_inner_node_cnt_);
        switch (buffer_idx) {
          case 0:
            t_bb0.set(slot_idx, fetch_n_bits(T, static_cast<i64>(scanner_states_[level].node_idx_) + delta, 8));
            l_bb0.set(slot_idx, fetch_n_bits(L, static_cast<i64>(scanner_states_[level].label_idx_), 8));
            break;
          case 1:
            t_bb1.set(slot_idx, fetch_n_bits(T, static_cast<i64>(scanner_states_[level].node_idx_) + delta, 8));
            l_bb1.set(slot_idx, fetch_n_bits(L, static_cast<i64>(scanner_states_[level].label_idx_), 8));
            break;
          case 2:
            t_bb2.set(slot_idx, fetch_n_bits(T, static_cast<i64>(scanner_states_[level].node_idx_) + delta, 8));
            l_bb2.set(slot_idx, fetch_n_bits(L, static_cast<i64>(scanner_states_[level].label_idx_), 8));
            break;
          case 3:
            t_bb3.set(slot_idx, fetch_n_bits(T, static_cast<i64>(scanner_states_[level].node_idx_) + delta, 8));
            l_bb3.set(slot_idx, fetch_n_bits(L, static_cast<i64>(scanner_states_[level].label_idx_), 8));
            break;
        }
      }

      // Initialize alpha. (which reads the first bits in the bit buffers)
      {
        u64 r0 = t_bb0.read();
        u64 r1 = t_bb1.read();
        u64 r2 = t_bb2.read();
        u64 r3 = t_bb3.read();
        alpha = (r3 << (8 * 3)) | (r2 << (8 * 2)) | (r1 << 8) | r0;
      }

      // Do multiple iteration to consume the bit buffers.
      for ($i32 b = 0; b < 7; ++b) {
        u64 length = n >> path_level;

        assert(pos < n);
        assert(pos + length <= n);
        assert(path_level >= 1);
        assert(length > 0);

        // Read the label of the current leaf node.
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
        results_[result_cnt].pos = pos;
        results_[result_cnt].length = length;
        result_cnt += label;

        // Increment the current position.
        pos += length;
        if (pos == n) goto done; // TODO eliminate branch

        // Walk upwards until a left child is found. (might be the current one)
        const auto advance_end = path_level + 1;
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
           level < teb_.encoded_tree_height_;
           ++level) {
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
    scan_path_ = path;
    scan_path_level_ = path_level;
    result_cnt_ = result_cnt;
    alpha_ = alpha;
  }

#ifdef __AVX2__
  void
  next_batch_avx2() noexcept __attribute__((flatten, hot, noinline)) {
    const auto h = tree_height_;
    const auto n = teb_.size();
    // Note: The path variable does NOT contain a sentinel bit. Instead, we
    //       keep track of the paths' level using a separate variable.
    register auto path = scan_path_;
    register auto path_level = scan_path_level_;
    auto pos = path << (h - path_level);
    auto length = n >> path_level;
    auto result_cnt = result_cnt_;
    auto alpha = 0u;
    auto alpha_la = 0u;
    auto beta = 0u;

    dtl::r256 tmp_t {
      .u8 = {
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0 }
    };
    dtl::r256 tmp_l {
      .u8 = {
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0,
          0, 0, 0, 0 }
    };

    const data_view<const word_type> T {
      // TODO use the bitmap view instead
      teb_.T_.data_.begin(),
      teb_.T_.data_.end(),
    };
    const data_view<const word_type> L {
      // TODO use the bitmap view instead
      teb_.L_.data_.begin(),
      teb_.L_.data_.end(),
    };

    // Bit buffers for the tree structure.
    dtl::bit_buffer_avx2<> tree_bit_buffer;

    // Bit buffers for the labels.
    dtl::bit_buffer_avx2<> label_bit_buffer;

    // Perfect levels. (all bits set to 1)
    for (std::size_t level = 0; level < perfect_levels_ - 1; ++level) {
      const auto slot_idx = level;
      tmp_t.u8[slot_idx] = ~u8(0);
      tmp_l.u8[slot_idx] = ~u8(0);
    }
    tree_bit_buffer.set_raw(tmp_t.i);
    label_bit_buffer.set_raw(tmp_l.i);

    assert(batch_size_ > 8);
    while (result_cnt < batch_size_ - 8
        && pos < n) {
      // Initialize the bit buffers.
      tree_bit_buffer.reset_read_mask();
      label_bit_buffer.reset_read_mask();
      for (std::size_t level = perfect_levels_ - 1;
           level < teb_.encoded_tree_height_;
           ++level) {
        const auto slot_idx = level;
        i64 delta = -static_cast<i64>(teb_.implicit_inner_node_cnt_);
        tmp_t.u8[slot_idx] = fetch_n_bits(T, static_cast<i64>(scanner_states_[level].node_idx_) + delta, 8);
        i64 delta_l = 0 - static_cast<i64>(teb_.implicit_leading_label_cnt_);
        tmp_l.u8[slot_idx] = fetch_n_bits(L, static_cast<i64>(scanner_states_[level].label_idx_) + delta_l, 8, false, false);
      }
      tree_bit_buffer.set_raw(tmp_t.i);
      label_bit_buffer.set_raw(tmp_l.i);

      // Initialize alpha. (which reads the first bits in the bit buffers)
      alpha = tree_bit_buffer.read();
      alpha_la = tree_bit_buffer.read_ahead();

      // Initialize beta.
      beta = label_bit_buffer.read();

      // Do multiple iteration to consume the bit buffers.
      $u1 is_not_done = true;
      for ($i32 b = 0; b < 6 && is_not_done; ++b) {
        u64 length = n >> path_level;

        assert(pos < n);
        assert(pos + length <= n);
        assert(path_level >= 1);
        assert(length > 0);

        // Produce output (a 1-fill).
        results_[result_cnt].pos = pos;
        results_[result_cnt].length = length;
        // Read the label of the current leaf node.
        u64 label = (beta >> path_level) & 1;
        result_cnt += label;

        // Increment the current position.
        pos += length;
        is_not_done = pos != n;

        // Advance the label scanner.
        label_bit_buffer.increment(u32(1) << path_level);
        beta = label_bit_buffer.read();

        // Walk upwards until a left child is found. (might be the current one)
        const auto advance_end = path_level + 1;
        u32 advance_mask_hi = (~u32(0)) >> (32 - advance_end);
        const auto up_steps = dtl::bits::tz_count(~path);
        path >>= up_steps;
        path_level -= up_steps;
        // Go to right sibling.
        path = path | 1;

        // Advance the scanners and update the alpha vector.
        const auto advance_begin = path_level;
        u32 advance_mask_lo = (~u32(0)) << advance_begin;
        u32 advance_mask = advance_mask_lo & advance_mask_hi;

        const auto alpha_next = (alpha_la & advance_mask) | (alpha & (~advance_mask));
        tree_bit_buffer.increment(advance_mask);

        // Walk downwards to the left-most leaf in that sub-tree.
        const auto down_steps = dtl::bits::tz_count(~(alpha_next >> path_level));
        path_level += down_steps;
        path <<= down_steps;

        // Update the alpha vector.
        alpha = alpha_next;
        alpha_la = tree_bit_buffer.read_ahead();
      }

      // Update the scanner states.
      dtl::r256 tmp_tm { .i = tree_bit_buffer.get_read_mask() };
      dtl::r256 tmp_lm { .i = label_bit_buffer.get_read_mask() };
      for ($u32 level = perfect_levels_ - 1;
           level < teb_.encoded_tree_height_;
           ++level) {
        scanner_states_[level].node_idx_ += dtl::bits::tz_count(tmp_tm.u8[level]);
        scanner_states_[level].label_idx_ += dtl::bits::tz_count(tmp_lm.u8[level]);
      }
    }

  done:
    if (result_cnt < batch_size_
        && pos == n) {
      results_[result_cnt].pos = n;
      results_[result_cnt].length = 0;
      ++result_cnt;
    }
    scan_path_ = path;
    scan_path_level_ = path_level;
    result_cnt_ = result_cnt;
    alpha_ = alpha;
  }
#endif // __AVX2__

#ifdef __AVX512BW__
  void //__teb_inline__
  next_batch_avx512() noexcept __attribute__((flatten, hot, noinline)) {
    const auto h = tree_height_;
    const auto eth = teb_.encoded_tree_height_;
    const auto u = perfect_levels_ - 1;
    const auto n = teb_.size();
    // Note: The path variable does NOT contain a sentinel bit. Instead, we
    //       keep track of the paths' level using a separate variable.
    register auto path = scan_path_;
    register auto path_level = scan_path_level_;
    auto pos = path << (h - path_level);
    auto length = n >> path_level;
    auto result_cnt = result_cnt_;
    auto alpha = 0u;
    auto alpha_la = 0u;
    auto beta = 0u;

    dtl::r512 tmp_t { .u16 = {
                          0, 0, 0, 0,
                          0, 0, 0, 0,
                          0, 0, 0, 0,
                          0, 0, 0, 0,
                          0, 0, 0, 0,
                          0, 0, 0, 0,
                          0, 0, 0, 0,
                          0, 0, 0, 0 } };
    dtl::r512 tmp_l { .u16 = {
                          0, 0, 0, 0,
                          0, 0, 0, 0,
                          0, 0, 0, 0,
                          0, 0, 0, 0,
                          0, 0, 0, 0,
                          0, 0, 0, 0,
                          0, 0, 0, 0,
                          0, 0, 0, 0 } };

    const data_view<const word_type> T {
      // TODO use the bitmap view instead
      teb_.T_.data_.begin(),
      teb_.T_.data_.end(),
    };
    const data_view<const word_type> L {
      // TODO use the bitmap view instead
      teb_.L_.data_.begin(),
      teb_.L_.data_.end(),
    };

    // Bit buffers for the tree structure.
    dtl::bit_buffer_avx512<> tree_bit_buffer;

    // Bit buffers for the labels.
    dtl::bit_buffer_avx512<> label_bit_buffer;

    // Perfect levels. (all bits set to 1)
    for (std::size_t level = 0; level < u; ++level) {
      const auto slot_idx = level;
      tmp_t.u16[slot_idx] = ~u16(0);
      tmp_l.u16[slot_idx] = ~u16(0);
    }
    tree_bit_buffer.set_raw(tmp_t.i);
    label_bit_buffer.set_raw(tmp_l.i);

    const auto batch_size = DEFAULT_BATCH_SIZE;
    while (result_cnt < batch_size - 16 && pos < n) {
      // Initialize the bit buffers.
      tree_bit_buffer.reset_read_mask();
      label_bit_buffer.reset_read_mask();
      for (std::size_t level = u; level < eth; ++level) {
        const auto slot_idx = level;
        i64 delta = 0 - static_cast<i64>(teb_.implicit_inner_node_cnt_);
        tmp_t.u16[slot_idx] = fetch_n_bits(T, static_cast<i64>(scanner_states_[level].node_idx_) + delta, 16);
        i64 delta_l = 0 - static_cast<i64>(teb_.implicit_leading_label_cnt_);
        tmp_l.u16[slot_idx] = fetch_n_bits(L, static_cast<i64>(scanner_states_[level].label_idx_) + delta_l, 16, false, false);
      }
      tree_bit_buffer.set_raw(tmp_t.i);
      label_bit_buffer.set_raw(tmp_l.i);

      // Initialize alpha. (which reads the first bits in the bit buffers)
      alpha = tree_bit_buffer.read();
      alpha_la = tree_bit_buffer.read_ahead();

      // Initialize beta.
      beta = label_bit_buffer.read();

      // Do multiple iteration to consume the bit buffers.
      $u1 is_not_done = true;
      for ($i32 b = 0; b < 15 && is_not_done; ++b) {
        u64 length = n >> path_level;

        assert(pos < n);
        assert(pos + length <= n);
        assert(path_level >= 1);
        assert(length > 0);

        // Produce output (a 1-fill).
        results_[result_cnt].pos = pos;
        results_[result_cnt].length = length;
        // Read the label of the current leaf node.
        u64 label = (beta >> path_level) & 1;
        result_cnt += label;

        // Increment the current position.
        pos += length;
        is_not_done = pos != n;

        const auto up_steps = dtl::bits::tz_count(~path);

        // Advance the label scanner.
        label_bit_buffer.increment(__mmask32(1) << path_level);
        beta = label_bit_buffer.read();

        // Walk upwards until a left child is found. (might be the current one)
        const auto advance_end = path_level + 1;
        const __mmask32 advance_mask_hi = (~__mmask32(0)) >> (32 - advance_end);
        path >>= up_steps;
        path_level -= up_steps;
        // Go to right sibling.
        path = path | 1;

        // Advance the scanners and update the alpha vector.
        const auto advance_begin = path_level;
        const __mmask32 advance_mask_lo = (~__mmask32(0)) << advance_begin;
        const __mmask32 advance_mask = advance_mask_lo & advance_mask_hi;

        const auto alpha_next = (alpha_la & advance_mask) | (alpha & (~advance_mask));
        tree_bit_buffer.increment(advance_mask);

        // Walk downwards to the left-most leaf in that sub-tree.
        const auto down_steps = dtl::bits::tz_count(~(alpha_next >> path_level));
        path_level += down_steps;
        path <<= down_steps;

        // Update the alpha vector.
        alpha = alpha_next;
        alpha_la = tree_bit_buffer.read_ahead();
      }

      // Update the scanner states.
      dtl::r512 tmp_tm { .i = tree_bit_buffer.get_read_mask() };
      dtl::r512 tmp_lm { .i = label_bit_buffer.get_read_mask() };
      for ($u32 level = u; level < eth; ++level) {
        scanner_states_[level].node_idx_ += dtl::bits::tz_count(tmp_tm.u16[level]);
        scanner_states_[level].label_idx_ += dtl::bits::tz_count(tmp_lm.u16[level]);
      }
    }

  done:
    if (result_cnt < batch_size && pos == n) {
      results_[result_cnt].pos = n;
      results_[result_cnt].length = 0;
      ++result_cnt;
    }
    scan_path_ = path;
    scan_path_level_ = path_level;
    result_cnt_ = result_cnt;
    alpha_ = alpha;
  }
#endif // __AVX512BW__

  /// A simple scalar implementation of the tree scan (without bit buffers).
  void //__teb_inline__
  next_batch_scalar() noexcept __attribute__((flatten, hot, noinline)) {
    const auto h = tree_height_;
    const auto n = teb_.size();
    register auto path = scan_path_;
    register auto path_level = scan_path_level_;
    auto pos = path << (h - path_level);
    auto length = n >> path_level;
    auto result_cnt = result_cnt_;
    $u32 alpha = 0u;

    auto advance_scanner = [&](std::size_t i) {
      u1 was_inner_node = dtl::bits::bit_test(alpha, i);
      scanner_states_[i].node_idx_ += 1;
      scanner_states_[i].label_idx_ += !was_inner_node;
      u1 is_inner_node = teb_.is_inner_node(scanner_states_[i].node_idx_);
      u1 x = was_inner_node ^ is_inner_node;
      alpha ^= u32(x) << i;
    };

    auto advance_scanners = [&](std::size_t advance_begin,
                                std::size_t advance_end) {
      for (std::size_t i = advance_begin; i < advance_end; ++i) {
        advance_scanner(i);
      }
    };

    // Initialize alpha.
    alpha = perfect_levels_ == 1 ? 0u : ~0u >> (sizeof(u32) * 8 - (perfect_levels_ - 1));
    for (std::size_t level = perfect_levels_ - 1;
         level < teb_.encoded_tree_height_;
         ++level) {
      u1 is_inner = teb_.is_inner_node(scanner_states_[level].node_idx_);
      alpha |= u32(is_inner) << level;
    }

    while (result_cnt < batch_size_ && pos < n) {
      assert(pos <= n);
      assert(pos + length <= n);
      assert(path_level >= 1);
      assert(length > 0);

      // Toggle sentinel bit (= highest bit set) and add offset.
      pos = path << (h - path_level);
      // The length of the 1-fill.
      length = n >> path_level;

      const auto label_idx = scanner_states_[path_level].label_idx_;
      if (label_idx >= first_1label_idx_) {
        const auto label = teb_.get_label_by_idx(scanner_states_[path_level].label_idx_);

        // Produce output (a 1-fill).
        results_[result_cnt].pos = pos;
        results_[result_cnt].length = length;
        result_cnt += label;
      }

      // Increment the current position.
      pos += length;
      if (pos == n) break; // TODO eliminate branch
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
      advance_scanners(advance_begin, advance_end);

      // Walk downwards to the left-most leaf in that sub-tree.
      const auto down_steps = dtl::bits::tz_count(~(alpha >> path_level));
      path_level += down_steps;
      path <<= down_steps;
    }

    if (result_cnt < batch_size_
        && pos == n) {
      results_[result_cnt].pos = n;
      results_[result_cnt].length = 0;
      ++result_cnt;
    }
    scan_path_ = path;
    scan_path_level_ = path_level;
    result_cnt_ = result_cnt;
  }

  /// A simple scalar implementation of the tree scan (without bit buffers).
  $u1 first_time_stack_ = true;
  void //__teb_inline__
  next_batch_scalar_stack() noexcept __attribute__((flatten, hot, noinline)) {
    if (first_time_stack_) {
      first_time_stack_ = false;
      // Advance the scanners.
      for (std::size_t level = 0; level < scan_path_level_; ++level) {
        scanner_states_[level].node_idx_++;
      }
    }
    const auto h = tree_height_;
    const auto n = teb_.size();
    register auto path = scan_path_;
    register auto path_level = scan_path_level_;
    auto pos = path << (h - path_level);
    auto result_cnt = result_cnt_;

    while (true) {
      auto node_idx = scanner_states_[path_level].node_idx_;
      auto node_is_inner = teb_.is_inner_node(node_idx);

      while (node_is_inner) {
        // ~~~ Push right child on the stack and go to left child.
        scanner_states_[path_level].node_idx_++;
        path_level++;
        path <<= 1;
        node_idx = scanner_states_[path_level].node_idx_;
        node_is_inner = teb_.is_inner_node(node_idx);
      }
      // Reached a leaf node.
      const auto label_idx = scanner_states_[path_level].label_idx_;
      const auto label = teb_.get_label_by_idx(label_idx);
      if (label) { // To branch or not to branch, that's the question.
        // Produce output (a 1-fill).
        results_[result_cnt].pos = path << (tree_height_ - path_level);
        results_[result_cnt].length = teb_.n_ >> path_level;
        result_cnt++;
      }
      scanner_states_[path_level].node_idx_++;
      scanner_states_[path_level].label_idx_++;

      // Walk upwards until a left child is found. (might be the current one)
      const auto up_steps = dtl::bits::tz_count(~path);
      path >>= up_steps;
      path_level -= up_steps;
      // Go to right sibling.
      path |= 1;
      if (result_cnt == batch_size_ || path_level == 0) break;
    }

    if (result_cnt < batch_size_
        && path_level == 0) {
      results_[result_cnt].pos = n;
      results_[result_cnt].length = 0;
      ++result_cnt;
    }

    scan_path_ = path;
    scan_path_level_ = path_level;
    result_cnt_ = result_cnt;
  }

  /// A simple scalar implementation of the tree scan (without bit buffers).
  /// In contrast to next_batch_scalar_stack(), this function respects the
  /// perfect levels.

  void //__teb_inline__
  next_batch_scalar_stack_perfect() noexcept __attribute__((flatten, hot, noinline)) {
    if (first_time_stack_) {
      first_time_stack_ = false;
      // Advance the scanners.
      for (std::size_t level = 0; level < scan_path_level_; ++level) {
        scanner_states_[level].node_idx_++;
      }
    }

    const auto h = tree_height_;
    const auto n = teb_.size();
    register auto path = scan_path_;
    register auto path_level = scan_path_level_;
    auto pos = path << (h - path_level);
    auto result_cnt = result_cnt_;
    while (true) {
      auto node_idx = scanner_states_[path_level].node_idx_;
      auto node_is_inner = teb_.is_inner_node(node_idx);

      while (node_is_inner) {
        // ~~~ Push right child on the stack and go to left child.
        scanner_states_[path_level].node_idx_++;
        path_level++;
        path <<= 1;
        node_idx = scanner_states_[path_level].node_idx_;
        node_is_inner = teb_.is_inner_node(node_idx);
      }
      // Reached a leaf node.
      const auto label_idx = scanner_states_[path_level].label_idx_;
      const auto label = teb_.get_label_by_idx(label_idx);
      if (label) { // To branch or not to branch, that's the question.
        // Produce output (a 1-fill).
        results_[result_cnt].pos = path << (tree_height_ - path_level);
        results_[result_cnt].length = teb_.n_ >> path_level;
        result_cnt++;
      }
      scanner_states_[path_level].node_idx_++;
      scanner_states_[path_level].label_idx_++;

      // Walk upwards until a left child is found. (might be the current one)
      const auto up_steps = dtl::bits::tz_count(~path);
      path >>= up_steps;
      path_level -= up_steps;
      // Go to right sibling.
      path |= 1;
      if (path_level < perfect_levels_) {
        ++top_node_idx_current_;
        path = path_t(top_node_idx_current_ - top_node_idx_begin_);
        path_level = perfect_levels_ - 1;
        if (top_node_idx_current_ == top_node_idx_end_) {
          path_level = 0;
          break;
        }
      }
      // Reserve one entry in the current batch, in case we reach the end.
      if (result_cnt == batch_size_ - 1 || path_level == 0) break;
    }

    if (result_cnt < batch_size_
        && path_level == 0) {
      results_[result_cnt].pos = n;
      results_[result_cnt].length = 0;
      ++result_cnt;
    }

    scan_path_ = path;
    scan_path_level_ = path_level;
    result_cnt_ = result_cnt;
  }

  /// Forwards the iterator to the given position. - Note: This functions
  /// simply calls next() until the desired position has been reached.
  void // __teb_inline__ TODO revert
  skip_to(const std::size_t to_pos) noexcept __attribute__((noinline)) {
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
      results_[result_read_pos_].length -= to_pos - results_[result_read_pos_].pos;
      results_[result_read_pos_].pos = to_pos;
    }
  }

  /// Returns true if the iterator reached the end, false otherwise.
  u1 __forceinline__ //__teb_inline__
  end() const noexcept {
    return results_[result_read_pos_].pos == teb_.n_;
  }

  /// Returns the starting position of the current 1-fill.
  u64 __forceinline__ //__teb_inline__
  pos() const noexcept {
    assert(result_read_pos_ < result_cnt_);
    return results_[result_read_pos_].pos;
  }

  /// Returns the length of the current 1-fill.
  u64 __forceinline__ //__teb_inline__
  length() const noexcept {
    assert(result_read_pos_ < result_cnt_);
    return results_[result_read_pos_].length;
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
