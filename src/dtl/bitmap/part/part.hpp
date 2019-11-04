#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/bitmap/iterator.hpp>
#include <dtl/dtl.hpp>
#include <dtl/math.hpp>

#include <boost/dynamic_bitset.hpp>
#include <memory>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// Applies a fixed size partitioning to the given bitmap.
template<
    /// The (compressed) bitmap type.
    typename B,
    /// The partition size in bits.
    std::size_t P>
class part {
  static_assert(dtl::is_power_of_two(P),
      "The partition size must be a power of two.");

protected:
  /// Pointer to the partitions.
  std::vector<std::unique_ptr<B>> parts_;
  /// The (total) length of the bitmap.
  std::size_t n_;

public:
  using bitmap_type = B;
  static constexpr std::size_t part_bitlength = P;

  /// C'tor (similar to all other implementations)
  explicit part(const boost::dynamic_bitset<$u32>& bitmap)
      : parts_(), n_(bitmap.size()) {
    const std::size_t part_cnt =
        (bitmap.size() + (part_bitlength - 1)) / part_bitlength;

    for (std::size_t p = 0; p < part_cnt; ++p) {
      // TODO avoid copy, or at least make it faster
      boost::dynamic_bitset<$u32> b(part_bitlength);
      const std::size_t p_begin = p * part_bitlength;
      const std::size_t p_end = std::min((p + 1) * part_bitlength, n_);
      auto i = (p_begin == 0)
          ? bitmap.find_first()
          : bitmap.find_next(p_begin - 1);
      while (i < p_end && i != boost::dynamic_bitset<$u32>::npos) {
        b[i % part_bitlength] = true;
        i = bitmap.find_next(i);
      }
      // Compress the current partition.
      auto compressed_part_ptr = std::make_unique<B>(b);
      parts_.push_back(std::move(compressed_part_ptr));
    }
  }

  part(const part& other) = delete;
  part(part&& other) noexcept = default;
  part& operator=(const part& other) = delete;
  part& operator=(part&& other) noexcept = default;

  /// Return the name of the implementation.
  static std::string
  name() noexcept {
    return std::string("part<")
        + B::name() + std::string(",")
        + std::to_string(dtl::log_2(part_bitlength))
        + std::string(">");
  }

  /// Returns the length of the original bitmap.
  std::size_t __forceinline__
  size() const noexcept {
    return n_;
  }

  /// Return the size in bytes.
  std::size_t __forceinline__
  size_in_byte() const noexcept {
    std::size_t s = 0;
    for (std::size_t p = 0; p < parts_.size(); ++p) {
      const auto& part_ptr = parts_[p];
      s += part_ptr->size_in_byte();
    }
    s += parts_.size() * sizeof(void*);
    return s;
  }

  /// Returns the name of the instance including the most important parameters
  /// in JSON.
  std::string
  info() const noexcept {
    return "{\"name\":\"" + name() + "\""
        + ",\"n\":" + std::to_string(n_)
        + ",\"part_cnt\":" + std::to_string(parts_.size())
        + ",\"size\":" + std::to_string(size_in_byte())
        + "}";
  }

  /// For debugging purposes.
  void
  print(std::ostream& os) const noexcept {
    for (std::size_t i = 0; i < parts_.size(); ++i) {
      os << std::setw(4) << i << ": " << *parts_[i] << std::endl;
    }
  }

  //===--------------------------------------------------------------------===//
  // Read related functions.
  //===--------------------------------------------------------------------===//
  /// Returns the value of the bit at the given position.
  u1 __forceinline__
  test(const std::size_t pos) const noexcept {
    assert(pos < n_);
    const auto part_idx = pos / part_bitlength;
    //    if (part_idx >= parts_.size()) {
    //      return false;
    //    }
    const auto& part_ptr = parts_[part_idx];
    return part_ptr->test(pos % part_bitlength);
  }

  //===--------------------------------------------------------------------===//
  /// 1-run iterator for partitioned bitmaps.
  template<run_iterator_type iter_type>
  class iter {
    /// Reference to the outer instance.
    const part& part_bitmap_;

    /// The iterator type of the nested bitmap.
    using nested_iter_type =
        typename obtain_run_iterator<bitmap_type, iter_type>::type;

    // Hack to place nested iterators on heap memory.
    struct heap_iter {
      nested_iter_type iter;

      heap_iter(B* bitmap)
          : iter(obtain_run_iterator<bitmap_type, iter_type>::from(*bitmap)) {}
      heap_iter(const heap_iter& other) = delete;
      heap_iter(heap_iter&& other) noexcept = delete;
      heap_iter& operator=(const heap_iter& other) = delete;
      heap_iter& operator=(heap_iter&& other) noexcept = delete;
      ~heap_iter() = default;
    };

    //===------------------------------------------------------------------===//
    // Iterator state
    //===------------------------------------------------------------------===//
    /// The current partition.
    std::size_t current_part_idx_;
    /// The nested iterator (of the current partition).
    std::unique_ptr<heap_iter> part_iter_;
    /// Points to the beginning of a 1-fill.
    $u64 pos_;
    /// The length of the current 1-fill
    $u64 length_;
    //===------------------------------------------------------------------===//

  public:
    explicit iter(const part& part)
        : part_bitmap_(part),
          current_part_idx_(0),
          part_iter_(nullptr) {
      // Instantiate the nested iterator.
      while (true) {
        auto& current_part = part_bitmap_.parts_[current_part_idx_];
        assert(current_part != nullptr);
        part_iter_ = std::make_unique<heap_iter>(current_part.get());

        if (part_iter_->iter.end()) {
          // The current partition is empty. Continue with the next partition,
          // if there is any.
          ++current_part_idx_;
          if (current_part_idx_ == part_bitmap_.parts_.size()) {
            // The entire bitmap is 0.
            pos_ = part_bitmap_.n_;
            length_ = 0;
            return;
          }
          continue;
        }
        else {
          // Update the iterator state.
          pos_ = (current_part_idx_ * part_bitlength) + part_iter_->iter.pos();
          length_ = part_iter_->iter.length();
          return;
        }
      }
    }

    iter(iter&&) = default;

    /// Advance to the next partition if there is any.
    void __forceinline__
    next_part() { // TODO should be private
      part_iter_ = nullptr;
      while (true) {
        ++current_part_idx_;
        if (current_part_idx_ >= part_bitmap_.parts_.size()) {
          // Reached the end.
          pos_ = part_bitmap_.n_;
          length_ = 0;
          return;
        }

        auto& current_part = part_bitmap_.parts_[current_part_idx_];

        // Instantiate the iterator of the current partition.
        assert(current_part != nullptr);
        part_iter_ = std::make_unique<heap_iter>(current_part.get());
        if (part_iter_->iter.end()) {
          // The current partition is empty. Continue with the next partition,
          // if there is any.
          continue;
        }
        else {
          // Update the iterator state.
          pos_ = (current_part_idx_ * part_bitlength) + part_iter_->iter.pos();
          length_ = part_iter_->iter.length();
          return;
        }
      }
    }

    void __forceinline__
    next() {
      assert(!end());
      assert(part_iter_.get() != nullptr);
      assert(!part_iter_->iter.end());

      part_iter_->iter.next();
      if (!part_iter_->iter.end()) {
        // Update the iterator state.
        pos_ = (current_part_idx_ * part_bitlength) + part_iter_->iter.pos();
        length_ = part_iter_->iter.length();
        return;
      }

      // Advance to the next partition if there is any.
      next_part();
    }

    void __forceinline__
    skip_to(const std::size_t to_pos) {
      if (to_pos < (pos_ + length_)) {
        length_ -= to_pos - pos_;
        pos_ = to_pos;
        return;
      }

      const auto dst_part = to_pos / part_bitlength;
      if (dst_part != current_part_idx_) {
        assert(dst_part > 0);
        // Skip to the destination partition.
        current_part_idx_ = dst_part - 1;
        next_part();
        // Check if we have reached the end or if we skipped over the
        // destination position.
        if (end() || pos_ > to_pos) {
          return;
        }
      }
      assert(current_part_idx_ == (to_pos / part_bitlength));

      // Skip within the current partition, if necessary.
      if (pos_ < to_pos) {
        part_iter_->iter.skip_to(to_pos % part_bitlength);
        if (!part_iter_->iter.end()) {
          // Update the iterator state.
          pos_ = (current_part_idx_ * part_bitlength) + part_iter_->iter.pos();
          length_ = part_iter_->iter.length();
          return;
        }
        else {
          // We skipped over the destination partition because there was no
          // 1-run. Thus, we simply forward the iterator to the next partition
          // and return.
          next_part();
          return;
        }
      }
    }

    u1 __forceinline__
    end() const noexcept {
      return length_ == 0;
    }

    u64 __forceinline__
    pos() const noexcept {
      return pos_;
    }

    u64 __forceinline__
    length() const noexcept {
      return length_;
    }
  };

  using skip_iter_type = iter<run_iterator_type::SKIP>;
  using scan_iter_type = iter<run_iterator_type::SCAN>;

  skip_iter_type __forceinline__
  it() const {
    return skip_iter_type(*this);
  }

  scan_iter_type __forceinline__
  scan_it() const {
    return scan_iter_type(*this);
  }
  //===--------------------------------------------------------------------===//
};
//===----------------------------------------------------------------------===//
} // namespace dtl
