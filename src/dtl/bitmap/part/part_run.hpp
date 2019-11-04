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
/// Applies a fixed size partitioning to the given bitmap. The implementation
/// is optimized for very sparse or very dense bitmaps. If a partition contains
/// only a single run (either 0's or 1's), then no (compressed) bitmap is created
/// for that partition.
template<
    /// The (compressed) bitmap type.
    typename B,
    /// The partition size in bits.
    std::size_t P>
class part_run {
  static_assert(dtl::is_power_of_two(P),
      "The partition size must be a power of two.");

  /// Encapsulate a partition pointer. Note that partitions which only consist
  /// of 0's or 1's do not have a bitmap instance on the heap memory. Instead
  /// the pointer is tagged to indicate whether the the partition is 'single-
  /// valued'. In that case the 2nd most significant bit of the tagged pointer
  /// contains the actual (partition-wide) value.
  struct part_ptr_t {
    static constexpr auto ptr_tag_mask =
        uintptr_t(1) << ((sizeof(uintptr_t) * 8) - 1);
    static constexpr auto value_mask =
        uintptr_t(1) << ((sizeof(uintptr_t) * 8) - 2);

    /// The tagged pointer. False by default.
    B* ptr;

    /// C'tor
    inline part_ptr_t() : ptr(reinterpret_cast<B*>(ptr_tag_mask)) {}
    inline part_ptr_t(const part_ptr_t& other) = delete;
    inline part_ptr_t(part_ptr_t&& other) noexcept
        : ptr(std::move(other.ptr)) {
      other.ptr = nullptr;
    }
    inline part_ptr_t& operator=(const part_ptr_t& other) = delete;
    inline part_ptr_t& operator=(part_ptr_t&& other) noexcept {
      ptr = other.ptr;
      other.ptr = nullptr;
    };

    /// D'tor
    inline ~part_ptr_t() {
      if (is_pointer() && ptr != nullptr) {
        delete ptr;
      }
    }

    /// Returns true if the tagged pointer actually is a pointer, false
    /// otherwise.
    inline u1
    is_pointer() const noexcept {
      return (reinterpret_cast<uintptr_t>(ptr) & ptr_tag_mask) == 0;
    }

    /// Returns the pointer to the partition.
    inline B*
    get_pointer() const noexcept {
      assert(is_pointer());
      return ptr;
    }

    /// Sets the pointer. (transfers ownership)
    inline void
    set_pointer(B* p) noexcept {
      ptr = p;
      assert(is_pointer());
    }

    /// Returns true if the tagged pointer is actually a boolean value, false
    /// otherwise.
    inline u1
    is_value() const noexcept {
      return !is_pointer();
    }

    /// Returns the boolean value (which is encoded in this tagged pointer).
    inline u1
    get_value() const noexcept {
      assert(is_value());
      return (reinterpret_cast<uintptr_t>(ptr) & value_mask) != 0;
    }

    /// Sets the boolean value.
    inline void
    set_value(u1 val) noexcept {
      assert(is_value() || ptr == nullptr); // otherwise we would leak memory
      if (val == true) {
        ptr = reinterpret_cast<B*>(ptr_tag_mask | value_mask);
      }
      else {
        ptr = reinterpret_cast<B*>(ptr_tag_mask);
      }
    }
  };

  /// Pointer to the partitions.
  std::vector<part_ptr_t> parts_;
  /// The (total) length of the bitmap.
  std::size_t n_;

public:
  using bitmap_type = B;
  static constexpr std::size_t part_bitlength = P;

  /// C'tor (similar to all other implementations)
  explicit part_run(const boost::dynamic_bitset<$u32>& bitmap)
      : parts_(), n_(bitmap.size()) {
    const std::size_t part_cnt =
        (bitmap.size() + (part_bitlength - 1)) / part_bitlength;

    for (std::size_t p = 0; p < part_cnt; ++p) {
      // TODO avoid copy
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
      parts_.emplace_back();
      auto& part_ptr = parts_.back();
      const std::size_t b_count = b.count();
      if (b_count == 0) {
        assert(part_ptr.is_value());
        assert(part_ptr.get_value() == false);
        // Nothing to do here, because the part_ptr is by default false.
      }
      else if (b_count == part_bitlength) {
        assert(part_ptr.is_value());
        part_ptr.set_value(true);
        assert(part_ptr.get_value() == true);
      }
      else {
        B* compressed_part = new B(b);
        part_ptr.set_pointer(compressed_part);
      }
    }
  }

  part_run(const part_run& other) = delete;
  part_run(part_run&& other) noexcept = default;
  part_run& operator=(const part_run& other) = delete;
  part_run& operator=(part_run&& other) noexcept = default;

  /// Return the name of the implementation.
  static std::string
  name() noexcept {
    return std::string("part_run<")
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
      if (part_ptr.is_pointer()) {
        s += part_ptr.get_pointer()->size_in_byte();
      }
    }
    s += parts_.size() * sizeof(part_ptr_t);
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
      os << std::setw(4) << i << ": ";
      if (parts_[i].is_pointer()) {
        os << *parts_[i].get_pointer() << std::endl;
      }
      else {
        if (parts_[i].get_value()) {
          os << "single run (1)" << std::endl;
        }
        else {
          os << "single run (0)" << std::endl;
        }
      }
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
    if (part_ptr.is_value()) {
      return part_ptr.get_value();
    }
    return part_ptr.get_pointer()->test(pos % part_bitlength);
  }

  //===--------------------------------------------------------------------===//
  /// 1-run iterator for partitioned bitmaps.
  template<run_iterator_type iter_type>
  class iter {
    /// Reference to the outer instance.
    const part_run& part_bitmap_;

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
    explicit iter(const part_run& part)
        : part_bitmap_(part),
          current_part_idx_(0),
          part_iter_(nullptr) {
      // Find the first partition that is not 0.
      while (current_part_idx_ < part_bitmap_.parts_.size()
          && part_bitmap_.parts_[current_part_idx_].is_value()
          && part_bitmap_.parts_[current_part_idx_].get_value() == false) {
        ++current_part_idx_;
      }
      // Is the entire bitmap 0?
      if (current_part_idx_ == part_bitmap_.parts_.size()) {
        pos_ = part_bitmap_.n_;
        length_ = 0;
        return;
      }
      // Is the current partition 1?
      if (part_bitmap_.parts_[current_part_idx_].is_value()) {
        assert(part_bitmap_.parts_[current_part_idx_].is_value());
        pos_ = (current_part_idx_ * part_bitlength);
        length_ = part_bitlength;
        return;
      }
      // Instantiate the iterator of the current partition.
      while (true) {
        auto& current_part = part_bitmap_.parts_[current_part_idx_];
        assert(current_part.is_pointer());
        assert(current_part.get_pointer() != nullptr);
        part_iter_ = std::make_unique<heap_iter>(current_part.get_pointer());

        if (part_iter_->iter.end()) {
          // The current partition contains a bitmap, but seems to be empty.
          // Continue with the next partition, if there is any.
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
        if (current_part.is_value()) {
          if (current_part.get_value() == true) {
            // The entire partition is 1.
            pos_ = current_part_idx_ * part_bitlength;
            length_ = part_bitlength;
            return;
          }
          else {
            // The entire partition is 0.
            continue; // with the next partition
          }
        }
        else {
          // Instantiate the iterator of the current partition.
          assert(current_part.is_pointer());
          assert(current_part.get_pointer() != nullptr);
          part_iter_ = std::make_unique<heap_iter>(current_part.get_pointer());
          if (part_iter_->iter.end()) {
            // The current partition contains a bitmap, but seems to be empty.
            // Continue with the next partition, if there is any.
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
    }

    void __forceinline__
    next() {
      assert(!end());
      assert(part_iter_.get() != nullptr
          || (part_bitmap_.parts_[current_part_idx_].is_value())
              && part_bitmap_.parts_[current_part_idx_].get_value() == true);
      if (part_iter_) {
        assert(!part_iter_->iter.end());
        part_iter_->iter.next();
        if (!part_iter_->iter.end()) {
          // Update the iterator state.
          pos_ = (current_part_idx_ * part_bitlength) + part_iter_->iter.pos();
          length_ = part_iter_->iter.length();
          return;
        }
      }
      else {
        // The current partition is a 1-run.
        assert(part_bitmap_.parts_[current_part_idx_].get_value() == true);
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
        if (part_iter_) {
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
        else {
          // The current partition is a 1-run.
          assert(part_bitmap_.parts_[current_part_idx_].is_value());
          assert(part_bitmap_.parts_[current_part_idx_].get_value() == true);
          pos_ = to_pos;
          length_ = part_bitlength - (to_pos % part_bitlength);
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

  /// Set the i-th bit to the given value.
  void __forceinline__
  set(std::size_t i, u1 val) noexcept {
    const auto part_idx = i / part_bitlength;
    auto& current_part = parts_[part_idx];
    // Helper function to decompress a partition.
    auto decompress = [&]() {
      if (current_part.is_pointer()) {
        B* b = current_part.get_pointer();
        auto dec = dtl::to_bitmap_using_iterator(*b);
        return std::move(dec);
      }
      else {
        boost::dynamic_bitset<$u32> dec(n_);
        if (current_part.get_value() == true) {
          dec.flip();
        }
        return std::move(dec);
      }

    };
    // Helper function to compress a partition.
    auto compress_and_install = [&](boost::dynamic_bitset<$u32>& b) {
      // Delete outdated bitmap (if exists).
      if (current_part.is_pointer()) {
        delete current_part.get_pointer();
        current_part.set_pointer(nullptr);
      }
      // Is it a single-run partition?
      const auto b_count = b.count();
      if (b_count == 0 || b_count == part_bitlength) {
        current_part.set_value(b_count == part_bitlength);
      }
      else {
        // Compress.
        current_part.set_pointer(new B(b));
      }
    };

    // Decompress.
    auto dec = decompress();
    // Apply the update.
    dec[i % part_bitlength] = val;
    // Re-compress and install.
    compress_and_install(dec);
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
