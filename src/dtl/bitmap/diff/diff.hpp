#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/bitmap/bitwise_operations.hpp>
#include <dtl/dtl.hpp>

#include <boost/dynamic_bitset.hpp>

#include <memory>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// Extends a bitmap by a differential data structure. All updates are written
/// to the differential data structure which is merged on request.
///
/// The differential data structure could be of any bitmap type that implements
/// the set(idx, bool) function and a constructor T(size_t) to create an empty
/// bitmap.
///
/// Note that read accesses become more costly, as the XOR of the bitmap and
/// diff needs to be computed on-the-fly.
template<
    /// The (compressed) bitmap type.
    typename B,
    /// The differential data structure to use.
    typename D>
class diff {
  /// The (compressed) bitmap.
  std::unique_ptr<B> bitmap_;
  /// Contains the pending updates.
  std::unique_ptr<D> diff_;

public:

  using bitmap_type = B;
  using diff_type = D;

  /// C'tor (similar to all other implementations)
  diff(const boost::dynamic_bitset<$u32>& bitmap)
      : bitmap_(std::make_unique<B>(bitmap)),
        diff_(std::make_unique<D>(bitmap.size())) {}

  /// Return the name of the implementation.
  static std::string
  name() noexcept {
    return std::string("diff<")
        + B::name() + std::string(",") + D::name()
        + std::string(">");
  }

  /// Returns the length of the original bitmap.
  std::size_t __forceinline__
  size() const noexcept {
    return bitmap_->size();
  }

  /// Return the size in bytes.
  std::size_t __forceinline__
  size_in_byte() const noexcept {
    return bitmap_->size_in_byte() + diff_->size_in_byte()
        + 16; // two pointers
  }

  /// Returns the name of the instance including the most important parameters
  /// in JSON.
  std::string
  info() const noexcept {
    return "{\"name\":\"" + name() + "\""
        + ",\"bitmap\":" + bitmap_->info()
        + ",\"diff\":" + diff_->info()
        + "}";
  }

  /// For debugging purposes.
  void
  print(std::ostream& os) const noexcept {
    os << *bitmap_ << ", diff = " << *diff_;
  }

  //===--------------------------------------------------------------------===//
  // Update related functions. - Updates are forwarded to the differential data
  // structure, and merged only on request.
  //===--------------------------------------------------------------------===//
  /// Set the i-th bit to the given value.
  void __forceinline__
  set(std::size_t i, u1 val) noexcept {
    diff_->set(i, bitmap_->test(i) ^ val);
  }

  /// Applies the pending updates and clears the diff.
  template<typename M>
  void
  merge() {
    // Apply the updates.
    M merge_strategy;
    merge_strategy.merge(bitmap_, diff_);

    // Clear the diff.
    std::unique_ptr<D> empty_diff = std::make_unique<D>();
    std::swap(diff_, empty_diff);
  }

  //===--------------------------------------------------------------------===//
  // Read related functions.
  //===--------------------------------------------------------------------===//
  /// Returns the value of the bit at the given position.
  u1 __forceinline__
  test(const std::size_t pos) const noexcept {
    return bitmap_->test(pos) ^ diff_->test(pos);
  }

  using skip_iter_type =
      decltype(dtl::bitwise_xor_it(bitmap_->it(), diff_->it()));
  // We always use the skip iterator for the diff, as the diffs are supposed
  // to be small.
  using scan_iter_type =
      decltype(dtl::bitwise_xor_it(bitmap_->scan_it(), diff_->it()));

  /// Returns a 1-run iterator, with efficient skip support.
  skip_iter_type __forceinline__
  it() const noexcept {
    return dtl::bitwise_xor_it(bitmap_->it(), diff_->it());
  }

  /// Returns a 1-run iterator, with WITHOUT efficient skip support.
  scan_iter_type __forceinline__
  scan_it() const noexcept {
    return dtl::bitwise_xor_it(bitmap_->scan_it(), diff_->it());
  }

  /// Returns a pointer to the compressed bitmap.
  B* get_bitmap() const noexcept {
    return bitmap_.get();
  }

  /// Returns a pointer to the diff.
  D* get_diff() const noexcept {
    return diff_.get();
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
