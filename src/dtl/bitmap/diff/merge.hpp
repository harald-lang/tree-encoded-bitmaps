#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/bitmap/util/bitmap_fun.hpp>
#include <dtl/bitmap/util/convert.hpp>
#include <dtl/dtl.hpp>

#include <boost/dynamic_bitset.hpp>

#include <memory>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
// Merge strategies for pending updates.
//===----------------------------------------------------------------------===//
/// Naive merge. - The bitmap is decompressed, the diff is applied on the
/// decompressed bitmap. Afterwards, the updated bitmap is re-compressed.
template<
    /// The (compressed) bitmap type.
    typename B,
    /// The differential data structure to use.
    typename D>
struct merge_naive {
  void
  merge(std::unique_ptr<B>& bitmap, std::unique_ptr<D>& diff) {
    // Decompress the bitmap.
    auto plain_bitmap = dtl::to_bitmap_using_iterator(*bitmap);
    auto diff_bitmap = dtl::to_bitmap_using_iterator(*diff);

    // Apply the diff.
    plain_bitmap ^= diff_bitmap;

    // Re-compress the bitmap.
    auto updated_bitmap = std::make_unique<B>(plain_bitmap);
    std::swap(bitmap, updated_bitmap);
  }

  static std::string
  name() {
    return "naive";
  }
};
//===----------------------------------------------------------------------===//
/// Naive merge using the iterator interface. - The bitmap is decompressed
/// whereas the diff is XORed on-the-fly. Afterwards, the updated bitmap is
/// re-compressed.
template<
    /// The (compressed) bitmap type.
    typename B,
    /// The differential data structure to use.
    typename D>
struct merge_naive_iter {
  void
  merge(std::unique_ptr<B>& bitmap, std::unique_ptr<D>& diff) {
    // Decompress the bitmap and XOR the diff on-the-fly.
    auto it = dtl::bitwise_xor_it(bitmap->scan_it(), diff->it());
    auto plain_bitmap = dtl::to_bitmap_from_iterator(it, bitmap->size());

    // Re-compress the bitmap.
    auto updated_bitmap = std::make_unique<B>(plain_bitmap);
    std::swap(bitmap, updated_bitmap);
  }

  static std::string
  name() {
    return "naive_iter";
  }
};
//===----------------------------------------------------------------------===//
/// In-place merge. Calls the (native) XOR function.
/// Note: This strategy requires B and D to be of the same type.
template<
    /// The (compressed) bitmap type.
    typename B,
    /// The differential data structure to use.
    typename D>
struct merge_inplace {
  static_assert(std::is_same<B, D>::value,
      "This strategy requires the bitmap and the diff to be of the same type.");
  void
  merge(std::unique_ptr<B>& bitmap, std::unique_ptr<D>& diff) {
    // TODO avoid the copy
    auto updated_bitmap = std::make_unique<B>(*bitmap);
    (*updated_bitmap) ^= *diff;
    // Try to save space.
    updated_bitmap->shrink();
    std::swap(bitmap, updated_bitmap);

  static std::string
  name() {
    return "inplace";
  }
};
//===----------------------------------------------------------------------===//
/// Mock. Do not merge at all.
template<
    /// The (compressed) bitmap type.
    typename B,
    /// The differential data structure to use.
    typename D>
struct merge_not {
  void
  merge(std::unique_ptr<B>& bitmap, std::unique_ptr<D>& diff) {
  }

  static std::string
  name() {
    return "not";
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
