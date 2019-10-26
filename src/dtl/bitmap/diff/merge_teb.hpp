#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/bitmap/bitwise_operations.hpp>
#include <dtl/bitmap/util/bitmap_fun.hpp>
#include <dtl/bitmap/util/convert.hpp>
#include <dtl/bitmap/util/mutable_bitmap_tree.hpp>
#include <dtl/dtl.hpp>

#include <boost/dynamic_bitset.hpp>

#include <memory>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
// Special merge strategies for TEBs.
//===----------------------------------------------------------------------===//
/// Transforms a TEB into a mutable bitmap tree (the intermediate
/// representation) and apply the updates there. Afterwards, the bitmap tree
/// is transformed back to a TEB.
template<
    /// The (compressed) bitmap type.
    typename B,
    /// The differential data structure to use.
    typename D>
struct merge_tree {
  void
  merge(std::unique_ptr<B>& bitmap, std::unique_ptr<D>& diff) {
    // Decompress the TEB into a bitmap tree.
    dtl::mutable_bitmap_tree<> mbt(*(bitmap->teb_));
    // Read the diff.
    auto it = diff->scan_it();
    while (!it.end()) {
      const auto b = it.pos();
      const auto e = it.length() + b;
      for (std::size_t i = b; i < e; ++i) {
        mbt.toggle(i); // TODO implement range-toggle
      }
      it.next();
    }

    // Re-compress the bitmap.
    auto updated_bitmap = std::make_unique<B>(std::move(mbt));
    std::swap(bitmap, updated_bitmap);
  }

  static std::string
  name() {
    return "tree_merge";
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
