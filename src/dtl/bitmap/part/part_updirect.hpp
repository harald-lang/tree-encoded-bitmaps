#pragma once
//===----------------------------------------------------------------------===//
#include "part.hpp"

#include <dtl/bitmap/iterator.hpp>
#include <dtl/dtl.hpp>
#include <dtl/math.hpp>

#include <boost/dynamic_bitset.hpp>
#include <memory>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// Applies a fixed size partitioning to the given bitmap.
/// Updates are handled by decompressing/re-compressing the target partition.
template<
    /// The (compressed) bitmap type.
    typename B,
    /// The partition size in bits.
    std::size_t P>
class part_updirect
    : public part<B,P> {
public:
  /// C'tor (similar to all other implementations)
  part_updirect(const boost::dynamic_bitset<$u32>& bitmap)
      : part<B,P>(bitmap) {}

  part_updirect(const part_updirect& other) = delete;
  part_updirect(part_updirect&& other) noexcept = default;
  part_updirect& operator=(const part_updirect& other) = delete;
  part_updirect& operator=(part_updirect&& other) noexcept = default;

  /// Return the name of the implementation.
  static std::string
  name() noexcept {
    return std::string("part_updirect<")
        + B::name() + std::string(",")
        + std::to_string(dtl::log_2(P))
        + std::string(">");
  }

  /// Set the i-th bit to the given value.
  void __forceinline__
  set(std::size_t i, u1 val) noexcept {
    const auto part_idx = i / P;
    // Decompress the partition.
    auto dec = dtl::to_bitmap_using_iterator(*this->parts_[part_idx]);
    // Apply the update.
    dec[i % P] = val;
    // Re-compress and install the partition.
    auto compressed_part_ptr = std::make_unique<B>(dec);
    std::swap(compressed_part_ptr, this->parts_[part_idx]);
  }

  /// Does nothing. Just for compatibility reasons.
  template<typename M>
  void __forceinline__
  merge() {}
};
//===----------------------------------------------------------------------===//
} // namespace dtl
