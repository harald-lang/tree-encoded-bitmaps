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
/// Updates are forwarded to the internally used bitmap type. This requires the
/// bitmap type to be updatable.
template<
    /// The (compressed) bitmap type.
    typename B,
    /// The partition size in bits.
    std::size_t P>
class part_upforward
    : public part<B,P> {
public:
  /// C'tor (similar to all other implementations)
  explicit part_upforward(const boost::dynamic_bitset<$u32>& bitmap)
      : part<B,P>(bitmap) {}

  part_upforward(const part_upforward& other) = delete;
  part_upforward(part_upforward&& other) noexcept = default;
  part_upforward& operator=(const part_upforward& other) = delete;
  part_upforward& operator=(part_upforward&& other) noexcept = default;

  /// Return the name of the implementation.
  static std::string
  name() noexcept {
    return std::string("part_upforward<")
        + B::name() + std::string(",")
        + std::to_string(dtl::log_2(P))
        + std::string(">");
  }

  /// Set the i-th bit to the given value.
  void __forceinline__
  set(std::size_t i, u1 val) noexcept {
    // Forward the call.
    const auto part_idx = i / P;
    this->parts_[part_idx]->set(i % P, val);
  }

  /// Apply the pending updates and clear the diff.
  template<typename M>
  void
  merge() {
    for (std::size_t part_idx = 0; part_idx < this->parts_.size(); ++part_idx) {
      this->parts_[part_idx]->template merge<M>();
    }
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
