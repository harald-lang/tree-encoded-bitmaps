#pragma once

#include <list>
#include <queue>
#include <vector>

#include "tree_mask_lo.hpp"

namespace dtl {


//===----------------------------------------------------------------------===//
template<std::size_t _N>
class partitioned_tree_mask {
public:

  static constexpr u64 N = _N;
  static constexpr u64 partition_cnt = 8; // must be a power of two
  static constexpr u64 part_n = N / partition_cnt;

  static_assert(dtl::is_power_of_two(_N), "Template parameter N must be a power of two.");
  static_assert(dtl::is_power_of_two(partition_cnt), "The number of partitions must be a power of two.");
  static_assert(partition_cnt <= _N, "The number of partitions must be less than or equal to N.");

  using tree_mask_t = tree_mask_lo<part_n>;

  std::vector<tree_mask_t> tree_masks_;

  /// C'tor
  explicit
  partitioned_tree_mask(const std::bitset<N>& bitmask) {
    for ($u64 pid = 0; pid < partition_cnt; pid++) {
      u64 offset = part_n * pid;
      std::bitset<part_n> part_bitmask;
      for ($u64 i = 0; i < part_n; i++) { // TODO optimize later
        part_bitmask[i] = bitmask[i + offset];
      }
      tree_masks_.emplace_back(part_bitmask);
      std::cout << "is true: " << tree_masks_.back().all()
                << ", is false: " << tree_masks_.back().none()
                << std::endl;
    }
  }


  /// Decodes the level-order encoding to a bitmap.
  std::bitset<N>
  to_bitset(){
    std::bitset<N> ret_val;
    for ($u64 pid = 0; pid < partition_cnt; pid++) {
      u64 offset = part_n * pid;
      std::bitset<part_n> part_bitmask = tree_masks_[pid].to_bitset();
      for ($u64 i = 0; i < part_n; i++) { // TODO optimize later
        ret_val[i + offset] = part_bitmask[i];
      }
    }
    return ret_val;
  }

  /// Return the size in bytes.
  std::size_t
  size_in_byte() {
    $u64 size = 0;
    for ($u64 pid = 0; pid < partition_cnt; pid++) {
      size += tree_masks_[pid].size_in_byte();
    }
    size += partition_cnt * 8; // pointer to tree masks
    return size;
  }

  bool operator!=(partitioned_tree_mask& other) const {
    // TODO
  }

  bool operator==(partitioned_tree_mask& other) const {
    // TODO
  }

  /// Bitwise XOR without compression of the resulting tree
  partitioned_tree_mask
  operator^(const partitioned_tree_mask& other) const{
    // TODO
  }

  /// Bitwise XOR (range encoding)
  partitioned_tree_mask
  xor_re(const partitioned_tree_mask& other) const {
    // TODO
  }

  void
  print(std::ostream& os) const {
    // TODO
  }

  /// Bitwise AND without compression of the resulting tree
  partitioned_tree_mask
  operator&(const partitioned_tree_mask& other) const {
    // TODO
  }

  /// Bitwise AND (range encoding)
  partitioned_tree_mask
  and_re(const partitioned_tree_mask& other) const {
    // TODO
  }

  /// Computes (a XOR b) & this
  /// Note: this, a and b must be different instances. Otherwise, the behavior is undefined.
  partitioned_tree_mask&
  fused_xor_and(const partitioned_tree_mask& a, const partitioned_tree_mask& b) const {
    // TODO
  }

  /// Bitwise XOR with compression of the resulting tree
  partitioned_tree_mask
  xor_compressed(const partitioned_tree_mask& other) const {
    // TODO
  }

  /// Bitwise AND with compression of the resulting tree
  partitioned_tree_mask
  and_compressed(const partitioned_tree_mask& other) const {
    // TODO
  }

  /// Return the name of the implementation.
  static std::string
  name() {
    return "partitioned_tree_mask";
  }

  /// Returns the value of the bit at the position pos.
  u1
  test(const std::size_t pos) const {
    u64 tree_mask_idx = pos / part_n;
    u64 in_part_pos = pos % part_n;
    return tree_masks_[tree_mask_idx].test(in_part_pos);
  }

};
//===----------------------------------------------------------------------===//

}; // namespace dtl