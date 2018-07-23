#pragma once

#include <list>
#include <queue>
#include <vector>

#include <boost/dynamic_bitset.hpp>

#include "dynamic_tree_mask_lo.hpp"

namespace dtl {


//===----------------------------------------------------------------------===//
class dynamic_partitioned_tree_mask {
public:

  u64 N;
  u64 partition_cnt; // must be a power of two
  u64 part_n;
  u64 part_n_log2;
  u64 part_n_mask;

  using tree_mask_t = dynamic_tree_mask_lo;

  std::vector<tree_mask_t> tree_masks_;

  /// C'tor
  explicit
  dynamic_partitioned_tree_mask(const boost::dynamic_bitset<$u32>& bitmask, u64 partition_cnt = 8)
      : N(bitmask.size()), partition_cnt(partition_cnt),
        part_n(N / partition_cnt), part_n_log2(dtl::log_2(part_n)), part_n_mask(u64(-1) >> (64 - dtl::log_2(part_n))) {

    if (!dtl::is_power_of_two(N)) {
      throw std::invalid_argument("The length of the bitmask must be a power of two.");
    }

    if (!dtl::is_power_of_two(partition_cnt)) {
      throw std::invalid_argument("The numbers of partitions must be a power of two.");
    }

    if (N < partition_cnt) {
      throw std::invalid_argument("The number of partitions must be less than or equal to N.");
    }

    for ($u64 pid = 0; pid < partition_cnt; pid++) {
      u64 offset = part_n * pid;
//      std::bitset<part_n> part_bitmask;
      boost::dynamic_bitset<$u32> part_bitmask(part_n);
      for ($u64 i = 0; i < part_n; i++) { // TODO optimize later
        part_bitmask[i] = bitmask[i + offset];
      }
      tree_masks_.emplace_back(part_bitmask);
      std::cout << "is true: " << tree_masks_.back().all()
                << ", is false: " << tree_masks_.back().none()
                << ", size: " << tree_masks_.back().size_in_byte()
                << std::endl;
    }
  }


  /// Decodes the level-order encoding to a bitmap.
  boost::dynamic_bitset<$u32>
  to_bitset(){
    boost::dynamic_bitset<$u32> ret_val(N);
    for ($u64 pid = 0; pid < partition_cnt; pid++) {
      u64 offset = part_n * pid;
      auto part_bitmask = tree_masks_[pid].to_bitset();
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

  bool operator!=(dynamic_partitioned_tree_mask& other) const {
    // TODO
  }

  bool operator==(dynamic_partitioned_tree_mask& other) const {
    // TODO
  }

  /// Bitwise XOR without compression of the resulting tree
  dynamic_partitioned_tree_mask
  operator^(const dynamic_partitioned_tree_mask& other) const{
    // TODO
  }

  /// Bitwise XOR (range encoding)
  dynamic_partitioned_tree_mask
  xor_re(const dynamic_partitioned_tree_mask& other) const {
    // TODO
  }

  void
  print(std::ostream& os) const {
    // TODO
  }

  /// Bitwise AND without compression of the resulting tree
  dynamic_partitioned_tree_mask
  operator&(const dynamic_partitioned_tree_mask& other) const {
    // TODO
  }

  /// Bitwise AND (range encoding)
  dynamic_partitioned_tree_mask
  and_re(const dynamic_partitioned_tree_mask& other) const {
    // TODO
  }

  /// Computes (a XOR b) & this
  /// Note: this, a and b must be different instances. Otherwise, the behavior is undefined.
  dynamic_partitioned_tree_mask&
  fused_xor_and(const dynamic_partitioned_tree_mask& a, const dynamic_partitioned_tree_mask& b) const {
    // TODO
  }

  /// Bitwise XOR with compression of the resulting tree
  dynamic_partitioned_tree_mask
  xor_compressed(const dynamic_partitioned_tree_mask& other) const {
    // TODO
  }

  /// Bitwise AND with compression of the resulting tree
  dynamic_partitioned_tree_mask
  and_compressed(const dynamic_partitioned_tree_mask& other) const {
    // TODO
  }

  /// Return the name of the implementation.
  static std::string
  name() {
    return "dynamic_partitioned_tree_mask";
  }

  /// Returns the value of the bit at the position pos.
  u1
  test(const std::size_t pos) const {
    u64 tree_mask_idx = pos >> part_n_log2;
    u64 in_part_pos = pos & part_n_mask;
    return tree_masks_[tree_mask_idx].test(in_part_pos);
  }

};
//===----------------------------------------------------------------------===//

}; // namespace dtl