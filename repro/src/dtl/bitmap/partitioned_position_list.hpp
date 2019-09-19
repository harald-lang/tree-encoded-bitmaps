#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/dtl.hpp>
#include <dtl/math.hpp>

#include <boost/dynamic_bitset.hpp>

#include <algorithm>
#include <cstddef>
#include <vector>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// Partitioned position list.
template<typename _block_type = $u32, typename _local_position_t = $u8>
struct partitioned_position_list {
  using position_t = uint32_t;
  using local_position_t = _local_position_t;

  /// Partition meta data.
  struct partition_info {
    /// The index of the first element in the partition.
    position_t begin;
    /// Offsets within the concatenated positions vector.
    position_t offset;

    void
    print(std::ostream& os) const {
      os << "[b:"
         << static_cast<u64>(begin)
         << ",o:"
         << static_cast<u64>(offset)
         << "]";
    }
  };

  /// The max number of entries per partition.
  static constexpr std::size_t partition_size =
      1ull << sizeof(local_position_t) * 8;

  /// The actual position list (all partitions concatenated).
  std::vector<local_position_t> positions_;
  /// Partition meta data.
  std::vector<partition_info> partitions_;

  /// The length of the range.
  $u64 n_;

  // TODO make private
  partitioned_position_list() = default;

  explicit partitioned_position_list(const boost::dynamic_bitset<_block_type>& in)
      : partitions_(), positions_(), n_(in.size()) {
    std::size_t current_pos = in.find_first();
    while (current_pos < in.size()) {
      push_back(static_cast<position_t>(current_pos));
      current_pos = in.find_next(current_pos);
    }
  }

  ~partitioned_position_list() = default;
  partitioned_position_list(const partitioned_position_list& other) = default;
  partitioned_position_list(partitioned_position_list&& other) noexcept = default;

  __forceinline__ partitioned_position_list&
  operator=(const partitioned_position_list& other) = default;

  __forceinline__ partitioned_position_list&
  operator=(partitioned_position_list&& other) noexcept = default;

  /// Return the size in bytes.
  std::size_t
  size_in_byte() const {
    return partitions_.size() * sizeof(partition_info) /* partitions */
        + sizeof(position_t) /* number of partitions */
        + positions_.size() * sizeof(local_position_t) /* positions */
        + sizeof(position_t) /* number of positions */
        + sizeof(n_) /* bit-length of the original bitmap */;
  }

  /// Return the length of the bitmap.
  std::size_t
  size() const {
    return n_;
  }

  /// Conversion to an boost::dynamic_bitset.
  boost::dynamic_bitset<_block_type>
  to_bitset() {
    boost::dynamic_bitset<_block_type> ret(n_);
    for (const partition_info& part : partitions_) {
      auto curr_local_pos = positions_[part.offset];
      ret[part.begin + curr_local_pos] = true;
      for (std::size_t i = part.offset + 1; i < positions_.size(); ++i) {
        if (positions_[i] <= curr_local_pos) {
          break;
        }
        curr_local_pos = positions_[i];
        ret[part.begin + curr_local_pos] = true;
      }
    }
    return ret;
  }

  // clang-format off
//  /// Bitwise AND
//  partitioned_position_list __forceinline__
//  operator&(const partitioned_position_list& other) const {
//    partitioned_position_list ret;
//    ret.n_ = n_;
//    auto this_it = positions_.begin();
//    const auto this_it_end = positions_.end();
//    auto other_it = other.positions_.begin();
//    const auto other_it_end = other.positions_.end();
//    while (!(this_it == this_it_end || other_it == other_it_end)) {
//      if (*this_it == *other_it) {
//        ret.positions_.push_back(*this_it);
//        ++this_it;
//        ++other_it;
//      }
//      else {
//        if (*this_it < *other_it) {
//          ++this_it;
//        }
//        else {
//          ++other_it;
//        }
//      }
//    }
//    return ret;
//  }
//
//  /// Bitwise AND (range encoding)
//  partitioned_position_list __forceinline__
//  and_re(const partitioned_position_list& other) const {
//    return *this & other; // nothing special here. fall back to AND
//  }
//
//  /// Bitwise XOR
//  partitioned_position_list __forceinline__
//  operator^(const partitioned_position_list& other) const {
//    partitioned_position_list ret;
//    ret.n_ = n_;
//
//    auto this_it = positions_.begin();
//    const auto this_it_end = positions_.end();
//    auto other_it = other.positions_.begin();
//    const auto other_it_end = other.positions_.end();
//    while (!(this_it == this_it_end || other_it == other_it_end)) {
//      if (*this_it < *other_it) {
//        ret.positions_.push_back(*this_it);
//        ++this_it;
//      }
//      else if (*other_it < *this_it) {
//        ret.positions_.push_back(*other_it);
//        ++other_it;
//      }
//      else {
//        ++this_it;
//        ++other_it;
//      }
//    }
//    if (this_it != this_it_end) {
//      while (this_it != this_it_end) {
//        ret.positions_.push_back(*this_it);
//        ++this_it;
//      }
//    }
//    if (other_it != other_it_end) {
//      while (other_it != other_it_end) {
//        ret.positions_.push_back(*other_it);
//        ++other_it;
//      }
//    }
//    return ret;
//  }
//
//  /// Bitwise XOR (range encoding)
//  partitioned_position_list __forceinline__
//  xor_re(const partitioned_position_list& other) const {
//    return *this ^ other; // nothing special here. fall back to XOR
//  }
//
//  /// Computes (a XOR b) & this
//  /// Note: this, a and b must be different instances. Otherwise, the behavior
//  /// is undefined.
//  __forceinline__ partitioned_position_list&
//  fused_xor_and(const partitioned_position_list& a, const partitioned_position_list& b) {
//    const auto x = a ^ b;
//    auto y = *this & x;
//    std::swap(positions_, y.positions_);
//    return *this;
//  }
  // clang-format on

  void
  print(std::ostream& os) const {
    os << "part: [";
    if (!partitions_.empty()) {
      os << partitions_[0];
      for (std::size_t i = 1; i < partitions_.size(); ++i) {
        os << "," << partitions_[i];
      }
    }
    os << "]";
    os << ", pos: [";
    if (!positions_.empty()) {
      os << static_cast<u64>(positions_[0]);
      for (std::size_t i = 1; i < positions_.size(); ++i) {
        os << "," << static_cast<u64>(positions_[i]);
      }
    }
    os << "]";
  }

  static std::string
  name() {
    return "partitioned_position_list_"
        + std::to_string(sizeof(local_position_t) * 8);
  }

  /// Returns the value of the bit at the position pos.
  u1 test(const std::size_t pos) const {
    // FIXME search the partitions first
    auto it = std::lower_bound(positions_.begin(), positions_.end(), pos);
    return *it == pos;
  }

  //===--------------------------------------------------------------------===//
  /// Iterator, with skip support.
  class iter {
    const partitioned_position_list& outer_;

    //===------------------------------------------------------------------===//
    // Iterator state
    //===------------------------------------------------------------------===//
    /// Read position within the list.
    $u64 partitions_read_pos_;
    /// Read position within the list.
    $u64 positions_read_pos_;
    /// Points to the beginning of the current range.
    $u64 range_begin_;
    /// The length of the current range.
    $u64 range_length_;
    //===------------------------------------------------------------------===//

  public:
    explicit __forceinline__
    iter(const partitioned_position_list& outer)
        : outer_(outer),
          partitions_read_pos_(0),
          positions_read_pos_(0),
          range_begin_(outer.positions_.empty()
                  ? outer_.n_
                  : outer_.positions_[0] + outer_.partitions_[0].begin),
          range_length_(outer.positions_.empty() ? 0 : 1) {
      const auto& parts = outer_.partitions_;
      const auto& pos = outer_.positions_;

      if (parts.empty()) {
        // Empty list.
        range_begin_ = outer_.n_;
        range_length_ = 0;
        return;
      }

      range_begin_ = parts[partitions_read_pos_].begin + pos[0];
      range_length_ = 1;
      ++positions_read_pos_;
      determine_length_of_current_run();
    }

    /// Determines the length of the current run.
    void
    determine_length_of_current_run() {
      const auto& parts = outer_.partitions_;
      const auto& pos = outer_.positions_;
      while (positions_read_pos_ < pos.size()
          && pos[positions_read_pos_ - 1] < pos[positions_read_pos_]
          && pos[positions_read_pos_] == range_begin_ + range_length_ - parts[partitions_read_pos_].begin) {
        ++positions_read_pos_;
        ++range_length_;
      }
      // Check, whether we reached the end of the partition.
      if (positions_read_pos_ < outer_.positions_.size()) {
        if (pos[positions_read_pos_] <= pos[positions_read_pos_ - 1]) {
          ++partitions_read_pos_;
        }
      }
    }

    void __forceinline__
    next() {
      const auto& parts = outer_.partitions_;
      const auto& pos = outer_.positions_;

      if (positions_read_pos_ < pos.size()) {
        range_begin_ =
            parts[partitions_read_pos_].begin + pos[positions_read_pos_];
        range_length_ = 1;
        ++positions_read_pos_;
        determine_length_of_current_run();
      }
      else {
        range_begin_ = outer_.n_;
        range_length_ = 0;
      }
    }

    void __forceinline__
    skip_to(const std::size_t to_pos) {
      assert(to_pos >= range_begin_);

      if (to_pos < (range_begin_ + range_length_)) {
        range_length_ -= to_pos - range_begin_;
        range_begin_ = to_pos;
        return;
      }

      const auto& parts = outer_.partitions_;
      const auto& pos = outer_.positions_;

      // Check if the destination position is within the current partition.
      if (to_pos > parts[partitions_read_pos_].begin + partition_size) {
        // Search the corresponding partition.
        auto part_search = std::lower_bound(
            parts.begin(), // + partitions_read_pos_,
            parts.end(),
            to_pos,
            [](const partition_info& part, const std::size_t pos) -> u1 {
              return (part.begin + partition_size) < pos;
            });
        partitions_read_pos_ = static_cast<u64>(std::distance(
            parts.begin(), part_search));
        if (part_search == parts.end()) {
          range_begin_ = outer_.n_;
          range_length_ = 0;
          return;
        }
      }

      const auto& part = parts[partitions_read_pos_];

      // Check partition boundary. - As multiple partitions are contiguous,
      // the target position may not be part of any partition.
      if (to_pos < part.begin) {
        positions_read_pos_ = part.offset;
        next();
        return;
      }

      // Search within the current partition.
      const auto part_offset_begin = part.offset;
      const auto part_offset_end =
          partitions_read_pos_ + 1 == outer_.partitions_.size()
          ? outer_.positions_.size()
          : (outer_.partitions_[partitions_read_pos_ + 1]).offset;
      const auto part_begin = pos.begin() + part_offset_begin;
      const auto part_end = pos.begin() + part_offset_end;
      auto pos_search = std::lower_bound(
          part_begin,
          part_end,
          to_pos - part.begin);
      if (pos_search == part_end) {
        // No match within the current partition. Forward the iterator to the
        // beginning of the next partition.
        ++partitions_read_pos_;
        positions_read_pos_ = part_offset_end;
        next();
      }
      else {
        range_begin_ = *pos_search + part.begin;
        positions_read_pos_ = std::distance(outer_.positions_.begin(), pos_search) + 1ull;
        range_length_ = 1;
        determine_length_of_current_run();
      }
      assert(range_begin_ <= outer_.n_);
      assert(range_begin_ + range_length_ <= outer_.n_);
    }

    u1 __forceinline__
    end() const noexcept {
      return range_length_ == 0;
    }

    u64 __forceinline__
    pos() const noexcept {
      return range_begin_;
    }

    u64 __forceinline__
    length() const noexcept {
      return range_length_;
    }
  };
  //===--------------------------------------------------------------------===//

  iter __forceinline__
  it() const {
    return std::move(iter(*this));
  }

  iter __forceinline__
  scan_it() const {
    return std::move(iter(*this));
  }

  /// Returns the name of the instance including the most important parameters
  /// in JSON.
  std::string
  info() const {
    return "{\"name\":\"" + name() + "\""
        + ",\"n\":" + std::to_string(n_)
        + ",\"size\":" + std::to_string(size_in_byte())
        + ",\"size_type_bytes\":" + std::to_string(sizeof(local_position_t))
        + ",\"partitions\":" + std::to_string(partitions_.size())
        + ",\"positions\":" + std::to_string(positions_.size())
        + "}";
  }

private:
  //===--------------------------------------------------------------------===//
  // Helper functions.
  //===--------------------------------------------------------------------===//
  /// Appends a new position. - Used during construction.
  inline void
  push_back(const position_t pos) {
    if (partitions_.empty()) {
      create_partition(pos);
    }
    if (pos - partitions_.back().begin >= partition_size) {
      create_partition(pos);
    }
    partition_info& part_info = partitions_.back();
    // Append the position to the position list.
    positions_.push_back(pos - part_info.begin);
  }

  /// Creates a new partition.
  inline void
  create_partition(const position_t pos) {
    partitions_.emplace_back();
    partition_info& part_info = partitions_.back();
    part_info.begin = pos;
    part_info.offset = static_cast<position_t>(positions_.size());
  }
  //===--------------------------------------------------------------------===//
};
//===----------------------------------------------------------------------===//
} // namespace dtl
