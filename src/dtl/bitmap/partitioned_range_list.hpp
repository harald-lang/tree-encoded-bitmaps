#pragma once

#include <cstddef>
#include <vector>

#include <dtl/dtl.hpp>
#include <dtl/math.hpp>

#include <boost/dynamic_bitset.hpp>

namespace dtl {

//===----------------------------------------------------------------------===//
/// Partitioned range list.
template<typename _block_type = $u32, typename _local_position_t = $u8>
struct partitioned_range_list {

  using position_t = uint32_t;
  using local_position_t = _local_position_t;

  struct range {
    local_position_t begin;
    local_position_t length;

    inline bool
    operator<(const range& other) {
      return begin < other.begin;
    }

    void
    print(std::ostream& os) const {
      os << "["
         << static_cast<u64>(begin) << ","
         << static_cast<u64>(length)
         << "]";
    }
  };

  /// Partition meta data.
  struct partition_info {
    /// The index of the first element in the partition.
    position_t begin;
    /// Offset within the concatenated range vector.
    position_t offset;

    inline bool operator<(const range& other) { return begin < other.begin; }

    void
    print(std::ostream& os) const {
      os << "["
         << static_cast<u64>(begin) << ","
         << static_cast<u64>(offset)
         << "]";
    }

  };

  /// The number of entries per partition.
  static constexpr std::size_t partition_size =
      1ull << (sizeof(local_position_t) * 8 - 1);

  /// The actual range list (all partitions concatenated).
  std::vector<range> ranges_;
  /// Partition meta data.
  std::vector<partition_info> partitions_;

  /// The length of the range.
  $u64 n_;

  // TODO make private
  partitioned_range_list() = default;

  explicit
  partitioned_range_list(const boost::dynamic_bitset<_block_type>& in)
    : partitions_(), ranges_(), n_(in.size()) {
    std::size_t current_begin = in.find_first();
    std::size_t current_length = 1;
    while (current_begin < n_) {
      while (current_begin + current_length < n_
          && in[current_begin + current_length]) {
        ++current_length;
      }
      push_back(current_begin ,current_length);
      current_begin = in.find_next(current_begin + current_length);
      current_length = 1;
    }
  }

  ~partitioned_range_list() = default;

  partitioned_range_list(const partitioned_range_list& other) = default;

  partitioned_range_list(partitioned_range_list&& other) noexcept = default;

  __forceinline__ partitioned_range_list&
  operator=(const partitioned_range_list& other) = default;

  __forceinline__ partitioned_range_list&
  operator=(partitioned_range_list&& other) noexcept = default;

  /// Return the size in bytes.
  std::size_t
  size_in_byte() const {
    return partitions_.size() * sizeof(partition_info) /* partitions */
        + sizeof(position_t) /* number of partitions */
        + ranges_.size() * sizeof(range) /* ranges */
        + sizeof(position_t) /* number of ranges */
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
      auto curr_local_range = ranges_[part.offset];
      set_bits(ret,
               part.begin + curr_local_range.begin,
               part.begin + curr_local_range.begin + curr_local_range.length);
      for (std::size_t i = part.offset + 1; i < ranges_.size(); ++i) {
        if (ranges_[i].begin  <= curr_local_range.begin + curr_local_range.length) {
          break;
        }
        curr_local_range = ranges_[i];
        set_bits(ret,
                 part.begin + curr_local_range.begin,
                 part.begin + curr_local_range.begin + curr_local_range.length);
      }
    }
    return ret;
  }

  /// Bitwise AND
  partitioned_range_list __forceinline__
  operator&(const partitioned_range_list& other) const {
    partitioned_range_list ret;
    ret.n_ = n_;
    auto this_it = ranges_.begin();
    const auto this_it_end = ranges_.end();
    auto other_it = other.ranges_.begin();
    const auto other_it_end = other.ranges_.end();
    while (!(this_it == this_it_end || other_it == other_it_end)) {
      if (*this_it == *other_it) {
        ret.ranges_.push_back(*this_it);
        ++this_it;
        ++other_it;
      }
      else {
        if (*this_it < *other_it) {
          ++this_it;
        }
        else {
          ++other_it;
        }
      }
    }
    return ret;
  }

  /// Bitwise AND (range encoding)
  partitioned_range_list __forceinline__
  and_re(const partitioned_range_list& other) const {
    return *this & other; // nothing special here. fall back to AND
  }

  /// Bitwise XOR
  partitioned_range_list __forceinline__
  operator^(const partitioned_range_list& other) const {
    partitioned_range_list ret;
    ret.n_ = n_;

    auto this_it = ranges_.begin();
    const auto this_it_end = ranges_.end();
    auto other_it = other.ranges_.begin();
    const auto other_it_end = other.ranges_.end();
    while (!(this_it == this_it_end || other_it == other_it_end)) {
      if (*this_it < *other_it) {
        ret.ranges_.push_back(*this_it);
        ++this_it;
      }
      else if (*other_it < *this_it) {
        ret.ranges_.push_back(*other_it);
        ++other_it;
      }
      else {
        ++this_it;
        ++other_it;
      }
    }
    if (this_it != this_it_end) {
      while (this_it != this_it_end) {
        ret.ranges_.push_back(*this_it);
        ++this_it;
      }
    }
    if (other_it != other_it_end) {
      while (other_it != other_it_end) {
        ret.ranges_.push_back(*other_it);
        ++other_it;
      }
    }
    return ret;
  }

  /// Bitwise XOR (range encoding)
  partitioned_range_list __forceinline__
  xor_re(const partitioned_range_list& other) const {
    return *this ^ other; // nothing special here. fall back to XOR
  }

  /// Computes (a XOR b) & this
  /// Note: this, a and b must be different instances. Otherwise, the behavior
  /// is undefined.
  __forceinline__ partitioned_range_list&
  fused_xor_and(const partitioned_range_list& a, const partitioned_range_list& b) {
    const auto x = a ^ b;
    auto y = *this & x;
    std::swap(ranges_, y.ranges_);
    return *this;
  }

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
    os << ", range: [";
    if (!ranges_.empty()) {
      os << ranges_[0];
      for (std::size_t i = 1; i < ranges_.size(); ++i) {
        os << "," << ranges_[i];
      }
    }
    os << "]";
  }

  static std::string
  name() {
    return "partitioned_range_list_"
        + std::to_string(sizeof(local_position_t) * 8);
  }

  /// Returns the value of the bit at the position pos.
  u1
  test(const std::size_t pos) const {
    auto it = std::lower_bound(ranges_.begin(), ranges_.end(), pos);
    return *it == pos;
  }


  //===--------------------------------------------------------------------===//
  /// Iterator, with skip support.
  class iter {

    const partitioned_range_list& outer_;

    //===------------------------------------------------------------------===//
    // Iterator state
    //===------------------------------------------------------------------===//
    /// Read position within the list.
    $u64 partitions_read_pos_;
    /// Read position within the list.
    $u64 ranges_read_pos_;
    /// Points to the beginning of the current range.
    $u64 range_begin_;
    /// The length of the current range.
    $u64 range_length_;
    //===------------------------------------------------------------------===//

  public:

    explicit __forceinline__
    iter(const partitioned_range_list& outer)
        : outer_(outer),
          partitions_read_pos_(0),
          ranges_read_pos_(0),
          range_begin_(outer.ranges_.empty()
                       ? outer_.n_
                       : outer_.partitions_[0].begin + outer_.ranges_[0].begin),
          range_length_(outer.ranges_.empty()
                        ? 0ul : outer_.ranges_[0].length) {
    }

    void __forceinline__
    next() {
      ++ranges_read_pos_;
      if (ranges_read_pos_ < outer_.ranges_.size()) {
        if (outer_.ranges_[ranges_read_pos_].begin
                <= outer_.ranges_[ranges_read_pos_ - 1].begin) {
          ++partitions_read_pos_;
        }
        range_begin_ = outer_.partitions_[partitions_read_pos_].begin
            + outer_.ranges_[ranges_read_pos_].begin;
        range_length_ = outer_.ranges_[ranges_read_pos_].length;
      }
      else {
        range_begin_ = outer_.n_;
        range_length_ = 0;
        ranges_read_pos_ = outer_.ranges_.size();
      }
    }

    void __forceinline__
    skip_to(const std::size_t to_pos) {
      constexpr position_t dont_care = 0;
      // Find the partition.
      auto part_search = std::lower_bound(
          outer_.partitions_.begin(), outer_.partitions_.end(),
          partition_info {to_pos, dont_care},
          [](const partition_info& lhs, const partition_info& rhs) -> bool {
            return lhs.begin + partition_size < rhs.begin;
          });
      if (part_search != outer_.partitions_.end()) {
        // Find the range (within the partition).
        const auto to_pos_local = to_pos - (*part_search).begin;

        const auto part_begin = (*part_search).offset;
        const auto part_end = part_search + 1 == outer_.partitions_.end()
            ? outer_.ranges_.size()
            : (*(part_search + 1)).offset;
        auto range_search = std::lower_bound(
            outer_.ranges_.begin() + part_begin,
            outer_.ranges_.begin() + part_end, range{to_pos_local, dont_care},
            [](const range& lhs, const range& rhs) -> bool {
              return lhs.begin + lhs.length <= rhs.begin;
            });
        if (range_search != outer_.ranges_.begin() + part_end) {
          ranges_read_pos_ = std::distance(outer_.ranges_.begin(), range_search) + 1ull;
          if (to_pos_local < (*range_search).begin) {
            // No match at the given position. Forward to the next match.
            range_begin_ = (*part_search).begin + (*range_search).begin;
            range_length_ = (*range_search).length;
          }
          else {
            // The given skip-to position matches with a range, but may not
            // necessarily match exactly with the beginning. Thus, the length
            // of the matching range may need to be adjusted.
            range_begin_ = to_pos;
            range_length_ = (*range_search).length - to_pos_local
                - (*range_search).begin;
          }
        }
        else {
          range_begin_ = outer_.n_;
          range_length_ = 0;
        }
      }
      else {
        range_begin_ = outer_.n_;
        range_length_ = 0;
      }
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
    return iter(*this);
  }

  /// Returns the name of the instance including the most important parameters
  /// in JSON.
  std::string
  info() {
    return "{\"name\":\"" + name() + "\""
        + ",\"n\":" + std::to_string(n_)
        + ",\"size\":" + std::to_string(size_in_byte())
        + ",\"size_type_bytes\":" + std::to_string(sizeof(local_position_t))
        + ",\"partitions\":" + std::to_string(partitions_.size())
        + ",\"ranges\":" + std::to_string(ranges_.size())
        + "}";
  }

private:

  //===--------------------------------------------------------------------===//
  // Helper functions.
  //===--------------------------------------------------------------------===//
  /// Appends a new range. - Used during construction.
  /// The given range is split into smaller ranges if necessary.
  inline void
  push_back(const position_t begin, const position_t length) {
    if (partitions_.empty()) {
      create_partition(begin);
    }
    if (begin - partitions_.back().begin >= partition_size) {
      create_partition(begin);
    }
    partition_info& part_info = partitions_.back();
    u1 exceeds_partition = (begin + length > part_info.begin + partition_size);
    range r;
    r.begin = static_cast<local_position_t>(begin - part_info.begin);
    r.length = exceeds_partition
        ? static_cast<local_position_t>(partition_size - r.begin)
        : static_cast<local_position_t>(length);
    ranges_.push_back(r);
    if (! exceeds_partition) return;
    push_back(begin + r.length, length - r.length);
  }

  /// Creates a new partition.
  inline void
  create_partition(const position_t pos) {
    partitions_.emplace_back();
    partition_info& part_info = partitions_.back();
    part_info.begin = pos;
    part_info.offset = static_cast<position_t>(ranges_.size());
  }

  /// Sets the bits [from, to) in the given bitset.
  inline void
  set_bits(boost::dynamic_bitset<_block_type>& bitset,
           const std::size_t from, const std::size_t to) {
    for (std::size_t i = from; i < to; ++i) {
      bitset[i] = true;
    }
  }
  //===--------------------------------------------------------------------===//

};
//===----------------------------------------------------------------------===//


} // namespace dtl
