#pragma once

#include <cstddef>
#include <vector>

#include <dtl/dtl.hpp>
#include <dtl/math.hpp>

#include <boost/dynamic_bitset.hpp>

namespace dtl {

//===----------------------------------------------------------------------===//
/// Position list.
template<typename _block_type = $u32>
struct range_list {

  using position_t = uint32_t;

  struct range {
    position_t begin;
    position_t length;

    inline bool
    operator<(const range& other) {
      return begin < other.begin;
    }

    void
    print(std::ostream& os) const {
      os << "[" << begin << "," << length << "]";
    }
  };

  std::vector<range> ranges_;

  /// The number of bits.
  $u64 n_;

  // TODO make private
  range_list() = default;

  explicit
  range_list(const boost::dynamic_bitset<_block_type>& in)
    : ranges_(), n_(in.size()) {
    std::size_t current_begin = in.find_first();
    std::size_t current_length = 1;
    while (current_begin < n_) {
      while (current_begin + current_length < n_
          && in[current_begin + current_length]) {
        ++current_length;
      }
      ranges_.emplace_back(
          range{static_cast<position_t>(current_begin),
                static_cast<position_t>(current_length)});
      current_begin = in.find_next(current_begin + current_length);
      current_length = 1;
    }
  }

  ~range_list() = default;

  range_list(const range_list& other) = default;

  range_list(range_list&& other) noexcept = default;

  __forceinline__ range_list&
  operator=(const range_list& other) = default;

  __forceinline__ range_list&
  operator=(range_list&& other) noexcept = default;

  /// Return the size in bytes.
  std::size_t
  size_in_byte() const {
    return ranges_.size() * sizeof(range) /* ranges */
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
    for (range& r : ranges_) {
      for (std::size_t j = r.begin; j < (r.begin + r.length); ++j) {
        ret[j] = true;
      }
    }
    return ret;
  }

  /// Bitwise AND
  range_list __forceinline__
  operator&(const range_list& other) const {
    range_list ret;
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
  range_list __forceinline__
  and_re(const range_list& other) const {
    return *this & other; // nothing special here. fall back to AND
  }

  /// Bitwise XOR
  range_list __forceinline__
  operator^(const range_list& other) const {
    range_list ret;
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
  range_list __forceinline__
  xor_re(const range_list& other) const {
    return *this ^ other; // nothing special here. fall back to XOR
  }

  /// Computes (a XOR b) & this
  /// Note: this, a and b must be different instances. Otherwise, the behavior
  /// is undefined.
  __forceinline__ range_list&
  fused_xor_and(const range_list& a, const range_list& b) {
    const auto x = a ^ b;
    auto y = *this & x;
    std::swap(ranges_, y.ranges_);
    return *this;
  }

  void
  print(std::ostream& os) const {
    os << "[";
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
    return "range_list_"
        + std::to_string(sizeof(position_t) * 8);;
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

    const range_list& outer_;

    //===------------------------------------------------------------------===//
    // Iterator state
    //===------------------------------------------------------------------===//
    /// Read position within the list.
    $u64 read_pos_;
    /// Points to the beginning of the current range.
    $u64 range_begin_;
    /// The length of the current range.
    $u64 range_length_;
    //===------------------------------------------------------------------===//

  public:

    explicit __forceinline__
    iter(const range_list& outer)
        : outer_(outer),
          read_pos_(0),
          range_begin_(read_pos_ < outer.ranges_.size()
                       ? outer_.ranges_[0].begin
                       : outer_.n_),
          range_length_(read_pos_ < outer.ranges_.size()
                        ? outer_.ranges_[0].length
                        : 0) {
    }

    iter(iter&&) noexcept = default;
    iter&
    operator=(iter&& other) noexcept = default;

    void __forceinline__
    next() {
      ++read_pos_;
      if (read_pos_ < outer_.ranges_.size()) {
        range_begin_ = outer_.ranges_[read_pos_].begin;
        range_length_ = outer_.ranges_[read_pos_].length;
      }
      else {
        range_begin_ = outer_.n_;
        range_length_ = 0;
        read_pos_ = outer_.ranges_.size();
      }
    }

    void __forceinline__
    skip_to(const std::size_t to_pos) {
      auto search = std::lower_bound(
          outer_.ranges_.begin(), outer_.ranges_.end(),
          range{to_pos, 1},
          [](const range& lhs, const range& rhs) -> bool {
            return lhs.begin + lhs.length < rhs.begin + rhs.length;
          });
      if (search != outer_.ranges_.end()) {
        read_pos_ = std::distance(outer_.ranges_.begin(), search);
        if (to_pos < (*search).begin) {
          range_begin_ = (*search).begin;
          range_length_ = (*search).length;
        }
        else {
          range_begin_ = to_pos;
          range_length_ = (*search).length - (to_pos - (*search).begin);
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
    return std::move(iter(*this));
  }

  /// Returns the name of the instance including the most important parameters
  /// in JSON.
  std::string
  info() const {
    return "{\"name\":\"" + name() + "\""
        + ",\"n\":" + std::to_string(n_)
        + ",\"size\":" + std::to_string(size_in_byte())
        + ",\"ranges\":" + std::to_string(ranges_.size())
        + "}";
  }

};
//===----------------------------------------------------------------------===//


} // namespace dtl
