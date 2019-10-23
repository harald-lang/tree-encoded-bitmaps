#pragma once
//===----------------------------------------------------------------------===//
#include "bitmap_fun.hpp"
#include "bitmap_writer.hpp"

#include <dtl/dtl.hpp>
#include <dtl/iterator.hpp>

#include <vector>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
template<typename _word_type>
struct bitmap_view {
  using word_type = _word_type;
  static constexpr std::size_t word_bitlength = sizeof(word_type) * 8;

  using fn = dtl::bitmap_fun<word_type>;
  data_view<word_type> data_;

  explicit bitmap_view(data_view<word_type> data)
      : data_(data) {
    assert(data_.size() > 0);
  }

  explicit bitmap_view(word_type* begin, word_type* end)
      : data_(dtl::data_view<word_type> { begin, end }) {
  }

  bitmap_view()
      : data_(data_view<word_type> { nullptr, nullptr }) {
  }

  void
  init(data_view<word_type> data) {
    data_ = data;
  }

  u1 __forceinline__
  operator[](std::size_t pos) const {
    assert(pos < (data_.size() * word_bitlength));
    return fn::test(data_.begin(), pos);
  }

  /// Test the i-th bit.
  u1 __forceinline__
  test(std::size_t i) const {
    assert((i / word_bitlength) < data_.size());
    return fn::test(data_.begin(), i);
  }

  std::size_t __forceinline__
  find_first() const {
    return fn::find_first(data_.begin(), data_.end());
  }

  std::size_t __forceinline__
  find_last() const {
    return fn::find_last(data_.begin(), data_.end());
  }

  void __forceinline__
  set(std::size_t i, u1 val) {
    assert(i < (data_.size() * word_bitlength));
    fn::set(data_.begin(), i, val);
  }

  /// Set the bits in [b,e) to 0.
  void __forceinline__
  clear(std::size_t b, std::size_t e) {
    if (e <= b) return;
    assert((b / word_bitlength) < data_.size());
//    assert(b <= n_);
//    assert(e <= n_);
    assert(b <= e);
    fn::clear(data_.begin(), b, e);
  }

  void
  print(std::ostream& os) const noexcept {
    const auto word_cnt = data_.size();
    for (std::size_t i = 0; i < word_cnt * word_bitlength; ++i) {
      os << (test(i) ? "1" : "0");
    }
  }

  bitmap_writer<word_type>
  writer(std::size_t start_idx) {
    return bitmap_writer<word_type>(data_.begin(), start_idx);
  }

  bitmap_limit_writer<word_type>
  writer(std::size_t start_idx, std::size_t end_idx /* non-inclusive */) {
    return bitmap_limit_writer<word_type>(data_.begin(), start_idx, end_idx);
  }
};
//===----------------------------------------------------------------------===//

} // namespace dtl
