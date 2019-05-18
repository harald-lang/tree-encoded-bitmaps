#pragma once

#include <vector>

#include <dtl/dtl.hpp>
#include <dtl/iterator.hpp>
#include <dtl/bitmap/util/bitmap_fun.hpp>

namespace dtl {

//===----------------------------------------------------------------------===//
template<typename _word_type>
struct bitmap_view {

  using word_type = _word_type;
  static constexpr std::size_t word_bitlength = sizeof(word_type) * 8;

  using fn = dtl::bitmap_fun<word_type>;
  data_view<const word_type> data_;

  explicit
  bitmap_view(data_view<const word_type> data)
      : data_(data) {
    assert(data_.size() > 0);
  }

  bitmap_view()
      : data_(data_view<const word_type> {nullptr, nullptr}) {
  }

  void
  init(data_view<const word_type> data) {
    data_ = data;
  }

  u1 __forceinline__
  operator[](std::size_t pos) const {
    return fn::test(data_.begin(), pos);
  }

  std::size_t
  find_first() const {
    return fn::find_first(data_.begin(), data_.end());
  }

  std::size_t
  find_last() const {
    return fn::find_last(data_.begin(), data_.end());
  }


};
//===----------------------------------------------------------------------===//


} // namespace dtl

