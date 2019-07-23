#pragma once

#include "teb.hpp"

namespace dtl {
//===----------------------------------------------------------------------===//
class teb_wrapper {

  using builder_t = teb<3>;

  std::vector<builder_t::word_type> data_;

public:

  explicit
  teb_wrapper(const boost::dynamic_bitset<$u32>& bitmap, f64 fpr = 0.0) {
    builder_t builder(bitmap, fpr);
    const auto word_cnt = (builder.serialized_size_in_bytes()
        + (sizeof(builder_t::word_type) - 1)) / sizeof(builder_t::word_type);
    data_.resize(word_cnt);
    builder.serialize(data_.data());
  }

};
//===----------------------------------------------------------------------===//
} // namespace dtl