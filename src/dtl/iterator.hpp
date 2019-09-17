#pragma once
//===----------------------------------------------------------------------===//
#include <iterator>
//===----------------------------------------------------------------------===//
// DTL experimental features
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// Range for arbitrary memory regions.
template<typename T>
struct data_view {
  T* p_begin;
  T* p_end;
  T* begin() const { return p_begin; }
  T* end() const { return p_end; }
  std::size_t size() const {
    return std::distance(std::begin(*this), std::end(*this));
  }
  inline bool
  operator==(const data_view& other){
    if (begin() == other.begin() && end() == other.end()) return true;
    const auto len = end() - begin();
    const auto other_len = other.end() - other.begin();
    if (len != other_len) return false;
    for (std::size_t i = 0; i < len; ++i) {
      if ((begin() + i) != (other.begin() + i)) return false;
    }
    return true;
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
//===----------------------------------------------------------------------===//
