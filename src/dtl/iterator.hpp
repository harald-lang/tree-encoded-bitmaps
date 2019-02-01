#pragma once

#include <iterator>

//===----------------------------------------------------------------------===//
// DTL experimental features
//===----------------------------------------------------------------------===//
namespace dtl {

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
};

}
//===----------------------------------------------------------------------===//
