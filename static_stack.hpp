#pragma once

#include <array>
#include <cassert>
#include <cstddef>

#include <dtl/dtl.hpp>

//===----------------------------------------------------------------------===//
/// A statically sized stack.
template<typename T, std::size_t N>
class static_stack {

private:

  std::array<T, N> stack_;
  std::size_t cnt_;

public:

  __forceinline__ explicit
  static_stack() : cnt_(0) { }

  __forceinline__
  ~static_stack() = default;

  __forceinline__
  static_stack(const static_stack& other) = default;

  __forceinline__
  static_stack(static_stack&& other) noexcept = default;

  __forceinline__ static_stack&
  operator=(const static_stack& other) = default;

  __forceinline__ static_stack&
  operator=(static_stack&& other) noexcept = default;

  __forceinline__ bool
  empty() const {
    return cnt_ == 0;
  }

  __forceinline__ void
  push_back(const T& item) {
    stack_[cnt_] = item;
    cnt_++;
    assert(cnt_ <= N);
  }

  __forceinline__ void
  push_back(const T&& item) {
    stack_[cnt_] = item;
    cnt_++;
    assert(cnt_ <= N);
  }

  __forceinline__ T
  back() const {
    return stack_[cnt_ - 1];
    assert(cnt_ >= 0 && cnt_ <= N);
  }

  __forceinline__ void
  pop_back() {
    cnt_--;
    assert(cnt_ >= 0 && cnt_ <= N);
  }

};
//===----------------------------------------------------------------------===//
