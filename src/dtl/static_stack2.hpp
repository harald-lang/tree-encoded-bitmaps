#pragma once

#include <array>
#include <cassert>
#include <cstddef>

#include <dtl/dtl.hpp>

namespace dtl {

//===----------------------------------------------------------------------===//
/// A statically sized stack.
template<typename T, std::size_t N>
class static_stack2 {

private:

  std::array<T, N> stack_;
  std::size_t cnt_;

  T top_element;

public:

  __forceinline__
  static_stack2() : cnt_(0) { }

  __forceinline__
  ~static_stack2() = default;

  __forceinline__
  static_stack2(const static_stack2& other) = default;

  __forceinline__
  static_stack2(static_stack2&& other) noexcept = default;

  __forceinline__ static_stack2&
  operator=(const static_stack2& other) = default;

  __forceinline__ static_stack2&
  operator=(static_stack2&& other) noexcept = default;

  __forceinline__ u1
  empty() const {
    return cnt_ == 0;
  }

  __forceinline__ void
  clear() {
    cnt_ = 0;
  }

  __forceinline__ void
  push_back(const T& item) {
    stack_[cnt_] = item;
    top_element = item;
    ++cnt_;
    assert(cnt_ <= N);
  }

  __forceinline__ void
  push(const T& item) {
    push_back(item);
  }

  __forceinline__ void
  push_back(const T&& item) {
    stack_[cnt_] = item;
    top_element = item;
    ++cnt_;
    assert(cnt_ <= N);
  }

  __forceinline__ void
  push(const T&& item) {
    push_back(item);
  }

  __forceinline__ T
  back() const {
//    return stack_[cnt_ - 1];
    return top_element;
    assert(cnt_ >= 0 && cnt_ <= N);
  }

  __forceinline__ T
  top() const {
    return back();
  }

  __forceinline__ void
  pop_back() {
    --cnt_;
    top_element = stack_[cnt_ - 1 > N ? 0 : cnt_ - 1];
    assert(cnt_ >= 0 && cnt_ <= N);
  }

  __forceinline__ void
  pop() {
    pop_back();
  }

};
//===----------------------------------------------------------------------===//
} // namespace dtl