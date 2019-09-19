#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/dtl.hpp>

#include <array>
#include <cassert>
#include <cstddef>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// A statically sized stack.
template<typename T, std::size_t N>
class alignas(64) static_stack {
public:
  std::array<T, N> stack_;
  std::size_t cnt_;

  __forceinline__
  static_stack() : cnt_(0) {}

  __forceinline__ ~static_stack() = default;

  __forceinline__
  static_stack(const static_stack& other) = default;

  __forceinline__
  static_stack(static_stack&& other) noexcept = default;

  __forceinline__ static_stack&
  operator=(const static_stack& other) = default;

  __forceinline__ static_stack&
  operator=(static_stack&& other) noexcept = default;

  __forceinline__ u1
  empty() const noexcept {
    return cnt_ == 0;
  }

  __forceinline__ void
  clear() noexcept {
    cnt_ = 0;
  }

  __forceinline__ void
  push_back(const T& item) noexcept {
    stack_[cnt_] = item;
    ++cnt_;
    assert(cnt_ <= N);
  }

  __forceinline__ void
  push(const T& item) noexcept {
    push_back(item);
  }

  __forceinline__ void
  push_back(const T&& item) noexcept {
    stack_[cnt_] = item;
    ++cnt_;
    assert(cnt_ <= N);
  }

  __forceinline__ void
  push(const T&& item) noexcept {
    push_back(item);
  }

  __forceinline__ T
  back() const noexcept {
    return stack_[cnt_ - 1];
    assert(cnt_ >= 0 && cnt_ <= N);
  }

  __forceinline__ T&
  back() noexcept {
    return stack_[cnt_ - 1];
    assert(cnt_ >= 0 && cnt_ <= N);
  }

  __forceinline__ T
  top() const noexcept {
    return back();
  }

  __forceinline__ T&
  top() noexcept {
    return back();
  }

  __forceinline__ void
  pop_back() noexcept {
    --cnt_;
    assert(cnt_ >= 0 && cnt_ <= N);
  }

  __forceinline__ void
  pop() noexcept {
    pop_back();
  }

  __forceinline__ T&
  push() noexcept {
    T& ret_val = stack_[cnt_];
    ++cnt_;
    return ret_val;
  }

  __forceinline__ const T&
  operator[](std::size_t i) const noexcept {
    return stack_[i];
  }

  __forceinline__ std::size_t
  size() const noexcept {
    return cnt_;
  }

  __forceinline__ void
  rewind(std::size_t i) noexcept {
    cnt_ = i;
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl