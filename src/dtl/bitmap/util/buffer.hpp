#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/dtl.hpp>

#include <cassert>
#include <cstddef>
#include <cstring>
#include <memory>
#include <type_traits>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
/// A simple memory buffer which allows allocations without initialization.
template<
    typename _word_type,
    typename _alloc = std::allocator<_word_type>>
class buffer {

  using word_type = _word_type;
  static constexpr std::size_t word_bitlength = sizeof(word_type) * 8;

  /// The number of elements.
  std::size_t size_;
  /// Pointer to the buffer.
  word_type* ptr_;
  /// The allocator.
  _alloc allocator_;

public:
  /// C'tor
  explicit
  buffer(std::size_t size, u1 init = true, const _alloc& alloc = _alloc())
      : size_(size), ptr_(nullptr), allocator_(alloc) {
    ptr_ = allocator_.allocate(size_);
    if (init) {
      std::memset(ptr_, 0, size_ * sizeof(word_type));
    }
  }

  buffer(const buffer& other)
      : size_(other.size_), ptr_(nullptr), allocator_(other.allocator_) {
    ptr_ = allocator_.allocate(size_);
    std::memcpy(ptr_, other.ptr_, size_ * sizeof(word_type));
  }

  buffer(buffer&& other) noexcept
      : size_(other.size_), ptr_(other.ptr_),
        allocator_(std::move(other.allocator_)) {
    other.ptr_ = nullptr;
  }

  buffer& operator=(const buffer& other) {
    assert(ptr_ != other.ptr_);
    if (size_ == other.size_) {
      std::memcpy(ptr_, other.ptr_, size_ * sizeof(word_type));
    }
    else {
      if (ptr_ != nullptr) {
        allocator_.deallocate(ptr_, size_);
      }
      size_ = other.size_;
      allocator_ = other.allocator_;
      ptr_ = allocator_.allocate(size_);
      std::memcpy(ptr_, other.ptr_, size_ * sizeof(word_type));
    }
    return *this;
  }

  buffer& operator=(buffer&& other) noexcept {
    if (this != &other) {
      if (ptr_ != nullptr) {
        allocator_.deallocate(ptr_, size_);
      }
      size_ = other.size_;
      allocator_ = std::move(other.allocator_);
      ptr_ = other.ptr_;
      other.ptr_ = nullptr;
    }
    return *this;
  }

  ~buffer() noexcept {
    if (ptr_ != nullptr) {
      allocator_.deallocate(ptr_, size_);
    }
  }

  inline std::size_t
  size() const noexcept {
    return size_;
  }

  inline word_type*
  data() noexcept {
    return ptr_;
  }

  inline const word_type*
  data() const noexcept {
    return ptr_;
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
