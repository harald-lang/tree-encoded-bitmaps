#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/dtl.hpp>

#include <cassert>
#include <cstddef>
#include <cstdint>
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
  static constexpr uintptr_t cache_line_size = 64; // bytes
  static constexpr auto elements_per_cache_line =
      cache_line_size / sizeof(word_type);

  /// The number of elements.
  std::size_t size_;
  /// Pointer to the buffer.
  word_type* ptr_;
  /// Cache line aligned pointer to the buffer.
  word_type* aligned_ptr_;
  /// The allocator.
  _alloc allocator_;

  static inline word_type*
  get_aligned_ptr(word_type* ptr) {
    const uintptr_t adjust_by =
        (reinterpret_cast<uintptr_t>(ptr) % cache_line_size) == 0
        ? 0
        : cache_line_size - (reinterpret_cast<uintptr_t>(ptr) % cache_line_size);
    auto p = reinterpret_cast<word_type*>(
        reinterpret_cast<uintptr_t>(ptr) + adjust_by);
    return reinterpret_cast<word_type*>(
        __builtin_assume_aligned(p, cache_line_size));
  }

  inline void
  do_allocate() {
    ptr_ = allocator_.allocate(size_ + elements_per_cache_line);
    aligned_ptr_ = get_aligned_ptr(ptr_);
  }

  inline void
  do_deallocate() {
    if (ptr_ != nullptr) {
      allocator_.deallocate(ptr_, size_ + elements_per_cache_line);
    }
  }

public:
  /// C'tor
  explicit buffer(std::size_t size, u1 init = true, const _alloc& alloc = _alloc())
      : size_(size), ptr_(nullptr), allocator_(alloc) {
    do_allocate();
    if (init) {
      std::memset(aligned_ptr_, 0, size_ * sizeof(word_type));
    }
  }

  buffer(const buffer& other)
      : size_(other.size_), ptr_(nullptr), aligned_ptr_(nullptr),
        allocator_(other.allocator_) {
    do_allocate();
    std::memcpy(aligned_ptr_, other.aligned_ptr_, size_ * sizeof(word_type));
  }

  buffer(buffer&& other) noexcept
      : size_(other.size_), ptr_(other.ptr_), aligned_ptr_(other.aligned_ptr_),
        allocator_(std::move(other.allocator_)) {
    other.ptr_ = nullptr;
    other.aligned_ptr_ = nullptr;
  }

  buffer& operator=(const buffer& other) {
    assert(ptr_ != other.ptr_);
    if (size_ == other.size_) {
      std::memcpy(aligned_ptr_, other.aligned_ptr_, size_ * sizeof(word_type));
    }
    else {
      do_deallocate();
      size_ = other.size_;
      allocator_ = other.allocator_;
      do_allocate();
      std::memcpy(aligned_ptr_, other.aligned_ptr_, size_ * sizeof(word_type));
    }
    return *this;
  }

  buffer& operator=(buffer&& other) noexcept {
    if (this != &other) {
      do_deallocate();
      size_ = other.size_;
      allocator_ = std::move(other.allocator_);
      ptr_ = other.ptr_;
      aligned_ptr_ = other.aligned_ptr_;
      other.ptr_ = nullptr;
      other.aligned_ptr_ = nullptr;
    }
    return *this;
  }

  ~buffer() noexcept {
    do_deallocate();
  }

  inline std::size_t
  size() const noexcept {
    return size_;
  }

  inline word_type*
  data() noexcept {
    return reinterpret_cast<word_type*>(
        __builtin_assume_aligned(aligned_ptr_, cache_line_size));
  }

  inline const word_type*
  data() const noexcept {
    return reinterpret_cast<const word_type*>(
        __builtin_assume_aligned(aligned_ptr_, cache_line_size));
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
