#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/dtl.hpp>
#include <dtl/bits.hpp>
#include <immintrin.h>
//===----------------------------------------------------------------------===//
namespace dtl {
//===----------------------------------------------------------------------===//
namespace bit_buffer_internal {
template<u64 n, u64 p>
struct every_nth_bit_set_helper {
  static constexpr u64 value = ((p % n) == 0 ? 1ul : 0ul) << p
      | every_nth_bit_set_helper<n, p - 1>::value;
};

template<u64 n>
struct every_nth_bit_set_helper<n, 0ull> {
  static constexpr u64 value = 1ul;
};

template<u64 n>
struct every_nth_bit_set {
  static constexpr u64 value = every_nth_bit_set_helper<n, 63>::value;
};
} // namespace bit_buffer_internal
//===----------------------------------------------------------------------===//
template<u64 _slot_bitlength = 8>
class bit_buffer {
  static constexpr u64 buffer_bitlength = 64;
  static constexpr u64 slot_bitlength = _slot_bitlength;
  static constexpr u64 slot_mask = (1ul << slot_bitlength) - 1;
  static constexpr u64 slot_cnt = buffer_bitlength / slot_bitlength ;
  static constexpr u64 initial_read_mask =
      bit_buffer_internal::every_nth_bit_set<slot_bitlength>::value;

  $u64 buf_;
  $u64 read_mask_;

public:

  bit_buffer() : buf_(0), read_mask_(initial_read_mask) {}
  explicit bit_buffer(u64 val) : buf_(val), read_mask_(initial_read_mask) {}
  bit_buffer(const bit_buffer& other) = default;
  bit_buffer(bit_buffer&& other) noexcept = default;
  bit_buffer& operator=(const bit_buffer& other) = default;
  bit_buffer& operator=(bit_buffer&& other) = default;
  ~bit_buffer() = default;

  inline u64
  get_raw() const noexcept {
    return buf_;
  }

  inline void
  set_raw(u64 val) noexcept {
    buf_ = val;
  }

  inline u64
  get(u64 slot_idx) const {
    return (buf_ >> (slot_idx * slot_bitlength)) & slot_mask;
  }

  inline void
  set(u64 slot_idx, u64 value) noexcept {
    buf_ ^= (get(slot_idx) ^ value) << (slot_idx * slot_bitlength);
  }

  inline void
  broadcast(u64 value) noexcept {
    buf_ = initial_read_mask * value;
  }

  inline void
  reset_read_mask() noexcept {
    read_mask_ = initial_read_mask;
  }

  inline u64
  get_read_mask() const noexcept {
    return read_mask_;
  }

  inline u64
  get_read_pos(u64 slot_idx) const noexcept {
    u64 t = (read_mask_ >> (slot_idx * slot_bitlength)) & slot_mask;
    return dtl::bits::tz_count(t);
  }

  inline void
  increment(u64 mask) noexcept {
    u64 slot_mask =
        _pdep_u64(mask, initial_read_mask) * ((1ul << slot_bitlength) - 1);
    read_mask_ =
        ((read_mask_ & slot_mask) << 1) | (read_mask_ & (~slot_mask));
  }

  inline u64
  read(u64 mask) const noexcept {
    return _pext_u64(buf_, read_mask_) & mask;
  }

  inline u64
  read() const noexcept {
    return _pext_u64(buf_, read_mask_);
  }

};
//===----------------------------------------------------------------------===//
} // namespace dtl