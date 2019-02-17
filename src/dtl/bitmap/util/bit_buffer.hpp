#pragma once

#include <dtl/dtl.hpp>
#include <immintrin.h>

namespace dtl {
//===----------------------------------------------------------------------===//
class bit_buffer {
  static constexpr u64 buffer_bitlength = 64;
  static constexpr u64 slot_bitlength = 8;
  static constexpr u64 slot_mask = (1ul << slot_bitlength) - 1;
  static constexpr u64 slot_cnt = buffer_bitlength / slot_bitlength ;

  $u64 buf_;
  $u64 read_mask_;

public:

  bit_buffer() : buf_(0), read_mask_(0x0101010101010101) {}
  explicit bit_buffer(u64 val) : buf_(val), read_mask_(0x0101010101010101) {}
  bit_buffer(const bit_buffer& other) = default;
  bit_buffer(bit_buffer&& other) = default;
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
    u64 t = ((value << slot_bitlength) | value)
        | (((value << slot_bitlength) | value) << (slot_bitlength * 2));
    buf_ = t | (t << (slot_bitlength * 4));
  }

  inline void
  reset_read_mask() noexcept {
    read_mask_ = 0x0101010101010101;
  }

  inline u64
  get_read_mask() noexcept {
    return read_mask_;
  }

  inline void
  increment(u64 mask) noexcept {
    u64 slot_mask =
        _pdep_u64(mask, 0x0101010101010101) * ((1ul << slot_bitlength) - 1);
    read_mask_ =
        ((read_mask_ & slot_mask) << 1) | (read_mask_ & (~slot_mask));
  }

  inline u64
  read(u64 mask) noexcept {
    return _pext_u64(buf_, read_mask_) & mask;
  }

  inline u64
  read() noexcept {
    return _pext_u64(buf_, read_mask_);
  }

};
//===----------------------------------------------------------------------===//
} // namespace dtl