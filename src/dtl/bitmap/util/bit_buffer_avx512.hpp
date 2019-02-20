#pragma once

#include <dtl/dtl.hpp>
#include <dtl/bits.hpp>
#include <dtl/simd.hpp>
#include <immintrin.h>

#ifdef __AVX512BW__
namespace dtl {
//===----------------------------------------------------------------------===//
namespace bit_buffer_avx512_internal {
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
} // namespace bit_buffer_avx512_internal
//===----------------------------------------------------------------------===//
template<u64 _slot_bitlength = 16 /* FIXME 16 is the only allowed value for now */>
class bit_buffer_avx512 {
  using buf_t = __m512i;

  static constexpr u64 buffer_bitlength = 512;
  static constexpr u64 slot_bitlength = _slot_bitlength;
  static constexpr u64 slot_mask = (1ul << slot_bitlength) - 1;
  static constexpr u64 slot_cnt = buffer_bitlength / slot_bitlength ;
//  static constexpr __m512i initial_read_mask = (r512 {
//      .u16 = {1u,1u,1u,1u,1u,1u,1u,1u,1u,1u,1u,1u,1u,1u,1u,1u,1u,1u,1u,1u}
//  }).i;

  const __m512i initial_read_mask_;
  __m512i buf_;
  __m512i read_mask_;

public:

  bit_buffer_avx512() :
      initial_read_mask_(_mm512_set1_epi16(1)),
      buf_(_mm512_setzero_si512()),
      read_mask_(initial_read_mask_) {}
  explicit bit_buffer_avx512(__m512i val) :
      initial_read_mask_(_mm512_set1_epi16(1)),
      buf_(val),
      read_mask_(initial_read_mask_) {}
  bit_buffer_avx512(const bit_buffer_avx512& other) = default;
  bit_buffer_avx512(bit_buffer_avx512&& other) noexcept = default;
  bit_buffer_avx512& operator=(const bit_buffer_avx512& other) = default;
  bit_buffer_avx512& operator=(bit_buffer_avx512&& other) = default;
  ~bit_buffer_avx512() = default;

  inline __m512i
  get_raw() const noexcept {
    return buf_;
  }

  inline void
  set_raw(const __m512i val) noexcept {
    buf_ = val;
  }

  inline u64
  get(u64 slot_idx) const {
    r512 copy {.i = buf_};
    return copy.u16[slot_idx];
  }

  inline void
  set(u64 slot_idx, u64 value) noexcept {
    r512 copy {.i = buf_};
    copy.u16[slot_idx] = value;
    buf_ = copy.i;
  }

  inline void
  broadcast(u64 value) noexcept {
//    u16 v = static_cast<u16>(value);
//    r512 t {.u16 = {v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v}};
//    buf_ = t.i;
    buf_ = _mm512_set1_epi16(value);
  }

  inline void
  reset_read_mask() noexcept {
    read_mask_ = initial_read_mask_;
  }

  inline __m512i
  get_read_mask() const noexcept {
    return read_mask_;
  }

  inline u64
  get_read_pos(u64 slot_idx) const noexcept {
    r512 copy {.i = read_mask_};
    return dtl::bits::tz_count(copy.u16[slot_idx]);
  }

  inline void
  increment(u64 mask) noexcept {
    read_mask_ = _mm512_mask_slli_epi16(read_mask_, mask, read_mask_, 1u);
  }

  inline u64
  read(u64 mask) const noexcept {
    return _mm512_mask_test_epi16_mask(mask, buf_, read_mask_);
  }

  inline u64
  read() const noexcept {
    return _mm512_test_epi16_mask(buf_, read_mask_);
  }

};
//===----------------------------------------------------------------------===//
} // namespace dtl
#endif // __AVX512BW__
