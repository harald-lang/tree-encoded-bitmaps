#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/bits.hpp>
#include <dtl/dtl.hpp>
#include <dtl/simd.hpp>

#ifdef __AVX2__
#include <immintrin.h>
#endif
//===----------------------------------------------------------------------===//
#ifdef __AVX2__
namespace dtl {
//===----------------------------------------------------------------------===//
namespace bit_buffer_avx2_internal {
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
} // namespace bit_buffer_avx2_internal
//===----------------------------------------------------------------------===//
template<u64 _slot_bitlength = 8 /* FIXME 8 is the only allowed value for now */>
class bit_buffer_avx2 {
  using buf_t = __m256i;

  static constexpr u64 buffer_bitlength = 256;
  static constexpr u64 slot_bitlength = _slot_bitlength;
  static constexpr u64 slot_mask = (1ul << slot_bitlength) - 1;
  static constexpr u64 slot_cnt = buffer_bitlength / slot_bitlength;

  const __m256i initial_read_mask_;
  __m256i buf_;
  __m256i read_mask_;

  // Taken from: https://stackoverflow.com/questions/21622212/how-to-perform-the-inverse-of-mm256-movemask-epi8-vpmovmskb
  inline __m256i
  get_mask3(const uint32_t mask) const {
    __m256i vmask(_mm256_set1_epi32(mask));
    const __m256i shuffle(_mm256_setr_epi64x(0x0000000000000000,
        0x0101010101010101, 0x0202020202020202, 0x0303030303030303));
    vmask = _mm256_shuffle_epi8(vmask, shuffle);
    const __m256i bit_mask(_mm256_set1_epi64x(0x7fbfdfeff7fbfdfe));
    vmask = _mm256_or_si256(vmask, bit_mask);
    return _mm256_cmpeq_epi8(vmask, _mm256_set1_epi64x(-1));
  }

public:
  bit_buffer_avx2()
      : initial_read_mask_(_mm256_set1_epi8(1)),
        buf_(_mm256_setzero_si256()),
        read_mask_(initial_read_mask_) {}
  explicit bit_buffer_avx2(__m256i val)
      : initial_read_mask_(_mm256_set1_epi8(1)),
        buf_(val),
        read_mask_(initial_read_mask_) {}
  bit_buffer_avx2(const bit_buffer_avx2& other) = default;
  bit_buffer_avx2(bit_buffer_avx2&& other) noexcept = default;
  bit_buffer_avx2& operator=(const bit_buffer_avx2& other) = default;
  bit_buffer_avx2& operator=(bit_buffer_avx2&& other) noexcept = default;
  ~bit_buffer_avx2() = default;

  inline __m256i
  get_raw() const noexcept {
    return buf_;
  }

  inline void
  set_raw(const __m256i val) noexcept {
    buf_ = val;
  }

  inline u64
  get(u64 slot_idx) const {
    r256 copy { .i = buf_ };
    return copy.u8[slot_idx];
  }

  inline void
  set(u64 slot_idx, u64 value) noexcept {
    r256 copy { .i = buf_ };
    copy.u8[slot_idx] = static_cast<u8>(value);
    buf_ = copy.i;
  }

  inline void
  broadcast(u64 value) noexcept {
    buf_ = _mm256_set1_epi8(static_cast<u8>(value));
  }

  inline void
  reset_read_mask() noexcept {
    read_mask_ = initial_read_mask_;
  }

  inline __m256i
  get_read_mask() const noexcept {
    return read_mask_;
  }

  inline u64
  get_read_pos(u64 slot_idx) const noexcept {
    r256 copy { .i = read_mask_ };
    return dtl::bits::tz_count(copy.u8[slot_idx]);
  }

  inline void
  increment(u32 mask) noexcept {
    const auto avx2_mask = get_mask3(mask);
    const auto inc = _mm256_slli_epi64(read_mask_, 1);
    read_mask_ = _mm256_blendv_epi8(read_mask_, inc, avx2_mask);
  }

  inline u32
  read(u32 mask) const noexcept {
    return _mm256_movemask_epi8(
               _mm256_cmpgt_epi8(_mm256_and_si256(buf_, read_mask_),
                   _mm256_setzero_si256()))
        & mask;
  }

  inline u32
  read() const noexcept {
    return static_cast<u32>(
        _mm256_movemask_epi8(
            _mm256_cmpgt_epi8(_mm256_and_si256(buf_, read_mask_),
                _mm256_setzero_si256())));
  }

  inline u32
  read_ahead() const noexcept {
    return static_cast<u32>(_mm256_movemask_epi8(
        _mm256_cmpgt_epi8(_mm256_and_si256(buf_,
                              _mm256_slli_epi64(read_mask_, 1u)),
            _mm256_setzero_si256())));
  }
};
//===----------------------------------------------------------------------===//
} // namespace dtl
#endif // __AVX2__
