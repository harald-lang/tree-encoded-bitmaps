#ifdef __AVX512BW__
#include "gtest/gtest.h"

#include <dtl/dtl.hpp>
#include <dtl/simd.hpp>
#include <dtl/bitmap/util/bit_buffer_avx512.hpp>

u1 _mm512_equals(__m512i a, __m512i b) {
  return _mm512_cmpneq_epi64_mask(a, b) == 0;
}

void print_u16(__m512i v) {
  dtl::r512 c { .i = v };
  std::cout << "[" << c.u16[0];
  for (std::size_t i = 1; i < 32; ++i) {
    std::cout << "," << c.u16[i];
  }
  std::cout << "]" << std::endl;
}

//===----------------------------------------------------------------------===//
TEST(bit_buffer_avx512,
     null_increment) {
  dtl::bit_buffer_avx512<> bb;
  auto a = bb.get_read_mask();
  bb.increment(0);
  auto b = bb.get_read_mask();
  ASSERT_TRUE(_mm512_equals(a, b));
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer_avx512,
     increment_all) {
  dtl::bit_buffer_avx512<> bb;
  bb.increment(~0ul);
  auto b = bb.get_read_mask();
  ASSERT_TRUE(_mm512_equals(b, _mm512_set1_epi16(2)));
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer_avx512,
     masked_increment) {
  dtl::bit_buffer_avx512<> bb;
  bb.increment((1ul << 0) | (1ul << 3) | (1ul << 7));
  auto a = bb.get_read_mask();
//  dtl::r512 b = {.u16 = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,2,1,1,2}};
  dtl::r512 b = {.u16 =   {2,1,1,2,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}};
  std::cout << "exp:"; print_u16(b.i);
  std::cout << "got:"; print_u16(a);
  ASSERT_TRUE(_mm512_equals(a, b.i));
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer_avx512,
     broadcast) {
  dtl::bit_buffer_avx512<> bb;
  bb.broadcast(42);
  for (std::size_t i = 0; i < 32; ++i) {
    ASSERT_EQ(bb.get(i), 42) << " i = " << i;
  }
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer_avx512,
     read) {
  dtl::bit_buffer_avx512<> bb0(_mm512_setzero_si512());
  ASSERT_EQ(0, bb0.read());
  dtl::bit_buffer_avx512<> bb1(_mm512_set1_epi32(-1));
  ASSERT_EQ(0xFFFFFFFF, bb1.read());
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer_avx512,
     masked_read) {
  dtl::bit_buffer_avx512<> bb(_mm512_set1_epi32(-1));
  ASSERT_EQ(0b00100101, bb.read(0b00100101));

  dtl::r512 b = {.u16 = {
      1,2,1,2,
      1,2,1,2,
      1,1,1,1,
      1,1,1,1,
      1,1,1,1,
      1,1,1,1,
      1,1,1,1,
      1,1,1,1}};
  bb.set_raw(b.i);
  ASSERT_EQ(0b01010101, bb.read(0b11111111));
  ASSERT_EQ(0b00000000, bb.read(0b10101010));
  bb.increment(0b11111111);
  ASSERT_EQ(0b10101010, bb.read(0b10101010));
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer_avx512,
     get_read_position) {
  dtl::bit_buffer_avx512<> bb;
  bb.increment((1ul << 0) | (1ul << 3) | (1ul << 7));
  bb.increment((1ul << 7));
  bb.increment((1ul << 20));
  bb.increment((1ul << 20));
  bb.increment((1ul << 20));
  ASSERT_EQ(bb.get_read_pos(0), 1);
  ASSERT_EQ(bb.get_read_pos(1), 0);
  ASSERT_EQ(bb.get_read_pos(2), 0);
  ASSERT_EQ(bb.get_read_pos(3), 1);
  ASSERT_EQ(bb.get_read_pos(4), 0);
  ASSERT_EQ(bb.get_read_pos(5), 0);
  ASSERT_EQ(bb.get_read_pos(6), 0);
  ASSERT_EQ(bb.get_read_pos(7), 2);
  ASSERT_EQ(bb.get_read_pos(20), 3);
}
//===----------------------------------------------------------------------===//
#endif // __AVX512BW__
