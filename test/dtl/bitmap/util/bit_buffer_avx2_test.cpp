#ifdef __AVX2__
//===----------------------------------------------------------------------===//
#include "gtest/gtest.h"

#include <dtl/dtl.hpp>
#include <dtl/simd.hpp>
#include <dtl/bitmap/util/bit_buffer_avx2.hpp>
//===----------------------------------------------------------------------===//
// Helper functions
u1 _mm256_equals(__m256i a, __m256i b) {
  const auto cmp_res = _mm256_cmpeq_epi8(a, b);
  return (~_mm256_movemask_epi8(cmp_res)) == 0;
}

void print_u8(__m256i v) {
  dtl::r256 c { .i = v };
  std::cout << "[" << static_cast<u32>(c.u8[0]);
  for (std::size_t i = 1; i < 32; ++i) {
    std::cout << "," << static_cast<u32>(c.u8[i]);
  }
  std::cout << "]" << std::endl;
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer_avx2,
     null_increment) {
  dtl::bit_buffer_avx2<> bb;
  auto a = bb.get_read_mask();
  bb.increment(0);
  auto b = bb.get_read_mask();
  print_u8(b);
  ASSERT_TRUE(_mm256_equals(a, b));
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer_avx2,
     increment_all) {
  dtl::bit_buffer_avx2<> bb;
  bb.increment(~0ul);
  auto b = bb.get_read_mask();
  print_u8(b);
  ASSERT_TRUE(_mm256_equals(b, _mm256_set1_epi8(2)));
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer_avx2,
     masked_increment) {
  dtl::bit_buffer_avx2<> bb;
  bb.increment((1ul << 0) | (1ul << 3) | (1ul << 7));
  auto a = bb.get_read_mask();
  dtl::r256 b = {.u8 =   {2,1,1,2,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}};
  std::cout << "exp:"; print_u8(b.i);
  std::cout << "got:"; print_u8(a);
  ASSERT_TRUE(_mm256_equals(a, b.i));
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer_avx2,
     broadcast) {
  dtl::bit_buffer_avx2<> bb;
  bb.broadcast(42);
  for (std::size_t i = 0; i < 32; ++i) {
    ASSERT_EQ(bb.get(i), 42) << " i = " << i;
  }
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer_avx2,
     read) {
  dtl::bit_buffer_avx2<> bb0(_mm256_setzero_si256());
  ASSERT_EQ(0, bb0.read());
  dtl::bit_buffer_avx2<> bb1(_mm256_set1_epi32(-1));
  ASSERT_EQ(0xFFFFFFFF, bb1.read());
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer_avx2,
     read_ahead) {
  dtl::r256 b = {.u8 =   {2,1,1,2,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2}};
  dtl::bit_buffer_avx2<> bb(b.i);
  $u64 ra = 0b10000000000000000000000010001001;
  ASSERT_EQ(ra, bb.read_ahead());
  bb.increment(~u32(0));
  ASSERT_EQ(ra, bb.read());
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer_avx2,
     masked_read) {
  dtl::bit_buffer_avx2<> bb(_mm256_set1_epi32(-1));
  ASSERT_EQ(0b00100101, bb.read(0b00100101));

  dtl::r256 b = {.u8 = {
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
TEST(bit_buffer_avx2,
     get_read_position) {
  dtl::bit_buffer_avx2<> bb;
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
#endif // __AVX2__
