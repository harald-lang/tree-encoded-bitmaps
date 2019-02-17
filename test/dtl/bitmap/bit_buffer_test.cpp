#include "gtest/gtest.h"

#include <dtl/dtl.hpp>
#include <dtl/bitmap/util/bit_buffer.hpp>

//===----------------------------------------------------------------------===//
TEST(bit_buffer,
     null_increment) {
  dtl::bit_buffer bb;
  auto a = bb.get_read_mask();
  bb.increment(0b00000000);
  auto b = bb.get_read_mask();
  ASSERT_EQ(a, b);
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer,
     increment_all) {
  dtl::bit_buffer bb;
  auto a = bb.get_read_mask();
  bb.increment(~0ul);
  auto b = bb.get_read_mask();
  ASSERT_EQ(a, b >> 1);
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer,
     masked_increment) {
  dtl::bit_buffer bb;
  bb.increment((1ul << 0) | (1ul << 3) | (1ul << 7));
  auto a = bb.get_read_mask();
  ASSERT_EQ(a, 0x201010102010102);
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer,
     broadcast) {
  dtl::bit_buffer bb;
  bb.broadcast(42);
  ASSERT_EQ(bb.get(0), 42);
  ASSERT_EQ(bb.get(1), 42);
  ASSERT_EQ(bb.get(2), 42);
  ASSERT_EQ(bb.get(3), 42);
  ASSERT_EQ(bb.get(4), 42);
  ASSERT_EQ(bb.get(5), 42);
  ASSERT_EQ(bb.get(6), 42);
  ASSERT_EQ(bb.get(7), 42);
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer,
     read) {
  dtl::bit_buffer bb0(0ul);
  ASSERT_EQ(0, bb0.read());
  dtl::bit_buffer bb1(~0ul);
  ASSERT_EQ(0xFF, bb1.read());
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer,
     masked_read) {
  dtl::bit_buffer bb(~0ul);
  ASSERT_EQ(0b00100101, bb.read(0b00100101));

  bb.set_raw(0x0201020102010201);
  ASSERT_EQ(0b01010101, bb.read(0b11111111));
  ASSERT_EQ(0b00000000, bb.read(0b10101010));
  bb.increment(0b11111111);
  ASSERT_EQ(0b10101010, bb.read(0b10101010));
}
//===----------------------------------------------------------------------===//
