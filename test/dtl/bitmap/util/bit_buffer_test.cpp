#include "gtest/gtest.h"

#include <dtl/dtl.hpp>
#include <dtl/bitmap/util/bit_buffer.hpp>
//===----------------------------------------------------------------------===//
TEST(bit_buffer,
     null_increment) {
  dtl::bit_buffer<> bb;
  auto a = bb.get_read_mask();
  bb.increment(0b00000000);
  auto b = bb.get_read_mask();
  ASSERT_EQ(a, b);
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer,
     increment_all) {
  dtl::bit_buffer<> bb;
  auto a = bb.get_read_mask();
  bb.increment(~0ul);
  auto b = bb.get_read_mask();
  ASSERT_EQ(a, b >> 1);
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer,
     masked_increment) {
  dtl::bit_buffer<> bb;
  bb.increment((1ul << 0) | (1ul << 3) | (1ul << 7));
  auto a = bb.get_read_mask();
  ASSERT_EQ(a, 0x201010102010102);
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer,
     broadcast) {
  dtl::bit_buffer<16> bb;
  bb.broadcast(42);
  for (std::size_t i = 0; i < 16; ++i) {
    ASSERT_EQ(bb.get(i), 42) << " i = " << i;
  }
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer,
     read) {
  dtl::bit_buffer<> bb0(0ul);
  ASSERT_EQ(0, bb0.read());
  dtl::bit_buffer<> bb1(~0ul);
  ASSERT_EQ(0xFF, bb1.read());
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer,
     masked_read) {
  dtl::bit_buffer<> bb(~0ul);
  ASSERT_EQ(0b00100101, bb.read(0b00100101));

  bb.set_raw(0x0201020102010201);
  ASSERT_EQ(0b01010101, bb.read(0b11111111));
  ASSERT_EQ(0b00000000, bb.read(0b10101010));
  bb.increment(0b11111111);
  ASSERT_EQ(0b10101010, bb.read(0b10101010));
}
//===----------------------------------------------------------------------===//
TEST(bit_buffer,
     get_read_position) {
  dtl::bit_buffer<> bb;
  bb.increment((1ul << 0) | (1ul << 3) | (1ul << 7));
  bb.increment((1ul << 7));
  ASSERT_EQ(bb.get_read_pos(0), 1);
  ASSERT_EQ(bb.get_read_pos(1), 0);
  ASSERT_EQ(bb.get_read_pos(2), 0);
  ASSERT_EQ(bb.get_read_pos(3), 1);
  ASSERT_EQ(bb.get_read_pos(4), 0);
  ASSERT_EQ(bb.get_read_pos(5), 0);
  ASSERT_EQ(bb.get_read_pos(6), 0);
  ASSERT_EQ(bb.get_read_pos(7), 2);
}
//===----------------------------------------------------------------------===//
