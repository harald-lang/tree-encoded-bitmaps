#include "gtest/gtest.h"

#include <dtl/bitmap/teb_types.hpp>
#include <dtl/bitmap/util/bitmap_fun.hpp>
#include <dtl/bitmap/util/bitmap_seq_reader.hpp>
#include <dtl/bitmap/util/buffer.hpp>
#include <dtl/bitmap/util/random.hpp>
#include <dtl/dtl.hpp>
#include <dtl/iterator.hpp>

#include <chrono>
#include <random>
//===----------------------------------------------------------------------===//
TEST(bitmap_seq_reader,
    no_bits_set) {
  using word_type = dtl::teb_word_type;
  using fn = dtl::bitmap_fun<word_type>;
  const auto word_bitlength = sizeof(word_type) * 8;
  const auto word_cnt = 10;
  const auto len = word_bitlength * word_cnt;

  // Repeat the test several times.
  for (std::size_t rep = 0; rep < 10; ++rep) {
    // Empty bitmap.
    dtl::buffer<word_type> bitmap(word_cnt);

    dtl::bitmap_seq_reader<word_type> reader(bitmap.data(), 0);
    for (std::size_t i = 0; i < len; ++i) {
      u1 actual_val = reader.next();
      u1 expected_val = fn::test(bitmap.data(), i);
      ASSERT_EQ(expected_val, actual_val);
    }
  }
}
//===----------------------------------------------------------------------===//
TEST(bitmap_seq_reader,
    all_bits_set) {
  using word_type = dtl::teb_word_type;
  using fn = dtl::bitmap_fun<word_type>;
  const auto word_bitlength = sizeof(word_type) * 8;
  const auto word_cnt = 10;
  const auto len = word_bitlength * word_cnt;

  // Repeat the test several times.
  for (std::size_t rep = 0; rep < 10; ++rep) {
    // Populate a bitmap.
    dtl::buffer<word_type> bitmap(word_cnt);
    fn::set(bitmap.data(), 0, len);

    dtl::bitmap_seq_reader<word_type> reader(bitmap.data(), 0);
    for (std::size_t i = 0; i < len; ++i) {
      u1 actual_val = reader.next();
      u1 expected_val = fn::test(bitmap.data(), i);
      ASSERT_EQ(expected_val, actual_val);
    }
  }
}
//===----------------------------------------------------------------------===//
TEST(bitmap_seq_reader,
    random_test) {
  using word_type = dtl::teb_word_type;
  using fn = dtl::bitmap_fun<word_type>;
  const auto word_bitlength = sizeof(word_type) * 8;
  const auto word_cnt = 10;
  const auto len = word_bitlength * word_cnt;

  // Repeat the test several times.
  for (std::size_t rep = 0; rep < 10; ++rep) {
    // Populate a bitmap.
    dtl::buffer<word_type> bitmap(word_cnt);
    two_state_markov_process mp(2, 0.1);
    for (std::size_t i = 0; i < len; ++i) {
      if (mp.next()) fn::set(bitmap.data(), i);
    }

    dtl::bitmap_seq_reader<word_type> reader(bitmap.data(), 0);
    for (std::size_t i = 0; i < len; ++i) {
      u1 actual_val = reader.next();
      u1 expected_val = fn::test(bitmap.data(), i);
      ASSERT_EQ(expected_val, actual_val);
    }
  }
}
//===----------------------------------------------------------------------===//
TEST(bitmap_seq_reader,
    random_test_with_read_offset) {
  using word_type = dtl::teb_word_type;
  using fn = dtl::bitmap_fun<word_type>;
  const auto word_bitlength = sizeof(word_type) * 8;
  const auto word_cnt = 10;
  const auto len = word_bitlength * word_cnt;

  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<> dis(1, len - 1);

  // Repeat the test several times.
  for (std::size_t rep = 0; rep < 10; ++rep) {
    // Populate a bitmap.
    dtl::buffer<word_type> bitmap(word_cnt);
    two_state_markov_process mp(2, 0.1);
    for (std::size_t i = 0; i < len; ++i) {
      if (mp.next()) fn::set(bitmap.data(), i);
    }

    const auto start_idx = dis(gen);
    dtl::bitmap_seq_reader<word_type> reader(bitmap.data(), start_idx);
    for (std::size_t i = start_idx; i < len; ++i) {
      u1 actual_val = reader.next();
      u1 expected_val = fn::test(bitmap.data(), i);
      ASSERT_EQ(expected_val, actual_val);
    }
  }
}
//===----------------------------------------------------------------------===//
TEST(bitmap_seq_reader,
    limit) {
  using word_type = dtl::teb_word_type;
  using fn = dtl::bitmap_fun<word_type>;
  const auto word_bitlength = sizeof(word_type) * 8;
  const auto word_cnt = 10;
  const auto len = word_bitlength * word_cnt;

  // Repeat the test several times.
  for (std::size_t rep = 0; rep < word_bitlength; ++rep) {
    // Populate a bitmap.
    dtl::buffer<word_type> bitmap(word_cnt);
    fn::set(bitmap.data(), 0, len);

    const auto start_idx = rep;
    const auto end_idx = len - rep;
    dtl::bitmap_limit_seq_reader<word_type> reader(bitmap.data(),
        start_idx, end_idx);
    for (std::size_t i = start_idx; i < len; ++i) {
      u1 actual_val = reader.next();
      u1 expected_val = i < end_idx
          ? fn::test(bitmap.data(), i)
          : false;
      ASSERT_EQ(expected_val, actual_val);
    }
  }
}
//===----------------------------------------------------------------------===//
