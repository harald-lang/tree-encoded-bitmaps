#include "api_types.hpp"
#include "gtest/gtest.h"

#include <dtl/bitmap.hpp>
#include <dtl/bitmap/util/random.hpp>
//===----------------------------------------------------------------------===//
// Tests for the differential update feature.
//===----------------------------------------------------------------------===//
using updatable_types_under_test = ::testing::Types<
    part_8_teb,
    part_8_wah,
    part_8_upfwd_wah,
    part_8_upfwd_diff_teb,
    
    part_run_8_teb,
    roaring_bitmap,
    wah,
    position_list_32,
    part_position_list_8,
    part_position_list_16
>;
//===----------------------------------------------------------------------===//
// Fixture for the parameterized test case.
template<typename T>
class update_test : public ::testing::Test {};
TYPED_TEST_CASE(update_test, updatable_types_under_test);
//===----------------------------------------------------------------------===//
/// Sets all bits in 'dst' that are set in 'src'.
template<typename T>
void set(T& dst, const dtl::bitmap& src) {
  auto i = src.find_first();
  while (i != boost::dynamic_bitset<$u32>::npos) {
    dst.set(i, true);
    i = src.find_next(i);
  }
}
//===----------------------------------------------------------------------===//
/// Clears all bits in 'dst' that are set in 'src'.
template<typename T>
void clear(T& dst, const dtl::bitmap& src) {
  auto i = src.find_first();
  while (i != boost::dynamic_bitset<$u32>::npos) {
    dst.set(i, false);
    i = src.find_next(i);
  }
}
//===----------------------------------------------------------------------===//
/// Test inserts on an initially empty bitmap.
TYPED_TEST(update_test, insert_on_empty_bitmap) {
  using T = TypeParam;
  const auto LEN = 8;
  for (auto i = 0; i < (1u << LEN); ++i) {
    // Start with an empty bitmap.
    dtl::bitmap empty(LEN);
    T t(empty);

    // Apply some changes
    dtl::bitmap bs(LEN, i);
    std::cout << bs << " (" << i << ")" << std::endl;
    set(t, bs);

    // Decompress and validate it.
    dtl::bitmap dec = dtl::to_bitmap_using_iterator(t);
    ASSERT_EQ(bs, dec)
        << "Decoding failed for i=" << i
        << ". - '" << bs << "' -> '" << t
        << "' -> '" << dec << "'"
        << std::endl;
  }
}
//===----------------------------------------------------------------------===//
/// Apply changes (inserts and deletes) to a bitmap.
TYPED_TEST(update_test, insert_and_delete) {
  using T = TypeParam;
  const auto LEN = 8;
  for (auto i = 0; i < (1u << LEN); ++i) {
    // Start with an empty bitmap.
    dtl::bitmap bs_initial(LEN, i);
    T t(bs_initial);

    // Toggle all bits.
    dtl::bitmap bs(LEN, i);
    std::cout << bs << " (" << i << ")" << std::endl;
    set(t, ~bs);
    clear(t, bs);

    // Decompress and validate it.
    dtl::bitmap dec = dtl::to_bitmap_using_iterator(t);
    ASSERT_EQ(~bs, dec)
        << "Decoding failed for i=" << i
        << ". - '" << bs << "' -> '" << t
        << "' -> '" << dec << "'"
        << std::endl;
  }
}
//===----------------------------------------------------------------------===//
/// Randomized test for updates (using varying bitmap sizes).
TYPED_TEST(update_test, random_updates) {
  using T = TypeParam;

  for (auto len = 1024; len <= 8192 * 2; len *= 2) {
    std::cout << "testing size " << len << std::endl;
    for (std::size_t rep = 0; rep < 10; ++rep) {
      const auto bm_initial = dtl::gen_random_bitmap_markov(len, 8.0, 0.1);

      T enc(bm_initial);
      auto plain = bm_initial;

      // Prepare the updates.
      const auto update = dtl::gen_random_bitmap_uniform(len, 0.01);
      plain |= update; // for comparison
      // Apply the updates.
      auto pos = update.find_first();
      while (pos != boost::dynamic_bitset<$u32>::npos) {
        enc.set(pos, true);
        pos = update.find_next(pos);
      }

      // Validate
      dtl::bitmap dec = dtl::to_bitmap_using_iterator(enc);
      ASSERT_EQ(plain, dec)
          << "Failed to apply updates:\n"
          << "              :'" << update << "'\n"
          << std::endl;
    }
  }
}
//===----------------------------------------------------------------------===//
