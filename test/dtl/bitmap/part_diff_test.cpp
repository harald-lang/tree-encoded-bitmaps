#include "api_types.hpp"
#include "gtest/gtest.h"

#include <dtl/bitmap.hpp>
#include <dtl/bitmap/diff/diff.hpp>
#include <dtl/bitmap/diff/merge.hpp>
#include <dtl/bitmap/diff/merge_teb.hpp>
#include <dtl/bitmap/part/part_upforward.hpp>
#include <dtl/bitmap/util/random.hpp>
//===----------------------------------------------------------------------===//
// Tests for the differential update feature with partitioning.
//===----------------------------------------------------------------------===//
template<
    typename BitmapType,
    typename DiffType,
    template<typename B, typename D> typename MergeType>
struct TestSetting {
  using type = dtl::part_upforward<dtl::diff<BitmapType, DiffType>, 1ull << 8>;
  using merge_type = MergeType<BitmapType, DiffType>;
};
//===----------------------------------------------------------------------===//
// Partitioned bitmaps, where each partition has a diff.
using part_diff_types_under_test = ::testing::Types<
    // Naive merge.
    TestSetting<teb_v2, roaring_bitmap, dtl::merge_naive>,
    TestSetting<teb_v2, wah, dtl::merge_naive>,
    TestSetting<wah, roaring_bitmap, dtl::merge_naive>,
    TestSetting<wah, wah, dtl::merge_naive>,
    // Naive merge (using the run iterator).
    TestSetting<teb_v2, roaring_bitmap, dtl::merge_naive_iter>,
    TestSetting<teb_v2, wah, dtl::merge_naive_iter>,
    TestSetting<wah, roaring_bitmap, dtl::merge_naive_iter>,
    TestSetting<wah, wah, dtl::merge_naive_iter>,
    // In-place merge.
    TestSetting<wah, wah, dtl::merge_inplace>,
    // Tree-based merge (TEB only).
    TestSetting<teb_v2, roaring_bitmap, dtl::merge_tree>,
    TestSetting<teb_v2, wah, dtl::merge_tree>
    >;
//===----------------------------------------------------------------------===//
// Fixture for the parameterized test case.
template<typename T>
class part_diff_test : public ::testing::Test {};
TYPED_TEST_CASE(part_diff_test, part_diff_types_under_test);
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
TYPED_TEST(part_diff_test, insert_on_empty_bitmap) {
  using T = typename TypeParam::type;
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
TYPED_TEST(part_diff_test, insert_and_delete) {
  using T = typename TypeParam::type;
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
// Apply and merge changes.
TYPED_TEST(part_diff_test, merge) {
  using T = typename TypeParam::type;
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

    // Merge the pending updates.
    using merge_type = typename TypeParam::merge_type;
    t.template merge<merge_type>();

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
/// Randomized test for differential updates (using varying bitmap sizes).
TYPED_TEST(part_diff_test, random_updates) {
  using T = typename TypeParam::type;

  for (auto len = 1024; len <= 8192 * 2; len *= 2) {
    std::cout << "testing size " << len << std::endl;
    for (std::size_t rep = 0; rep < 10; ++rep) {
      const auto bm_initial = dtl::gen_random_bitmap_markov(len, 8.0, 0.1);

      T enc(bm_initial);
      auto plain = bm_initial;

      // Perform multiple batch updates.
      for (std::size_t i = 0; i < 5; ++i) {
        // prepare the updates
        const auto update = dtl::gen_random_bitmap_uniform(len, 0.01);
        plain |= update; // for comparison
        // Enqueue the updates.
        auto pos = update.find_first();
        while (pos != boost::dynamic_bitset<$u32>::npos) {
          enc.set(pos, true);
          pos = update.find_next(pos);
        }
        // Merge the pending updates.
        using merge_type = typename TypeParam::merge_type;
        enc.template merge<merge_type>();

        dtl::bitmap dec = dtl::to_bitmap_using_iterator(enc);
        ASSERT_EQ(plain, dec)
            << "Failed to apply updates:\n"
            << "              :'" << update << "'\n"
            << "i=" << i << std::endl;
      }
    }
  }
}
//===----------------------------------------------------------------------===//
