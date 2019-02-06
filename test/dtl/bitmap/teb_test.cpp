#include "gtest/gtest.h"

#include <chrono>

#include <dtl/dtl.hpp>
#include <dtl/bitmap.hpp>
#include <dtl/bitmap/dynamic_bitmap.hpp>
#include <dtl/bitmap/dynamic_roaring_bitmap.hpp>
#include <dtl/bitmap/dynamic_wah.hpp>
#include <dtl/bitmap/partitioned_position_list.hpp>
#include <dtl/bitmap/partitioned_range_list.hpp>
#include <dtl/bitmap/position_list.hpp>
#include <dtl/bitmap/range_list.hpp>
#include <dtl/bitmap/teb.hpp>
#include <dtl/bitmap/util/random.hpp>

#include <dtl/bitmap/util/markov_process.hpp>


u64
now_nanos() {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(
      std::chrono::steady_clock::now().time_since_epoch()).count();
}

//===----------------------------------------------------------------------===//
// Test the compression levels of TEBs.
TEST(teb,
     ensure_that_higher_compression_do_not_increase_the_size) {
  for (auto n = 1ull << 10; n <= 1ull << 10; n <<= 1) {
    for (auto f = 1; f <= 16; f *= 2) {
      for (auto d : {0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5} ) {
        for (std::size_t rep = 0; rep < 1; ++rep) {
          dtl::bitmap bs = dtl::gen_random_bitmap_markov(n, f, d);
          std::cout << "plain: " << bs << std::endl;
          dtl::teb<0> teb_o0(bs);
          std::cout << "o0:\n" << teb_o0 << std::endl;
          dtl::teb<1> teb_o1(bs);
          std::cout << "o1:\n" << teb_o1 << std::endl;
          dtl::teb<2> teb_o2(bs);
          std::cout << "o2:\n" << teb_o2 << std::endl;

          ASSERT_LE(teb_o1.size_in_byte(), teb_o0.size_in_byte())
              << "\n"
              << "plain:  " << bs
              << "\n"
              << "TEB-o0: " << teb_o0
              << "\n"
              << "TEB-o1: " << teb_o1
              << std::endl;
          ASSERT_LE(teb_o2.size_in_byte(), teb_o1.size_in_byte())
              << "\n"
              << "plain:  " << bs
              << "\n"
              << "TEB-o1: " << teb_o1
              << "\n"
              << "TEB-o2: " << teb_o2
              << std::endl;
        }
      }
    }
  }
}
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
TEST(teb,
     ensure_the_correct_number_of_perfect_levels) {
  using teb = dtl::teb<2>;
  //                      implicit_inner_node_cnt, expected # of perfect levels
  ASSERT_EQ(teb::determine_perfect_tree_levels( 0), 1);

  ASSERT_EQ(teb::determine_perfect_tree_levels( 1), 2);
  ASSERT_EQ(teb::determine_perfect_tree_levels( 2), 2);

  ASSERT_EQ(teb::determine_perfect_tree_levels( 3), 3);
  ASSERT_EQ(teb::determine_perfect_tree_levels( 4), 3);
  ASSERT_EQ(teb::determine_perfect_tree_levels( 5), 3);
  ASSERT_EQ(teb::determine_perfect_tree_levels( 6), 3);

  ASSERT_EQ(teb::determine_perfect_tree_levels( 7), 4);
  ASSERT_EQ(teb::determine_perfect_tree_levels( 8), 4);
  ASSERT_EQ(teb::determine_perfect_tree_levels( 9), 4);
  ASSERT_EQ(teb::determine_perfect_tree_levels(10), 4);
  ASSERT_EQ(teb::determine_perfect_tree_levels(11), 4);
  ASSERT_EQ(teb::determine_perfect_tree_levels(12), 4);
  ASSERT_EQ(teb::determine_perfect_tree_levels(13), 4);
  ASSERT_EQ(teb::determine_perfect_tree_levels(14), 4);

  ASSERT_EQ(teb::determine_perfect_tree_levels(15), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(16), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(17), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(18), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(19), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(20), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(21), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(22), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(23), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(24), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(25), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(26), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(27), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(28), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(29), 5);
  ASSERT_EQ(teb::determine_perfect_tree_levels(30), 5);

  ASSERT_EQ(teb::determine_perfect_tree_levels(31), 6);
  ASSERT_EQ(teb::determine_perfect_tree_levels(32), 6);
}
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
TEST(teb,
     ensure_the_correct_top_node_idxs) {
  using teb = dtl::teb<2>;

  auto top_node_idx_end =  [&](u64 perfect_levels) {
    return (1ull << perfect_levels) - 1;
  };
  auto top_node_idx_begin = [&](u64 perfect_levels) {
    return (1ll << (perfect_levels - 1)) - 1;
  };

  ASSERT_EQ(top_node_idx_begin( 1), 0);
  ASSERT_EQ(top_node_idx_end  ( 1), 1);

  ASSERT_EQ(top_node_idx_begin( 2), 1);
  ASSERT_EQ(top_node_idx_end  ( 2), 3);

  ASSERT_EQ(top_node_idx_begin( 3), 3);
  ASSERT_EQ(top_node_idx_end  ( 3), 7);

}
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
// Compare TEB and BitmapTree
TEST(teb,
     compare_explicit_node_cnt_of_teb_o1_and_bitmaptree) {
  for (auto n = 1ull << 10; n <= 1ull << 10; n <<= 1) {
    for (auto f = 1; f <= 16; f *= 2) {
      for (auto d : {0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5} ) {
        for (std::size_t rep = 0; rep < 1; ++rep) {
          dtl::bitmap bs = dtl::gen_random_bitmap_markov(n, f, d);
          dtl::teb<1> teb(bs);
          dtl::bitmap_tree<1> bt(bs);
          ASSERT_LE(bt.estimate_encoded_size_in_bytes(), teb.size())
              << " - bitmap=" << bs << std::endl;
          std::cout << std::endl;
        }
      }
    }
  }
}

TEST(teb,
     compare_size_of_teb_o2_and_bitmaptree) {
  for (auto n = 1ull << 10; n <= 1ull << 12; n <<= 1) {
    for (auto f = 1; f <= 16; f *= 2) {
      for (auto d : {0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5} ) {
        for (std::size_t rep = 0; rep < 1; ++rep) {
          std::cout << "n=" << n << std::flush;
          dtl::bitmap bs = dtl::gen_random_bitmap_markov(n, f, d);
          const auto teb_nanos_start = now_nanos();
          dtl::teb<2> teb(bs);
          std::cout << teb.info() << std::endl;
          const auto teb_nanos_end = now_nanos();

//          std::cout
//              << "[teb] total: " << (teb.implicit_inner_node_cnt_
//                + teb.structure_.size()+ teb.implicit_leaf_node_cnt_)
//              << ", leading: " << teb.implicit_inner_node_cnt_
//              << ", explicit: " << teb.structure_.size()
//              << ", trailing: " << teb.implicit_leaf_node_cnt_
//              << ", labels: " << teb.labels_.size()
//              << std::endl;
          const auto bt_nanos_start = now_nanos();
          dtl::bitmap_tree<2> bt(bs);
          const auto bt_nanos_end = now_nanos();
          ASSERT_LE(bt.estimate_encoded_size_in_bytes(), teb.size())
              << " - bitmap=" << bs << std::endl;
          std::cout
              << ", teb_millis=" << (teb_nanos_end - teb_nanos_start) / 1000 / 1000
              << ", bt_millis=" << (bt_nanos_end - bt_nanos_start) / 1000 / 1000
              << ", teb/bt=" << (teb_nanos_end - teb_nanos_start) / (bt_nanos_end - bt_nanos_start);
          std::cout << std::endl;
        }
      }
    }
  }
}
//===----------------------------------------------------------------------===//
TEST(teb,
     fooooo) {
//  for (auto n = 1ull << 10; n <= 1ull << 10; n <<= 1) {
//  for (auto n = 1ull << 8; n <= 1ull << 8; n <<= 1) {
  for (auto n = 1ull << 20; n <= 1ull << 20; n <<= 1) {
    for (auto f : {2, 4, 8} ) {
      std::cout << "--- f=" << f << std::endl;
      for (std::size_t rep = 0; rep < 1; ++rep) {
        const auto k = 3;
        std::cout << "f=" << f << ", ";
        std::vector<std::size_t> cnt(k, f);
        std::vector<$u32> vals(n, 0);
        markov_process mp(k, f);
        for (std::size_t i = 0; i < n; ++i) {
          const auto val = mp.next();
          cnt[val]++;
          vals[i] = val;
//          std::cout << val << ",";
        }
        std::cout << std::endl;
        std::cout << "card: ";
        for (std::size_t i = 0; i < k; ++i) {
          std::cout << i << ":" << cnt[i] << ",";
        }
        std::cout << std::endl;
        std::cout << "f: ";
        for (std::size_t i = 0; i < k; ++i) {
          std::cout << i << ":"
              << dtl::determine_clustering_factor(vals, i) << ", ";
        }
        std::cout << std::endl;

      }
    }
  }
}
