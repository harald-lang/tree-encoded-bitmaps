#include "gtest/gtest.h"

#include <dtl/dtl.hpp>
#include <dtl/bitmap.hpp>
#include <dtl/bitmap/util/random.hpp>

//===----------------------------------------------------------------------===//
// Test the different compression levels of TEBs.
//===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
// Test the compression levels of TEBs.
TEST(teb_compression_levels_test,
     ensure_that_higher_compression_do_not_increase_the_size) {
  for (auto n = 1ull << 10; n <= 1ull << 10; n <<= 1) {
    for (auto f = 1; f <= 16; f *= 2) {
      for (auto d : {0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5} ) {
        for (std::size_t rep = 0; rep < 1; ++rep) {
          dtl::bitmap bs = dtl::gen_random_bitmap(n, f, d);
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

          // Reconstruct the original bitmap.
          dtl::bitmap dec0 = to_bitmap_using_iterator(teb_o0);
//          dtl::bitmap dec0 = teb_o0.to_bitset();
          ASSERT_EQ(bs, dec0)
                        << "Decoding failed: "
                        << "'" << bs << "' -> '" << teb_o0
                        << "' -> '" << dec0 << "'"
                        << std::endl;
//          dtl::bitmap dec1 = to_bitmap_using_iterator(teb_o1);
          dtl::bitmap dec1 = teb_o1.to_bitset();
          ASSERT_EQ(bs, dec1)
                        << "Decoding failed: "
                        << "'" << bs << "' -> '" << teb_o1
                        << "' -> '" << dec1 << "'"
                        << std::endl;

//          dtl::bitmap dec2 = to_bitmap_using_iterator(teb_o2);
          dtl::bitmap dec2 = teb_o2.to_bitset();
          ASSERT_EQ(bs, dec2)
                        << "Decoding failed: "
                        << "'" << bs << "' -> '" << teb_o2
                        << "' -> '" << dec2 << "'"
                        << std::endl;
//          std::cout << std::endl;
        }
      }
    }
  }
}
//===----------------------------------------------------------------------===//

