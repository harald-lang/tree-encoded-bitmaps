#include "gtest/gtest.h"

#include "boost/dynamic_bitset.hpp"
#include "dynamic_tree_mask_lo.hpp"
#include <bitset>
#include "dynamic_tree_mask_lo.hpp"
#include <dtl/dtl.hpp>

//===----------------------------------------------------------------------===//
// Tests for iterators on TEB_LO.
//===----------------------------------------------------------------------===//

// Instantiate the type templates to a fixed length.
constexpr std::size_t LEN = 8;
using iter_simple = dtl::dynamic_tree_mask_lo::iter_and_simple;

// Fixture for the parameterized test case.
template<typename T>
class iter_test : public ::testing::Test {};

// Specify the types for which we want to run the API tests.
using types_under_test = ::testing::Types<
    iter_simple
>;
TYPED_TEST_CASE(iter_test, types_under_test);

//===----------------------------------------------------------------------===//
// Bitwise AND.
TYPED_TEST(iter_test, bitwise_and_correct) {
  using T = TypeParam;

  boost::dynamic_bitset<$u32> bs_1{LEN,0};
  boost::dynamic_bitset<$u32> bs_2{LEN,0};

  for (std::size_t a = 0; a < (1u << LEN); a++) {

    std::bitset<LEN> bs_a(a);
    for(auto i = 0; i < LEN; i++)
      bs_1[i] = bs_a[i];

    dtl::dynamic_tree_mask_lo tm_1(bs_1);

    std::cout << std::endl << "bs_a: " << bs_a << std::endl;

    for (std::size_t b = 0; b < (1u << LEN); b++) {
      std::bitset<LEN> bs_b(b);

      std::cout << "   bs_b: " << bs_b << std::endl;

      for(auto i = 0; i < LEN; i++)
        bs_2[i] = bs_b[i];

      dtl::dynamic_tree_mask_lo tm_2(bs_2);

      std::bitset<LEN> and_check(bs_a & bs_b); // used to check correctness
      std::bitset<LEN> and_it; // used for the result of the iterator

      dtl::dynamic_tree_mask_lo::iter_and_simple it_s(tm_1, tm_2);

      while(!it_s.end()){
        dtl::dynamic_tree_mask_lo::range r = it_s.matches();

        for (auto i = r.first; i < r.last; i++)
          and_it.set(i);

        if(!it_s.end())
          it_s.next();
      }

      ASSERT_EQ(and_check, and_it) << "Iterator was not correct! " << std::endl
                                   << "bs_a: " << bs_a << std::endl
                                   << "bs_b: " << bs_b << std::endl
                                   << "res : " << (bs_a & bs_b) << std::endl
                                   << "it_s: " << and_it << std::endl;

    }
  }
}
//===----------------------------------------------------------------------===//