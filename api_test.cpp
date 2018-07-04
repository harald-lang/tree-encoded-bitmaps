#include "gtest/gtest.h"

#include <bitset>

#include "roaring_bitmap.hpp"
#include "tree_mask_lo.hpp"
#include "tree_mask_po.hpp"
#include "wah.hpp"

//===----------------------------------------------------------------------===//
// Typed API tests for loss-less compressed bitmaps.
//===----------------------------------------------------------------------===//

// Instantiate the type templates to a fixed length.
constexpr std::size_t LEN = 16;
using roaring_bitmap = dtl::roaring_bitmap<LEN>;
using tree_mask_lo = dtl::tree_mask_lo<LEN>;
using tree_mask_po = dtl::tree_mask_po<LEN>;
using wah32 = dtl::wah32<LEN>;
using wah64 = dtl::wah64<LEN>;

// Fixture for the parameterized test case.
template<typename T>
class api_test : public ::testing::Test {};

// Specify the types for which we want to run the API tests.
using types_under_test = ::testing::Types<
    roaring_bitmap,
    tree_mask_lo,
    tree_mask_po,
    wah32
//    wah64 // FIXME: Position iterator does not work with 64-bit impl.
>;
TYPED_TEST_CASE(api_test, types_under_test);


//===----------------------------------------------------------------------===//
// Encode, decode, compare the results
TYPED_TEST(api_test, encode_decode) {
  using T = TypeParam;

  for (auto i = 0; i < (1u << LEN); ++i) {
    std::bitset<LEN> bs(i);
    T t(bs);
    std::bitset<LEN> dec = t.to_bitset();
    ASSERT_EQ(bs, dec) << "Decoding failed for i=" << i << ". - '" << bs << "' -> '" << t << "' -> '" << dec << "'"
                       << std::endl;
  }
}
//===----------------------------------------------------------------------===//

