#include "gtest/gtest.h"

#include <dtl/dtl.hpp>
#include <dtl/bitmap.hpp>
#include <experiments/compression/common.hpp>

TEST(common,
     encode_decode_pregenerated_bitmaps_test) {

  // The implementations under test.
  std::vector<bitmap_t> bitmap_types;
  for (auto bitmap_type = static_cast<int>(bitmap_t::_first);
       bitmap_type <= static_cast<int>(bitmap_t::_last);
       ++bitmap_type) {
    bitmap_types.push_back(static_cast<bitmap_t>(bitmap_type));
  }

  const auto bitmap_ids = db.ids();
  for (auto id : bitmap_ids) {
    std::cout << "id=" << id << std::endl;
    for (auto t : bitmap_types) {
      config c;
      c.bitmap_id = id;
      c.bitmap_type = t;
      run(c, std::cout);
    }
  }
}
