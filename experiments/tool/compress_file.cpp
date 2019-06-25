#include <fstream>
#include <vector>

#include <dtl/dtl.hpp>
#include <dtl/env.hpp>
#include <dtl/bitmap/util/convert.hpp>
#include <dtl/bitmap/teb.hpp>
#include <dtl/bitmap/dynamic_bitmap.hpp>
#include <dtl/bitmap/dynamic_roaring_bitmap.hpp>
#include <dtl/bitmap/dynamic_wah.hpp>

#include <boost/algorithm/string.hpp>
#include <dtl/bitmap/util/random.hpp>

#include "../util/bitmap_db.hpp"
#include "../util/threading.hpp"

//===----------------------------------------------------------------------===//
$i32 main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
  }
  std::ifstream input(argv[1], std::ios::binary);

  // Copy the file content into a buffer.
  std::vector<unsigned char> buffer(std::istreambuf_iterator<char>(input), {});

  std::cout << "Read " << buffer.size() << " bytes." << std::endl;

  while (! dtl::is_power_of_two(buffer.size())) {
    buffer.push_back(char(0));
  }

  dtl::data_view<u32> dv {
      reinterpret_cast<u32*>(&buffer[0]),
      reinterpret_cast<u32*>(&buffer[buffer.size() / 4])
  };

  dtl::bitmap_view<u32> bv(dv);


  const std::size_t n = buffer.size() * 8;
  dtl::bitmap bm;
  for (std::size_t i = 0; i < n; ++i) {
    bm.push_back(bv[i]);
  }
  std::cout << "d=" << (bm.count() * 1.0 / bm.size())
      << ", f=" << dtl::determine_clustering_factor(bm)
      << std::endl;
  dtl::teb<> teb(bm);
  dtl::dynamic_roaring_bitmap roaring(bm);
  std::cout << "TEB size =     " << teb.size_in_byte() << std::endl;
  std::cout << "TEB info =     " << teb.info() << std::endl;
  std::cout << "Roaring size = " << roaring.size_in_byte() << std::endl;
}
//===----------------------------------------------------------------------===//
