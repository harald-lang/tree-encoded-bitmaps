#include <atomic>
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
$i32 main() {
  const std::string DB_FILE =
      dtl::env<std::string>::get("DB_FILE", "./random_bitmaps.sqlite3");
  bitmap_db db(DB_FILE);

  const auto bitmap_id = dtl::env<$i64>::get("BITMAP_ID", -1);
  if (bitmap_id == -1) {
    std::cerr << "No bitmap id specified." << std::endl;
  }

  const auto implicit_nodes = dtl::env<$i64>::get("IMPLICIT", 0) == 1;

  const auto bm = db.load_bitmap(bitmap_id);
  std::cout << "% bitmap_id=" << bitmap_id
      << ", d=" << dtl::determine_bit_density(bm)
      << ", f=" << dtl::determine_clustering_factor(bm)
      << std::endl;

  const auto n = bm.size();
  const auto h = dtl::log_2(n) + 1;
  std::vector<$u32> histo(h, 0);
  dtl::bitmap_tree<3> bt(bm);
  std::size_t node_cntr = 0;

  if (!implicit_nodes) {
    // Count explicit tree nodes.
    for (auto it = bt.breadth_first_begin(); it != bt.breadth_first_end(); ++it) {

      u64 idx = (*it).idx;
      u64 level = (*it).level;

      if (idx < bt.get_first_explicit_node_idx()) {
        continue;
      }

      if (idx > bt.get_last_explicit_node_idx()) {
        continue;
      }

      assert(level < histo.size());
      node_cntr++;
      histo[level]++;
    }
  }
  else {
    // Count implicit tree nodes.
    for (auto it = bt.breadth_first_begin(); it != bt.breadth_first_end(); ++it) {

      u64 idx = (*it).idx;
      u64 level = (*it).level;

      if (idx >= bt.get_first_explicit_node_idx()
          && idx <= bt.get_last_explicit_node_idx()) {
        continue;
      }

      assert(level < histo.size());
      node_cntr++;
      histo[level]++;
    }
  }

  for (std::size_t i = 0; i < histo.size(); ++i) {
    std::cout << i << " " << histo[i] << std::endl;
  }
  std::cout << "% total number of nodes: " << node_cntr << std::endl;
}
//===----------------------------------------------------------------------===//
