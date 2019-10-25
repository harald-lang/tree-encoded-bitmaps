#include "common.hpp"
#include "experiments/util/bitmap_db.hpp"
#include "experiments/util/bitmap_types.hpp"
#include "experiments/util/differential_bitmap_types.hpp"
#include "experiments/util/gen.hpp"
#include "experiments/util/params.hpp"
#include "experiments/util/prep_data.hpp"

#include <dtl/dtl.hpp>

#include <dtl/bitmap/bitwise_operations.hpp>
#include <dtl/bitmap/diff/diff.hpp>
#include <dtl/bitmap/diff/merge.hpp>
#include <dtl/bitmap/diff/merge_teb.hpp>
#include <iostream>
#include <random>
#include <vector>
//===----------------------------------------------------------------------===//
// Experiment: Measures the memory consumption of compressed bitmaps with
//             differential data structures. In this experiment we perform
//             several updates and track the size of the diff structure.
//             When it exceeds a given threshold, e.g. 1% of the uncompressed
//             bitmap size, we perform a merge operation.
//===----------------------------------------------------------------------===//
/// Used to queue up pending updates.
struct update_entry {
  $u32 pos = -1;
  $u1 value = false;
  update_entry(u32 pos, u1 value) : pos(pos), value(value) {}
  bool operator<(const update_entry& o) const { return pos < o.pos; }
  void
  print(std::ostream& os) const noexcept {
    os << "update_entry[pos=" << pos
        << ",value=" << (value ? "true" : "false")
        << "]";
  }
};
struct range_update_entry {
  $u32 pos = -1;
  $u32 length = 0;
  $u1 value = false;
  bool operator<(const range_update_entry& o) const { return pos < o.pos; }

  void
  print(std::ostream& os) const noexcept {
    os << "range_update_entry[pos=" << pos
       << ",length=" << length
       << ",value=" << (value ? "true" : "false")
       << "]";
  }
};
struct task {
  diff_bitmap_t bitmap_type; // under test
  $u64 bitmap_id_a;
  $u64 bitmap_id_b;
  $f64 diff_limit_frac;
};
//===----------------------------------------------------------------------===//
std::vector<range_update_entry>
prepare_range_updates(dtl::bitmap bm_a, dtl::bitmap bm_b) {
  //  std::cout << "a: " << bm_a << std::endl;
  //  std::cout << "b: " << bm_b << std::endl;
  //  std::cout << "d: " << (bm_a ^ bm_b) << std::endl;
  dtl::dynamic_bitmap<$u32> a(bm_a);
  dtl::dynamic_bitmap<$u32> b(bm_b);
  auto xor_it = dtl::bitwise_xor_it(a.scan_it(), b.scan_it());
  std::vector<range_update_entry> updates;
  while (!xor_it.end()) {
    {
      // Start a new update range.
      updates.emplace_back();
      auto& update_entry = updates.back();
      update_entry.pos = xor_it.pos();
      update_entry.length = 1;
      update_entry.value = !a.test(xor_it.pos());
    }
    for (std::size_t i = 1; i < xor_it.length(); ++i) {
      auto val = !a.test(xor_it.pos() + i);
      if (val == updates.back().value) {
        // Extend the current update range.
        ++updates.back().length;
      }
      else {
        // Start a new update range.
        updates.emplace_back();
        auto& update_entry = updates.back();
        update_entry.pos = xor_it.pos() + i;
        update_entry.length = 1;
        update_entry.value = !a.test(xor_it.pos() + i);
      }
    }
    xor_it.next();
  }

  // Validation
  auto a_up = bm_a;
  for (std::size_t u = 0; u < updates.size(); ++u) {
    //    std::cout << std::setw(4) << u << ": " << updates[u] << std::endl;
    auto& entry = updates[u];
    for (std::size_t i = entry.pos; i < (entry.pos + entry.length); ++i) {
      a_up[i] = entry.value;
    }
  }
  //  std::cout << "r: " << a_up << std::endl;
  if (a_up != bm_b) {
    std::cerr << "Validation failed." << std::endl;
    std::exit(1);
  }

  return updates;
}
//===----------------------------------------------------------------------===//
template<typename T>
void do_measurement(task t, std::ostream& os) {
  // The bitmap that is used as starting point.
  const auto bm_a = db.load_bitmap(t.bitmap_id_a);
  const auto bm_a_density = dtl::determine_bit_density(bm_a);
  const auto bm_a_clustering_factor = dtl::determine_clustering_factor(bm_a);
  // The second bitmap determines the updates to perform.
  const auto bm_b = db.load_bitmap(t.bitmap_id_b);
  const auto bm_b_density = dtl::determine_bit_density(bm_b);
  const auto bm_b_clustering_factor = dtl::determine_clustering_factor(bm_b);

  const auto n = bm_a.size();
  assert(bm_a.size() == bm_b.size());

  // Prepare the updates. The updates are clustered, so that the characteristics
  // of the bitmap remains roughly the same during modifications.
  auto range_updates = prepare_range_updates(bm_a, bm_b);

  // Shuffle the update order.
  std::mt19937 gen(42); // for reproducible results
  std::shuffle(range_updates.begin(), range_updates.end(), gen);

  std::vector<update_entry> updates;
  updates.reserve(bm_b.count());
  for (auto& range : range_updates) {
    for (std::size_t i = range.pos; i < range.pos + range.length; ++i) {
      updates.emplace_back(i, range.value);
    }
  }

  const std::size_t total_update_cnt = updates.size();

  // Determine the size limit for the differential data structure. Used to
  // trigger merges.
  const auto diff_limit_frac = t.diff_limit_frac;
  const std::size_t diff_size_limit_bytes =
      static_cast<std::size_t>((bm_a.size() / 8) * diff_limit_frac);

  // The encoded bitmap.
  T b(bm_a);
  // A plain bitmap used for validation.
  auto bm_expected = bm_a;

  using merge_type = dtl::merge_naive<typename T::bitmap_type, typename T::diff_type>;

  // Perform the updates.
  std::size_t performed_updates_cnt = 0;
  std::size_t pending_updates_cnt = 0;
  std::size_t merge_cnt = 0;

  auto print = [&]() {
    std::string type_info = b.info();
    boost::replace_all(type_info, "\"", "\"\""); // Escape JSON for CSV output.

    std::cout << RUN_ID
        << ",\"" << BUILD_ID << "\""
        << "," << "\"" << b.name() << "\""
        << "," << t.bitmap_id_a
        << "," << t.bitmap_id_b
        << "," << diff_limit_frac
        << "," << performed_updates_cnt // the number of performed updates
        << "," << total_update_cnt // the total number of updates
        << "," << merge_cnt
        << "," << pending_updates_cnt // the number of updates not yet merged
        << "," << diff_size_limit_bytes
        << "," << b.diff_size_in_byte()
        << "," << b.size_in_byte()
        << "," << dtl::determine_bit_density(bm_expected)
        << "," << dtl::determine_clustering_factor(bm_expected)
        << "," << "\"" << type_info << "\""
        << std::endl;
  };

  print();
  for (std::size_t i = 0; i < total_update_cnt; ++i) {
    auto& u = updates[i];
    b.set(u.pos, u.value);
    bm_expected.set(u.pos, u.value);
    ++pending_updates_cnt;
    ++performed_updates_cnt;

    const auto diff_size = b.diff_size_in_byte();
    if (diff_size >= diff_size_limit_bytes) {
      b.template merge<merge_type>();
      pending_updates_cnt = 0;
      ++merge_cnt;
    }

    if (performed_updates_cnt % 100 == 0 || pending_updates_cnt == 0) {
      print();
    }
  }
}
//===----------------------------------------------------------------------===//
// Type switch.
#define __GEN_CASE(NAME)                                             \
  case diff_bitmap_t::NAME:                                          \
    do_measurement<diff_type_of<diff_bitmap_t::NAME>::type>(std::forward<_Params>(params)...); \
    break;
template<typename... _Params>
void do_measurement(diff_bitmap_t bitmap_type, _Params&&... params) {
  switch (bitmap_type) {
    __GEN_CASE(roaring_roaring)
    __GEN_CASE(roaring_wah)
    __GEN_CASE(teb_roaring)
    __GEN_CASE(teb_wah)
    __GEN_CASE(wah_roaring)
    __GEN_CASE(wah_wah)
  }
}
#undef __GEN_CASE
//===----------------------------------------------------------------------===//
$i32 main() {
  // Prepare benchmark settings.
  u64 n_min = 1ull << 20;
  u64 n_max = 1ull << 20;

  // Use the same setting as with the compression skyline.
  std::cerr << "run_id=" << RUN_ID << std::endl;
  std::cerr << "build_id=" << BUILD_ID << std::endl;

  std::vector<$f64> clustering_factors;
  clustering_factors.push_back(8.0);

  std::vector<$f64> bit_densities;
  bit_densities.push_back(0.1);

  std::vector<$u64> n_values;
  for ($u64 n = n_min; n <= n_max; n <<= 1) {
    n_values.push_back(n);
  }

  if (GEN_DATA) {
    std::vector<params_markov> params;
    for (auto f : clustering_factors) {
      for (auto d : bit_densities) {
        for (auto n : n_values) {
          if (!markov_parameters_are_valid(n, f, d)) continue;
          params_markov p;
          p.n = n;
          p.clustering_factor = f;
          p.density = d;
          params.push_back(p);
        }
      }
    }
    prep_data(params, RUNS, db);
    std::exit(0);
  }
  else {
    if (db.empty()) {
      std::cerr << "Bitmap database is empty. Use GEN_DATA=1 to populate the "
                   "database."
                << std::endl;
      std::exit(1);
    }
  }

  // The implementations under test.
  std::vector<diff_bitmap_t> bitmap_types;
  bitmap_types.push_back(diff_bitmap_t::roaring_wah);
  bitmap_types.push_back(diff_bitmap_t::roaring_roaring);
  bitmap_types.push_back(diff_bitmap_t::teb_roaring);
  bitmap_types.push_back(diff_bitmap_t::teb_wah);
  bitmap_types.push_back(diff_bitmap_t::wah_roaring);
  bitmap_types.push_back(diff_bitmap_t::wah_wah);

  std::vector<task> tasks;
  for (auto f : clustering_factors) {
    for (auto d : bit_densities) {
      for (auto n : n_values) {
        if (!markov_parameters_are_valid(n, f, d)) continue;

        auto bitmap_ids = db.find_bitmaps(n, f, d);
        if (bitmap_ids.empty()) {
          continue;
        }
        if (bitmap_ids.size() < RUNS) {
          std::cerr << "There are only " << bitmap_ids.size() << " prepared "
                    << "bitmaps for the parameters n=" << n << ", f=" << f
                    << ", d=" << d << ", but " << RUNS << " are required."
                    << std::endl;
        }
        if (bitmap_ids.size() >= 2) {
          task t;
          t.bitmap_id_a = bitmap_ids[0];
          t.bitmap_id_b = bitmap_ids[1];
          t.bitmap_type = diff_bitmap_t::teb_wah;
          for (auto limit : { 0.01, 0.05, 0.10 }) {
            t.diff_limit_frac = limit;
            for (auto b : bitmap_types) {
              t.bitmap_type = b;
              tasks.push_back(t);
            }
          }
        }
      }
    }
  }

  // Run the actual benchmark.
  std::function<void(const task&, std::ostream&)> fn =
      [](const task t, std::ostream& os) -> void {
        do_measurement(t.bitmap_type, t, os);
      };
  dispatch(tasks, fn);
}
//===----------------------------------------------------------------------===//
