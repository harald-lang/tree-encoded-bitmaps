#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/bitmap.hpp>
#include <dtl/bitmap/bitwise_operations.hpp>
#include <dtl/bitmap/dynamic_bitmap.hpp>
#include <dtl/dtl.hpp>

#include <cstddef>
#include <vector>
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
//===----------------------------------------------------------------------===//
/// Used to queue up pending updates.
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
//===----------------------------------------------------------------------===//
/// Computes the necessary updates to turn bitmap A into bitmap B.
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
