#include "gtest/gtest.h"

#include <dtl/dtl.hpp>
#include <dtl/bitmap.hpp>
#include <dtl/iterator.hpp>
#include <experiments/compression/common.hpp>

TEST(single,
     test) {

//  const auto ids = db.ids();
//  for (auto id : ids) {
//    using T = dtl::teb_scan<2>;
//    auto plain = db.load_bitmap(id);
//    if (plain.size() > 1ull << 10) continue;
//    T enc(plain);
//    auto dec = dtl::to_bitmap_using_iterator(enc);
//    if (plain != dec) {
//      std::cout << "id=" << id << ", n=" << plain.size() << std::endl;
//    }
//  }
// id = 71
  // id = 36905

//  using Ta = dtl::teb<2>;
  using Ta = dtl::teb_scan<1>;

//  using Tb = dtl::teb<2>;
  using Tb = dtl::teb_scan<2>;
  auto plain = db.load_bitmap(36905);
  auto valid_a = plain;
  valid_a.reset();
  auto valid_b = plain;
  valid_b.reset();

  Ta enc_a(plain);
  std::cout << "\na info: " << enc_a.info() << std::endl;
  std::cout << "a dump: " << enc_a << std::endl;
  Tb enc_b(plain);
  std::cout << "b info: " << enc_b.info() << std::endl;
  std::cout << "b dump: " << enc_b << std::endl;
  auto it_a = enc_a.it();
  auto it_b = enc_b.it();
  auto pos_a = it_a.pos();
  auto pos_b = it_b.pos();
  auto len_a = it_a.length();
  auto len_b = it_b.length();
  while (!it_a.end()) {
    pos_a = it_a.pos();
    pos_b = it_b.pos();
    len_a = it_a.length();
    len_b = it_b.length();

    auto end_a = pos_a + len_a;
    auto end_b = pos_b + len_b;

    for (std::size_t i = pos_a; i < end_a; ++i) { valid_a[i] = true; }
    for (std::size_t i = pos_b; i < end_b; ++i) { valid_b[i] = true; }

    std::cout << "a=[" << pos_a << "," << end_a << ")" << std::endl;
    std::cout << "b=[" << pos_b << "," << end_b << ")" << std::endl;
    std::cout << std::endl;
    if (end_a < end_b) {
      it_a.next();
    }
    else if (end_a > end_b) {
      it_b.next();
    }
    else {
      it_a.next();
      it_b.next();
    }
  }

  $u1 failed = false;
  if (plain != valid_a) {
    failed |= true;
    std::cout << "Validation failed for A" << std::endl;
    std::cout << "exp: " << plain << std::endl;
    std::cout << "got: " << valid_a << std::endl;
  }

  if (plain != valid_b) {
    failed |= true;
    std::cout << "Validation failed for B" << std::endl;
    std::cout << "exp: " << plain << std::endl;
    std::cout << "got: " << valid_b << std::endl;
  }

  ASSERT_FALSE(failed);
//  using T = dtl::teb<2>;
//
//
//  auto dec = dtl::to_bitmap_using_iterator(enc);
//  if (plain != dec) {
//    std::cout << "info: " << enc.info() << std::endl;
//    std::cout << "count: exp=" << plain.count()
//              << ", act=" << dec.count()
//              << std::endl;
//    std::cout << "exp: " << plain << std::endl;
//    std::cout << "act: " << dec << std::endl;
//  }
//
//
//  std::cout << "n=" << plain.size() << std::endl;
//  for (std::size_t i = 0; i < plain.size(); ++i) {
//    if (plain[i] != dec[i]) {
//      std::cout << "at pos " << i << ": exp=" << plain[i]
//                << ", act=" << dec[i] << std::endl;
//    }
//  }

////  std::cerr << plain << std::endl;
//  ASSERT_EQ(plain.count(), dec.count());
}


