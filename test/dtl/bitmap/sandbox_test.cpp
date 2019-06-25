#include "gtest/gtest.h"

#include <random>

#include <dtl/dtl.hpp>
#include <dtl/bitmap.hpp>
#include <dtl/iterator.hpp>
#include <dtl/bitmap/teb.hpp>
#include <navin/experiments/util/gen.hpp>
#include <dtl/env.hpp>
#include <navin/experiments/util/bitmap_db.hpp>
#include <dtl/bitmap/util/convert.hpp>

#include "thirdparty/perfevent/PerfEvent.hpp"

//TEST(single,
//     DISABLED_test) {
//
//  // Paste area.
////  std::string a(
////      "1000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000001"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000001");
////  std::string b(
////      "1000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000001"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000101");
//
////  std::string a(
////      "0000000000000000000000000000000000001000000000000000000000000000"
////      "0000000000000000000000000000010000000100000000000000000001000100"
////      "0001000001000000000000000000000000000000000000000000000000000000"
////      "0000001000000001000000001000000000001000001000000000100000000000"
////      "0000000010000000000000000000000000000001010000000000000000000000"
////      "0000000000000000000010000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0100000000000010000000000000000000000000100000000000000000000000");
////  std::string b(
////      "0000001000000001000000000000000000000000011011000000000000000000"
////      "0000000000000000000000000100000000000000000010000001000000011000"
////      "0000000000000000000000010010000110000000010010000000000000100000"
////      "0000000000000001000000010010000010000000000000000000000100000100"
////      "0000000000100001000000001000100000000010001000000100000000000000"
////      "0000000011000000000000000001000010000000000000000000000100000010"
////      "0000000001000001100000000000010100000100000000000000000000000000"
////      "0000000000000000000000010000000110100001000000000000000010001010");
//
////  std::string a(
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000001000000000000000000000000000000000000000000000010"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000001");
////  std::string b(
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000001111000000000000000000000000000000000000000000000100"
////      "0000000000100000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000000"
////      "0000000000000000000000000000000000000000000000000000000000000010");
//
////  std::string a("10001010000100010000101101000100101000000110000011001100000000100010100100001000110010000000000010000100001000111000010000111000100011011001000000000001100000001101000000001001000000000010000000000000100000001000000100100100011000000100000001001000010100111000010000010100000001010000000100011000010000000000000010011001000010000000010100110000000100100100100000000000011010011000001010000100110000100000011101000010001000000011010100000000010000110100001011001000001011000001010010000000000001010010101000000010");
////  std::string b("10111011110010101110101111101000011010110111111010111100110110001110111011111011111111001101111110101110011111000010111101110110101111111100111110010101110110000101011111101101100011100010110101011111011010111110101011111010110110000111000111001100110111101101100011100110101111101111110100110010010101110111111111001111100000111111101100011110001110111000000111010011110100111111111101011111100001111010100110010111101100001111010010001001010101111000111111110110100011000101110110111100101110011010100110001101");
////
////  std::cout << "len: " << a.length() << std::endl;
////  assert(a.length() == b.length());
////  dtl::bitmap bm_a(a.length());
////  dtl::bitmap bm_b(b.length());
////  for (std::size_t i = 0; i < a.length(); ++i) {
////    bm_a[a.length() - i - 1] = a[i] == '1';
////    bm_b[a.length() - i - 1] = b[i] == '1';
////  }
////  std::cout << "a:   " << bm_a << std::endl;
////  std::cout << "b:   " << bm_b << std::endl;
////
////
////  auto expected_result = bm_a ^ bm_b;
////
////  using enc_t = dtl::partitioned_position_list<$u32, $u8>;
////
////  using enc_iter_t = decltype(std::declval<enc_t>().it());
////  const std::type_info& ti = typeid(enc_iter_t);
////  std::cout << "iter type: " << ti.name() << std::endl;
////
////  enc_t enc_a(bm_a);
////  enc_t enc_b(bm_b);
////  std::cout << "enc_a: " << enc_a << std::endl;
////  std::cout << "enc_b: " << enc_b << std::endl;
////  auto result = bitwise_xor(enc_a, enc_b);
////
////  std::cout << "exp: " << expected_result << std::endl;
////  std::cout << "act: " << result << std::endl;
////
////  if (expected_result != result) {
////    std::cerr << "Validation failed!" << std::endl;
////  }
////  ASSERT_EQ(result, expected_result);
//
//
//  //===--------------------------------------------------------------------===//
//  auto now_nanos = []() {
//    return std::chrono::duration_cast<std::chrono::nanoseconds>(
//        std::chrono::steady_clock::now().time_since_epoch()).count();
//  };
//  //===--------------------------------------------------------------------===//
//  const std::string DB_FILE =
//      dtl::env<std::string>::get("DB_FILE", "./random_bitmaps.sqlite3");
//  bitmap_db db(DB_FILE);
//
//  const auto n = 1ull << dtl::env<std::size_t>::get("N_LOG2", 22);
//  const auto f = dtl::env<$f64>::get("F", 8);
//  const auto d = dtl::env<$f64>::get("D", 0.01);
//////  const auto d = 0.01;
////  auto d2 = types::Numeric<5,4>::castString("0.01");
////  for (std::size_t i = 0; i < 100; ++i) {
////    const auto t =
////        (d2 * types::Numeric<5, 4>::castString("1.25")).cast_precision<4>();
////    std::cout << "numeric: " << t << std::endl;
////    d2 = t;
////    auto t2 = t.cast_precision<8>();
////    std::cout << t2 << std::endl;
////    std::cout << static_cast<double>(t2) << std::endl;
////  }
//////  const auto d4 = (d3 * types::Numeric<5,4>::castString("1.25")).cast_precision<4>();
////
//////  std::cout << "numeric: " << d2 << std::endl;
//////  std::cout << "numeric: " << d3 << std::endl;
//////  std::cout << "numeric: " << d4 << std::endl;
//
//  const auto bitmap_ids = db.find_bitmaps(n, f, d);
//  if (bitmap_ids.empty()) {
//    const auto bitmap = gen_random_bitmap_markov(n, f, d);
//    db.store_bitmap(n,f,d,bitmap);
//    std::cout << "generated a random bitmap. please restart" << std::endl;
//    std::exit(0);
//  }
//  for (auto id : bitmap_ids) {
//    std::cout << id << " ";
//  }
//  std::cout << std::endl;
//
////  const auto skip_distance = 100;
//
//  const auto bid = bitmap_ids[0];
//  const auto plain_bitmap = db.load_bitmap(bid);
//  std::cout << "loaded bitmap: n=" << plain_bitmap.size()
//      << ", f=" << dtl::determine_clustering_factor(plain_bitmap)
//      << ", d=" << dtl::determine_bit_density(plain_bitmap)
//      << std::endl;
//
////  using T = dtl::dynamic_wah32;
////  using T = dtl::teb<0>;
//  using T = dtl::teb<2>;
////  using T = dtl::dynamic_roaring_bitmap;
//  T enc_bitmap(plain_bitmap);
//
//  // Warm up
//  auto dec_bitmap = plain_bitmap;
////  {
////    dec_bitmap.reset();
////    for(auto it = enc_bitmap.it(); !it.end(); it.next()) {
////      for (std::size_t i = it.pos(); i < (it.pos() + it.length()); ++i) {
////        dec_bitmap[i] = true;
////      }
////    }
////  }
//
//
//  $f64 nanos_per_next_call = 0.0;
////  {
////    auto it = enc_bitmap.it();
////    std::size_t next_cntr = 0;
////    const auto nanos_begin = now_nanos();
////    while (!it.end()) {
////      it.next();
////      next_cntr++;
////    }
////    const auto nanos_end = now_nanos();
////    nanos_per_next_call = (nanos_end - nanos_begin) / (next_cntr * 1.0);
////  }
//
//
////  for (std::size_t skip_distance = 10; skip_distance < n/32; skip_distance += 100){
//  std::cout << enc_bitmap.info() << std::endl;
//  for (std::size_t skip_distance = 1; skip_distance <= 1000000; skip_distance *= 10) {
//    if (skip_distance <= f) continue;
//    std::size_t skip_cntr = 0;
//    PerfEvent e;
//    e.startCounters();
//    const auto nanos_begin = now_nanos();
//    do {
//      auto it = enc_bitmap.it();
//      while (!it.end()) {
//        it.skip_to(it.pos() + skip_distance);
//        skip_cntr++;
//      }
//    } while ((now_nanos() - nanos_begin) < 2000000000);
//    const auto nanos_end = now_nanos();
//    e.stopCounters();
//    std::cout << skip_distance
//        << "," << skip_cntr
//        << "," << nanos_per_next_call
//        << "," << (nanos_end - nanos_begin) / (skip_cntr * 1.0)
//        << std::endl;
//    e.printReport(std::cout, skip_cntr);
//    std::cout << std::endl;
//  }
//
//  if (plain_bitmap != dec_bitmap) {
//    std::cerr << "Encoding/decoding failed." << std::endl;
//    std::exit(1);
//  }
//  std::cout << std::endl;
//  //===--------------------------------------------------------------------===//
//
//}
//
//TEST(single,
//    DISABLED_test_bitmap_tree) {
//  std::string lit = "1111001000001100";
//  dtl::bitmap bm(lit);
//  std::cout << bm << std::endl;
//  constexpr auto opt_level = 0;
//  f64 fpr = 0.2; //0.5;
//  dtl::bitmap_tree<opt_level> bmt(bm, fpr);
////  dtl::bitmap_tree<opt_level> bmt(bm);
//  const auto explicit_node_cnt = (bmt.max_node_cnt()
//      - bmt.get_leading_inner_node_cnt() - bmt.get_trailing_leaf_node_cnt());
//
//
//  std::size_t node_cntr = 0;
//  dtl::bitmap e;
//  dtl::bitmap s;
//  dtl::bitmap l;
//  $u1 prev_was_explicit = false;
//  for (auto it = bmt.breadth_first_begin(); it != bmt.breadth_first_end(); ++it) {
////    ++node_cntr;
//
//    u64 idx = (*it).idx;
//    u64 level = (*it).level;
//    u1 is_explicit = opt_level == 0
//        || ((idx >= bmt.get_first_explicit_node_idx())
//           && (idx <= bmt.get_last_explicit_node_idx()));
//    if (is_explicit != prev_was_explicit) {
//      std::cout << (prev_was_explicit ? "]" : "[");
//      prev_was_explicit = !prev_was_explicit;
//    }
//    std::cout << (bmt.is_inner_node(idx) ? "1" : "0");
//    e.push_back(is_explicit);
//    if (bmt.is_inner_node(idx)) {
//      s.push_back(true);
//    }
//    else {
//      s.push_back(false);
//      l.push_back(bmt.label_of_node(idx));
//    }
//  }
//  if (prev_was_explicit) std::cout << "]";
//  std::cout << std::endl;
//
//  std::cout
//      << "\nbmt.get_leading_inner_node_cnt(): " << bmt.get_leading_inner_node_cnt()
//      << "\nbmt.get_trailing_inner_node_cnt(): " << bmt.get_trailing_leaf_node_cnt()
//      << "\nexplicit node cnt: " << explicit_node_cnt
//      << "\nfirst explicit node idx: " << bmt.get_first_explicit_node_idx()
//      << "\nlast explicit node idx: " << bmt.get_last_explicit_node_idx()
//      << std::endl;
//  std::cout << "e:" << e << std::endl;
//  std::cout << "s:" << s << std::endl;
//  std::cout << "l:" << l << std::endl;
//  std::cout << "size = " << e.count() << " * 1.0625 + " << l.size()
//      << " = " << (e.count() * 1.0625 + l.size()) << std::endl;
////  std::cout << dtl::to_bitmap_using_iterator()
//  std::cout << bmt;
//  std::cout << "      ";
//  for (std::size_t i = 0; i < bm.size(); ++i) {
//    std::cout << (bm[i] ? "1" : "0");
//  }
//  std::cout << std::endl;
//
//  dtl::teb<opt_level> teb(bm, fpr);
//  std::cout << teb << std::endl;
//  const auto dec_bm = dtl::to_bitmap_using_iterator(teb);
//  std::cout << "false positive cnt: " << (dec_bm.count() - bm.count())
//      << std::endl;
//
//  std::vector<$f64> fprs = {
//      0.0,
//      0.0001,
//      0.00025,
//      0.0005,
//      0.00075,
//      0.001,
//      0.0025,
//      0.005,
//      0.0075,
//      0.01,
//      0.025,
//      0.05,
//      0.075,
//      0.1
//  };
//  const auto bitmap_ids = db.ids();
//  for (auto bid : bitmap_ids) {
//    const auto bm = db.load_bitmap(bid);
//    const auto n = bm.size();
//    if (n < 1ull << 20) continue;
//    const auto c = bm.count();
//
//    dtl::teb<opt_level> teb_lossless(bm);
//
//    for (auto fpr : fprs) {
//      dtl::teb<opt_level> teb_lossy(bm, fpr);
//      const auto dec_lossy = dtl::to_bitmap_using_iterator(teb_lossy);
//      const auto fp_cnt = dec_lossy.count() - c;
//      std::cout << bid
//          << "," << fpr
//          << "," << teb_lossless.size_in_byte()
//          << "," << teb_lossy.size_in_byte()
//          << "," << fp_cnt
//          << "," << (teb_lossy.size_in_byte() * 1.0 / teb_lossless.size_in_byte())
//          << std::endl;
//
//      // Validation
//      {
//        if ((bm & dec_lossy) != bm) {
//          std::cerr << "Validation failed. Aborting." << std::endl;
//          std::exit(1);
//        }
//      }
//    }
//  }
//}

////===----------------------------------------------------------------------===//
//extern std::size_t h = 22;
//TEST(single,
//    test) {
//
//  std::vector<std::size_t> src_paths;
//  std::vector<std::size_t> dst_positions;
//  std::random_device rd;
//  std::mt19937 gen(rd());
//  std::uniform_int_distribution<uint32_t> dis(0, 1u << 20);
//  u64 CNT = 10000000;
//  for (std::size_t i = 0; i < CNT; ++i) {
//    src_paths.push_back(dis(gen));
//    dst_positions.push_back(dis(gen));
//  }
//
//  std::size_t sink = 0;
//  PerfEvent e;
//  e.startCounters();
//  const auto tsc_begin = _rdtsc();
//  for (std::size_t i = 0; i < CNT; ++i) {
//    const auto src_path = src_paths[i];
//    const auto dst_pos = src_paths[i];
//
//    const auto dst_path = dst_pos | (1u << h);
//
//    $u64 common_ancestor_path;
//    $u64 common_ancestor_level;
//    dtl::teb<>::determine_common_ancestor_path(
//        src_path, dst_path, common_ancestor_path, common_ancestor_level);
//
//    const auto right_child_of_common_ancestor_path =
//        (common_ancestor_path << 1) | 1ull;
//    const auto right_child_of_common_ancestor_level = common_ancestor_level + 1;
//
//    sink += right_child_of_common_ancestor_path
//        + right_child_of_common_ancestor_level;
//  }
//  const auto tsc_end = _rdtsc();
//  e.stopCounters();
//
//  std::cout << (tsc_end - tsc_begin) * 1.0 / CNT
//      << "," << sink
//      << std::endl;
//  e.printReport(std::cout, CNT);
//
//}
////===----------------------------------------------------------------------===//

//===----------------------------------------------------------------------===//
TEST(single,
    DISABLED_scan_test) {
//
//  const std::string DB_FILE = "./sandbox_bitmaps.sqlite3";
//  bitmap_db db(DB_FILE);
//
//  const auto n_log2 = 16u;
//  u64 n = 1ull << n_log2;
//  const auto f = 7.0;
//  const auto d = 0.47;
//
//  auto scan_levels = 0ull;
//  dtl::bitmap bs;
//  while (scan_levels != 2) {
//    bs = gen_random_bitmap_markov(n, f, d);
//    dtl::teb<> enc_bs(bs);
//    const auto u =
//        enc_bs.determine_perfect_tree_levels(enc_bs.implicit_inner_node_cnt_);
//    scan_levels = n_log2 - u;
//  }

//  dtl::bitmap bs(8, 0b10101001);
//  dtl::bitmap bs(8, 0b00111100);
//  dtl::bitmap bs(8, 0b11001100);
//  dtl::bitmap bs(8, 0b11110001);
//  dtl::bitmap bs(8, 0b00000001);
//  dtl::bitmap bs(8, 0b00011111);
//  dtl::bitmap bs(8, 0b00000001);

//  dtl::bitmap bs;
//  bs.append(0b00111100001111000011110000111100);
//  bs.append(0b00111100001111000011110000111100);
//  bs.append(0b00111100001111000011110000111100);
//  bs.append(0b00111100001111000011110000111100);
//  bs.append(0b10101001101010011010100110101001);
//  bs.append(0b10101001101010011010100110101001);
//  bs.append(0b10101001101010011010100110101001);
//  bs.append(0b10101001101010011010100110101011);

  dtl::bitmap bs;
  bs.append(0b00011111000111110001111100011111);
  bs.append(0b00011111000111110001111100011111);
  bs.append(0b00011111000111110001111100011111);
  bs.append(0b00011111000111110001111100011111);
  bs.append(0b00011111000111110001111100011111);
  bs.append(0b00011111000111110001111100011111);
  bs.append(0b00011111000111110001111100011111);
  bs.append(0b00011111000111110001111100011111);
  bs.append(0b00011111000111110001111100011111);
  bs.append(0b00011111000111110001111100011111);
  bs.append(0b00011111000111110001111100011111);
  bs.append(0b00011111000111110001111100011111);
  bs.append(0b00011111000111110001111100011111);
  bs.append(0b00011111000111110001111100011111);
  bs.append(0b00011111000111110001111100011111);
  bs.append(0b00011111000111110001111100011111);
  std::cout << bs << std::endl;

  dtl::teb<2> enc_bs(bs);
  std::cout << enc_bs << std::endl;
  std::cout << enc_bs.info() << std::endl;

  auto it = enc_bs.scan_it();
  dtl::bitmap val = bs;
  while (!it.end()) {
    std::cout << "[" << it.pos() << "," << it.length() << ")" << std::endl;
    for (std::size_t i = it.pos(); i < it.pos() + it.length(); ++i) {
      assert(val[i] == true);
      val[i] = false;

    }
    it.next();
  }
  std::cout << bs << std::endl;
  std::cout << val << std::endl;
  assert(val.count() == 0);
}
//===----------------------------------------------------------------------===//
TEST(single,
    DISABLED_lossy_test) {

  // lossy

//  dtl::bitmap bs(16, 0b0000000011111101);
//  dtl::bitmap bs(16, 0b1110111111111101);
//  dtl::bitmap bs(16, 0b1110111111111101);
//  dtl::bitmap bs(16, 0b1111111111011101);
//  dtl::bitmap bs(8, 0b00000011);
//  dtl::bitmap bs(8, 0b00001101);
//  dtl::bitmap bs(8, 0b00000111);

//  const auto bitmap_id = 12;
//  const auto bitmap_id = 3901;
  const auto bitmap_id = 271;
  const std::string DB_FILE = "./random_bitmaps.sqlite3";
  bitmap_db db(DB_FILE);
  dtl::bitmap bs = db.load_bitmap(bitmap_id);

  const auto fpr = 0.2;
  dtl::teb<> enc_bs(bs, fpr);
  std::cout << enc_bs.info() << std::endl;

  auto it = enc_bs.scan_it();
//  auto it = enc_bs.it();
  dtl::bitmap val = bs;
  val.reset();

  for (std::size_t i = 0; i < std::min(std::size_t(64), bs.size()); ++i) {
    std::cout << (bs[i] ? "1" : "0");
  }
  std::cout << std::endl;

  while (!it.end()) {
//    std::cout << "[" << it.pos() << "," << it.length() << ")" << std::endl;
    for (std::size_t i = it.pos(); i < it.pos() + it.length(); ++i) {
      assert(fpr > 0.0 || bs[i] == true);
      assert(val[i] == false);
      val[i] = true;
    }
    it.next();
  }
  for (std::size_t i = 0; i < std::min(std::size_t(64), bs.size()); ++i) {
    std::cout << (val[i] ? "1" : "0");
  }
  std::cout << std::endl;

  u64 max_fp_cnt = static_cast<u64>(bs.size() * fpr);
  const auto fp_cnt = (bs ^ val).count();
  std::cout << fp_cnt << " / " << max_fp_cnt << std::endl;
  ASSERT_TRUE((bs & val) == bs);
  ASSERT_TRUE(fp_cnt <= max_fp_cnt);
  std::cout << bs << std::endl;
  std::cout << val << std::endl;
//  ASSERT_TRUE(val.count() == 0);
}
//===----------------------------------------------------------------------===//
TEST(single,
    test) {

  // lossy

  // sparse
  dtl::bitmap bs(64, 0b0000000010000001000000000000000000000000000000000000000000000000);
  // dense
//  dtl::bitmap bs(64,   0b0100000111000001111100000000001101000000001100000111111000100000);
//  dtl::bitmap bs(16, 0b1110111111111101);
//  dtl::bitmap bs(16, 0b1110111111111101);
//  dtl::bitmap bs(16, 0b1111111111011101);
//  dtl::bitmap bs(8, 0b00000011);
//  dtl::bitmap bs(8, 0b00001101);
//  dtl::bitmap bs(8, 0b00000111);

//  const auto bitmap_id = 12;
//  const auto bitmap_id = 3901;
  const auto bitmap_id = 271;
  const std::string DB_FILE = "./random_bitmaps.sqlite3";
  bitmap_db db(DB_FILE);
//  dtl::bitmap bs = db.load_bitmap(bitmap_id);



  const auto fpr = 0.0;

  std::cout << "basic TEB:" << std::endl;
  std::cout << dtl::teb<0>(bs).info() << std::endl;
  std::cout << dtl::bitmap_tree<0>(bs) << std::endl;

  std::cout << "optimized TEB:" << std::endl;
  std::cout << dtl::teb<>(bs).info() << std::endl;
  std::cout << dtl::bitmap_tree<>(bs) << std::endl;

  std::exit(0);

  dtl::teb<2> enc_bs(bs, fpr);
  std::cout << enc_bs.info() << std::endl;
  dtl::bitmap_tree<2> tree(bs);
  std::cout << dtl::bitmap_tree<0>(bs) << std::endl;
  std::cout << dtl::bitmap_tree<2>(bs) << std::endl;


  auto it = enc_bs.scan_it();
//  auto it = enc_bs.it();
  dtl::bitmap val = bs;
  val.reset();

  for (std::size_t i = 0; i < std::min(std::size_t(64), bs.size()); ++i) {
    std::cout << (bs[i] ? "1" : "0");
  }
  std::cout << std::endl;

  while (!it.end()) {
//    std::cout << "[" << it.pos() << "," << it.length() << ")" << std::endl;
    for (std::size_t i = it.pos(); i < it.pos() + it.length(); ++i) {
      assert(fpr > 0.0 || bs[i] == true);
      assert(val[i] == false);
      val[i] = true;
    }
    it.next();
  }
  for (std::size_t i = 0; i < std::min(std::size_t(64), bs.size()); ++i) {
    std::cout << (val[i] ? "1" : "0");
  }
  std::cout << std::endl;

  u64 max_fp_cnt = static_cast<u64>(bs.size() * fpr);
  const auto fp_cnt = (bs ^ val).count();
  std::cout << fp_cnt << " / " << max_fp_cnt << std::endl;
  ASSERT_TRUE((bs & val) == bs);
  ASSERT_TRUE(fp_cnt <= max_fp_cnt);
  std::cout << bs << std::endl;
  std::cout << val << std::endl;
//  ASSERT_TRUE(val.count() == 0);
}
//===----------------------------------------------------------------------===//
