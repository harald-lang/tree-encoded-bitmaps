#include <iostream>

#include <bitset>

#include <boost/dynamic_bitset.hpp>

#include <dtl/dtl.hpp>
#include <dtl/bits.hpp>
#include <dtl/bitset.hpp>
#include <dtl_storage/tree_mask.hpp>

#include <fastbit/bitvector.h>
#include <fastbit/bitvector64.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/bp_support_gg.hpp>

#include "static_stack.hpp"
#include "two_state_markov_process.hpp"
#include "tree_mask_po.hpp"
#include "tree_mask_util.hpp"
#include "index.hpp"
#include "util.hpp"
#include <thread>
#include <dtl/thread.hpp>

#include "tree_mask_lo.hpp"
#include "dynamic_tree_mask_lo.hpp"


void a() {
  //                     0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9
  //                     ( ( ( ) ( ) ) ( ) ( ( ( ) ( ) ) ( ) ) )
  sdsl::bit_vector bv = {1,1,1,0,1,0,0,1,0,1,1,1,0,1,0,0,1,0,0,0};
  std::cout << bv << std::endl;
  std::cout << "rank:" << std::endl;
  sdsl::bit_vector::rank_1_type bv_rank1(&bv);
  std::cout << bv_rank1(1) << std::endl;
  std::cout << bv_rank1(3) << std::endl;
  std::cout << bv_rank1(5) << std::endl;

  std::cout << "select:" << std::endl;
  sdsl::bit_vector::select_1_type bv_select1(&bv);
  std::cout << bv_select1(1) << std::endl;
  std::cout << bv_select1(2) << std::endl;

  std::cout << "balanced parentheses:" << std::endl;
  sdsl::bp_support_gg<> bps(&bv);
  std::cout << bps.enclose(1) << std::endl;
  std::cout << bps.enclose(4) << std::endl;
  std::cout << bps.enclose(5) << std::endl;
  std::cout << bps.find_close(1) << std::endl;

  std::cout << "bv size:  " << bv.size() << std::endl;
  std::cout << "bps size: " << bps.size() << std::endl;

  boost::dynamic_bitset<uint8_t> dbs;
  dbs.push_back(1);
  dbs.push_back(1);
  dbs.push_back(1);
  dbs.push_back(0);
  dbs.push_back(1);
  std::cout << "dbs:      " << dbs << std::endl;
  std::cout << "dbs size: " << dbs.size() << std::endl;
}

template<std::size_t N>
struct tree {

  dtl::match_tree<N> m;
  static constexpr u64 height = dtl::match_tree<N>::height;

  explicit
  tree(const std::bitset<N>& b) : m(b) {}

  // Point access - O(log N), where N is the length of the bitmap.
  inline constexpr u1
  operator[](std::size_t bit_idx) const noexcept {
    const std::size_t path = bit_idx;
    std::size_t path_pos = m.height - 1;
    std::size_t node_idx = 0; // root node
    // Navigate
    while (m.tree.is_inner_node(node_idx)) {
      node_idx = 2 * node_idx + 1 + dtl::bits::bit_test(path, path_pos);
      path_pos--;
    }
    return m.labels[node_idx];
  }

};

void test_tree_mask_po_xor_re();
void test_tree_mask_po_and();
void test_tree_mask_po_fused_xor();

void test_treemask_lo_encoding_decoding();
void test_treemask_lo_size();
void test_tree_mask_lo_xor_re();
void test_tree_mask_lo_and_re();
void test_tree_mask_lo_and();
void test_tree_mask_lo_fused_xor();
void test_tree_mask_lo_xor_re_compressed();

template<u64 N>
__forceinline__ tree<N>
operator&(const tree<N>& lhs, const tree<N>& rhs) noexcept {
  tree<N> ret_val(lhs);
  static_stack<std::size_t, tree<N>::height * 2> stack;
  stack.push_back(0); // root node
  // Navigate
  while (!stack.empty()) {
    std::size_t node_idx = stack.back();
    stack.pop_back();
    const auto lhs_is_inner = lhs.m.tree.is_inner_node(node_idx);
    const auto rhs_is_inner = rhs.m.tree.is_inner_node(node_idx);
    if (lhs_is_inner & rhs_is_inner) {
      stack.push_back(2 * node_idx + 2); // right child
      stack.push_back(2 * node_idx + 1); // left child
    }
    else {
      // Hit a leaf node
      const auto lhs_is_leaf = !lhs_is_inner;
      const auto rhs_is_leaf = !rhs_is_inner;
      if (lhs_is_leaf & rhs_is_leaf) {
        // Trivial case.
        ret_val.m.labels.set(node_idx, lhs.m.labels[node_idx] & rhs.m.labels[node_idx]);
      }
      else {
        const tree<N>& l = lhs_is_leaf ? lhs : rhs;
        const tree<N>& t = lhs_is_leaf ? rhs : lhs;
        const auto label = l.m.labels[node_idx];
        if (label == false) {
          // Simple case: The entire sub-tree evaluates to 0.
          ret_val.m.tree.set_leaf(node_idx);
          ret_val.m.labels[node_idx] = 0;
        }
        else {
          // Walk the sub-tree.
          static_stack<std::size_t, tree<N>::height * 2> s;
          s.push_back(node_idx);
          while (!s.empty()) {
            std::size_t i = s.back();
            s.pop_back();
            if (t.m.tree.is_inner_node(i)) {
              ret_val.m.tree.set_inner(i);
              s.push_back(2 * node_idx + 1); // left child
              s.push_back(2 * node_idx + 2); // right child
            }
            else {
              ret_val.m.tree.set_leaf(i);
              ret_val.m.labels[i] = t.m.labels[i];
            }
          }
        }
      }
    }
  }
  return ret_val;
};


__forceinline__ bool
bittest(const int64_t word, const uint32_t idx) {
  asm goto("btl	%1, (%0)\n\t"
           "jnc	%l[not_set]"
           :: "r"(&word), "r"(idx) : "cc" : not_set);
  return true;
not_set:
  return false;
}

int t_main() {


//  {
//    constexpr std::size_t skip_width_log2 = 5;
//    constexpr std::size_t skip_width = 1ull << skip_width_log2; // [bits]
//
////    constexpr std::size_t lut_size_log2 = 8;
////    constexpr std::size_t lut_size = 1ull << lut_size_log2;
////    constexpr std::size_t lut_mask = lut_size - 1;
////    std::vector<std::size_t> lut(lut_size, 0);
//    std::vector<uint64_t> histo(skip_width, 0);
//    for (std::size_t i = 0; i < (1ull << skip_width) - 1; i++) {
//      std::bitset<skip_width> s(i);
//      int32_t c = 1;
//      for (std::size_t j = skip_width - 2; j < skip_width; j--) {
//        c = s[j] ? c - 1 : c + 1;
//      }
//      c = c < 0 ? 0 : c;
//      histo[c]++;
//    }
//
//    std::cout << "Histo:" << std::endl;
//    for (std::size_t i = 0; i < histo.size(); i++) {
//      std::cout << i << ";" << histo[i] << std::endl;
//    }
//    std::cout << "---" << std::endl;
//
//
//  }

//  std::exit(0);




  constexpr std::size_t N = 1u << 17;
  std::bitset<N> bitmask_a;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, N-1);
  for (std::size_t i = 0; i < (N/2); i++) {
    bitmask_a[dis(gen)] = true;
  }
  dtl::match_tree<N> match_tree(bitmask_a);

  std::cout << "---" << std::endl;

  tree<N> tree_a(bitmask_a);
  std::bitset<N> bitmask_a_reconstructed;
  for (std::size_t i = 0; i < N; i++) {
    bitmask_a_reconstructed[i] = tree_a[i];
  }
  if (bitmask_a == bitmask_a_reconstructed) {
    std::cout << "passed" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
  }


  std::bitset<N> bitmask_b;
  bitmask_b[8] = true;
  bitmask_b[9] = true;
  bitmask_b[14] = true;
  std::bitset<N> bitmask_c(bitmask_a & bitmask_b);

  std::cout << "a: " << bitmask_a << std::endl;
  std::cout << "b: " << bitmask_b << std::endl;
  std::cout << "c: " << bitmask_c << std::endl;

  tree<N> tree_b(bitmask_b);
  tree<N> tree_c = (tree_a & tree_b);
  std::bitset<N> bitmask_c_reconstructed;
  for (std::size_t i = 0; i < N; i++) {
    bitmask_c_reconstructed[i] = tree_c[i];
  }
  std::cout << "-  " << bitmask_c_reconstructed << std::endl;

  if ((bitmask_a & bitmask_b) == bitmask_c_reconstructed) {
    std::cout << "passed" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
  }

  static_stack<std::size_t, 8> stack;

  two_state_markov_process mp(10, 0.01);
//  dtl::bitset<N> mp_out;
  std::bitset<N> mp_out;
//  boost::bitmap_view<> mp_out(N);
//  sdsl::bit_vector mp_out(N); // does not have a count() function
  for (std::size_t i = 0; i < N; i++) {
    mp_out[i] = mp.next();
  }
  std::cout << mp_out << std::endl;
  std::cout << "\npop count:       " << mp_out.count() << " (" << (mp_out.count()*1.0/N) << ")" << std::endl;
  auto zero_fills_cnt = count_0fills(mp_out);
  auto one_fills_cnt = count_1fills(mp_out);
  std::cout << "# of 0-fills:      " << zero_fills_cnt << std::endl;
  std::cout << "# of 1-fills:      " << one_fills_cnt << std::endl;
  std::cout << "clustering factor: " << (mp_out.count() * 1.0 / one_fills_cnt) << std::endl;
  std::cout << "bit density:       " << (mp_out.count() * 1.0 / mp_out.size()) << std::endl;


//  std::cout << "WAH32" << std::endl;
//  ibis::fileManager::instance(); // initialize file manager, which is responsible for memory management in IBIS.
//  ibis::bitvector wah32;
//  ibis::bitvector64 wah64;
////  wah32.reserve()
//  for (std::size_t i = 0; i < N; i++) {
//    wah32.setBit(i,mp_out[i]);
//    wah64.setBit(i,mp_out[i]);
//  }
////  std::cout << "wah32: " << wah32 << std::endl;
////  std::cout << "size: " << wah32.bytes() << std::endl;
//  wah32.compress();
//  wah64.compress();
//  std::cout << "wah32 size: " << wah32.bytes() << std::endl;
//  std::cout << "wah64 size: " << wah64.bytes() << std::endl;
//
//  ibis::bitvector::pit pit(wah32);
//  while (*pit < N) {
//    std::cout << *pit << std::endl;
//    pit.next();
//  }
//  std::cout << *pit << std::endl;

  std::bitset<N> mp2_out;
  std::bitset<N> mp3_out;
  for (std::size_t i = 0; i < N; i++) {
    mp2_out[i] = mp.next();
  }
  for (std::size_t i = 0; i < N; i++) {
    mp3_out[i] = mp.next();
  }
  dtl::roaring_bitmap<N> a(mp_out);
  dtl::roaring_bitmap<N> b(mp2_out);
  dtl::roaring_bitmap<N> c(mp3_out);
  c.fused_xor_and(a, b);
  std::cout << a.size_in_byte() << std::endl;
  std::cout << b.size_in_byte() << std::endl;
  std::cout << c.size_in_byte() << std::endl;

  auto val = c.to_bitset();
  if (val == ((mp_out^mp2_out)&mp3_out)) {
    std::cout << "passed" << std::endl;
  }
  else {
    std::cout << "FAILED" << std::endl;
  }

  std::cout << mp_out << std::endl;
  dtl::tree_mask_po<N> tm(mp_out);
  std::cout << tm << std::endl;
  std::cout << tm.size_in_byte() << std::endl;
  std::cout << ~mp_out << std::endl;
  dtl::tree_mask_po<N> tmn(~mp_out);
  std::cout << tmn << std::endl;
  std::cout << tmn.size_in_byte() << std::endl;

  auto word = ~0ll >> 1;
  for (uint32_t i = 0; i < 64; i++) {
    assert(bittest(word, i) == dtl::bits::bit_test64(word, i));
  }

  {
    constexpr std::size_t LEN = 8;
    for (std::size_t i = 0; i < (1u << LEN); i++) {
//    for (std::size_t i = 20; i < 21; i++) {
      std::bitset<LEN> bm(i);
      std::cout << "\n-----testing: " << bm << " (" << i << ")" << std::endl;
//      bm[1] = 1;
//      bm[2] = 1;
//      bm[3] = 1;
//      std::cout << "bm: " << bm << std::endl;
      dtl::tree_mask_po<LEN> tm(bm);
//      std::cout << "tm: " << tm << std::endl;
      auto res = tm.to_bitset();
//      std::cout << res << std::endl;
      assert(res == bm);

    }

    std::cout << "navi" << std::endl;
    {

      dtl::tree_mask_po<N> tm(mp_out);
      assert(tm.to_bitset() == mp_out);
    }
//    dtl::tree_mask_po<4>::traversal tmt(tm);
//    std::cout << tmt.is_inner_node() << std::endl;
//    tmt.goto_left_child();
//    std::cout << tmt.is_inner_node() << std::endl;
//    tmt.goto_left_child();
//    std::cout << tmt.is_inner_node() << std::endl;
//    std::cout << tmt.get_label() << std::endl;
//    std::cout << tmt.is_inner_node() << std::endl;
//    tmt.goto_parent();


    {
      std::size_t it = 1;
      // test navigation
      constexpr std::size_t LEN = 1ull << 17;
      // create a bitmap with every second bit set (to prevent tree compression)
      std::bitset<LEN> bm;
      for (std::size_t i = 0; i < LEN; i+=2) {
        bm[i] = true;
      }
      dtl::tree_mask_po<LEN> tm(bm);
      dtl::tree_mask_po<LEN>::traversal tmt(tm);
      do {
        if (!tmt.is_inner_node()) continue;
        const auto current_node_idx = tmt.get_node_idx();
        const auto expected_node_idx = dtl::binary_tree_structure<LEN>::right_child_of(current_node_idx);

        dtl::tree_mask_po<LEN>::traversal tmt2_naive(tmt); // for validation
        tmt2_naive.goto_right_child_naive();
        dtl::tree_mask_po<LEN>::traversal tmt2(tmt);
        tmt2.goto_right_child();
        assert(tmt2.s_pos_ == tmt2_naive.s_pos_);
        assert(tmt2.l_pos_ == tmt2_naive.l_pos_);
        const auto actual_node_idx = tmt2.get_node_idx();
//        std::cout << "path: " << std::bitset<64>(tmt.path_)
//                  << ", current node idx: " << current_node_idx
//                  << ", expected node idx: " << expected_node_idx
//                  << ", actual node idx: " << actual_node_idx
//                  << std::endl;
        assert(actual_node_idx == expected_node_idx);
        if (it == 0) break;
        it--;
      } while (tmt.next());
    }

  }

  std::cout << "xor_re:" << std::endl;
  test_tree_mask_po_xor_re();

}

static void test_tree_mask_lo_and(u32 begin, u32 end){

    constexpr std::size_t LEN = 16;
    for (std::size_t a = begin; a < end; a++) {
      std::bitset<LEN> bm_a(a);
      dtl::tree_mask_lo<LEN> tm_a(bm_a);
      std::cout << "a:" << bm_a << std::endl;

      for (std::size_t b = 0; b < (1u << LEN); b++) {

        //std::cout << "-------------------" << std::endl;

        // make sure, that all bit that are set in a are also set in b (as guaranteed in range encoding)
        //if ((a & b) != a) continue;
        std::bitset<LEN> bm_b(b);
        std::bitset<LEN> bm_expected = bm_a & bm_b;

        dtl::tree_mask_lo<LEN> tm_b(bm_b);

        dtl::tree_mask_lo<LEN> tm_c = tm_a & tm_b;

        std::bitset<LEN> bm_actual = tm_c.to_bitset();

        if (bm_actual != bm_expected) {
          std::cout << "a:" << bm_a << " -> " << tm_a << std::endl;
          std::cout << "b:" << bm_b << " -> " << tm_b << std::endl;
          std::cout << "c:" << bm_actual << " -> " << tm_c << std::endl;
          std::cout << std::endl;

          std::cout << "test: " << bm_a << " (" << a << ") AND " << bm_b << " (" << b << ")" << std::endl;
          std::cout << "Validation failed: expected=" << bm_expected << ", actual=" << bm_actual << std::endl;
          std::exit(1);
        }

        // TODO compression
        /*
        if (!dtl::is_compressed(tm_c)) {
          std::cout << "Validation failed: Resulting tree mask is not compressed." << std::endl;
          dtl::tree_mask_po<LEN> tm_res(bm_expected);
          std::cout << "compressed treemask: " << tm_res << std::endl;
          std::exit(1);
        }
        */

        assert(bm_actual == bm_expected);
      }
    }
}

void run_parallel(){

  u32 max = (1ul << 16);
  u32 step_size = max/std::thread::hardware_concurrency();
  std::mutex m;
  $u32 begin = 0;
  $u32 end = step_size;

  auto func = [&](){

      m.lock();
      u32 this_begin = begin;
      u32 this_end = end;
      begin += step_size;
      end += step_size;
      m.unlock();

      constexpr std::size_t LEN = 16;
      for (std::size_t a = this_begin; a < this_end; a++) {
        std::bitset<LEN> bm_a(a);
        dtl::tree_mask_lo<LEN> tm_a(bm_a);

        std::cout << "a: " << bm_a << std::endl;

        for (std::size_t b = 0; b < (1u << LEN); b++) {

          //std::cout << "-------------------" << std::endl;

          // make sure, that all bit that are set in a are also set in b (as guaranteed in range encoding)
          //if ((a & b) != a) continue;
          std::bitset<LEN> bm_b(b);
          std::bitset<LEN> bm_expected = bm_a & bm_b;

          dtl::tree_mask_lo<LEN> tm_b(bm_b);

          dtl::tree_mask_lo<LEN> tm_c = tm_a & tm_b;

          std::bitset<LEN> bm_actual = tm_c.to_bitset();

          if (bm_actual != bm_expected) {
            std::cout << "a:" << bm_a << " -> " << tm_a << std::endl;
            std::cout << "b:" << bm_b << " -> " << tm_b << std::endl;
            std::cout << "c:" << bm_actual << " -> " << tm_c << std::endl;
            std::cout << std::endl;

            std::cout << "test: " << bm_a << " (" << a << ") AND " << bm_b << " (" << b << ")" << std::endl;
            std::cout << "Validation failed: expected=" << bm_expected << ", actual=" << bm_actual << std::endl;
            std::exit(1);
          }

          // TODO compression
          /*
          if (!dtl::is_compressed(tm_c)) {
            std::cout << "Validation failed: Resulting tree mask is not compressed." << std::endl;
            dtl::tree_mask_po<LEN> tm_res(bm_expected);
            std::cout << "compressed treemask: " << tm_res << std::endl;
            std::exit(1);
          }
          */

          assert(bm_actual == bm_expected);
        }
      }
  };

  std::vector<std::thread> workers;

  dtl::run_in_parallel_async(func, workers);
  std::cout << "Number of threads: " << workers.size() << std::endl;
  dtl::wait_for_threads(workers);

  /*
  std::cout << "waiting for threads to finish..." << std::endl;
  t1.join();
  std::cout << "t1 started" << std::endl;
  t2.join();
  std::cout << "t2 started" << std::endl;
  t3.join();
  std::cout << "t3 started" << std::endl;
  t4.join();
  std::cout << "t4 started" << std::endl;
  std::cout << "done!" << std::endl;
  */
}

template<u64 S>
void test_iterator(double f, double d){

  boost::dynamic_bitset<$u32> bitset{S,0};
  /*two_state_markov_process p(f, d);

  for(auto i = 0; i < S; i++)
    bitset[i] = p.next();
  */

  std::bitset<S> bs_1(0b1111100110111111);

  for(auto i = 0; i < S; i++){
    bitset[i] = bs_1[i];
  }

  std::cout << "Printing the bitset" << std::endl;
  for(auto i = 0; i < S; i++)
    std::cout << bitset[i];
  std::cout << std::endl;

  dtl::dynamic_tree_mask_lo tm_test(bitset);
  tm_test.print(std::cout);
  std::cout << "\n" << std::endl;

  dtl::dynamic_tree_mask_lo::iter it(tm_test);

  dtl::dynamic_tree_mask_lo::iter it_ranges(tm_test);

  while(!it_ranges.end()){
    std::cout << " ["<< it_ranges.pos() << ";" << it_ranges.pos() + it_ranges.length() << "), ";
    it_ranges.next();
  }
  std::cout << std::endl << std::endl;

  for(auto i = it.pos(); i < S; i++){
    std::cout << "Test skipping to: " << i << std::endl;
    dtl::dynamic_tree_mask_lo::iter it_tmp(tm_test);

    std::cout << "Current range: [" << it_tmp.pos() << ";" << it_tmp.length() << ")" << std::endl;

    it_tmp.skip_to(i);
    std::cout << "pos: " << it_tmp.pos() << " length: " << it_tmp.length() << std::endl;
    boost::dynamic_bitset<$u32> tmp{S,0};

    for(auto j = it_tmp.pos(); j < it_tmp.pos() + it_tmp.length(); j++){
      tmp.set(j);
    }

    std::cout << "Resulting bitmap: ";
    for(auto e = 0; e < S; e++){
      std::cout << tmp[e];
    }
    std::cout << std::endl;
    std::cout << "Input bitmap:     ";
    for(auto e = 0; e < S; e++){
      std::cout << bs_1[e];
    }
    std::cout << std::endl;

    std::cout << "Result after next(): ";
    it_tmp.next();
    std::cout << "[" << it_tmp.pos() << ";" << it_tmp.pos()+it_tmp.length() << ")" << std::endl;

    if((tmp | bitset) != bitset){
      std::cout << "Skip to is incorrect" << std::endl;
    }

    std::cout << std::endl << std::endl;
  }

  std::cout << "Done" << std::endl;
}

template<u64 S>
void test_iterator(double f_1, double d_1, double f_2, double d_2){

  boost::dynamic_bitset<$u32> bs_1{S,0};
  boost::dynamic_bitset<$u32> bs_2{S,0};

  two_state_markov_process p_1(f_1, d_1);
  two_state_markov_process p_2(f_2, d_2);

  for(auto i = 0; i < S; i++)
    bs_1[i] = p_1.next();

  for(auto i = 0; i < S; i++)
    bs_2[i] = p_2.next();


  std::cout << "Printing the bitset_1: ";
  for(auto i = 0; i < S; i++)
    std::cout << bs_1[i];
  std::cout << std::endl;

  std::cout << "Printing the bitset_2: ";
  for(auto i = 0; i < S; i++)
    std::cout << bs_2[i];
  std::cout << std::endl;

  boost::dynamic_bitset<$u32> bs_and = bs_1 & bs_2;
  std::cout << "Printing the bitset_and: ";
  for(auto i = 0; i < S; i++)
    std::cout << bs_and[i];
  std::cout << std::endl;

  dtl::dynamic_tree_mask_lo tm_1(bs_1);
  tm_1.print(std::cout);
  std::cout << std::endl;

  dtl::dynamic_tree_mask_lo tm_2(bs_2);
  tm_2.print(std::cout);
  std::cout << "\n" << std::endl;


  std::cout << "Iterator 1" << std::endl;
  dtl::dynamic_tree_mask_lo::iter it_1(tm_1);

  while(!it_1.end()){
    std::cout << "pos: " << it_1.pos() << " length: " << it_1.length() << std::endl;
    it_1.next();
  }

  std::cout << std::endl << "Iterator 2" << std::endl;
  dtl::dynamic_tree_mask_lo::iter it_2(tm_2);

  while(!it_2.end()){
    std::cout << "pos: " << it_2.pos() << " length: " << it_2.length() << std::endl;
    it_2.next();
  }

  dtl::dynamic_tree_mask_lo::iter_and_simple it_s(tm_1, tm_2);

  std::cout << std::endl << "Matches" << std::endl;
  while(!it_s.end()){
    std::cout << "interval: ";
    it_s.matches().print();
    it_s.next();
  }

}

int main(){
  test_iterator<16>(5, 0.5);
  //test_iterator<16>(5, 0.5, 4, 0.4);
  //test_tree_mask_lo_and();
  return 0;
}

/// treemask_po tests

void test_tree_mask_po_xor_re() {
  constexpr std::size_t LEN = 16;
  for (std::size_t a = 0; a < (1u << LEN); a++) {
    for (std::size_t b = 0; b < (1u << LEN); b++) {
//  for (std::size_t a = 0; a < (1u << LEN); a++) {
//    for (std::size_t b = 0; b < (1u << LEN); b++) {
      // make sure, that all bit that are set in a are also set in b (as guaranteed in range encoding)
      if ((a & b) != a) continue;
      std::bitset<LEN> bm_a(a);
      std::bitset<LEN> bm_b(b);
      std::bitset<LEN> bm_expected = bm_a ^ bm_b;
      dtl::tree_mask_po<LEN> tm_a(bm_a);
      dtl::tree_mask_po<LEN> tm_b(bm_b);
      std::cout << "a:" << bm_a << " -> " << tm_a << std::endl;
      std::cout << "b:" << bm_b << " -> " << tm_b << std::endl;


      dtl::tree_mask_po<LEN> tm_c = tm_a ^ tm_b;

      //dtl::tree_mask_po<LEN> tm_c = tm_a.xor_re(tm_b);
      std::bitset<LEN> bm_actual = tm_c.to_bitset();
      std::cout << "c:" << bm_actual << " -> " << tm_c << std::endl;
      std::cout << std::endl;
      if (bm_actual != bm_expected) {
        std::cout << "test: " << bm_a << " (" << a << ") XOR " << bm_b << " (" << b << ")" << std::endl;
        std::cout << "Validation failed: expected=" << bm_expected << ", actual=" << bm_actual << std::endl;
        std::exit(1);
      }
      if (!dtl::is_compressed(tm_c)) {
        std::cout << "Validation failed: Resulting tree mask is not compressed." << std::endl;
        dtl::tree_mask_po<LEN> tm_res(bm_expected);
        std::cout << "compressed treemask: " << tm_res << std::endl;
        std::exit(1);
      }
      assert(bm_actual == bm_expected);
    }
  }
}

void test_tree_mask_po_and() {
  constexpr std::size_t LEN = 16;
  for (std::size_t a = 0; a < (1u << LEN); a++) {
    for (std::size_t b = 0; b < (1u << LEN); b++) {

      // make sure, that all bit that are set in a are also set in b (as guaranteed in range encoding)
      if ((a & b) != a) continue;
      std::bitset<LEN> bm_a(a);
      std::bitset<LEN> bm_b(b);
      std::bitset<LEN> bm_expected = bm_a & bm_b;
      dtl::tree_mask_po<LEN> tm_a(bm_a);
      dtl::tree_mask_po<LEN> tm_b(bm_b);
      std::cout << "a:" << bm_a << " -> " << tm_a << std::endl;
      std::cout << "b:" << bm_b << " -> " << tm_b << std::endl;

      dtl::tree_mask_po<LEN> tm_c = tm_a & tm_b;

      std::bitset<LEN> bm_actual = tm_c.to_bitset();
      std::cout << "c:" << bm_actual << " -> " << tm_c << std::endl;
      std::cout << std::endl;
      if (bm_actual != bm_expected) {
        std::cout << "test: " << bm_a << " (" << a << ") XOR " << bm_b << " (" << b << ")" << std::endl;
        std::cout << "Validation failed: expected=" << bm_expected << ", actual=" << bm_actual << std::endl;
        std::exit(1);
      }
      if (!dtl::is_compressed(tm_c)) {
        std::cout << "Validation failed: Resulting tree mask is not compressed." << std::endl;
        dtl::tree_mask_po<LEN> tm_res(bm_expected);
        std::cout << "compressed treemask: " << tm_res << std::endl;
        std::exit(1);
      }
      assert(bm_actual == bm_expected);
    }
  }
}

void test_tree_mask_po_fused_xor() {
  // test for correctness of: (a XOR b) AND c

  constexpr std::size_t LEN = 8;
  for (std::size_t a = 0; a < (1u << LEN); a++) {
    for (std::size_t b = 0; b < (1u << LEN); b++) {
      // (a XOR b)
      // make sure, that all bit that are set in a are also set in b (as guaranteed in range encoding)
      if ((a & b) != a) continue;
      std::bitset<LEN> bm_a(a);
      std::bitset<LEN> bm_b(b);
      std::bitset<LEN> bm_xor_expected = bm_a ^ bm_b;
      dtl::tree_mask_po<LEN> tm_a(bm_a);
      dtl::tree_mask_po<LEN> tm_b(bm_b);
      std::cout << "a:" << bm_a << " -> " << tm_a << std::endl;
      std::cout << "b:" << bm_b << " -> " << tm_b << std::endl;

      dtl::tree_mask_po<LEN> tm_xor = tm_a ^ tm_b;
      std::bitset<LEN> bm_xor_actual = tm_xor.to_bitset();

      if (bm_xor_actual != bm_xor_expected) {
        std::cout << "test: " << bm_a << " (" << a << ") XOR " << bm_b << " (" << b << ")" << std::endl;
        std::cout << "Validation failed: expected=" << bm_xor_expected << ", actual=" << bm_xor_actual << std::endl;
        std::exit(1);
      }

      std::cout << "xor:" << bm_xor_actual << " -> " << tm_xor << std::endl << "----- all values for c -----" << std::endl;
      std::cout << std::endl;

      if (!dtl::is_compressed(tm_xor)) {
        std::cout << "Validation failed: Resulting tree mask is not compressed." << std::endl;
        dtl::tree_mask_po<LEN> tm_xor_res(bm_xor_expected);
        std::cout << "compressed treemask: " << tm_xor_res << std::endl;
        std::exit(1);
      }
      assert(bm_xor_actual == bm_xor_expected);

      for(std::size_t c = 0; c < (1u << LEN); c++){

        std::bitset<LEN> bm_c(c);
        std::bitset<LEN> bm_fused_expected = bm_xor_expected & bm_c;
        dtl::tree_mask_po<LEN> tm_c(bm_c);
        std::cout << "c:" << bm_c << " -> " << tm_c << std::endl;

        dtl::tree_mask_po<LEN> tm_fused = tm_c & tm_xor;
        std::bitset<LEN> bm_fused_actual = tm_fused.to_bitset();

        if (bm_fused_actual != bm_fused_expected) {
          std::cout << "test: " << bm_c << " (" << c << ") AND " << bm_xor_expected << std::endl;
          std::cout << "Validation failed: expected=" << bm_fused_expected << ", actual=" << bm_fused_actual << std::endl;
          std::exit(1);
        }

        std::cout << "fused:" << bm_fused_actual << " -> " << tm_fused << std::endl;
        std::cout << std::endl;

        if (!dtl::is_compressed(tm_fused)) {
          std::cout << "Validation failed: Resulting tree mask is not compressed." << std::endl;
          dtl::tree_mask_po<LEN> tm_fused_res(bm_xor_expected);
          std::cout << "compressed treemask: " << tm_fused_res << std::endl;
          std::exit(1);
        }
        assert(bm_fused_actual == bm_fused_expected);
      }
    }
  }
}

/// treemask_lo tests

void test_treemask_lo_encoding_decoding(){

  constexpr std::size_t LEN = 16;

  for(auto i = 0; i < (1u << LEN); ++i) {
    std::cout << "i: " << i << std::endl;
    std::bitset<LEN> bs(i);
    dtl::tree_mask_lo<LEN> t(bs);
    std::bitset<LEN> dec = t.to_bitset();

    if(dec != bs){
      std::cout << "Wrong for:" << std::endl;
      std::cout << "Bitset:  " << bs << std::endl;
      std::cout << "Decoded: " << dec << std::endl;
    }

    if(LEN == 64){
      if(i == std::numeric_limits<u64>::max()){
        std::cout << "Reached the max" << std::endl;
        break;
      }
    }
  }
}

template<u64 N>
void test_treemask_lo_size(){
  std::bitset<N> bs;

  for(auto i = 0; i < bs.size(); i++){
    bs[i] = std::rand() % 2;
  }

  dtl::tree_mask_lo<N> t(bs);

  std::cout << "Bitset: " << bs << std::endl;
  t.size_in_byte();
}

void test_tree_mask_lo_xor_re(){

  constexpr std::size_t LEN = 16;
  for (std::size_t a = 0; a < (1u << LEN); a++) {
    std::bitset<LEN> bm_a(a);
    dtl::tree_mask_lo<LEN> tm_a(bm_a);
    std::cout << "a:" << bm_a << std::endl;

    for (std::size_t b = 0; b < (1u << LEN); b++) {

      //std::cout << "-------------------" << std::endl;

      // make sure, that all bit that are set in a are also set in b (as guaranteed in range encoding)
      if ((a & b) != a) continue;
      std::bitset<LEN> bm_b(b);
      std::bitset<LEN> bm_expected = bm_a ^ bm_b;

      dtl::tree_mask_lo<LEN> tm_b(bm_b);

      dtl::tree_mask_lo<LEN> tm_c = tm_a ^ tm_b;

      std::bitset<LEN> bm_actual = tm_c.to_bitset();

      if (bm_actual != bm_expected) {
        std::cout << "a:" << bm_a << " -> " << tm_a << std::endl;
        std::cout << "b:" << bm_b << " -> " << tm_b << std::endl;
        std::cout << "c:" << bm_actual << " -> " << tm_c << std::endl;
        std::cout << std::endl;

        std::cout << "test: " << bm_a << " (" << a << ") XOR " << bm_b << " (" << b << ")" << std::endl;
        std::cout << "Validation failed: expected=" << bm_expected << ", actual=" << bm_actual << std::endl;
        std::exit(1);
      }

      // TODO compression
      /*
      if (!dtl::is_compressed(tm_c)) {
        std::cout << "Validation failed: Resulting tree mask is not compressed." << std::endl;
        dtl::tree_mask_po<LEN> tm_res(bm_expected);
        std::cout << "compressed treemask: " << tm_res << std::endl;
        std::exit(1);
      }
      */

      assert(bm_actual == bm_expected);
    }
  }

}

void test_tree_mask_lo_and_re(){

  constexpr std::size_t LEN = 16;
  for (std::size_t a = 0; a < (1u << LEN); a++) {
    std::bitset<LEN> bm_a(a);
    dtl::tree_mask_lo<LEN> tm_a(bm_a);
    std::cout << "a:" << bm_a << std::endl;

    for (std::size_t b = 0; b < (1u << LEN); b++) {

      //std::cout << "-------------------" << std::endl;

      // make sure, that all bit that are set in a are also set in b (as guaranteed in range encoding)
      if ((a & b) != a) continue;
      std::bitset<LEN> bm_b(b);
      std::bitset<LEN> bm_expected = bm_a & bm_b;

      dtl::tree_mask_lo<LEN> tm_b(bm_b);

      dtl::tree_mask_lo<LEN> tm_c = tm_a & tm_b;

      std::bitset<LEN> bm_actual = tm_c.to_bitset();

      if (bm_actual != bm_expected) {
        std::cout << "a:" << bm_a << " -> " << tm_a << std::endl;
        std::cout << "b:" << bm_b << " -> " << tm_b << std::endl;
        std::cout << "c:" << bm_actual << " -> " << tm_c << std::endl;
        std::cout << std::endl;

        std::cout << "test: " << bm_a << " (" << a << ") AND " << bm_b << " (" << b << ")" << std::endl;
        std::cout << "Validation failed: expected=" << bm_expected << ", actual=" << bm_actual << std::endl;
        std::exit(1);
      }

      // TODO compression
      /*
      if (!dtl::is_compressed(tm_c)) {
        std::cout << "Validation failed: Resulting tree mask is not compressed." << std::endl;
        dtl::tree_mask_po<LEN> tm_res(bm_expected);
        std::cout << "compressed treemask: " << tm_res << std::endl;
        std::exit(1);
      }
      */

      assert(bm_actual == bm_expected);
    }
  }

}

void test_tree_mask_lo_and(){

  constexpr std::size_t LEN = 16;
  for (std::size_t a = 0; a < (1u << LEN); a++) {
    std::bitset<LEN> bm_a(a);
    dtl::tree_mask_lo<LEN> tm_a(bm_a);
    std::cout << "a:" << bm_a << std::endl;

    for (std::size_t b = 0; b < (1u << LEN); b++) {

      //std::cout << "-------------------" << std::endl;

      // make sure, that all bit that are set in a are also set in b (as guaranteed in range encoding)
      //if ((a & b) != a) continue;
      std::bitset<LEN> bm_b(b);
      std::bitset<LEN> bm_expected = bm_a & bm_b;

      dtl::tree_mask_lo<LEN> tm_b(bm_b);

      dtl::tree_mask_lo<LEN> tm_c = tm_a & tm_b;

      std::bitset<LEN> bm_actual = tm_c.to_bitset();

      if (bm_actual != bm_expected) {
        std::cout << "a:" << bm_a << " -> " << tm_a << std::endl;
        std::cout << "b:" << bm_b << " -> " << tm_b << std::endl;
        std::cout << "c:" << bm_actual << " -> " << tm_c << std::endl;
        std::cout << std::endl;

        std::cout << "test: " << bm_a << " (" << a << ") AND " << bm_b << " (" << b << ")" << std::endl;
        std::cout << "Validation failed: expected=" << bm_expected << ", actual=" << bm_actual << std::endl;
        std::exit(1);
      }

      // TODO compression
      /*
      if (!dtl::is_compressed(tm_c)) {
        std::cout << "Validation failed: Resulting tree mask is not compressed." << std::endl;
        dtl::tree_mask_po<LEN> tm_res(bm_expected);
        std::cout << "compressed treemask: " << tm_res << std::endl;
        std::exit(1);
      }
      */

      assert(bm_actual == bm_expected);
    }
  }

}

void test_tree_mask_lo_fused_xor() {
  // test for correctness of: (a XOR b) AND c

  constexpr std::size_t LEN = 8;
  for (std::size_t a = 0; a < (1u << LEN); a++) {
    for (std::size_t b = 0; b < (1u << LEN); b++) {
      // (a XOR b)
      // make sure, that all bit that are set in a are also set in b (as guaranteed in range encoding)
      if ((a & b) != a) continue;
      std::bitset<LEN> bm_a(a);
      std::bitset<LEN> bm_b(b);
      std::bitset<LEN> bm_xor_expected = bm_a ^ bm_b;
      dtl::tree_mask_lo<LEN> tm_a(bm_a);
      dtl::tree_mask_lo<LEN> tm_b(bm_b);
      std::cout << "a:" << bm_a << " -> " << tm_a << std::endl;
      std::cout << "b:" << bm_b << " -> " << tm_b << std::endl;

      dtl::tree_mask_lo<LEN> tm_xor = tm_a ^ tm_b;
      std::bitset<LEN> bm_xor_actual = tm_xor.to_bitset();

      if (bm_xor_actual != bm_xor_expected) {
        std::cout << "test: " << bm_a << " (" << a << ") XOR " << bm_b << " (" << b << ")" << std::endl;
        std::cout << "Validation failed: expected=" << bm_xor_expected << ", actual=" << bm_xor_actual << std::endl;
        std::exit(1);
      }

      std::cout << "xor:" << bm_xor_actual << " -> " << tm_xor << std::endl << "----- all values for c -----" << std::endl;
      std::cout << std::endl;
      /*
      if (!dtl::is_compressed(tm_xor)) {
        std::cout << "Validation failed: Resulting tree mask is not compressed." << std::endl;
        dtl::tree_mask_po<LEN> tm_xor_res(bm_xor_expected);
        std::cout << "compressed treemask: " << tm_xor_res << std::endl;
        std::exit(1);
      }
      */
      assert(bm_xor_actual == bm_xor_expected);

      for(std::size_t c = 0; c < (1u << LEN); c++){

        std::bitset<LEN> bm_c(c);
        std::bitset<LEN> bm_fused_expected = bm_xor_expected & bm_c;
        dtl::tree_mask_lo<LEN> tm_c(bm_c);
        std::cout << "c:" << bm_c << " -> " << tm_c << std::endl;

        dtl::tree_mask_lo<LEN> tm_fused = tm_c & tm_xor;
        std::bitset<LEN> bm_fused_actual = tm_fused.to_bitset();

        if (bm_fused_actual != bm_fused_expected) {
          std::cout << "test: " << bm_c << " (" << c << ") AND " << bm_xor_expected << std::endl;
          std::cout << "Validation failed: expected=" << bm_fused_expected << ", actual=" << bm_fused_actual << std::endl;
          std::exit(1);
        }

        std::cout << "fused:" << bm_fused_actual << " -> " << tm_fused << std::endl;
        std::cout << std::endl;
        /*
        if (!dtl::is_compressed(tm_fused)) {
          std::cout << "Validation failed: Resulting tree mask is not compressed." << std::endl;
          dtl::tree_mask_po<LEN> tm_fused_res(bm_xor_expected);
          std::cout << "compressed treemask: " << tm_fused_res << std::endl;
          std::exit(1);
        }
        */
        assert(bm_fused_actual == bm_fused_expected);
      }
    }
  }
}

void test_tree_mask_lo_xor_re_compressed(){

  constexpr std::size_t LEN = 16;
  for (std::size_t a = 0; a < (1u << LEN); a++) {
    std::bitset<LEN> bm_a(a);
    dtl::tree_mask_lo<LEN> tm_a(bm_a);
    std::cout << bm_a << std::endl;

    for (std::size_t b = 0; b < (1u << LEN); b++) {
      //std::cout << "-------------------" << std::endl;

      // make sure, that all bit that are set in a are also set in b (as guaranteed in range encoding)
      if ((a & b) != a) continue;
      std::bitset<LEN> bm_b(b);
      std::bitset<LEN> bm_expected = bm_a ^ bm_b;

      dtl::tree_mask_lo<LEN> tm_b(bm_b);

      //std::cout << "a:" << bm_a << " -> " << tm_a << std::endl;
      //std::cout << "b:" << bm_b << " -> " << tm_b << std::endl;
      dtl::tree_mask_lo<LEN> tm_c = tm_a.xor_compressed(tm_b);
      //std::cout << "c:" << bm_expected << " -> " << tm_c << std::endl;
      std::bitset<LEN> bm_actual = tm_c.to_bitset();

      if (bm_actual != bm_expected) {
        std::cout << "a:" << bm_a << " -> " << tm_a << std::endl;
        std::cout << "b:" << bm_b << " -> " << tm_b << std::endl;
        std::cout << "c:" << bm_actual << " -> " << tm_c << std::endl;
        std::cout << std::endl;

        std::cout << "test: " << bm_a << " (" << a << ") XOR " << bm_b << " (" << b << ")" << std::endl;
        std::cout << "Validation failed: expected=" << bm_expected << ", actual=" << bm_actual << std::endl;
        std::exit(1);
      }

      dtl::tree_mask_lo<LEN> tm_comp(bm_expected);
      if (tm_comp != tm_c) {
        std::cout << "Validation failed: Resulting tree mask is not compressed." << std::endl;
        std::cout << "compressed treemask: " << tm_comp << std::endl;
        std::exit(1);
      }

      assert(bm_actual == bm_expected);
    }
  }

}

