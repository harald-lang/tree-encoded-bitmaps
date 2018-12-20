#include <chrono>
#include <iostream>
#include <set>

#include <dtl/dtl.hpp>

#include <dtl/bitmap/util/two_state_markov_process.hpp>
#include <dtl/bitmap/dynamic_roaring_bitmap.hpp>
#include <dtl/bitmap/static/tree_mask_po.hpp>
#include <dtl/bitmap/static/wah.hpp>
#include <dtl/bitmap/static/roaring_bitmap.hpp>
#include <dtl/bitmap/static/tree_mask_lo.hpp>
#include <dtl/bitmap/static/partitioned_tree_mask.hpp>
#include <dtl/bitmap/dynamic_partitioned_tree_mask.hpp>
#include <dtl/bitmap/dynamic_tree_mask_lo.hpp>
#include <navin/util.hpp>
#include <dtl/bitmap/dynamic_bitmap.hpp>
#include <dtl/thread.hpp>

// The number of independent runs.
static constexpr u64 RUNS = 100;
//static constexpr u64 RUNS = 1;
static constexpr u64 N = 1u << 20;
//static constexpr u64 N = 1u << 22;

static constexpr u64 M = 1024;

static const i64 RUN_ID = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();


u64
now_nanos() {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count();
}


std::bitset<N>
gen_bitmap(f64 f, f64 d) {
  // init bitset
  f64 f_min = d >= 1.0 ? N : d/(1-d);
  f64 f_actual = std::max(f, f_min);

  two_state_markov_process mp(f_actual, d);
  std::bitset<N> bs;

  const std::size_t max_tries = 5;
  for (std::size_t retry = 1; retry < max_tries; retry++) {
    for ($u64 i = 0; i < N; i++) {
      bs[i] = mp.next();
    }

    $f64 d_actual = (bs.count() * 1.0) / N;
    if (std::abs(d - d_actual) > 1
        || std::abs(f - f_actual) > 0.25) {
      if (max_tries == retry - 1) {
        throw std::runtime_error("Failed to generate a bitmap with f=" + std::to_string(f) + " and d=" + std::to_string(f) + ".");
      }
    }
    else {
      break;
    }
  }
  return bs;
}

template<typename T, typename B>
void run(f64 f, f64 d, const B& bitmap, const std::vector<$u32>& lookup_pos, std::ostream& os) {

  T enc_bitmap(bitmap);

  const auto size_in_bytes = enc_bitmap.size_in_byte();

  // performance measurement
  $u32 match_pos[1024];

  $u32* match_writer = match_pos;

  u64 nanos_begin = now_nanos();
  $u64 repetition_cntr = 0;
  while (now_nanos() - nanos_begin < 2000000000ull) {
    for (std::size_t r = 0; r < 10; r++) {
      match_writer = match_pos;
      for ($u32 i = 0; i < lookup_pos.size(); i++) {
        u32 pos = lookup_pos[i];
        u1 bit = enc_bitmap.test(pos);
        *match_writer = pos;
        match_writer += bit;
      }
    }
    repetition_cntr += 10;
  }
  u64 match_cnt = match_writer - match_pos;
  u64 nanos_end = now_nanos();
  u64 nanos_cnt = nanos_end - nanos_begin;

  // validation code
  $u32* match_reader = match_pos;
  for ($u32 i = 0; i < lookup_pos.size(); i++) {
    u32 pos = lookup_pos[i];
    u1 bit = bitmap.test(pos);
    if (bit && *match_reader != pos) {
      std::cerr << "validation failed" << std::endl;
      std::exit(1);
    }
    match_reader += bit;
  }

  const auto d_actual = (bitmap.count() * 1.0) / bitmap.size();
  const auto f_actual = (bitmap.count() * 1.0) / count_1fills(bitmap);
  os << RUN_ID
     << "," << N
     << "," << enc_bitmap.name()
//     << "," << T::name()
     << "," << d
     << "," << d_actual
     << "," << f
     << "," << f_actual
     << "," << size_in_bytes
     << "," << (nanos_cnt * 1.0) / (lookup_pos.size() * repetition_cntr)
     << std::endl;
}

void run(f64 f, f64 d, const std::vector<$u32>& lookup_pos, std::ostream& os) {
  auto bitmap = gen_bitmap(f, d);
  auto dynamic_bitmap = to_dynamic_bitset(bitmap);
  run<dtl::roaring_bitmap<N>>(f, d, bitmap, lookup_pos, std::cout);
  run<dtl::wah32<N>>(f, d, bitmap, lookup_pos, std::cout);
  run<dtl::tree_mask_po<N>>(f, d, bitmap, lookup_pos, std::cout);
  run<dtl::tree_mask_lo<N>>(f, d, bitmap, lookup_pos, std::cout);
  run<dtl::dynamic_tree_mask_lo>(f, d, dynamic_bitmap, lookup_pos, std::cout);
  run<dtl::partitioned_tree_mask<N>>(f, d, bitmap, lookup_pos, std::cout);
  run<dtl::dynamic_partitioned_tree_mask>(f, d, dynamic_bitmap, lookup_pos, std::cout);
}

std::vector<$u32>
gen_uniform_unique(u32 from, u32 to, u32 cnt) {
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<$u32> dis(from, to - 1);

  std::set<$u32> s;
  while (s.size() < cnt) {
    s.insert(dis(gen));
  }
  std::vector<$u32> r(s.begin(), s.end());
  return r;
}

std::vector<$u32>
gen_uniform(u32 from, u32 to, u32 cnt) {
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_int_distribution<$u32> dis(from, to - 1);

  std::vector<$u32> r;
  r.reserve(cnt);
  std::size_t i = 0;
  while (r.size() < cnt) {
//    r.push_back(dis(gen));
    r.push_back(i%N);
    i++;
  }
  return r;
}


void point_queries() {
  const auto lookup_pos = M < N ? gen_uniform_unique(0, N, 1024)
                                : gen_uniform(0, N, 1024);

  auto f = 64;
//  auto d = 0.005;
  auto d = 0.01;
  run(f, d, lookup_pos, std::cout);
}


//===----------------------------------------------------------------------===//
template<typename T, typename F>
void __forceinline__
bitwise_and_using_skip_iterator(T bitmap_a, T bitmap_b, F consumer) {
  auto it_a = bitmap_a.it();
  auto it_b = bitmap_b.it();
  while (!(it_a.end() || it_b.end())) { // TODO should be ||
    const auto a_begin = it_a.pos();
    const auto a_end   = it_a.pos() + it_a.length();
    const auto b_begin = it_b.pos();
    const auto b_end   = it_b.pos() + it_b.length();

    const auto begin_max = std::max(a_begin, b_begin);
    const auto end_min   = std::min(a_end,   b_end);

    u1 overlap = begin_max < end_min;

    if (overlap) {
      if (a_end == b_end) {
        it_a.next();
        it_b.next();
      }
      else if (a_end <= b_end) {
        it_a.next();
      }
      else {
        it_b.next();
      }
    }
    else {
      // no overlap
      if (a_end == b_end) {
        it_a.next();
        it_b.next();
      }
      else if (a_end < b_end) {
        it_a.nav_to(b_begin);
      }
      else {
        it_b.nav_to(a_begin);
      }
    }

    if (overlap) {
      consumer(begin_max, end_min);
    }
  }
}
//===----------------------------------------------------------------------===//


template<typename T, typename ...Params>
void run_and(f64 f_a, f64 d_a, f64 f_b, f64 d_b,
             const boost::dynamic_bitset<$u32>& bitmap_a,
             const boost::dynamic_bitset<$u32>& bitmap_b,
             std::ostream& os,
             Params&&... constructor_params) {

  const T enc_bitmap_a(bitmap_a, std::forward<Params>(constructor_params)...);
  const T enc_bitmap_b(bitmap_b, std::forward<Params>(constructor_params)...);

  const auto size_in_bytes_a = enc_bitmap_a.size_in_byte();
  const auto size_in_bytes_b = enc_bitmap_b.size_in_byte();


  // performance measurement
  u64 nanos_begin = now_nanos();
  $u64 repetition_cntr = 0;
  $u64 match_cnt = 0;
  while (now_nanos() - nanos_begin < 200000000ull) {
    static constexpr std::size_t MIN_REP = 50;
    for (std::size_t r = 0; r < MIN_REP; r++) {
      bitwise_and_using_skip_iterator(enc_bitmap_a, enc_bitmap_b, [&](u64 begin, u64 end) { match_cnt += end - begin; });
    }
    repetition_cntr += MIN_REP;
  }
  u64 nanos_end = now_nanos();
  u64 nanos_cnt = nanos_end - nanos_begin;

  {
    // Validation code
    const auto expected = bitmap_a & bitmap_b;
    auto actual = expected;
    actual.reset();
    bitwise_and_using_skip_iterator(enc_bitmap_a, enc_bitmap_b,
        [&](u64 begin, u64 end) {
          for (std::size_t i = begin; i < end; ++i) {
            actual[i] = true;
          }
        });
    if (actual != expected) {
      std::cerr << "expected: " << expected << std::endl;
      std::cerr << "actual:   " << actual << std::endl;
      std::cerr << "validation failed" << std::endl;
      std::exit(1);
    }
  }

  const auto d_actual_a = (bitmap_a.count() * 1.0) / bitmap_a.size();
  const auto d_actual_b = (bitmap_b.count() * 1.0) / bitmap_b.size();
  const auto f_actual_a = (bitmap_a.count() * 1.0) / count_1fills(bitmap_a);
  const auto f_actual_b = (bitmap_b.count() * 1.0) / count_1fills(bitmap_b);

  os << RUN_ID
     << "," << N
     << "," << enc_bitmap_a.name()
//     << "," << T::name()
     << "," << d_a
     << "," << d_actual_a
     << "," << f_a
     << "," << f_actual_a
     << "," << size_in_bytes_a
     << "," << d_b
     << "," << d_actual_b
     << "," << f_b
     << "," << f_actual_b
     << "," << size_in_bytes_b
     << "," << repetition_cntr
     << "," << (nanos_cnt * 1.0) / ((bitmap_a.size() * 1.0) * repetition_cntr)
     << std::endl;
}




template<typename T, typename B>
void run_and_iter(f64 f_a, f64 d_a, f64 f_b, f64 d_b,
             const B& bitmap_a, const B& bitmap_b,
             std::ostream& os) {

  const T enc_bitmap_a(bitmap_a);
  const T enc_bitmap_b(bitmap_b);

  const auto size_in_bytes_a = enc_bitmap_a.size_in_byte();
  const auto size_in_bytes_b = enc_bitmap_b.size_in_byte();

  std::vector<dtl::dynamic_tree_mask_lo::range> res;

  // performance measurement
  u64 nanos_begin = now_nanos();
  $u64 repetition_cntr = 0;
  $u64 dep = 0;

  while (now_nanos() - nanos_begin < 2000000000ull) {
    for (std::size_t r = 0; r < 10; r++) {

      dtl::dynamic_tree_mask_lo::iter_and it_and(enc_bitmap_a, enc_bitmap_b);
      while(!it_and.end()){
        //res.push_back(it_and.matches());
        it_and.next();
      }

      dep += res.size();
      //res.clear();
    }
    repetition_cntr += 10;
  }
  u64 nanos_end = now_nanos();
  u64 nanos_cnt = nanos_end - nanos_begin;

  // validation code
  {
    const T res = enc_bitmap_a & enc_bitmap_b;
    const B res_valid = bitmap_a & bitmap_b;
    for ($u32 i = 0; i < res_valid.size(); i++) {
      //TODO
      //if (res_valid.test(i) != res.test(i)) {
      //  std::cerr << "validation failed" << std::endl;
      //  std::exit(1);
      //}
    }
  }

  const auto d_actual_a = (bitmap_a.count() * 1.0) / bitmap_a.size();
  const auto d_actual_b = (bitmap_b.count() * 1.0) / bitmap_b.size();
  const auto f_actual_a = (bitmap_a.count() * 1.0) / count_1fills(bitmap_a);
  const auto f_actual_b = (bitmap_b.count() * 1.0) / count_1fills(bitmap_b);

  os << RUN_ID
     << "," << N
     << "," << T::name() << " via Iterator"
     //<< "," << d_a
     //<< "," << d_actual_a
     //<< "," << f_a
     //<< "," << f_actual_a
     //<< "," << size_in_bytes_a
     //<< "," << d_b
     //<< "," << d_actual_b
     //<< "," << f_b
     //<< "," << f_actual_b
     //<< "," << size_in_bytes_b
     //<< "," << dep / repetition_cntr
     << "," << (nanos_cnt * 1.0) / ((bitmap_a.size() * 1.0) * repetition_cntr)
     << std::endl;
}

template<typename _bitmap_t = dtl::dynamic_bitmap<$u32>, typename ...params>
void
iterate_1fills(const boost::dynamic_bitset<$u32>& input, params&&... constructor_params) {

  // encode the input bitmap
  _bitmap_t tm(input, std::forward<params>(constructor_params)...);

  // iterate over 1fills
  const auto nanos_begin = now_nanos();
  $u64 pos = 0;
  $u64 len = 0;
  for (std::size_t r = 0; r < RUNS; r++) {
    auto it = tm.it();
    while (!it.end()) {
      pos += it.pos();
      len += it.length();
      it.next();
    }
  }
  const auto nanos_end = now_nanos();

  const auto one_fill_cnt = count_1fills(input);
  const auto pop_count = input.count();
  const auto density = (input.count() * 1.0) / input.size();
  std::cout << "size=" << tm.size_in_byte() << " (" << (input.size() + 7 / 8) << ")" <<  std::endl;
  std::cout << "1fill count=" << one_fill_cnt << std::endl;
  std::cout << pos << ", " << len << std::endl;
  std::cout << ((nanos_end - nanos_begin) * 1.0) / (N * RUNS) << " [ns/elem]" << std::endl;
  std::cout << ((nanos_end - nanos_begin) * 1.0) / (N * RUNS * density) << " [ns/1bit]" << std::endl;
  std::cout << ((nanos_end - nanos_begin) * 1.0) / (one_fill_cnt * RUNS) << " [ns/1fill]" << std::endl;
}

$i32 main() {

  dtl::this_thread::set_cpu_affinity(0);

  std::cerr << "run_id=" << std::endl;

  //point_queries();

  //static constexpr u64 N = 1u << 20;
  std::vector<std::size_t> n_values { 1u << 10, 1u << 12, 1u << 14, 1u << 16, 1u << 18, 1u << 20 };
  std::vector<double> densities {0.0001, 0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 0.9, 1};

  for(auto f_a = 1; f_a < 400; f_a += 10) {
    for(auto f_b = 1; f_b < 400; f_b += 10) {
      for (auto d_a : densities) {
        for (auto d_b : densities) {

          const auto bitmap_a = gen_bitmap(f_a, d_a);
          const auto dynamic_bitmap_a = to_dynamic_bitset(bitmap_a);
          const auto bitmap_b = gen_bitmap(f_b, d_b);
          const auto dynamic_bitmap_b = to_dynamic_bitset(bitmap_b);
          run_and<dtl::dynamic_roaring_bitmap>(f_a, d_a, f_b, d_b, dynamic_bitmap_a, dynamic_bitmap_b, std::cout);
          run_and<dtl::dynamic_tree_mask_lo>  (f_a, d_a, f_b, d_b, dynamic_bitmap_a, dynamic_bitmap_b, std::cout);
          run_and<dtl::dynamic_partitioned_tree_mask>(f_a, d_a, f_b, d_b, dynamic_bitmap_a, dynamic_bitmap_b, std::cout,  64);
        }
      }
    }
  }
//      const auto f_a = i;
//      const auto d_a = 0.05;
//      const auto f_b = j;
//      const auto d_b = 0.27;
//      const auto bitmap_a = gen_bitmap(f_a, d_a);
//      const auto dynamic_bitmap_a = to_dynamic_bitset(bitmap_a);
//      const auto bitmap_b = gen_bitmap(f_b, d_b);
//      const auto dynamic_bitmap_b = to_dynamic_bitset(bitmap_b);

//      run_and<dtl::dynamic_roaring_bitmap>            (f_a, d_a, f_b, d_b, dynamic_bitmap_a, dynamic_bitmap_b, std::cout);
  //    run_and<dtl::wah32<N>>                     (f_a, d_a, f_b, d_b, bitmap_a, bitmap_b, std::cout);
//      run_and<dtl::dynamic_tree_mask_lo>         (f_a, d_a, f_b, d_b, dynamic_bitmap_a, dynamic_bitmap_b, std::cout);
  //    run_and<dtl::dynamic_partitioned_tree_mask>(f_a, d_a, f_b, d_b, dynamic_bitmap_a, dynamic_bitmap_b, std::cout,   8);
  //    run_and<dtl::dynamic_partitioned_tree_mask>(f_a, d_a, f_b, d_b, dynamic_bitmap_a, dynamic_bitmap_b, std::cout,  16);
  //    run_and<dtl::dynamic_partitioned_tree_mask>(f_a, d_a, f_b, d_b, dynamic_bitmap_a, dynamic_bitmap_b, std::cout,  32);
  //    run_and<dtl::dynamic_partitioned_tree_mask>(f_a, d_a, f_b, d_b, dynamic_bitmap_a, dynamic_bitmap_b, std::cout,  64);
  //    run_and<dtl::dynamic_partitioned_tree_mask>(f_a, d_a, f_b, d_b, dynamic_bitmap_a, dynamic_bitmap_b, std::cout, 128);
//    }
//  }


  /*
  const auto f_a = 20;
  const auto d_a = 0.05;
  const auto f_b = 8;
  const auto d_b = 0.27;
  const auto bitmap_a = gen_bitmap(f_a, d_a);
  const auto dynamic_bitmap_a = to_dynamic_bitset(bitmap_a);
  const auto bitmap_b = gen_bitmap(f_b, d_b);
  const auto dynamic_bitmap_b = to_dynamic_bitset(bitmap_b);

  run_and<dtl::roaring_bitmap<N>>                 (f_a, d_a, f_b, d_b, bitmap_a, bitmap_b, std::cout);
  run_and<dtl::wah32<N>>                          (f_a, d_a, f_b, d_b, bitmap_a, bitmap_b, std::cout);
  run_and<dtl::dynamic_tree_mask_lo>              (f_a, d_a, f_b, d_b, dynamic_bitmap_a, dynamic_bitmap_b, std::cout);
  //run_and_iter<dtl::dynamic_tree_mask_lo>         (f_a, d_a, f_b, d_b, dynamic_bitmap_a, dynamic_bitmap_b, std::cout);
   */
}

/*
$i32 main() {

  std::cerr << "run_id=" << std::endl;

//  point_queries();
//
//  const auto f_a = 200;
//  const auto d_a = 0.15;
//  const auto f_b = 8;
//  const auto d_b = 0.27;
//  const auto bitmap_a = gen_bitmap(f_a, d_a);
//  const auto dynamic_bitmap_a = to_dynamic_bitset(bitmap_a);
//  const auto bitmap_b = gen_bitmap(f_b, d_b);
//  const auto dynamic_bitmap_b = to_dynamic_bitset(bitmap_b);
//
//  run_and<dtl::roaring_bitmap<N>>            (f_a, d_a, f_b, d_b, bitmap_a, bitmap_b, std::cout);
//  run_and<dtl::wah32<N>>                     (f_a, d_a, f_b, d_b, bitmap_a, bitmap_b, std::cout);
////  run_and<dtl::tree_mask_o<N>>              (f_a, d_a, f_b, d_b, bitmap_a, bitmap_b, std::cout);
//  run_and<dtl::dynamic_tree_mask_lo>         (f_a, d_a, f_b, d_b, dynamic_bitmap_a, dynamic_bitmap_b, std::cout);
//  run_and<dtl::dynamic_partitioned_tree_mask>(f_a, d_a, f_b, d_b, dynamic_bitmap_a, dynamic_bitmap_b, std::cout);

  const auto f_a = 18;
  const auto d_a = 0.05;
  const auto bitmap_a = gen_bitmap(f_a, d_a);
  const auto dynamic_bitmap_a = to_dynamic_bitset(bitmap_a);
  const auto one_fill_cnt = count_1fills(dynamic_bitmap_a);

  std::cout << "input:"
            << " popcnt=" << bitmap_a.count()
            << " 1-fill cnt=" << one_fill_cnt
            << std::endl;


  iterate_1fills<dtl::dynamic_tree_mask_lo>(dynamic_bitmap_a);
  std::cout << std::endl;
  iterate_1fills<dtl::dynamic_bitmap<$u32>>(dynamic_bitmap_a);
  std::cout << std::endl;
  iterate_1fills<dtl::dynamic_partitioned_tree_mask>(dynamic_bitmap_a, 64);
  std::cout << std::endl;

}
*/