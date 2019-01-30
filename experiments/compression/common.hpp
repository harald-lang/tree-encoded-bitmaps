#pragma once

#include <bitset>
#include <chrono>
#include <cmath>
#include <limits>
#include <iostream>

#include <boost/algorithm/string.hpp>

#include <dtl/dtl.hpp>
#include <dtl/env.hpp>
#include <dtl/thread.hpp>

#include <dtl/bitmap/dynamic_bitmap.hpp>
#include <dtl/bitmap/dynamic_roaring_bitmap.hpp>
#include <dtl/bitmap/teb.hpp>
#include <dtl/bitmap/dynamic_wah.hpp>
#include <dtl/bitmap/util/two_state_markov_process.hpp>
#include <dtl/bitmap/position_list.hpp>
#include <dtl/bitmap/partitioned_position_list.hpp>
#include <dtl/bitmap/range_list.hpp>
#include <dtl/bitmap/partitioned_range_list.hpp>

// The number of independent runs.
static constexpr u64 RUNS = 10;
static constexpr u64 N = 1u << 20;

static const i64 RUN_ID = std::chrono::duration_cast<std::chrono::seconds>(
    std::chrono::system_clock::now().time_since_epoch()).count();

// Read the CPU affinity for the process.
static const auto cpu_mask = dtl::this_thread::get_cpu_affinity();


enum class bitmap_t {
  bitmap,
  roaring,
  teb,
  wah,
  position_list,
  partitioned_position_list_u8,
  partitioned_position_list_u16,
  range_list,
  partitioned_range_list_u8,
  partitioned_range_list_u16
};

struct config {
  bitmap_t bitmap_type;
  $f64 density;
  $f64 clustering_factor;
};


template<typename T>
void run($f64 f, $f64 d, std::ostream& os) {

  // Construct a random bitmap.
  f64 f_min = d >= 1.0 ? N : d/(1-d);
  f64 f_actual = std::max(f, f_min);
  two_state_markov_process mp(f_actual, d);
  boost::dynamic_bitset<$u32> bs(N);
  for ($u64 i = 0; i < N; i++) {
    bs[i] = mp.next();
  }

  $f64 d_actual = (bs.count() * 1.0) / N;
  if (std::abs(d - d_actual) > 1
      || std::abs(f - f_actual) > 0.25) {
    return;
  }

  // Encode the bitmap.
  T enc_bs(bs);

  const auto size_in_bytes = enc_bs.size_in_byte();

  std::string type_info = enc_bs.info();
  boost::replace_all(type_info, "\"", "\"\""); // Escape JSON for CSV output.

  os << RUN_ID
     << "," << N
     << "," << T::name()
     << "," << d
     << "," << d_actual
     << "," << f
     << "," << f_actual
     << "," << size_in_bytes
     << "," << "\"" << type_info << "\""
     << std::endl;
}

void run(config c, std::ostream& os) {
  switch (c.bitmap_type) {
    case bitmap_t::bitmap:
      run<dtl::dynamic_bitmap<$u32>>(c.clustering_factor, c.density, os);
      break;
    case bitmap_t::roaring:
      run<dtl::dynamic_roaring_bitmap>(c.clustering_factor, c.density, os);
      break;
    case bitmap_t::teb:
      run<dtl::teb<>>(c.clustering_factor, c.density, os);
      break;
    case bitmap_t::wah:
      run<dtl::dynamic_wah32>(c.clustering_factor, c.density, os);
      break;
    case bitmap_t::position_list:
      run<dtl::position_list<$u32>>(c.clustering_factor, c.density, os);
      break;
    case bitmap_t::partitioned_position_list_u8:
      run<dtl::partitioned_position_list<$u32, $u8>>(
          c.clustering_factor, c.density, os);
      break;
    case bitmap_t::partitioned_position_list_u16:
      run<dtl::partitioned_position_list<$u32, $u16>>(
          c.clustering_factor, c.density, os);
      break;
    case bitmap_t::range_list:
      run<dtl::range_list<$u32>>(c.clustering_factor, c.density, os);
      break;
    case bitmap_t::partitioned_range_list_u8:
      run<dtl::partitioned_range_list<$u32, $u8>>(
          c.clustering_factor, c.density, os);
      break;
    case bitmap_t::partitioned_range_list_u16:
      run<dtl::partitioned_range_list<$u32, $u16>>(
          c.clustering_factor, c.density, os);
      break;
  }
}


template<typename T>
void dispatch(const std::vector<T>& tasks,
              std::function<void(const T&, std::ostream&)> fn) {

  i64 thread_cnt = dtl::env<$u64>::get("THREAD_CNT", cpu_mask.count());
  i64 config_cnt = tasks.size();
  i64 min_batch_size = 1;
  i64 max_batch_size = 16;

  const auto time_start = std::chrono::system_clock::now();
  std::atomic<$i64> cntr { 0 };
  auto thread_fn = [&](u32 thread_id) {
    while (true) {
      // Grab work.
      const auto inc = std::min(std::max(min_batch_size, (config_cnt - cntr) / thread_cnt), max_batch_size);
      std::stringstream s;
      s << "thread " << thread_id << " got " << inc << " task(s)" << std::endl;
      std::cerr << s.str();
      const auto config_idx_begin = cntr.fetch_add(inc);
      const auto config_idx_end = std::min(config_idx_begin + inc, config_cnt);
      if (config_idx_begin >= config_cnt) break;

      std::stringstream str;
      for ($i64 ci = config_idx_begin; ci < config_idx_end; ci++) {
        fn(tasks[ci], str);
      }
      std::cout << str.str();

      if (thread_id == 0) {
        i64 i = cntr;
        i64 r = std::min(config_cnt, config_cnt - i);
        // Estimate time until completion.
        const auto now = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = now - time_start;
        f64 avg_sec_per_config = elapsed_seconds.count() / i;
        u64 remaining_sec = avg_sec_per_config * r;
        u64 h = (remaining_sec / 3600);
        u64 m = (remaining_sec % 3600) / 60;
        std::stringstream str;
        str << "Progress: [" << (i + 1) << "/" << config_cnt << "]";
        str << " - estimated time until completion: " << h << "h " << m << "m" << std::endl;
        std::cerr << str.str();
      }
    }
  };
  dtl::run_in_parallel(thread_fn, cpu_mask, cpu_mask.count());
}

