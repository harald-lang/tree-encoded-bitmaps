#pragma once

#include <iostream>
#include <fstream>

#include <dtl/adept.hpp>
#include <dtl/storage.hpp>
#include <dtl/thread.hpp>
#include <dtl/env.hpp>
#include <dtl/psma.hpp>

#include <tbb/concurrent_queue.h>
#include <tbb/concurrent_vector.h>
#include <sys/resource.h>

#include "dtl/storage/file_parser.hpp"
#include "dtl/storage/h_psma.hpp"


namespace tests{

static dtl::schema ontime_schema {
    {"year", dtl::rtt::i16},
    {"month", dtl::rtt::i8},
    {"day_of_month", dtl::rtt::i8},
    {"day_of_week", dtl::rtt::i8},
    {"dep_time", dtl::rtt::i16},
    {"scheduled_dep_time", dtl::rtt::i16},
    {"arr_time", dtl::rtt::i16},
    {"scheduled_arr_time", dtl::rtt::i16},
    {"carrier", dtl::rtt::str},
    {"flight_num", dtl::rtt::i16},
    {"tail_num", dtl::rtt::str},
    {"actual_elapsed_time", dtl::rtt::i16},
    {"csr_elapsed_time", dtl::rtt::i16},
    {"arr_time", dtl::rtt::i16},
    {"arr_delay", dtl::rtt::i16},
    {"dep_delay", dtl::rtt::i16},
    {"origin", dtl::rtt::str},
    {"dest", dtl::rtt::str},
    {"distance", dtl::rtt::u16},
    {"taxi_in", dtl::rtt::u16},
    {"taxi_out", dtl::rtt::u16},
    {"cancelled", dtl::rtt::u8},
    {"cancellation_code", dtl::rtt::str},
    {"diverted", dtl::rtt::u8},
    {"carrier_delay", dtl::rtt::i16},
    {"weather_delay", dtl::rtt::i16},
    {"nas_delay", dtl::rtt::i16},
    {"security_delay", dtl::rtt::i16},
    {"late_aircraft_delay", dtl::rtt::i16}
};


struct result{
  std::vector<dtl::mem_info> mem_cons_treemasks;

  result(){
  }

  result(u8 n){
    mem_cons_treemasks.resize(n);
  }

  void set($u32 pos, const dtl::mem_info& m) {
    if(pos < mem_cons_treemasks.size())
      mem_cons_treemasks[pos] = m;
    else
      std::cout << "Elem is to large! Elem: " << pos << "Size vector: " << mem_cons_treemasks.size() << std::endl;
  };

  dtl::mem_info agg_result(const std::vector<dtl::mem_info>& r){
    std::string idx = "";
    $u64 mask_size = 0;
    $u64 size_pointer = 0;
    $u64 table_bitvectors = 0;
    $u64 size_local_var = 0;

    for(auto e : r){
      size_pointer += e.size_pointer;
      table_bitvectors += e.table_bitvectors;
      size_local_var += e.size_local_var;
    }

    if(r.size() > 0){
      idx = r[0].idx;
      mask_size = r[0].mask_size;
    }

    return dtl::mem_info(idx, mask_size, size_pointer, table_bitvectors, size_local_var);
  }

  void agg_result(const std::vector<result>& v){

    if(mem_cons_treemasks.size() != v[0].mem_cons_treemasks.size())
      std::cout << "Number of elements of the aggregated result and the blocks doesn't match! #agg: "
                << mem_cons_treemasks.size() << "#in blocks: " << v[0].mem_cons_treemasks.size() << std::endl;

    for(auto e: mem_cons_treemasks)
      e.reset();

    for(result vec : v){
      for(auto i = 0; i < vec.mem_cons_treemasks.size(); i++) {
        mem_cons_treemasks[i] = mem_cons_treemasks[i] + vec.mem_cons_treemasks[i];
      }
    }
  }
};

template<u64 block_size>
struct task {
  std::vector<std::string> csv;
};

template<typename T, u64 N>
histo::samples<T> sample_of_column(const auto& column){

  using block_type = dtl::column_block<T, N>;
  const block_type* column_block = reinterpret_cast<const block_type*>(column);
  auto block = column_block->data.data();

  return histo::samples<T>(block, column_block->size());
};


template<typename T, u64 N>
std::vector<histo::samples<T>> col_distr(const auto& blocks){

  std::vector<histo::samples<T>> distr;

  for(auto i = 0; i < ontime_schema.size(); i++) {
    std::cout << "Scanning Column: " << ontime_schema[i].name << std::endl;
    const auto& column = blocks[i];
    auto type = ontime_schema[i].type;

    switch(type){
      case dtl::rtt::u8:
        sample_of_column<u8,N>(column);
        break;
      /*
      case dtl::rtt::u16:
        sample_of_column<u16,N>(column);
        break;
      case dtl::rtt::u32:
        sample_of_column<u32,N>(column);
        break;
      case dtl::rtt::u64:
        sample_of_column<u64,N>(column);
        break;
      case dtl::rtt::i8:
        sample_of_column<i8,N>(column);
        break;
      case dtl::rtt::i16:
        sample_of_column<i16,N>(column);
        break;
      case dtl::rtt::i32:
        sample_of_column<i32,N>(column);
        break;
      case dtl::rtt::i64:
        sample_of_column<i64,N>(column);
        break;
        */
    }


  }

  return distr;
}

/// N = the number of tuples of the data_block
/// M = the number of bits per table entry.
/// L = true: loss_less compression false: lossy compression
template<typename T, u64 N, u64 M, bool L>
result process_block(task<N>* t, u16 column_id) {
  // parse input and create column blocks
  auto blocks = fp::parse_input_tuples<N>(ontime_schema, t->csv);

  const auto& column = blocks[column_id];

  using block_type = dtl::column_block<T, N>;
  const block_type* column_block = reinterpret_cast<const block_type*>(column);
  auto block = column_block->data.data();
  auto type = ontime_schema[column_id].type;

  auto build_index = [&](const auto* block, auto&& idx) {
    auto idx_build = idx.builder();
    idx_build(block->data.data(),
              block->size(),
              [&](u64 i) { return block->is_null(i); });
    idx_build.done();
    idx.print();

    return idx.memory_footprint();
  };

  //col_distr<T,N>(blocks); // TODO

  result r;
/*
  //Treemask-experiments loss_less, here M does not matter
  r.mem_cons_treemasks.push_back(
      build_index(column_block, dtl::psma_tree_mask<T, N, M, true>()));

  //Treemask-experiments lossy
  r.mem_cons_treemasks.push_back(
      build_index(column_block, dtl::psma_tree_mask<T, N, 64>()));
  r.mem_cons_treemasks.push_back(
      build_index(column_block, dtl::psma_tree_mask<T, N, 128>()));
  r.mem_cons_treemasks.push_back(
      build_index(column_block, dtl::psma_tree_mask<T, N, 256>()));
  r.mem_cons_treemasks.push_back(
      build_index(column_block, dtl::psma_tree_mask<T, N, 512>()));
  r.mem_cons_treemasks.push_back(
      build_index(column_block, dtl::psma_tree_mask<T, N, 1024>()));
  r.mem_cons_treemasks.push_back(
      build_index(column_block, dtl::psma_tree_mask<T, N, N/(64/sizeof(T))>())); // cacheline granularity

  r.mem_cons_treemasks.push_back(
      build_index(column_block, h_psmas::h_psma<T, 16, true, true>()));
  r.mem_cons_treemasks.push_back(
      build_index(column_block, h_psmas::h_psma_zone_mask<T,N,M,16,true, true>()));

  r.mem_cons_treemasks.push_back(
      build_index(column_block, h_psmas::h_psma_tree_mask<T,N,M,true,16,true, true>()));
*/
  // delete blocks
  std::cout << "." << std::flush;
  for (auto b : blocks) {
    delete b;
  }

  return r;
};

template<typename T, u64 N, u64 M, bool L>
result construct_treemasks(std::string input_src, u16 column_id){

  std::cout << std::endl << "Start Treemask-Tests" << std::endl;

  using namespace std::chrono_literals;

  // multi threading
  u32 worker_cnt = std::thread::hardware_concurrency() - 2; // main thread does I/O and dispatching
  u64 max_queue_size = worker_cnt * 2;
  tbb::concurrent_queue<task<N>*> task_queue;
  std::atomic<$u1> termination_flag(false);
  std::atomic<$u64> task_cnt_dispatched(0);
  std::atomic<$u64> task_cnt_processed(0);

  tbb::concurrent_vector<result> block_result;

  // worker function
  auto worker_fn = [&]() {
    $u64 local_task_cntr = 0;
    $u64 local_idle_cntr = 0;
    while(!termination_flag.load()) {
      task<N>* t;
      u1 got_work = task_queue.try_pop(t);
      if (!got_work) {
        std::this_thread::sleep_for(1us);
        local_idle_cntr++;
        continue;
      }

      block_result.push_back(process_block<T,N,M,L>(t, column_id));
      local_task_cntr++;
      task_cnt_processed++;
      delete t;
    }

    std::stringstream str;
    str << "thread id: " << std::this_thread::get_id()
        << ", processed: " << local_task_cntr << ", idle: " << local_idle_cntr << " [us])" << std::endl;
    std::cout << str.str();
  };

  // spawn workers
  std::vector<std::thread> workers;

  dtl::run_in_parallel_async(worker_fn, workers, worker_cnt);


  //---------------------------------------- read the input file ----------------------------------------
  std::string input_file = dtl::env<std::string>::get("IN", input_src);
  std::ifstream is(input_file);

  std::string line;
  std::string header("Year,Month,DayofMonth");
  task<N>* t = new task<N>;

  while(std::getline(is, line)) {
    if (line.compare(0, header.length(), header) == 0) continue; // ignore header(s)
    t->csv.push_back(line);
    if (t->csv.size() == N) {
      if (task_queue.unsafe_size() >= max_queue_size) {
        // do it yourself
        std::cout << "*" << std::flush;
        block_result.push_back(process_block<T,N,M,L>(t, column_id));
        task_cnt_processed++;
        delete t;
      }
      else {
        task_queue.push(t);
      }
      task_cnt_dispatched++;
      t = new task<N>;
    }
  }
  if (t->csv.size() >= N) {
    std::cout << "*" << std::flush;
    block_result.push_back(process_block<T,N,M,L>(t, column_id));
    task_cnt_dispatched++;
    task_cnt_processed++;
  }

  // wait for all threads to finish
  std::cout << "dispatched " << task_cnt_dispatched << " tasks." << std::endl;
  while (task_cnt_dispatched.load() != task_cnt_processed.load()) {
    std::this_thread::sleep_for(100ms);
  }
  termination_flag = true;
  std::cout << "processed " << task_cnt_processed << " tasks. done." << std::endl;

  dtl::wait_for_threads(workers);

  std::cout << "Aggregate results" << std::endl;

  std::vector<result> block_res;

  for(auto e: block_result)
    block_res.push_back(e);

  result mem(block_res[0].mem_cons_treemasks.size());
  mem.agg_result(block_res);

  std::cout << std::endl;
  std::cout << "Block_size : " << N << " tuple" << std::endl;
  std::cout << "Column size: " << block_result.size() * N * sizeof(T) << " bytes" << std::endl;


  for(auto e: mem.mem_cons_treemasks)
    e.print();

  return mem;
};
} // namespace tests