#pragma once

#ifndef _DTL_STORAGE_INCLUDED
#error "Never use <dtl/storage/relation.hpp> directly; include <dtl/storage.hpp> instead."
#endif

#include <array>
#include <atomic>
#include <mutex>

#include <dtl/dtl.hpp>
#include <dtl/lockstripe.hpp>
#include <dtl/mem.hpp>
#include <dtl/storage/column_block.hpp>
#include <dtl/storage/schema.hpp>


namespace dtl {


template<typename ColumnBlock>
struct partition {

  using allocator_t = std::allocator<ColumnBlock>; //Alloc;
  using column_block_t = ColumnBlock;
  using column_block_ptr_t = ColumnBlock*;
  using block_t = std::vector<column_block_ptr_t>;

  static constexpr u64 block_size = column_block_t::max_size;

  const allocator_t allocator;
  const dtl::schema& schema;

  /// locks for concurrent writers
  static constexpr u64 max_concurrent_writer = 8;
  std::unique_ptr<dtl::lockstripe<max_concurrent_writer>> lockstripe;

  /// short-term lock (for the entire partition)
  std::mutex partition_latch;

  /// store column blocks PAX-alike
  std::vector<block_t> blocks;

  /// points to the last block
  $u64 tail_block_id;

  /// points to the first underutilized block (writer should use this block for space efficiency reasons)
  $u64 writeable_block_id;

  /// creates a new data block consisting of multiple column blocks
  block_t
  create_block() const {
    block_t column_blocks;
    column_blocks.reserve(schema.size());
    for (auto& attr : schema) {
      column_blocks.push_back(dtl::make_column_block<block_size>(attr.type)); // TODO use allocator
    }
    return column_blocks;
  }

  /// c'tor
  partition(const dtl::schema& schema,
            const dtl::mem::allocator_config allocator_config)
      : allocator(allocator),
        schema(schema),
        lockstripe(std::make_unique<dtl::lockstripe<max_concurrent_writer>>()),
        tail_block_id(0),
        writeable_block_id(0) {
    blocks.push_back(create_block());
  }


  std::vector<column_block_ptr_t>
  acquire_block() {
    std::vector<column_block_ptr_t> blocks;
    blocks.reserve(schema.size());
    for (auto& attr : schema) {
      blocks.push_back(dtl::make_column_block<block_size>(attr.type));
    }
    return blocks;
  }


  void
  release_block(std::vector<column_block_ptr_t> blocks) {

  };

};

struct relation {

  static constexpr u64 block_size = 1ull << 16;
  using column_block_t = dtl::column_block_base<block_size>;
  using column_block_ptr_t = column_block_t*;
  using partition_t = partition<column_block_t>;

  const dtl::schema schema;

  /// the number of partitions. should be somewhere between #numa_nodes and #hardware_threads
  u64 partition_cnt;

  /// references to all blocks of this column
  std::vector<partition_t> partitions;

  /// c'tor
  explicit
  relation(const dtl::schema& schema)
      : schema(schema),
        partition_cnt(dtl::mem::get_cpu_nodes().size() * 8) {
    // create partitions
    partitions.reserve(partition_cnt);
    const auto numa_nodes = dtl::mem::get_cpu_nodes();
    const auto numa_node_cnt = numa_nodes.size();
    for ($u64 i = 0; i < partition_cnt; i++) {
      const auto alloc = dtl::mem::allocator_config::on_node(numa_nodes[i % numa_node_cnt]);
      partitions.emplace_back(schema, alloc);
    }
  }


  /// exclusively acquires a data block for writing
  std::vector<column_block_ptr_t>
  acquire_block() {
    std::vector<column_block_ptr_t> blocks;
    blocks.reserve(schema.size());
    for (auto& attr : schema) {
      blocks.push_back(dtl::make_column_block<block_size>(attr.type));
    }
    return blocks;
  }


  void
  release_block(std::vector<column_block_ptr_t> blocks) {

  };

};


} // namespace dtl