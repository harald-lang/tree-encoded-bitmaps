#pragma once

#include "dtl/index.hpp"
#include "dtl/sma.hpp"
#include "dtl/tree_mask.hpp"
#include "dtl/zone_mask.hpp"
#include <algorithm>
#include <cstring>
#include <type_traits>
#include <vector>
#include "memory"
#include "dtl/storage.hpp"

#include <immintrin.h>
#include "histograms.hpp"
#include "dtl/mem_info.hpp"

namespace h_psmas {

/// A PSMA lookup table for the type T. Each table entry consists of an instance of V (e.g., a range in the
/// default implementation).
template<typename T, typename V, u32 BINS_CNT, bool EQ_W, bool BLOCK_WISE>
class h_psma_table {

public:

  histo::histogram<T,BINS_CNT,EQ_W,BLOCK_WISE> h;
  V entries[BINS_CNT];

  // compute the PSMA slot for a given value
  inline u64
  get_slot(const T value) const noexcept {
    return h.get_slot(value);
  }

};

/// PSMA implementation based on the paper of Lang et al. 'Data Blocks: Hybrid OLTP and OLAP on Compressed Storage
/// using both Vectorization and Compilation'
template<typename T, u64 N, u32 BINS_CNT, bool EQ_W, bool BLOCK_WISE>
class h_psma {
public:
  using value_t = typename std::remove_cv<T>::type;

  using sma_t = dtl::sma<value_t>;
  using table_t = h_psma_table<T, dtl::range<N>, BINS_CNT, EQ_W, BLOCK_WISE>;

  static constexpr u32 size = BINS_CNT;

  sma_t _sma;
  table_t table;

  struct psma_builder {
    h_psma& ref;

    inline void
    operator()(const T *const values, const size_t n, std::function<u1(u64)> is_null) noexcept {

      // requires three passes
      // 1. build the SMA
      auto sma_build = ref._sma.builder();
      sma_build(values, n, is_null);
      sma_build.done();

      // 2. create histogram
      histo::samples<T> smpl(values, n);
      histo::histogram<T,BINS_CNT,EQ_W,BLOCK_WISE> hist(smpl);
      ref.table.h = hist;

      // 3. populate the PSMA table
      for ($u32 i = 0; i != n; i++) {
        u64 slot_id = ref.table.get_slot(values[i]);
        auto& range = ref.table.entries[slot_id];

        if (range.is_empty())
          range = {i, i + 1};
        else
          range.end = i + 1;
      }

    }

    inline void
    done() {
      // nothing to do here.
      // builder performs in-place updates
    }

  };

  inline psma_builder
  builder() {
    // reset table entries
    for ($u32 i = 0; i != size; i++) {
      table.entries[i].reset();
    }
    // return a builder instance
    return psma_builder { *this };
  }

  // c'tor
  h_psma() noexcept {
    // initialize the lookup table with empty range
    for ($u32 i = 0; i != size; i++) {
      table.entries[i].reset();
    }
    table.h.reset();
  }

  inline dtl::range<N>
  lookup(const dtl::predicate& p) const noexcept {

    value_t value = *reinterpret_cast<value_t*>(p.value_ptr);
    value_t second_value; // in case of between predicates

    $u32 s = table.get_slot(value);
    dtl::range<N> r;
    r.reset();

    if(value < _sma.min_value){
      switch (p.comparison_operator){
        case dtl::op::EQ:
        case dtl::op::LE:
        case dtl::op::LT:
          return r;
        case dtl::op::GT:
        case dtl::op::GE:
          r.begin = 0;
          r.end = size;
          return r;
      }
    }

    if(value > _sma.max_value){
      switch (p.comparison_operator){
        case dtl::op::EQ:
        case dtl::op::GT:
        case dtl::op::GE:
          return r;
        case dtl::op::LE:
        case dtl::op::LT:
          r.begin = 0;
          r.end = size;
          return r;
      }
    }

    $u32 b = 0;
    $u32 e = 0;
    switch (p.comparison_operator) {
      case dtl::op::EQ:
        return table.entries[s];

      case dtl::op::LT:
        if(value > _sma.min_value) // handle the case where value == min of the histogram
          e = table.get_slot(value - 1); // get slot of the next smaller value
        else
          return r;
        break;

      case dtl::op::LE:
        e = s;
        break;

      case dtl::op::GT:
        if(value >= _sma.max_value) // the searched value is the largest possible value => no larger values
          return r;

        b = table.get_slot(value + 1); // get slot of the next larger value
        e = size - 1;
        break;

      case dtl::op::GE:
        if(value < _sma.min_value )
          b = 0;
        else
          b = s;

        e = size - 1;
        break;

      case dtl::op::BETWEEN:
        if(value > _sma.max_value)
          return r;

        second_value = *reinterpret_cast<value_t*>(p.second_value_ptr);
        if(second_value < _sma.min_value)
          return r;

        if(value < _sma.min_value)
          b = 0;
        else
          b = s;

        e = table.get_slot(second_value);
        break;

      case dtl::op::BETWEEN_LO: // (x,y]
        if(value >= _sma.max_value) // the searched value is >= the largest possible value => no larger values
          return r;

        second_value = *reinterpret_cast<value_t*>(p.second_value_ptr);
        if(second_value < _sma.min_value) // upper bound smaller than smallest value
          return r;

        if(value < _sma.min_value) // first value smaller than minimum
          b = 0;
        else
          b = table.get_slot(value + 1); //exclude x -> [x+1;y)

        e = table.get_slot(second_value);
        break;

      case dtl::op::BETWEEN_RO: // [x,y)
        if(value > _sma.max_value)
          return r;

        second_value = *reinterpret_cast<value_t*>(p.second_value_ptr);
        if(second_value <= _sma.min_value) // upper bound smaller than smallest value
          return r;

        if(value < _sma.min_value) // first value smaller than minimum
          b = 0;
        else
          b = s;

        e = table.get_slot(second_value - 1); // exclude y -> [x;y-1]
        break;

      case dtl::op::BETWEEN_O: // (x,y)
        if(value >= _sma.max_value)
          return r;

        second_value = *reinterpret_cast<value_t*>(p.second_value_ptr);
        if(second_value <= _sma.min_value)
          return r;

        if(value < _sma.min_value)
          b = 0;
        else
          b = table.get_slot(value + 1); // exclude x -> [x+1;y)

        e = table.get_slot(second_value - 1); // exclude y -> [x+1;y-1]
        break;
    }

    for (size_t i = b; i <= e; i++) {
      r = r | table.entries[i];
    }
    return r;
  }

  /// Memory footprint of the index structure in bytes
  inline dtl::mem_info memory_footprint(){
    auto bytes_pointer = 0;
    auto bytes_ranges = sizeof(table.entries);
    auto bytes_sma = sizeof(T) * 2 + 1; // boundaries + has_null_values
    auto local_var = sizeof(size) + bytes_sma; //size-variable + sma
    auto hist = table.h.memory_footprint();

    struct dtl::mem_info memory("Histogram_PSMA",0,bytes_pointer,bytes_ranges, local_var, hist);

    return memory;
  }

  inline void
  print() noexcept {
    for(auto i = 0; i < BINS_CNT; i++)
      std::cout << "[" << table.h.bin_min[i] << ";" << table.h.bin_max[i] << ") "
                << ": [" << table.entries[i].begin << ";" << table.entries[i].end << ")" << std::endl;
  }
};


/// A combination of a PSMA lookup table and a Zone Mask.
/// M = the number of bits per table entry.
template<typename T, u64 N, u64 M, u32 BINS_CNT, bool EQ_W, bool BLOCK_WISE>
class h_psma_zone_mask {
public:

  using value_t = typename std::remove_cv<T>::type;

  using mask_t = dtl::zone_mask<N, M>;

  using sma_t = dtl::sma<value_t>;
  using table_t = h_psma_table<T, mask_t, BINS_CNT, EQ_W, BLOCK_WISE>;

  static constexpr u32 size = BINS_CNT;

  sma_t _sma;
  table_t table;

  struct psma_builder {
    h_psma_zone_mask& ref;

    inline void
    operator()(const T *const values, const size_t n, std::function<u1(u64)> is_null) noexcept {
      // requires three passes
      // 1. build the SMA
      auto sma_build = ref._sma.builder();
      sma_build(values, n, is_null);
      sma_build.done();

      // 2. create histogram
      histo::samples<T> smpl(values, n);
      histo::histogram<T,BINS_CNT,EQ_W,BLOCK_WISE> hist(smpl);
      ref.table.h = hist;



      // 3. build populate PSMA table
      for ($u32 i = 0; i != n; i++) {
        u64 slot_id = ref.table.get_slot(values[i]);
        auto& entry = ref.table.entries[slot_id];
        entry.set(i);
        //cout << i << " Value: " << to_string(values[i]) << " Slot: " << slot_id << " Entry: " << entry.get_data() << endl;
      }
    }

    inline void
    done() {
      // nothing to do here.
      // builder performs in-place updates
    }

  };

  inline psma_builder
  builder() {
    // reset table entries
    for (uint32_t i = 0; i != size; i++) {
      table.entries[i].reset();
    }
    // return a builder instance
    return psma_builder { *this };
  }


//  inline void
//  update(const T* const values, const size_t n) noexcept {
//    for (uint32_t i = 0; i != n; i++) {
//      auto& entry = table.entries[table.get_slot(values[i])];
//      entry.set(i);
//    }
//  }

  inline mask_t
  lookup(const dtl::predicate& p) const noexcept {

    value_t value = *reinterpret_cast<value_t*>(p.value_ptr);
    value_t second_value; // in case of between predicates

    $u32 s = table.get_slot(value);
    mask_t r;
    r.reset();

    if(value < _sma.min_value){
      switch (p.comparison_operator) {
        case dtl::op::EQ:
        case dtl::op::LE:
        case dtl::op::LT:
          return r;
        case dtl::op::GT:
        case dtl::op::GE:
          r.set();
          return r;
      }
    }

    if(value > _sma.max_value){
      switch (p.comparison_operator){
        case dtl::op::EQ:
        case dtl::op::GT:
        case dtl::op::GE:
          return r;
        case dtl::op::LE:
        case dtl::op::LT:
          r.set();
          return r;
      }
    }

    $u32 b = 0;
    $u32 e = 0;
    switch (p.comparison_operator) {
      case dtl::op::EQ:
        return table.entries[s];

      case dtl::op::LT:
        if(value > _sma.min_value) // handle the case where value == min of the histogram
          e = table.get_slot(value - 1); //get slot of the next smaller value
        else
          return r;
        break;

      case dtl::op::LE:
        e = s;
        break;

      case dtl::op::GT:
        if(value >= _sma.max_value)  // the searched value is the largest value one -> no larger values
          return r;

        b = table.get_slot(value + 1); // get slot of the next larger value
        e = size - 1;
        break;

      case dtl::op::GE:
        if(value < _sma.min_value )
          b = 0;
        else
          b = s;

        e = size - 1;
        break;

      case dtl::op::BETWEEN:
        if(value > _sma.max_value)
          return r;

        second_value = *reinterpret_cast<value_t*>(p.second_value_ptr);
        if(second_value < _sma.min_value)
          return r;

        if(value < _sma.min_value)
          b = 0;
        else
          b = s;

        e = table.get_slot(second_value);
        break;

      case dtl::op::BETWEEN_LO: // (x,y]
        if(value >= _sma.max_value)  // the searched value is the largest possible value => no larger values
          return r;

        second_value = *reinterpret_cast<value_t*>(p.second_value_ptr);
        if(second_value < _sma.min_value)       // upper bound smaller than smallest value
          return r;

        if(value < _sma.min_value)              // first value smaller than minimum
          b = 0;
        else
          b = table.get_slot(value + 1);        // exclude x -> [x+1;y)

        e = table.get_slot(second_value);
        break;

      case dtl::op::BETWEEN_RO: // [x,y)
        if(value > _sma.max_value) // upper bound smaller than smallest value
          return r;

        second_value = *reinterpret_cast<value_t*>(p.second_value_ptr);
        if(second_value <= _sma.min_value) // upper bound smaller than smallest value
          return r;

        if(value < _sma.min_value) // first value smaller than minimum
          b = 0;
        else
          b = s;

        e = table.get_slot(second_value - 1); // exclude y -> [x;y-1]
        break;

      case dtl::op::BETWEEN_O: // (x,y)
        if(value >= _sma.max_value)
          return r;

        second_value = *reinterpret_cast<value_t*>(p.second_value_ptr);
        if(second_value <= _sma.min_value)
          return r;

        if(value < _sma.min_value)
          b = 0;
        else
          b = table.get_slot(value + 1); // exclude x -> [x+1;y)

        e = table.get_slot(second_value - 1); // exclude y -> [x+1;y-1]
        break;
    }

    for (size_t i = b; i <= e; i++) {
      r = r | table.entries[i];
    }
    return r;
  }

  /// Memory footprint of the index structure in bytes
  inline dtl::mem_info memory_footprint(){
    auto bytes_pointer = 0;
    auto bytes_zonemasks = size * sizeof(mask_t);
    auto bytes_sma = sizeof(T) * 2 + 1; // boundaries + has_null_values
    auto local_var = sizeof(size) + bytes_sma; // size-variable + sma
    auto hist = table.h.memory_footprint();

    std::string idx;
    if(N != M)
      idx = "Histogram Zonemask";
    else
      idx = "Histogram Bitmask";

    struct dtl::mem_info memory(idx, M, bytes_pointer,bytes_zonemasks, local_var, hist);

    return memory;
  }

  inline void
  print() noexcept {
    for(auto i = 0; i < BINS_CNT; i++)
      std::cout << "[" << table.h.bin_min[i] << ";" << table.h.bin_max[i] << ") "
                << ": " << table.entries[i].data << std::endl;
  }
};


template<typename T, u64 N, u32 BINS_CNT, bool EQ_W, bool BLOCK_WISE>
using  h_psma_bitmask = h_psma_zone_mask<T, N, N, BINS_CNT, EQ_W, BLOCK_WISE>;


/// A combination of a PSMA lookup table and a Zone Mask.
/// N = the number of tuples of the data_block
/// M = the number of bits per table entry.
/// L = true: loss_less compression false: lossy compression
template<typename T, u64 N, u64 M, bool L, u32 BINS_CNT, bool EQ_W, bool BLOCK_WISE>
class h_psma_tree_mask {
public:
  using value_t = typename std::remove_cv<T>::type;

  using sma_t = dtl::sma<value_t>;
  using mask_t = dtl::tree_mask<N, M>;
  using table_t = h_psma_table<T, mask_t, BINS_CNT, EQ_W, BLOCK_WISE>;

  static constexpr u32 size = BINS_CNT;

  sma_t _sma;
  table_t table;

  struct psma_builder {
    h_psma_tree_mask& ref;

    inline void
    operator()(const value_t *const values, const size_t n, std::function<u1(u64)> is_null) noexcept {
      // built on top of psma_bitmap
      auto bm = std::make_unique<h_psma_bitmask<value_t, N, BINS_CNT, EQ_W, BLOCK_WISE>>();
      auto bm_build = bm->builder();
      bm_build(values, n, is_null);
      ref.table.h = bm_build.ref.table.h;
      bm_build.done();

      // copy SMA
      ref._sma = bm->_sma;

      // encode bitmaps to tree masks
      for ($u32 i = 0; i != size; i++) {
        if (L){
          *(ref.table.entries[i].data) =
              dtl::tree_mask<N, M>::encode(bm->table.entries[i].data);
        } else
          *(ref.table.entries[i].data) = dtl::tree_mask<N,M>::compress(bm->table.entries[i].data);

      }

    }

    inline void
    done() {
      // nothing to do here.
      // builder performs in-place updates
    }

  };

  inline psma_builder
  builder() {
    // reset table entries
    for ($u32 i = 0; i != size; i++) {
      table.entries[i].reset();
    }
    // return a builder instance
    return psma_builder { *this };
  }

//  inline void
//  update(const psma_bitmask<T, N>& src) noexcept {
//    for (uint32_t i = 0; i != table_t::size; i++) {
//      table.entries[i].set(src.table.entries[i].data);
//    }
//  }

  inline std::bitset<N>
  lookup(const dtl::predicate& p) const noexcept {

    value_t value = *reinterpret_cast<value_t*>(p.value_ptr);
    std::bitset<N> r;
    value_t second_value; // in case of between predicates
    $u32 s = table.get_slot(value);

    if(value < _sma.min_value){
      switch (p.comparison_operator) {
        case dtl::op::EQ:
        case dtl::op::LE:
        case dtl::op::LT:
          return r;
        case dtl::op::GT:
        case dtl::op::GE:
          r.set();
          return r;
      }
    }

    if(value > _sma.max_value){
      switch (p.comparison_operator){
        case dtl::op::EQ:
        case dtl::op::GT:
        case dtl::op::GE:
          return r;
        case dtl::op::LE:
        case dtl::op::LT:
          r.set();
          return r;
      }
    }

    $u32 b = 0;
    $u32 e = 0;
    switch (p.comparison_operator) {
      case dtl::op::EQ:
        return table.entries[s].get();

      case dtl::op::LT:
        if(value > _sma.min_value)
          e = table.get_slot(value - 1); // get slot of the next smaller value
        else
          return r;
        break;

      case dtl::op::LE:
        e = s;
        break;

      case dtl::op::GT:
        if(value >= _sma.max_value) // the searched value is the largest value one -> no larger values
          return r;

        b = table.get_slot(value + 1); // get slot of the next larger value
        e = size - 1;
        break;

      case dtl::op::GE:
        if(value < _sma.min_value )
          b = 0;
        else
          b = s;

        e = size - 1;
        break;

      case dtl::op::BETWEEN:
        if(value > _sma.max_value)
          return r;

        second_value = *reinterpret_cast<value_t*>(p.second_value_ptr);
        if(second_value < _sma.min_value)
          return r;

        if(value < _sma.min_value)
          b = 0;
        else
          b = s;

        e = table.get_slot(second_value);
        break;

      case dtl::op::BETWEEN_LO: // (x,y]
        if(value >= _sma.max_value) // the searched value is the largest possible value => no larger values
          return r;

        second_value = *reinterpret_cast<value_t*>(p.second_value_ptr);
        if(second_value < _sma.min_value) // upper bound smaller than smallest value
          return r;

        if(value < _sma.min_value) // first value smaller than minimum
          b = 0;
        else
          b = table.get_slot(value + 1); // exclude x -> [x+1;y)

        e = table.get_slot(second_value);
        break;

      case dtl::op::BETWEEN_RO: // [x,y)
        if(value > _sma.max_value) // upper bound smaller than smallest value
          return r;

        second_value = *reinterpret_cast<value_t*>(p.second_value_ptr);
        if(second_value <= _sma.min_value) // upper bound smaller than smallest value
          return r;

        if(value < _sma.min_value) // first value smaller than minimum
          b = 0;
        else
          b = s;

        e = table.get_slot(second_value - 1); // exclude y -> [x;y-1]
        break;

      case dtl::op::BETWEEN_O: // (x,y)
        if(value >= _sma.max_value)
          return r;

        second_value = *reinterpret_cast<value_t*>(p.second_value_ptr);
        if(second_value <= _sma.min_value)
          return r;

        if(value < _sma.min_value)
          b = 0;
        else
          b = table.get_slot(value + 1); // exclude x -> [x+1;y)

        e = table.get_slot(second_value - 1); // exclude y -> [x+1;y-1]
        break;
    }

    for (size_t i = b; i <= e; i++) {
      r = r | table.entries[i].get();
    }
    return r;
  }

  /// Memory footprint of the index structure in bytes
  inline dtl::mem_info memory_footprint(){
    auto bytes_pointer = size*sizeof(mask_t);
    auto bytes_sma = sizeof(T)*2 + 1; //boundaries + has_null_values
    auto local_var = sizeof(size) + bytes_sma; //size-variable + sma
    auto bytes_treemask = 0;
    auto hist = table.h.memory_footprint();
    $u32 mask_size;

    for(auto e: table.entries){
      bytes_treemask += e.size_in_bytes();
    }

    if(L)
      mask_size = N; // for loos-less treemasks the masksize is the blocksize
    else
      mask_size = M; // lossy compressed treemasks of size M

    struct dtl::mem_info memory("Histogram Treemask", mask_size, bytes_pointer,bytes_treemask, local_var, hist);

    return memory;
  }

  inline void
  print() noexcept {
    for(auto i = 0; i < BINS_CNT; i++){
      std::cout << "[" << table.h.bin_min[i] << ";" << table.h.bin_max[i] << ")" << " : ";
      table.entries[i].print_decoded();
      std::cout << std::endl;
    }
  }

};

//template<typename T, u64 N, u64 M, u32 BINS_CNT, bool EQ_W, bool BLOCK_WISE>
//using  psma_tree_mask = psma_tree_mask<T, N, M, false, BINS_CNT, EQ_W, BLOCK_WISE>;
} // namespace h_psmas