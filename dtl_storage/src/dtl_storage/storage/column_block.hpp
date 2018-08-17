#pragma once

#ifndef _DTL_STORAGE_INCLUDED
#error "Never use <dtl/storage/column_block.hpp> directly; include <dtl/storage.hpp> instead."
#endif

#include <memory>

#include <dtl/dtl.hpp>
#include <dtl/bitset.hpp>
#include <dtl/storage/types.hpp>


namespace dtl {

const u16 file_signature = 0x484c;


/// A naive column-block implementation which does not support deletions
template<u64 N>
struct column_block_base {
  static constexpr u64 max_size = N;

  $u64 write_pos;
  alignas(64) dtl::bitset<N> null;

  virtual ~column_block_base() {};

  inline u64
  size() const {
    return write_pos;
  }

  inline void
  push_back(const dtl::null /* val */) {
    null[write_pos] = true;
    write_pos++;
  }

  inline u1
  is_null(u64 n) const {
    return null[n];
  }

};


template<typename T, u64 N>
struct column_block : public column_block_base<N> {
  using column_block_base<N>::write_pos;
  using column_block_base<N>::is_null;
  using column_block_base<N>::push_back;

  alignas(64) std::array<T, N> data;

  inline void
  push_back(const T& val) {
    data[write_pos++] = val;
  }

  inline void
  push_back(T&& val) {
    data[write_pos++] = std::move(val);
  }

  inline T&
  operator[](u64 n) {
    return data[n];
  }

  inline const T&
  operator[](u64 n) const {
    return data[n];
  }

};


/// construct a block (high level function, intended to be used in interpreted code)
template<u64 N>
dtl::column_block_base<N>*
make_column_block(const dtl::rtt type) {
  dtl::column_block_base<N>* ptr;
  switch (type) {
#define DTL_GENERATE(T) \
    case dtl::rtt::T:  \
      ptr = new dtl::column_block<dtl::map<dtl::rtt::T>::type, N>(); \
      break;
    DTL_GENERATE(u8)
    DTL_GENERATE(i8)
    DTL_GENERATE(u16)
    DTL_GENERATE(i16)
    DTL_GENERATE(u32)
    DTL_GENERATE(i32)
    DTL_GENERATE(u64)
    DTL_GENERATE(i64)
    DTL_GENERATE(f32)
    DTL_GENERATE(f64)
    DTL_GENERATE(str)
#undef DTL_GENERATE
  }
  assert(ptr != nullptr);
  return ptr;
}


/// insert data into a block (high level function for ingestion, intended to be used in interpreted code)
template<u64 N>
inline void
column_block_insert(dtl::column_block_base<N>* block_ptr,
                    const dtl::rtt type,
                    const std::string& value,
                    const std::string& null_indicator) {
  switch (type) {
#define DTL_GENERATE(T) \
    case dtl::rtt::T: {                                                                \
      using block_type = dtl::column_block<dtl::map<dtl::rtt::T>::type, N>;            \
      block_type* b = static_cast<block_type*>(block_ptr);                             \
      if (value == null_indicator) {                                                   \
        b->push_back(dtl::null::value);                                                \
      }                                                                                \
      else {                                                                           \
        const auto parse = dtl::parse<dtl::rtt::T>();                                  \
        b->push_back(parse(value));                                                    \
      }                                                                                \
      break;                                                                           \
    }
    DTL_GENERATE(u8)
    DTL_GENERATE(i8)
    DTL_GENERATE(u16)
    DTL_GENERATE(i16)
    DTL_GENERATE(u32)
    DTL_GENERATE(i32)
    DTL_GENERATE(u64)
    DTL_GENERATE(i64)
    DTL_GENERATE(f32)
    DTL_GENERATE(f64)
    DTL_GENERATE(str)
#undef DTL_GENERATE
  }
}


} // namespace dtl