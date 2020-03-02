#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/bitmap.hpp>
#include <dtl/bitmap/bitwise_operations.hpp>
#include <dtl/dtl.hpp>
//===----------------------------------------------------------------------===//
template<typename T>
dtl::bitmap
bitwise_and(const T& bitmap_a, const T& bitmap_b) {
  dtl::bitmap ret_val(bitmap_a.size());
  auto it_a = bitmap_a.scan_it();
  auto it_b = bitmap_b.it();
  while (!(it_a.end() || it_b.end())) {
    const auto a_begin = it_a.pos();
    const auto a_end = it_a.pos() + it_a.length();
    const auto b_begin = it_b.pos();
    const auto b_end = it_b.pos() + it_b.length();

    const auto begin_max = (a_begin < b_begin) ? b_begin : a_begin;
    const auto end_min = (a_end < b_end) ? a_end : b_end;
    u1 overlap = begin_max < end_min;

    if (overlap) {
      // Produce an output.
      for (std::size_t i = begin_max; i < end_min; ++i) {
        // Make sure, no bits are set more than once.
        assert(ret_val[i] == false);
        ret_val[i] = true;
      }

      if (a_end <= b_end) {
        it_a.next();
      }
      if (b_end <= a_end) {
        it_b.next();
      }
    }
    else {
      if (a_end < b_end) {
        it_a.skip_to(b_begin);
      }
      else {
        it_b.skip_to(a_begin);
      }
    }
  }
  return ret_val;
}
//===----------------------------------------------------------------------===//
template<typename T>
dtl::bitmap
bitwise_and_iter(const T& bitmap_a, const T& bitmap_b) {
  dtl::bitmap ret_val(bitmap_a.size());
  auto and_it = dtl::bitwise_and_it(bitmap_a.scan_it(), bitmap_b.it());
//  auto it_a = bitmap_a.scan_it();
//  auto it_b = bitmap_b.it();
//  auto and_it = dtl::bitwise_and_it(it_a, it_b);
  while (!and_it.end()) {
    const auto begin = and_it.pos();
    const auto end = and_it.pos() + and_it.length();
    for (std::size_t i = begin; i < end; ++i) {
      // Make sure, no bits are set more than once.
      assert(ret_val[i] == false);
      ret_val[i] = true;
    }
    and_it.next();
  }
  return ret_val;
}
//===----------------------------------------------------------------------===//
template<typename T>
dtl::bitmap
bitwise_or(const T& bitmap_a, const T& bitmap_b) {
  dtl::bitmap ret_val(bitmap_a.size());
  auto it_a = bitmap_a.scan_it();
  auto it_b = bitmap_b.it();
  while (!(it_a.end() || it_b.end())) {
    const auto a_begin = it_a.pos();
    const auto a_end = it_a.pos() + it_a.length();
    const auto b_begin = it_b.pos();
    const auto b_end = it_b.pos() + it_b.length();

    const auto begin_min = (a_begin > b_begin) ? b_begin : a_begin;
    const auto begin_max = (a_begin < b_begin) ? b_begin : a_begin;
    const auto end_min = (a_end < b_end) ? a_end : b_end;
    const auto end_max = (a_end > b_end) ? a_end : b_end;
    u1 contiguous = begin_max <= end_min;

    if (contiguous) {
      // Produce an output.
      for (std::size_t i = begin_min; i < end_max; ++i) {
        // Make sure, no bits are set more than once.
        assert(ret_val[i] == false);
        ret_val[i] = true;
      }
      it_a.skip_to(end_max); // TODO not sure if that is the most efficient way
      it_b.skip_to(end_max);
    }
    else {
      if (a_end < b_begin) {
        // Produce an output.
        for (std::size_t i = a_begin; i < a_end; ++i) {
          // Make sure, no bits are set more than once.
          assert(ret_val[i] == false);
          ret_val[i] = true;
        }
        it_a.next();
      }
      else if (b_end < a_begin) {
        for (std::size_t i = b_begin; i < b_end; ++i) {
          // Make sure, no bits are set more than once.
          assert(ret_val[i] == false);
          ret_val[i] = true;
        }
        it_b.next();
      }
    }
  }
  while (!it_a.end()) {
    for (std::size_t i = it_a.pos(); i < (it_a.pos() + it_a.length()); ++i) {
      // Make sure, no bits are set more than once.
      assert(ret_val[i] == false);
      ret_val[i] = true;
    }
    it_a.next();
  }
  while (!it_b.end()) {
    for (std::size_t i = it_b.pos(); i < (it_b.pos() + it_b.length()); ++i) {
      // Make sure, no bits are set more than once.
      assert(ret_val[i] == false);
      ret_val[i] = true;
    }
    it_b.next();
  }
  return ret_val;
}
//===----------------------------------------------------------------------===//
template<typename T>
dtl::bitmap
bitwise_or_iter(const T& bitmap_a, const T& bitmap_b) {
  dtl::bitmap ret_val(bitmap_a.size());
  auto or_it = dtl::bitwise_or_it(bitmap_a.scan_it(), bitmap_b.it());
//  auto it_a = bitmap_a.scan_it();
//  auto it_b = bitmap_b.it();
//  auto or_it = dtl::bitwise_or_it(it_a, it_b);
  while (!or_it.end()) {
    const auto begin = or_it.pos();
    const auto end = or_it.pos() + or_it.length();
    for (std::size_t i = begin; i < end; ++i) {
      // Make sure, no bits are set more than once.
      assert(ret_val[i] == false);
      ret_val[i] = true;
    }
    or_it.next();
  }
  return ret_val;
}
//===----------------------------------------------------------------------===//
template<typename T>
dtl::bitmap
bitwise_xor(const T& bitmap_a, const T& bitmap_b) {
  dtl::bitmap ret_val(bitmap_a.size());
  auto it_a = bitmap_a.scan_it();
  auto it_b = bitmap_b.it();
  while (!(it_a.end() || it_b.end())) {
    const auto a_begin = it_a.pos();
    const auto a_end = it_a.pos() + it_a.length();
    const auto b_begin = it_b.pos();
    const auto b_end = it_b.pos() + it_b.length();

    const auto begin_min = (a_begin > b_begin) ? b_begin : a_begin;
    const auto begin_max = (a_begin < b_begin) ? b_begin : a_begin;
    const auto end_min = (a_end < b_end) ? a_end : b_end;
    const auto end_max = (a_end > b_end) ? a_end : b_end;
    u1 overlapping = begin_max <= end_min;

    if (overlapping) {
      // Produce an output.
      for (std::size_t i = begin_min; i < begin_max; ++i) {
        assert(i < ret_val.size());
        // Make sure, no bits are set more than once.
        assert(ret_val[i] == false);
        ret_val[i] = true;
      }
      if (a_end < b_end) {
        it_a.next();
        it_b.skip_to(a_end);
      }
      else if (b_end < a_end) {
        it_b.next();
        it_a.skip_to(b_end);
      }
      else {
        it_a.next();
        it_b.next();
      }
    }
    else {
      if (a_end < b_begin) {
        // Produce an output.
        for (std::size_t i = a_begin; i < a_end; ++i) {
          assert(i < ret_val.size());
          // Make sure, no bits are set more than once.
          assert(ret_val[i] == false);
          ret_val[i] = true;
        }
        it_a.next();
      }
      else /*if (b_end < a_begin)*/ {
        for (std::size_t i = b_begin; i < b_end; ++i) {
          assert(i < ret_val.size());
          // Make sure, no bits are set more than once.
          assert(ret_val[i] == false);
          ret_val[i] = true;
        }
        it_b.next();
      }
    }
  }
  while (!it_a.end()) {
    for (std::size_t i = it_a.pos(); i < (it_a.pos() + it_a.length()); ++i) {
      // Make sure, no bits are set more than once.
      assert(i < ret_val.size());
      assert(ret_val[i] == false);
      ret_val[i] = true;
    }
    it_a.next();
  }
  while (!it_b.end()) {
    for (std::size_t i = it_b.pos(); i < (it_b.pos() + it_b.length()); ++i) {
      // Make sure, no bits are set more than once.
      assert(i < ret_val.size());
      assert(ret_val[i] == false);
      ret_val[i] = true;
    }
    it_b.next();
  }
  return std::move(ret_val);
}
//===----------------------------------------------------------------------===//
template<typename T>
dtl::bitmap
bitwise_xor_iter(const T& bitmap_a, const T& bitmap_b) {
  dtl::bitmap ret_val(bitmap_a.size());
  auto xor_it = dtl::bitwise_xor_it(bitmap_a.scan_it(), bitmap_b.it());
//  auto it_a = bitmap_a.scan_it();
//  auto it_b = bitmap_b.it();
//  auto xor_it = dtl::bitwise_xor_it(it_a, it_b);
  while (!xor_it.end()) {
    const auto begin = xor_it.pos();
    const auto end = xor_it.pos() + xor_it.length();
    for (std::size_t i = begin; i < end; ++i) {
      // Make sure, no bits are set more than once.
      assert(ret_val[i] == false);
      ret_val[i] = true;
    }
    xor_it.next();
  }
  return ret_val;
}
//===----------------------------------------------------------------------===//
template<typename T>
dtl::bitmap
bitwise_xor_re(const T& bitmap_a, const T& bitmap_b) {
  dtl::bitmap ret_val(bitmap_a.size());
  auto it_a = bitmap_a.scan_it();
  auto it_b = bitmap_b.it();
  while (!(it_a.end() || it_b.end())) {
    const auto a_begin = it_a.pos();
    const auto a_end = it_a.pos() + it_a.length();
    const auto b_begin = it_b.pos();
    const auto b_end = it_b.pos() + it_b.length();

    if (a_begin < b_begin) {
      const auto end = std::min(a_end, b_begin);
      for (std::size_t i = a_begin; i < b_begin; ++i) {
        // Make sure, no bits are set more than once.
        assert(ret_val[i] == false);
        ret_val[i] = true;
      }
      if (a_end < b_begin) {
        it_a.next();
      }
      else {
        it_a.skip_to(b_end);
        it_b.next();
      }
    }
    else if (b_begin < a_begin) {
      const auto end = std::min(b_end, a_begin);
      for (std::size_t i = b_begin; i < end; ++i) {
        // Make sure, no bits are set more than once.
        assert(ret_val[i] == false);
        ret_val[i] = true;
      }
      if (b_end < a_begin) {
        it_b.next();
      }
      else {
        it_b.skip_to(a_end);
        it_a.next();
      }
    }
    else /* a_begin == b_begin */ {
      const auto end_min = (a_end < b_end) ? a_end : b_end;
      it_a.skip_to(end_min);
      it_b.skip_to(end_min);
    }
  }
  while (!it_a.end()) {
    for (std::size_t i = it_a.pos(); i < (it_a.pos() + it_a.length()); ++i) {
      // Make sure, no bits are set more than once.
      assert(i < ret_val.size());
      assert(ret_val[i] == false);
      ret_val[i] = true;
    }
    it_a.next();
  }
  while (!it_b.end()) {
    for (std::size_t i = it_b.pos(); i < (it_b.pos() + it_b.length()); ++i) {
      // Make sure, no bits are set more than once.
      assert(i < ret_val.size());
      assert(ret_val[i] == false);
      ret_val[i] = true;
    }
    it_b.next();
  }
  return std::move(ret_val);
}
//===----------------------------------------------------------------------===//
template<typename T>
dtl::bitmap
bitwise_xor_re_iter(const T& bitmap_a, const T& bitmap_b) {
  dtl::bitmap ret_val(bitmap_a.size());
  auto xor_it = dtl::bitwise_xor_it(bitmap_a.scan_it(), bitmap_b.it());
//  auto it_a = bitmap_a.scan_it();
//  auto it_b = bitmap_b.it();
//  auto xor_it = dtl::bitwise_xor_re_it(it_a, it_b);
  while (!xor_it.end()) {
    const auto begin = xor_it.pos();
    const auto end = xor_it.pos() + xor_it.length();
    for (std::size_t i = begin; i < end; ++i) {
      // Make sure, no bits are set more than once.
      assert(ret_val[i] == false);
      ret_val[i] = true;
    }
    xor_it.next();
  }
  return ret_val;
}
//===----------------------------------------------------------------------===//
