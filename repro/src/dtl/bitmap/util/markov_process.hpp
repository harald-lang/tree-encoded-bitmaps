#pragma once
//===----------------------------------------------------------------------===//
#include <random>
#include <dtl/dtl.hpp>
//===----------------------------------------------------------------------===//
/// Implementation of a k-state Markov process to generate (clustered)
/// attribute values as defined in the TODS paper "Optimizing Bitmap Indices
/// With Efficient Compression" of Wu et al.
/// (http://www.yajun.info/LBNL-49626-tods.pdf).
class markov_process {

  $f64 q;
  $u32 k;
  $u32 state = 0;

  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<> dis;
  std::uniform_int_distribution<> state_dis;

 public:

  /// C'tor
  explicit
  markov_process(u32 k, f64 f)
      : q(1 / f), k(k), rd(), gen(rd()), dis(0.0, 1.0), state_dis(1, k - 1) {
    state = static_cast<u32>(dis(gen) * k);
    assert(f >= 1.0);
    assert(k >= 2);
  }

  ~markov_process() = default;
  markov_process(const markov_process& other) = delete;
  markov_process(markov_process&& other) = delete;
  markov_process& operator=(const markov_process& other) = delete;
  markov_process& operator=(markov_process&& other) = delete;

  /// Generate and return the next value.
  $u32
  next() {
    const double r = dis(gen);
    if (r > q) return state;
    const auto j = state_dis(gen);
    state = (state + j) % k;
    return state;
  }

};
//===----------------------------------------------------------------------===//
