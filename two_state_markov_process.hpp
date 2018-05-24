#pragma once

#include <random>

#include <dtl/dtl.hpp>

//===----------------------------------------------------------------------===//
/// Implementation of a two-state Markov process to generate (clustered) bit
/// sequences as defined in the TODS paper "Optimizing Bitmap Indices With
/// Efficient Compression" of Wu et al. (http://www.yajun.info/LBNL-49626-tods.pdf)
class two_state_markov_process {

  double p;
  double q;
  bool state = true;

  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<> dis;

 public:

  /// C'tor
  /// f: Clustering factor, the average number of bits in a 1-fill.
  /// d: Bit density, the fraction of bits that are 1.
  __forceinline__ explicit
  two_state_markov_process(const double f, const double d)
      : p(d / ((1 - d) * f)), q(1 / f), rd(), gen(rd()), dis(0.0, 1.0) {
    // In the paper of Wu et al., the authors decided to always start in state '1'.
    // Thus the generated bit sequence always starts with a 1-fill.
    // The given implementation differs by choosing the initial state randomly.
    state = dis(gen) < 0.5;
    assert(f >= 1.0);
    assert(d >= 0.0 && d <= 1.0);
  }

  __forceinline__
  ~two_state_markov_process() = default;

  __forceinline__
  two_state_markov_process(const two_state_markov_process& other) = default;

  __forceinline__
  two_state_markov_process(two_state_markov_process&& other) = default;

  __forceinline__ two_state_markov_process&
  operator=(const two_state_markov_process& other) = default;

  __forceinline__ two_state_markov_process&
  operator=(two_state_markov_process&& other) = default;

  /// Generate and return the next bit.
  __forceinline__
  bool
  next() {
    const double r = dis(gen);
    const double t = state ? q : p;
    state = r < t == !state;
    return state;
  }

};
//===----------------------------------------------------------------------===//
