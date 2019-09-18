#pragma once
//===----------------------------------------------------------------------===//
#include <dtl/dtl.hpp>

#include <cmath>
#include <random>
#include <tuple>
//===----------------------------------------------------------------------===//
/// Implementation of a two-state Markov process to generate (clustered) bit
/// sequences as defined in the TODS paper "Optimizing Bitmap Indices With
/// Efficient Compression" of Wu et al.
/// (http://www.yajun.info/LBNL-49626-tods.pdf).
class two_state_markov_process {

  double p;
  double q;
  double d;
  double f;
  bool state = true;

  std::random_device rd;
  std::mt19937 gen;
  std::uniform_real_distribution<> dis;

 public:

  /// C'tor
  /// f: Clustering factor, the average number of bits in a 1-fill.
  /// d: Bit density, the fraction of bits that are 1.
  explicit
  two_state_markov_process(const double f, const double d)
      : p(d / ((1 - d) * f)), q(1 / f), d(d), f(f), rd(), gen(rd()), dis(0.0, 1.0) {
    // In the paper of Wu et al., the authors decided to always start in state '1'.
    // Thus the generated bit sequence always starts with a 1-fill.
    // The given implementation differs by choosing the initial state randomly.
//    state = dis(gen) < 0.5;
    state = dis(gen) < (1-d) ? false : true;
    assert(f >= 1.0);
    assert(d >= 0.0 && d <= 1.0);
  }

  ~two_state_markov_process() = default;

  two_state_markov_process(const two_state_markov_process& other) = delete;
  two_state_markov_process(two_state_markov_process&& other) = delete;
  two_state_markov_process&
  operator=(const two_state_markov_process& other) = delete;

  __forceinline__ two_state_markov_process&
  operator=(two_state_markov_process&& other) = delete;

  /// Generate and return the next bit.
  __forceinline__
  bool
  next() {
    const double r = dis(gen);
    const double t = state ? q : p;
    state = r < t == !state;
    return state;
  }

  std::tuple<double, double>
  state_distribution(const std::tuple<double, double>& init,
      const std::size_t steps = 1) {
    double i0 = std::get<0>(init);
    double i1 = std::get<1>(init);
    double r0 = (1 - p) * i0 + q * i1;
    double r1 = p * i0 + (1 - q) * i1;
//    double r0 = (1 - p) * i0 + p * i1;
//    double r1 = q * i0 + (1 - q) * i1;
    for (std::size_t s = 1; s < steps; ++s) {
      i0 = r0;
      i1 = r1;
      r0 = (1 - p) * i0 + q * i1;
      r1 = p * i0 + (1 - q) * i1;
//      r0 = (1 - p) * i0 + p * i1;
//      r1 = q * i0 + (1 - q) * i1;
    }
    return std::make_tuple(r0, r1);
  };

  double
  entropy_rate() {
    // The asymptotic distribution of the Markov chain.
    const auto mu0 = 1-d ; //q / (p + q);
    const auto mu1 = d;    //p / (p + q);
    // The transition matrix P.
    const auto P00 = 1 - p;
    const auto P01 = p;
    const auto P10 = q;
    const auto P11 = 1 - q;
//    const auto P00 = 1 - p;
//    const auto P01 = q;
//    const auto P10 = p;
//    const auto P11 = 1 - q;

    double r = 0.0;
    // i = 0, j = 0
//    r += mu0 * P00 * std::log(P00); // TODO verify
    r += mu0 * P00 * (std::log(P00) / std::log(2));
    // i = 0, j = 1
//    r += mu0 * P01 * std::log(P01);
    r += mu0 * P01 * (std::log(P01) / std::log(2));
    // i = 1, j = 0
//    r += mu1 * P10 * std::log(P10);
    r += mu1 * P10 * (std::log(P10) / std::log(2));
    // i = 1, j = 1
//    r += mu1 * P11 * std::log(P11);
    r += mu1 * P11 * (std::log(P11) / std::log(2));

    r *= -1.0;
    return r;
  }

};
//===----------------------------------------------------------------------===//
