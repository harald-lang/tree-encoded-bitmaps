#pragma once
#include <bitset>
#include <cstring>
#include <iostream>
#include <limits>
#include "dtl/adept.hpp"

namespace cvt_float {

  $u32 offset_norm_neg = 1;
  $u32 offset_denorm_neg = offset_norm_neg + (static_cast<$u32>(254) * (1<<23)); // = 1 + 254 * 2^23
  $u32 offset_zero_neg = offset_denorm_neg + (1<<23) - 1; // = 255 * 2^23
  $u32 offset_zero_pos = offset_zero_neg + 1; // = 255 * 2^23 + 1
  $u32 offset_denorm_pos = offset_zero_pos + 1; // = 255 * 2^23 + 2
  $u32 offset_norm_pos = offset_denorm_pos + (1<<23) - 1; // = 256 * 2^23 +1
  $u32 offset_inf_pos = offset_norm_pos + static_cast<$u32>(254) * (1<<23); // = 510 * 2^23 +1
  $u32 offset_NaN = offset_inf_pos + 1; // = 510 * 2^23 + 2 @Harald die zwei NaN Faelle unterscheiden?

  void print_bits(float n){
    std::bitset<32> bit;
    std::memcpy(&bit, &n, sizeof(n));

    std::cout << bit[31] << "|";

    auto i = 30;
    for(; i > 22; i--)
      std::cout << bit[i];

    std::cout << "|";

    for(;i>0; i--)
      std::cout << bit[i];

    std::cout << bit[0] << std::endl;
  }

  inline $u32
  rank_norm($u8 exp, $u32 frac, $u32 offset) { //-inf = rank 0, max: 1 + (256-1 -1) + 2^23
    return offset + (exp - 1) * (1 << 23) + frac;
  }

  inline $u32
  rank_denorm($u32 frac, $u32 offset) {
    return offset + frac - 1;
  }

  inline float
  float_norm($u32 rank, $u32 offset) {
    float r = 0;
    $u32 *f = reinterpret_cast<$u32 *>(&r);

    $u32 temp = rank - offset;
    $u32 exp = (temp & 0x7F800000) +(1<<23);
    $u32 frac = temp & 0x007FFFFF;
    *f = (exp | frac);
    if(offset == offset_norm_neg)
      *f |= 0x80000000;

    return r;
  }

  inline float
  float_denorm($u32 rank, $u32 offset) {
    float r = 0;
    $u32 *f = reinterpret_cast<$u32 *>(&r);
    *f = rank - offset + 1; // frac can't be zero for denorm_numbers, its mapped the beginning of the offset

    if(offset == offset_denorm_neg)
      *f |= 0x80000000;

    return r;
  }

  $u32 float_to_rank(float n) {

    $u32 frac;
    std::memcpy(&frac, &n, sizeof(n));
    bool sign = (frac & 0x80000000) > 0;
    $u8 exp = (frac & 0x7F800000) >> 23;
    frac = frac & 0x007FFFFF;

    // calc rank
    if (sign) {
      if (exp != 0 && exp != 255)
        return rank_norm(exp, frac, offset_norm_neg); // Normalized negative
      else {
        if (exp == 0) {
          if (frac == 0)
            return offset_zero_neg; // 0-
          else
            return rank_denorm(frac, offset_denorm_neg); // Denormalized negative
        } else {
          if (frac == 0)
            return 0; // -inf
          else
            return offset_NaN; // NaN
        }
      }
    } else {
      if (exp != 0 && exp != 255)
        return rank_norm(exp, frac, offset_norm_pos); // Normalized positive
      else {
        if (exp == 0) {
          if (frac == 0)
            return offset_zero_pos; // 0+
          else
            return rank_denorm(frac, offset_denorm_pos); // Denormalized positive
        } else {
          if (frac == 0)
            return offset_inf_pos; // inf+
          else
            return offset_NaN; // NaN
        }
      }
    }
  }

  float rank_to_float($u32 r){
    float ret;
    $u32 *tmp = reinterpret_cast<$u32 *>(&ret);

    if(r < offset_zero_pos) { // negative
      if(r < offset_denorm_neg){
        if(r != 0) // normal_negative
          return float_norm(r, offset_norm_neg);
        else // -infinity
          *tmp = 0xFF800000;
      } else {
        if(r != offset_zero_neg) // denormal_negative
          return float_denorm(r, offset_denorm_neg);
        else // 0-
          *tmp = 0x80000000;
      }
    } else { // positive
      if(r < offset_norm_pos){ // denormal or zero
        if(r != offset_zero_pos) // denormal_positive
          return float_denorm(r, offset_denorm_pos);
        else // 0+
          *tmp = 0x00000000;
      } else { // normal, +infinity, NAN
        if(r < offset_inf_pos) // normal_positive
          return float_norm(r, offset_norm_pos);
        else {
          if(r == offset_inf_pos) // +infinity
            *tmp = 0x7F800000;
          else // NaN
            *tmp = 0x7F800001;
        }
      }
    }
    return ret;
  }

};

namespace cvt_double{

  $u64 offset_norm_neg = 1;
  $u64 offset_denorm_neg = offset_norm_neg + (2046 * (static_cast<$u64>(1)<<52)); // = 2046 * (2 ^ 52) + 1
  $u64 offset_zero_neg = offset_denorm_neg + (static_cast<$u64>(1)<<52) - 1; // = 2047 * (2 ^ 52)
  $u64 offset_zero_pos = offset_zero_neg + 1; // = 2047 * (2 ^ 52) +1
  $u64 offset_denorm_pos = offset_zero_pos + 1; // = 2047 * (2 ^ 52) + 2
  $u64 offset_norm_pos = offset_denorm_pos + (static_cast<$u64>(1)<<52) - 1; // = 2048 * (2 ^ 52) + 1
  $u64 offset_inf_pos = offset_norm_pos + (2046 * (static_cast<$u64>(1)<<52)); // = 4094 * (2 ^ 52) + 1
  $u64 offset_NaN = offset_inf_pos + 1; // = 4094 * (2 ^ 52) + 2

  void print_bits(double n){
    std::bitset<64> bit;
    std::memcpy(&bit, &n, sizeof(n));

    std::cout << bit[63] << "|";

    auto i = 62;
    for(; i > 51; i--)
      std::cout << bit[i];

    std::cout << "|";

    for(; i > 0; i--)
      std::cout << bit[i];

    std::cout << bit[0] << std::endl;
  }

  inline $u64
  rank_norm($u16 exp, $u64 frac, $u64 offset) {
    return offset + (exp - 1) * (static_cast<$u64>(1) << 52) + frac;
  }

  inline $u64
  rank_denorm($u64 frac, $u64 offset) {
    return offset + frac - 1;
  }

  inline double
  double_norm($u64 rank, $u64 offset) {
    double r = 0;
    $u64 *f = reinterpret_cast<$u64 *>(&r);

    $u64 temp = rank - offset;
    $u64 exp = (temp & 0x7FF0000000000000) + (static_cast<$u64>(1)<<52);
    $u64 frac = temp & 0x000FFFFFFFFFFFFF;
    *f = (exp | frac);
    if(offset == offset_norm_neg)
      *f |= 0x8000000000000000;

    return r;
  }

  inline double
  double_denorm($u64 rank, $u64 offset) {
    double r = 0;
    $u64 *f = reinterpret_cast<$u64 *>(&r);
    *f = rank - offset + 1; // frac can't be zero for denorm_numbers, its mapped the beginning of the offset

    if(offset == offset_denorm_neg)
      *f |= 0x8000000000000000;

    return r;
  }

  $u64 double_to_rank(double n) {

    $u64 frac;
    std::memcpy(&frac, &n, sizeof(n));
    bool sign = (frac & 0x8000000000000000) > 0;
    $u16 exp = (frac & 0x7FF0000000000000) >> 52;
    frac = frac & 0x000FFFFFFFFFFFFF;

    // calc rank
    if (sign) {
      if (exp != 0 && exp != 2047)
        return rank_norm(exp, frac, offset_norm_neg); // Normalized negative
      else {
        if (exp == 0) {
          if (frac == 0)
            return offset_zero_neg; // 0-
          else
            return rank_denorm(frac, offset_denorm_neg); // Denormalized negative
        } else {
          if (frac == 0)
            return 0; // -inf
          else
            return offset_NaN; // NaN
        }
      }
    } else {
      if (exp != 0 && exp != 2047)
        return rank_norm(exp, frac, offset_norm_pos); // Normalized positive
      else {
        if (exp == 0) {
          if (frac == 0)
            return offset_zero_pos; // 0+
          else
            return rank_denorm(frac, offset_denorm_pos); // Denormalized positive
        } else {
          if (frac == 0)
            return offset_inf_pos; // inf+
          else
            return offset_NaN; // NaN
        }
      }
    }
  }

  double rank_to_double($u64 r){

    double ret = 0;
    $u64 *f = reinterpret_cast<$u64 *>(&ret);

    if(r < offset_zero_pos) { // negative
      if(r < offset_denorm_neg){
        if(r != 0) // normal_negative
          return double_norm(r, offset_norm_neg);
        else // -infinity
          *f = 0xFFF0000000000000;
      } else {
        if(r != offset_zero_neg) // denormal_negative
          return double_denorm(r, offset_denorm_neg);
        else // 0-
          *f = 0x8000000000000000;
      }
    } else { // positive
      if(r < offset_norm_pos){ // denormal or zero
        if(r != offset_zero_pos) // denormal_positive
          return double_denorm(r, offset_denorm_pos);
        else // 0+
          *f = 0x0000000000000000;
      } else { // normal, +infinity, NAN
        if(r < offset_inf_pos) // normal_positive
          return double_norm(r, offset_norm_pos);
        else {
          if(r == offset_inf_pos) // +infinity
            *f = 0x7FF0000000000000;
          else // NaN
            *f = 0x7FF0000000000001;
        }
      }
    }
    return ret;
  }
};