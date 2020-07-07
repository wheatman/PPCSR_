// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef PLATFORM_ATOMICS_H_
#define PLATFORM_ATOMICS_H_


/*
GAP Benchmark Suite
File:   Platform Atomics
Author: Scott Beamer

Wrappers for compiler intrinsics for atomic memory operations (AMOs)
 - If not using OpenMP (serial), provides serial fallbacks
*/



// gcc/clang/icc instrinsics

template<typename T, typename U>
T fetch_and_add(T &x, U inc) {
  return __sync_fetch_and_add(&x, inc);
}

template<typename T>
bool compare_and_swap(T &x, const T &old_val, const T &new_val) {
  return __sync_bool_compare_and_swap(&x, old_val, new_val);
}

template<>
bool compare_and_swap(float &x, const float &old_val, const float &new_val) {
  return __sync_bool_compare_and_swap(reinterpret_cast<uint32_t*>(&x),
                                      reinterpret_cast<const uint32_t&>(old_val),
                                      reinterpret_cast<const uint32_t&>(new_val));
}

template<>
bool compare_and_swap(double &x, const double &old_val, const double &new_val) {
  return __sync_bool_compare_and_swap(reinterpret_cast<uint64_t*>(&x),
                                      reinterpret_cast<const uint64_t&>(old_val),
                                      reinterpret_cast<const uint64_t&>(new_val));
}
#endif  // PLATFORM_ATOMICS_H_
