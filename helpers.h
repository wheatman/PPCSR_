#ifndef HELPERS_H
#define HELPERS_H
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include "parallel.h"
#ifdef debug
#define dprintf(fmt, args...) fprintf(stderr, fmt, ##args)
#else
#define dprintf(fmt, args...) /* Don't do anything in release builds */
#endif

// find index of first 1-bit (least significant bit)
static inline int bsf_word(int word) {
  int result;
  __asm__ volatile("bsf %1, %0" : "=r"(result) : "r"(word));
  return result;
}

static inline int bsr_word(int word) {
  int result;
  __asm__ volatile("bsr %1, %0" : "=r"(result) : "r"(word));
  return result;
}

typedef struct _pair_int {
  uint32_t x; // length in array
  uint32_t y; // depth
} pair_int;

typedef struct _edge {
  uint32_t value; // length in array
  uint32_t dest;
} edge_t;

typedef struct _pair_double {
  double x;
  double y;
} pair_double;

int isPowerOfTwo(int x) { return ((x != 0) && !(x & (x - 1))); }

// same as find_leaf, but does it for any level in the tree
// index: index in array
// len: length of sub-level.
int find_node(int index, int len) { return (index / len) * len; }

// I'm assuming that this is between 0 and 2^31
uint64_t get_worker_num() {
  return __cilkrts_get_worker_number() + 1;
  // return 0;
}

std::vector<uint32_t> parent_to_depth2(int32_t *parents, uint32_t size) {
  std::vector<uint32_t> depths(size, UINT32_MAX);
  cilk_for (uint32_t j = 0; j < size; j++) {
    uint32_t current_depth = 0;
    int32_t current_parent = j;
    if (parents[j] >= 0) {
      while (current_parent != parents[current_parent]) {
        current_depth += 1;
        current_parent = parents[current_parent];
      }
      depths[j] = current_depth;
    }
  }
  return depths;
}

template <typename ET>
inline bool atomic_compare_and_swap(ET *ptr, ET oldv, ET newv) {
//    return __sync_bool_compare_and_swap(ptr, oldv, newv);
  if (sizeof(ET) == 1) {
    return __sync_bool_compare_and_swap((bool*)ptr, *((bool*)&oldv), *((bool*)&newv));
  } else if (sizeof(ET) == 4) {
    return __sync_bool_compare_and_swap((int*)ptr, *((int*)&oldv), *((int*)&newv));
  } else if (sizeof(ET) == 8) {
    return __sync_bool_compare_and_swap((long*)ptr, *((long*)&oldv), *((long*)&newv));
  } else if (sizeof(ET) == 16) {
    return __sync_bool_compare_and_swap_16((__int128*)ptr,*((__int128*)&oldv),*((__int128*)&newv));
  }
  else {
    std::cout << "CAS bad length : " << sizeof(ET) << std::endl;
    abort();
  }
}


#endif
