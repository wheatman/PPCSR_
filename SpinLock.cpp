#pragma once
#include <unistd.h>
#include <stdio.h>

// basic spin lock 
class SpinLock {
public:
  SpinLock();

  void lock();
  void unlock();
  uint32_t x;
};

//IMPORTANT!!!!!!!!! locks need to be initialized as all zeros
SpinLock::SpinLock() {
  x = 0;
}

// exclusive lock
void SpinLock::lock() {
  bool success = false;

  while(!success) {
    uint32_t old_val = __sync_fetch_and_add(&x, 1);
    success = (old_val == 0);
    if(!success) {
      __sync_fetch_and_add(&x, -1);
    } else {
      break;
    }
  }
  assert(x > 0);
}

void SpinLock::unlock() {
  assert(x > 0);
  __sync_fetch_and_add(&x, -1);
}
