#include <unistd.h>
#include <stdio.h>

#define num_tries 3

// small basic recursive lock
class MiniLock {
public:
  MiniLock();

  void lock();
  void unlock();
  uint32_t x;
};

#define MAX_WORKERS 46
static int primes[MAX_WORKERS] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199};

uint64_t get_worker_num2() {
  return __cilkrts_get_worker_number() + 1;
}

uint32_t get_id() {
  uint32_t worker_num = __cilkrts_get_worker_number();
  assert(worker_num < MAX_WORKERS);
  return primes[worker_num];
}

//IMPORTANT!!!!!!!!! locks need to be initialized as all zeros
MiniLock::MiniLock() {
  x = 0;
}

// exclusive lock
void MiniLock::lock() {
  printf("trying to grab minilock %p, by %lu\n", this, get_worker_num2());
  assert(num_tries > 0);
  int tries;
  bool success = false;

  while(!success) {
    tries = 0;
    while (tries < num_tries) {
      uint32_t old_val = __sync_fetch_and_add(&x, get_id());
      success = ((old_val % get_id()) == 0);
      if(!success) {
        __sync_fetch_and_add(&x, -get_id());
        tries++;
      } else {
        break;
      }
    }
    // sleep in ms
    if(!success) {
      usleep(wait_time);
    }
  }
  assert(x > 0);
}

void MiniLock::unlock() {
  printf("trying to release minilock %p, by %lu\n", this, get_worker_num2());
  assert(x > 0);
  assert(x % get_id() == 0);
  __sync_fetch_and_add(&x, -get_id());
}
