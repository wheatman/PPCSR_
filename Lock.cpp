#include <atomic>
#include <unistd.h>
#include <stdio.h>
#include "helpers.h"

#define num_tries 3
#define wait_time 1


//TODO simplify since they don't need to be recursive anymore
//TODO get rid of locks being deleted

typedef enum reasons {GENERAL, REBALANCE, DOUBLE, SAME} REASONS;
class Lock {
public:
  Lock();

  REASONS reason;
  uint32_t reason_set_by;
  bool i_own_lock(uint32_t task_id);
  bool lock(uint32_t task_id, REASONS r = GENERAL);
  bool try_lock(uint32_t task_id, REASONS r = GENERAL);
  void unlock(uint32_t task_id, REASONS r = GENERAL);
  bool lock_shared(uint32_t task_id = 0);
  void unlock_shared(uint32_t task_id = 0);
  void name(char const about[10]);
  void number(int i);
  void make_shared(uint32_t task_id);
  void make_exclusive(uint32_t task_id);
  bool check_unlocked();
//private:
   uint32_t x;
#ifdef debug
  char info[10];
   int num = 0;
#endif
#ifndef NDEBUG
  uint32_t owner;
#endif
};

//IMPORTANT!!!!!!!!! locks need to be initialized as all zeros
Lock::Lock() {
  x = 0;
  reason = GENERAL;
  reason_set_by = 0;
#ifndef NDEBUG
  owner = 0;
#endif
}

void Lock::name(char const about[10]){
#ifdef debug
  for (int i = 0; i < 10; i++) {
    info[i] = about[i];
  }
#endif
}


void Lock::number(int i){
#ifdef debug
  num = i;
#endif
}

// should only be used in sequntial areas
bool Lock::check_unlocked() {
  return x == 0;
}


bool cas( uint32_t *old, uint32_t old_val, uint32_t new_val) {
  return __sync_bool_compare_and_swap(old, old_val, new_val);
} 


// exclusive lock
bool Lock::lock(uint32_t task_id, REASONS r) {
lock_start:
  assert(task_id > 0);
  assert(r != SAME);
  assert(num_tries > 0);
  int tries;
  bool success = false;


  while(!success) {
    tries = 0;
    while (tries < num_tries) {

      success = cas(&x, 0, 1);
      
      
      if(!success) {
        tries++;
      } else {
        break;
      }
    }
    // sleep in ms
    if(!success) {
      //usleep(wait_time);
      //sched_yield();
    }
  }

  if(r < reason) {
    cas(&x, 1, 0);
    dprintf("worker %lu is trying to lock %p, %s,%d, with the wrong reason %d, needed reason %d\n",task_id,this,info,num, r, reason);
    //usleep(wait_time);
    //sched_yield();
    dprintf("###############################################################################\n");
    goto lock_start;
    //return lock(task_id, r);
  }
#ifndef NDEBUG
  assert(x == 1);
  assert(owner == 0);
  owner = task_id;
#endif
  return true;
}

bool Lock::try_lock(uint32_t task_id, REASONS r) {
  assert(task_id > 0);
  assert(r != SAME);
  assert(num_tries > 0);
  int tries;
  bool success = false;

  tries = 0;
  while (tries < num_tries) {

    if (x == 0) {
      success = cas(&x, 0, 1);
      if(success) { break; }
    }

    tries++;
  }

  if(success && r < reason) {    
    cas(&x, 1, 0);
    dprintf("worker %lu is trying to lock %p, %s,%d, with the wrong reason %d, needed reason %d\n",task_id,this,info,num, r, reason);
    return false;
  }
#ifndef NDEBUG
  if (success) {
    assert(owner == 0);
    owner = task_id;
  }
#endif
  return success;
}


void Lock::unlock(uint32_t task_id, REASONS r) {
  assert(task_id > 0);
#ifdef debug
  assert(num != -1);;
#endif
  assert(num_tries > 0);
  assert(owner == task_id);
  assert(x == 1);

  if (r != SAME && (reason_set_by == task_id || reason == GENERAL || reason_set_by == 0)) {
    reason = r;
    if (r != GENERAL) {
      reason_set_by = task_id;
    } else {
      reason_set_by = 0;
    }
  }
  


#ifdef debug
  assert(num != -1);;
#endif

#ifndef NDEBUG
  owner = 0;
  bool success = cas(&x, 1, 0);
  assert(success);
  
#else
  x = 0;
#endif

}

bool Lock::i_own_lock(uint32_t task_id) {
#ifdef debug
  assert(num != -1);;
#endif
#ifndef NDEBUG
  return owner == task_id;
#else
  printf("i_own_lock should never be called unless we are in debug mode");
  exit(0);
  return true;
#endif
}

void Lock::make_shared(uint32_t task_id) {
  assert(task_id > 0);
  assert(task_id == owner);
#ifdef debug
  assert(num != -1);;
#endif
#ifndef NDEBUG
  owner = 0;
#endif
  bool done = cas(&x, 1, 2);
  assert(done);
}

void Lock::make_exclusive(uint32_t task_id) {
  assert(task_id > 0);
#ifdef debug
  assert(num != -1);;
#endif
  assert(num_tries > 0);
  assert(owner == 0);
  bool success = false;
  int tries = 0;
  while(!success) {
    uint64_t new_x = x;
    dprintf("worker %lu is trying to make exclusive %p,%s,%d with current value %lx\n",task_id,this,info,num, new_x);
    tries = 0;
    if (new_x == 1) {
      while (tries < num_tries) {
        success = cas(&x, 1, 2);
        if (!success) {
          tries++;
          new_x = x;
        } else {
          break;
        }
      }
    }
    //usleep(wait_time);
    //sched_yield();
  }
#ifndef NDEBUG
  owner = task_id;
#endif
}

bool Lock::lock_shared(uint32_t task_id) {
  int tries;
  bool success = false;
  assert(owner == 0);
  while(!success) {
    tries = 0;
    uint64_t old_x = x;
    dprintf("SHARED worker %lu is trying to grab lock %p,%s,%d with current value %lx\n",task_id,this,info,num, old_x);
    if (x % 2 == 0) {
      while (tries < num_tries) {
        success = cas(&x, old_x, old_x + 2);
        if(!success) {
          tries++;
          old_x = x;
        } else {
          break;
        }
      }
    }
    // sleep in ms
    //usleep(wait_time);
    //sched_yield();
  }
  #ifdef debug
    uint64_t lock_val = x;
  #endif
  dprintf("SHARED worker %lu got lock %p with current value %lx\n",task_id,this, lock_val);
  return true;
}

void Lock::unlock_shared(uint32_t task_id) {
#ifdef debug
  assert(num != -1);;
#endif
  assert(num_tries > 0);
  #ifdef debug
    uint64_t old_x = x;
  #endif
  dprintf("SHARED worker %lu is trying to unlock %p,%s,%d with current value %lx\n",task_id,this,info,num, old_x);
  assert(x > 1);
  assert(owner == 0);
  assert(x%2==0);
  __sync_fetch_and_add(&x,-2);
}
