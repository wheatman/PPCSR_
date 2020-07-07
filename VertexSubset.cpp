#pragma once
#include "BitArray.cpp"
#include <mutex>
#include "sliding_queue.h"


class VertexSubset {
  public:
  bool all;
  bool is_sparse;
  uint64_t max_el;
  BitArray *curr_ba = NULL;
  BitArray *next_ba = NULL;
  SlidingQueue<int32_t> *queue = NULL;
  QueueBuffer<int32_t> *queue_array = NULL;


  bool has(uint64_t i) {
    if (all) {
      return true;
    }
    if (is_sparse) {
      printf("shouldn't be calling has, is currently sparse\n");
      exit(-1);
      return false;
    } else {
      return curr_ba->get(i);
    }
  }
  bool has_dense_no_all(uint64_t i) {
    return curr_ba->get(i);
  }
  void has_dense_no_all_prefetch(uint64_t i) {
    return curr_ba->prefetch(i);
  }


  uint64_t get_n() {
    //printf("get_n: is_sparse = %d, remove_duplicates = %d\n", is_sparse, remove_duplicates);
    if (all) {
      return max_el;
    } else if (is_sparse) {
      return queue->size();
    } else {
      //printf("count = %lu\n", curr_ba->count());
      return curr_ba->count();
    }
  }
  void print() {
    printf("is_sparse = %d\n", is_sparse);
    if (all) {
      printf("{0,...,%lu}\n", max_el);
    } else if (is_sparse) {
      const uint32_t start = queue->shared_out_start;
      const uint32_t end = queue->shared_out_end;
      printf("{");
      for(uint32_t i = start; i < end; i++) {
        printf("%d, ",queue->shared[i]); 
      }
      printf("}\n");
    } else {
      printf("{");
      for(uint32_t i = 0; i < max_el; i++) {
        if (curr_ba->get(i)) {
          printf("%d, ",i); 
        }
      }
      printf("}\n");
    }
  }
  void insert(uint64_t i) {
    if (is_sparse) {
      queue_array[4*getWorkerNum()].push_back(i);
      return;
    } else {
      return next_ba->set(i);
    }
  }
  void insert_sparse(uint64_t i) {
    queue_array[4*getWorkerNum()].push_back(i);
    return;
  }
  void insert_dense(uint64_t i) {
    return next_ba->set(i);
  }
  void move_next_to_current() {
    all = false;
    //printf("moving next to current, %d, %d\n", is_sparse, remove_duplicates);
    if (is_sparse) {
      //printf("should be calling flush soon, workers = %u\n", getWorkers());
      queue->reset();
      parallel_for (int i = 0; i < getWorkers(); i++) {
        queue_array[i*4].flush();
      }
      queue->slide_window();
    } else {
      if (curr_ba) {
        delete curr_ba;
      }
      curr_ba = next_ba;
      next_ba = new BitArray(max_el);
    }
  }
  template <class F> void map(F &f) {
    if (all) {
      parallel_for (uint64_t i = 0; i < max_el; i++) {
        f.update(i);
      }
      return;
    }
    if (is_sparse) {
      const uint32_t start = queue->shared_out_start;
      const uint32_t end = queue->shared_out_end;
      parallel_for(uint32_t i = start; i < end; i++) {
        f.update(queue->shared[i]); 
      }
    } else {
      return curr_ba->map(f);
    }
  }
  VertexSubset() = delete;
  VertexSubset(uint32_t e, uint64_t max_el_, bool all_ = false, bool next = false) : 
        all(all_), is_sparse(true), max_el(max_el_) {
    if (all) {
      is_sparse = false;
      curr_ba = NULL;
      next_ba = NULL;
      queue = NULL;
      queue_array = NULL;
      if (next) {
        next_ba = new BitArray(max_el);
        queue = new SlidingQueue<int32_t>(max_el);
        queue_array = (QueueBuffer<int32_t> *) malloc(4*sizeof(QueueBuffer<int32_t>) * getWorkers());
        for (int i = 0; i < getWorkers(); i++) {
          new(&queue_array[i*4]) QueueBuffer<int32_t>(*queue, max_el);
        }
      }
      return;
    }
    queue = new SlidingQueue<int32_t>(max_el);
    queue_array = (QueueBuffer<int32_t> *) malloc(4*sizeof(QueueBuffer<int32_t>) * getWorkers());
    queue->push_back(e);
    queue->slide_window();
    for (int i = 0; i < getWorkers(); i++) {
      new(&queue_array[i*4]) QueueBuffer<int32_t>(*queue, max_el);
    }
  }
  VertexSubset(bool *els, uint64_t len) : 
        all(false), is_sparse(false), max_el(len) {
    queue = new SlidingQueue<int32_t>(max_el);
    queue_array = (QueueBuffer<int32_t> *) malloc(4*sizeof(QueueBuffer<int32_t>) * getWorkers());
    for (int i = 0; i < getWorkers(); i++) {
      new(&queue_array[i*4]) QueueBuffer<int32_t>(*queue, max_el);
    }
    curr_ba = new BitArray(max_el);
    next_ba = new BitArray(max_el);
    parallel_for(uint64_t j = 0; j < max_el; j+=512) {
      uint64_t end = std::min(j+512, max_el);
      for (uint64_t i = j; i < end; i++) {
        if (els[i]) {
          curr_ba->set(i);
        }
      }
    }
  }

  // can't add anything to these once they have been copied, 
  // just for keeping state like pushing past frontiers into a vector
  VertexSubset(const VertexSubset &other) : all(other.all), is_sparse(other.is_sparse), max_el(other.max_el) {
    if (all) {
      return;
    }
    curr_ba = NULL;
    next_ba = NULL;
    queue = NULL;
    queue_array = NULL;
    if (is_sparse) {
      if (other.queue) {
        queue = new SlidingQueue<int32_t>(*other.queue, max_el);
      }
      if (other.queue_array) {
        queue_array = (QueueBuffer<int32_t> *) malloc(4*sizeof(QueueBuffer<int32_t>) * getWorkers());
        for (int i = 0; i < getWorkers(); i++) {
          new(&queue_array[i*4]) QueueBuffer<int32_t>(*queue, other.queue_array[i*4].local_size);
          queue_array[i*4].in = other.queue_array[i*4].in;
          memcpy(queue_array[i*4].local_queue, other.queue_array[i*4].local_queue, queue_array[i*4].in*sizeof(int32_t));
        }
      }

    } else {
      if (other.curr_ba) {
        curr_ba = new BitArray(*other.curr_ba);
      }
      if(other.next_ba) {
        next_ba = new BitArray(*other.next_ba);
      }
      queue = new SlidingQueue<int32_t>(max_el);
      queue_array = (QueueBuffer<int32_t> *) malloc(4*sizeof(QueueBuffer<int32_t>) * getWorkers());
      for (int i = 0; i < getWorkers(); i++) {
        new(&queue_array[i*4]) QueueBuffer<int32_t>(*queue, max_el);
      }
    }
  }
  void del() {
    this->~VertexSubset();
  }
  ~VertexSubset() {
    if (curr_ba != NULL) {
      delete curr_ba;
    }
    if (next_ba != NULL) {
      delete next_ba;
    }
    if (queue != NULL) {
      delete queue;
      for (int i = 0; i < getWorkers(); i++) {
        queue_array[i*4].~QueueBuffer();
      }
      free(queue_array);
    }
  }

  void convert_to_dense() {
    if (all || !is_sparse) {
      return;
    }
    //printf("converting sparse to dense\n");
    is_sparse = false;
    curr_ba = new BitArray(max_el);
    next_ba = new BitArray(max_el);
    //need an atomic setter to work in parallel
    for (uint32_t i = queue->shared_out_start; i < queue->shared_out_end; i++) {
      curr_ba->set(queue->shared[i]);
    }
    queue->slide_window();

  }
  void convert_to_sparse() {
    if (all || is_sparse) {
      return;
    }
    //printf("converting dense to sparse\n");
    is_sparse = true;
    parallel_for(uint32_t i = 0; i < max_el; i++) {
      if (curr_ba->get(i)) {
        queue_array[4*getWorkerNum()].push_back(i);
      }
    }
    parallel_for (int i = 0; i < getWorkers(); i++) {
      queue_array[i*4].flush();
    }
    queue->slide_window();

    delete curr_ba;
    delete next_ba;
    curr_ba = NULL;
    next_ba = NULL;
  }
};
