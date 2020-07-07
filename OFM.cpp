#pragma once
#include "Graph.hpp"
#include "helpers.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <immintrin.h>
#include <atomic>
#include "tbassert.h"
#include <cilk/cilk.h>
#include <cilk/reducer_opadd.h>
#include "Lock.cpp"


// potential TODO: explicitly store the implicit tree (index, len) as a struct
// potential problem: there's a real graph as well as the implicit one.

//TODO order node lock to avoid deadlock

typedef struct _node {
  // beginning and end of the associated region in the edge list
  uint32_t beginning;     // deleted = max int
  uint32_t end;           // end pointer is exclusive
  uint32_t num_neighbors; // number of edgess with this node as source
  Lock lock;
} node_t;

// each node has an associated sentinel (max_int, offset) that gets back to its
// offset into the node array
// UINT32_MAX
//
// if value == 0, read it as null.
//typedef struct _edge { // moved to helpers for general use, left here for helpful comment
//  uint32_t value;
//  uint32_t dest; // destination of this edge in the graph, MAX_INT if this is a
//                 // sentinel
//} edge_t;

typedef union _edgeu {
  edge_t e;
  uint64_t i;
} edge_u;

typedef struct csr_ {
  uint32_t *nodes;
  uint32_t *dests;
  uint32_t *values;
} CSR_simple;


#define NULL_VAL (UINT32_MAX)
#define SENT_VAL (UINT32_MAX -1)
//#define REDISTRIBUTE_PAR_SIZE (UINT32_MAX)
#define REDISTRIBUTE_PAR_SIZE (1 << 15)


typedef struct edge_list {
  volatile uint64_t N;
  volatile uint32_t H;
  volatile uint32_t logN;
  volatile uint32_t loglogN;
  volatile uint32_t mask_for_leaf;
  uint32_t volatile  * volatile vals;
  uint32_t volatile  * volatile dests;

  // Lock list_lock;

  volatile double density_limit;
} edge_list_t;

class OFM : public Graph {
public:
  // data members
  edge_list_t edges;
  std::vector<node_t> nodes;
  Lock node_lock;
  uint64_t next_task_id;

  double upper_density_bound[32];
  double lower_density_bound[32];
  // graph_t g;

  // function headings
  OFM(uint32_t init_n);
  OFM(OFM &other);
  ~OFM();
  void double_list(uint64_t task_id, std::vector<uint64_t> &sub_counts, uint64_t num_elements);
  void half_list(uint64_t task_id, std::vector<uint64_t> &sub_counts, uint64_t num_elements);
  //void half_list();
  void slide_right(uint64_t index, uint32_t *vals, uint32_t *dests);
  void slide_left(uint64_t index, uint32_t *vals, uint32_t *dests);
  void redistribute(uint64_t index, uint64_t len);
  void redistribute_par(uint64_t index, uint64_t len, std::vector<uint64_t> &sub_counts, uint64_t num_elements, bool for_double = false);
  void fix_sentinel(uint32_t node_index, uint64_t in);
  void print_array(uint64_t worker_num = 0);
  uint32_t find_value(uint32_t src, uint32_t dest);
  void print_graph();
  void add_node();
  void add_node_fast(uint64_t task_id);
  void add_edge(uint32_t src, uint32_t dest, uint32_t value);
  void remove_edge(uint32_t src, uint32_t dest);
  void remove_edge_fast(uint32_t src, uint32_t dest, uint64_t task_id);
  void remove_edge_batch(uint32_t *srcs, uint32_t *dest, uint32_t edge_count);
  void add_edge_update(uint32_t src, uint32_t dest, uint32_t value);

  void add_edge_update_fast(uint32_t src, uint32_t dest, uint32_t value, uint64_t task_id);
  
  void add_edge_batch_update(uint32_t *srcs, uint32_t *dests, uint32_t *values, uint32_t edge_count);
  void add_edge_batch_update_no_val(uint32_t *srcs, uint32_t *dests, uint32_t edge_count);


  void insert(uint64_t task_id, uint64_t index, uint32_t elem_dest, uint32_t elem_value, uint32_t src, pair_int held_locks);
  void remove(uint64_t task_id, uint64_t index, uint32_t elem_dest, uint32_t src, pair_int held_locks);
  uint64_t get_size();
  uint64_t get_n();
  void convert(Graph *g);
  void add_file3(string filename);
  vector<tuple<uint32_t, uint32_t, uint32_t>> get_edges();
  void clear();
  bool check_no_locks();
  bool check_no_locks_for_me(uint64_t task_id);
  bool check_no_node_locks_for_me(uint64_t task_id);
  bool check_no_node_locks_held_by_me(uint64_t task_id);
  bool grab_all_locks(uint64_t task_id, bool exclusive, REASONS reason = GENERAL);
  void release_all_locks(uint64_t task_id, bool exclusive, REASONS reason = GENERAL);
  pair_int grab_locks_in_range_with_resets(uint64_t task_id, uint64_t index, uint64_t len, REASONS reason, uint64_t guess);
  pair_int grab_locks_for_leaf_with_resets(uint64_t task_id, uint32_t src, REASONS reason = GENERAL);
  void release_locks_in_range(uint64_t task_id, pair_int locks, REASONS reason = GENERAL);
  uint32_t find_contaning_node(uint64_t index);

  pair_int which_locks_in_range(uint64_t index, uint64_t len, uint64_t guess);
  pair_int which_locks_for_leaf(uint32_t src);
  bool check_every_lock_in_leaf(uint64_t task_id, uint64_t index);
  bool check_every_lock_in_node(uint64_t task_id, uint64_t index, uint64_t len);



  
  uint32_t num_neighbors(uint32_t node) {
    return nodes[node].num_neighbors;
  }
  uint64_t num_edges() {
    uint64_t num = 0;
    for (uint32_t i = 0; i < get_n(); i++) {
      num += num_neighbors(i);
    }
    return num;
  }


  class iterator {
  public:
    uint64_t place;
    uint64_t end;
    uint32_t * vals;
    uint32_t * dests;
    uint8_t loglogN;
    iterator(OFM *G, uint32_t node, bool start) {
      if (!start) {
        place = G->nodes[node].end;
        return;
      }
      place = G->nodes[node].beginning + 1;
      end = G->nodes[node].end;
      vals = (uint32_t *)G->edges.vals;
      dests = (uint32_t *)G->edges.dests;
      loglogN = G->edges.loglogN;
      while ((place < end) && (dests[place] == NULL_VAL)) {
        place = ((place >> loglogN) + 1) << (loglogN);
      }
      if (place > end) {
        place = end;
      }
      return;
    }
    bool operator==(const iterator& other) const {
      return (place == other.place);
    }
    bool operator!=(const iterator& other) const {
      return (place != other.place);
    }
    iterator& operator++() {
      place += 1;
      while ((place < end) && (dests[place] == NULL_VAL)) {
        place = ((place >> loglogN) + 1) << (loglogN);
      }
      if (place > end) {
        place = end;
      }
      return *this;
    }
    edge_t operator*() const {
      return {vals[place], dests[place]};
    }
  };
  iterator begin(uint32_t node) {
    return iterator(this, node, true);
  }
  iterator end(uint32_t node) {
    return iterator(this, node, false);
  }

  void make_symetric() {
    vector<tuple<uint32_t, uint32_t, uint32_t> > edges_to_add = get_edges();
    for(uint64_t i = 0; i < edges_to_add.size(); i++) {
      uint32_t src = get<0>(edges_to_add[i]);
      uint32_t dest = get<1>(edges_to_add[i]);
      uint32_t val = get<2>(edges_to_add[i]);
      add_edge_update(src,dest,val); 
      add_edge_update(dest,src,val);
    }
    return; 
  }
  BFS
  PAGERANK
  SPMV
  TRIANGLE_COUNT_SORTED
  pvector<int32_t> parallel_bfs(int32_t source, int32_t total_edges, int alpha, int beta) {
    int n_workers = __cilkrts_get_nworkers();
    pvector<int32_t> parent(get_n());
    cilk_for (int32_t n = 0; n < get_n(); n++) {
      parent[n] = num_neighbors(source) != 0 ? - num_neighbors(source) : -1;
    }
    parent[source] = source;
    SlidingQueue<int32_t> queue(get_n());
    queue.push_back(source);
    queue.slide_window();
    Bitmap curr(get_n());
    curr.reset();
    Bitmap front(get_n());
    front.reset();
    int64_t edges_to_check = total_edges;
    int64_t scout_count = num_neighbors(source);
    uint32_t *dests = (uint32_t *) edges.dests;
    while (!queue.empty()) {
      if (scout_count > edges_to_check / alpha) {
        int64_t awake_count, old_awake_count;
        cilk_for (int32_t i = queue.shared_out_start; i < queue.shared_out_end; i++) {
          int32_t u = queue.shared[i];
          front.set_bit_atomic(u);
        }
        awake_count = queue.size();
        queue.slide_window();
        do {
          old_awake_count = awake_count;
          awake_count = 0;
          curr.reset();
          vector<int64_t> awake_count_vector(n_workers, 0);
          parallel_for (int32_t u=0; u < get_n(); u++) {
            uint32_t worker_num = __cilkrts_get_worker_number();
            if (parent[u] < 0) {
              int64_t awake_count_local = awake_count_vector[worker_num];
              uint64_t start = nodes[u].beginning + 1;
              uint64_t end = nodes[u].end;
              for (uint32_t i = start; i < end; i++) {
                int32_t v = dests[i];
                if (v < NULL_VAL) {
                  if (front.get_bit(v)) {
                    parent[u] = v;
                    awake_count_local+=1;
                    curr.set_bit_atomic(u);
                    break;
                  }
                }
              }
              /*
              for (iterator it = begin(u); it != end(u); ++it) {
                edge_t edge = *it;
                int32_t v = edge.dest;
                if (front.get_bit(v)) {
                  parent[u] = v;
                  awake_count_local+=1;
                  curr.set_bit(u);
                  break;
                }
              }
              */
              awake_count_vector[worker_num] = awake_count_local;
            }
          }
          for (auto &item : awake_count_vector) {
            awake_count+=item;
          }
          front.swap(curr);
        } while ((awake_count >= old_awake_count) || (awake_count > get_n() / beta));
        QueueBuffer<int32_t> *queue_array = (QueueBuffer<int32_t> *)malloc(4*sizeof(QueueBuffer<int32_t>) * n_workers);
        if (queue_array == NULL) {
		printf("bad malloc in bfs\n");
		exit(-1);
	}
	for (int i = 0; i < n_workers; i++) {
          new(&queue_array[i*4]) QueueBuffer(queue);
        }
        cilk_for (int32_t n=0; n < get_n(); n++) {
          if (front.get_bit(n)) {
            queue_array[__cilkrts_get_worker_number()*4].push_back(n);
          }
        }
        cilk_for (int i = 0; i < n_workers; i++) {
          queue_array[i*4].flush();
	  delete &queue_array[i*4];
        }
	free(queue_array);
        queue.slide_window();
        scout_count = 1;
      } else {
        edges_to_check -= scout_count;
        scout_count = 0;
        vector<int64_t> scout_count_vector(n_workers, 0);
        QueueBuffer<int32_t> *queue_array = (QueueBuffer<int32_t> *)malloc(4*sizeof(QueueBuffer<int32_t>) * n_workers);
        if (queue_array == NULL) {
		printf("bad malloc in bfs\n");
		exit(-1);
	}
        for (int i = 0; i < n_workers; i++) {
          new(&queue_array[i*4]) QueueBuffer(queue);
        }
        cilk_for (int32_t i = queue.shared_out_start; i < queue.shared_out_end; i++) {
          int32_t u = queue.shared[i];
          uint32_t worker_num = __cilkrts_get_worker_number();
          for (iterator it = begin(u); it != end(u); ++it) {
            edge_t edge = *it;
            int32_t v = edge.dest;
            int32_t curr_val = parent[v];
            if (curr_val < 0) {
              if (compare_and_swap(parent[v], curr_val, u)) {
                queue_array[4*worker_num].push_back(v);
                scout_count_vector[worker_num] += -curr_val;
              }
            }
          }
        }
        cilk_for (int i = 0; i < n_workers; i++) {
          queue_array[i*4].flush();
	  delete &queue_array[i*4];
        }
	free(queue_array);
        for (auto &item : scout_count_vector) {
          scout_count+=item;
        }
        queue.slide_window();
      }
    }
    for (int32_t n = 0; n < get_n(); n++) {
      if (parent[n] < -1) {
        parent[n] = -1;
      }
    }
    return parent;
  }
  template <class F, typename VS, bool output> 
  void map_sparse(F &f, VS &vs, uint64_t self_index) {
    uint32_t idx = nodes[self_index].beginning + 1;
    uint32_t idx_end = nodes[self_index].end;
    while ( idx < idx_end) {
      uint32_t v = edges.dests[idx];
      if ( v != NULL_VAL) {
        if (f.cond(v) == 1 && f.updateAtomic(self_index, v) == 1) {
          if constexpr (output) {
            vs.insert_sparse(v);
          }
        }
        idx++;
      } else {
        idx = ((idx >> edges.loglogN) +1 ) << (edges.loglogN);
      }
    }
  }
  template <class F, typename VS, bool output, bool vs_all> 
  void map_dense(F &f, VS &vs, uint64_t self_index) {
    if constexpr (!vs_all) {
      uint32_t idx = nodes[self_index].beginning + 1;
      uint32_t idx_end = nodes[self_index].end;
      while ( idx < idx_end) {
        uint32_t v = edges.dests[idx];
        if ( v != NULL_VAL) {
          if (vs.has_dense_no_all(v) && f.update(v, self_index) == 1) {
            if constexpr(output) {
              vs.insert_dense(self_index);
            }
          }
          if (f.cond(self_index) == 0) {
            return;
          }
          idx++;
        } else {
          idx = ((idx >> edges.loglogN) +1 ) << (edges.loglogN);
        }
      }
    } else {
      uint32_t idx = nodes[self_index].beginning + 1;
      uint32_t idx_end = nodes[self_index].end;
      while ( idx < idx_end) {
        uint32_t v = edges.dests[idx];
        if ( v != NULL_VAL) {
          if (f.update(v, self_index) == 1) {
            if constexpr(output) {
              vs.insert_dense(self_index);
            }
          }
          if (f.cond(self_index) == 0) {
            return;
          }
          idx++;
        } else {
          idx = ((idx >> edges.loglogN) +1 ) << (edges.loglogN);
        }
      }
    }
  }


  CSR_simple __attribute__ ((noinline)) convert_to_csr() {
    uint32_t N = nodes.size();
    uint32_t M = num_edges();
    CSR_simple csr = {NULL, NULL};
    csr.nodes = (uint32_t *) malloc(N * sizeof(uint32_t));
    uint32_t start = 0;
    // parallel_prefix_sum
    for (uint32_t i = 0; i < N; i++) {
      csr.nodes[i] = start;
      start += num_neighbors(i);
    }
    csr.dests = (uint32_t *) malloc(M * sizeof(uint32_t));
    csr.values = (uint32_t *) malloc(M * sizeof(uint32_t));
    cilk_for(uint32_t i = 0; i < N; i++) {
      uint32_t * dest_array = &csr.dests[csr.nodes[i]];
      uint32_t * value_array = &csr.values[csr.nodes[i]];
      iterator it_end = end(i);
      uint32_t count = 0;
      for (iterator it = begin(i); it != it_end; ++it) {
        edge_t edge = *it;
        dest_array[count] = edge.dest;
        value_array[count] = edge.value;
        count++;
      }
    }
    return csr;
  }

};


// given index, return the starting index of the leaf it is in
//TODO this could be aster if we store a mask and just do a single and
uint64_t find_leaf(edge_list_t *list, uint64_t index) {
  return index & list->mask_for_leaf;
}


uint64_t find_prev_valid(uint32_t volatile  * volatile dests, uint64_t start) {
  while (dests[start] == NULL_VAL) {
    start--;
  }
  return start;
}

bool OFM::check_no_locks() {
  bool ret = true;
  ret = ret && node_lock.check_unlocked();
  assert(node_lock.check_unlocked());
  // ret = ret && edges.list_lock.check_unlocked();
  // assert(edges.list_lock.check_unlocked());
  for (uint32_t i = 0; i < nodes.size(); i++) {
    ret = ret && nodes[i].lock.check_unlocked();
    if (!ret) {
      printf("found lock on iter %d\n", i);
    }
    ret = ret && (nodes[i].lock.reason == GENERAL);
    if (!ret) {
      printf("found lock on iter %d with %d reason\n", i, nodes[i].lock.reason);
    }
    assert(nodes[i].lock.check_unlocked());
    assert(nodes[i].lock.reason == GENERAL);
  }
  return ret;
}

bool OFM::check_no_locks_for_me(uint64_t task_id) {
  bool ret = true;
  ret = ret && !node_lock.i_own_lock(task_id);
  assert(!node_lock.i_own_lock(task_id));
  // ret = ret && edges.list_lock.check_unlocked();
  // assert(edges.list_lock.check_unlocked());
  for (uint32_t i = 0; i < nodes.size(); i++) {
    ret = ret && !nodes[i].lock.i_own_lock(task_id);
    if (!ret) {
      printf("found lock on iter %d\n", i);
    }
    ret = ret && (nodes[i].lock.reason_set_by != task_id);
    if (!ret) {
      printf("found lock on iter %d that worker %lu set the reason with %d\n", i,task_id, nodes[i].lock.reason);
    }
    assert(!nodes[i].lock.i_own_lock(task_id));
    assert(nodes[i].lock.reason_set_by != task_id);
  }
  return ret;
}  

bool OFM::check_no_node_locks_for_me(uint64_t task_id) {
  bool ret = true;
  for (uint32_t i = 0; i < nodes.size(); i++) {
    ret = ret && !nodes[i].lock.i_own_lock(task_id);
    if (!ret) {
      printf("found lock on iter %d\n", i);
    }
    ret = ret && (nodes[i].lock.reason_set_by != task_id);
    if (!ret) {
      printf("found lock on iter %d that worker %lu set the reason with %d\n", i,task_id, nodes[i].lock.reason);
    }
    assert(!nodes[i].lock.i_own_lock(task_id));
    assert(nodes[i].lock.reason_set_by != task_id);
  }
  return ret;
}  

bool OFM::check_no_node_locks_held_by_me(uint64_t task_id) {
  bool ret = true;
  for (uint32_t i = 0; i < nodes.size(); i++) {
    ret = ret && !nodes[i].lock.i_own_lock(task_id);
    if (!ret) {
      printf("found lock on iter %d\n", i);
    }
    assert(!nodes[i].lock.i_own_lock(task_id));
  }
  return ret;
}  


uint64_t next_leaf(uint64_t index, uint32_t loglogN) {
  return ((index >> loglogN) + 1) << (loglogN);
}

bool OFM::grab_all_locks(uint64_t task_id, bool exclusive, REASONS reason) {
  for (uint32_t i = 0; i < nodes.size(); i++) {
    if (exclusive) {
      if (!nodes[i].lock.lock(task_id, reason)) {
        return false;
      }
    } else  {
      if (!nodes[i].lock.lock_shared(task_id)) {
        return false;
      }
    }
  }
  return true;
}
void OFM::release_all_locks(uint64_t task_id, bool exclusive, REASONS reason) {
  cilk_for (uint32_t i = 0; i < nodes.size(); i++) {
    if (exclusive) {
      nodes[i].lock.unlock(task_id, reason);
    } else  {
      nodes[i].lock.unlock_shared(task_id);
    }
  }
}

void OFM::clear() {
  printf("clear called\n");
  uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2);
  grab_all_locks(task_id, true, GENERAL);
  
  free((void*)edges.vals);
  free((void*)edges.dests);
  edges.N = 2 << bsr_word(0);
  // printf("%d\n", bsf_word(list->N));
  edges.loglogN = bsr_word(bsr_word(edges.N) + 1);
  edges.logN = (1 << edges.loglogN);
  edges.H = bsr_word(edges.N / edges.logN);
}
// for soc-
// starting at 1
void OFM::add_file3(string filename) {
  ifstream myfile(filename.c_str());
  string line;
  if (myfile.is_open()) {
    while (getline(myfile, line)) {
      vector<string> elems = split(line, '\t');
      int src = atoi(elems[0].c_str()) - 1;

      while (src >= get_n()) {
        add_node();
      }
      int dest = atoi(elems[1].c_str()) - 1;

      while (dest >= get_n()) {
        add_node();
      }
      add_edge(src, dest, 1);
      // if (line_num++ > 400000000) {
      //  break;
      // }
    }
    myfile.close();
    // return 0;
  } else {
    printf("file was not opened\n");
  }
}

// TODO jump to next leaf
vector<tuple<uint32_t, uint32_t, uint32_t>> OFM::get_edges() {
  // TODO grab locks in the lock list
  // for now, grabs the global lock
  node_lock.lock_shared(); // lock node array
  // edges.list_lock.lock();
  uint64_t n = get_n();
  vector<tuple<uint32_t, uint32_t, uint32_t>> output;

  for (uint64_t i = 0; i < n; i++) {
    uint64_t start = nodes[i].beginning;
    uint64_t end = nodes[i].end;
    nodes[i].lock.lock_shared();
    for (uint64_t j = start + 1; j < end; j++) {
      if (edges.dests[j]!=NULL_VAL) {
        output.push_back(
            make_tuple(i, edges.dests[j], edges.vals[j]));
      }
    }
    nodes[i].lock.unlock_shared();
  }
  // edges.list_lock.unlock();
  node_lock.unlock_shared(); // lock node array
  return output;
}


//TODO this might only work if you start empty
void OFM::convert(Graph *g) {
  uint64_t n = g->get_n();
  free((void*)edges.vals);
  free((void*)edges.dests);
  edges.N = 2UL << bsr_word(n);
  // printf("%d\n", bsf_word(list->N));
  edges.loglogN = bsr_word(bsr_word(edges.N) + 1);
  edges.logN = (1 << edges.loglogN);
  assert(edges.logN > 0);
  edges.H = bsr_word(edges.N / edges.logN);
  // printf("N = %d, logN = %d, loglogN = %d, H = %d\n", list->N, list->logN,
  // list->loglogN, list->H);

  edges.dests = (uint32_t *)malloc(edges.N * sizeof(*(edges.dests)));
	if (edges.dests == NULL) {
		printf("bad malloc in convert\n");
		exit(-1);
	}
  edges.vals = (uint32_t *)malloc(edges.N * sizeof(*(edges.vals)));
	if (edges.vals == NULL) {
		printf("bad malloc in convert\n");
		exit(-1);
	}
  for (uint64_t i = 0; i < edges.N; i++) {
    edges.dests[i] = NULL_VAL;
    edges.vals[i] = 0;
  }

  for (uint64_t i = 0; i < n; i++) {
    add_node();
  }
  //TODO why does making these cilk_for crash the compilier 
  for (uint64_t i = 0; i < n; i++) {
    for (uint64_t j = 0; j < n; j++) {
      // find_value returns 0 if not found.
      uint32_t value = g->find_value(i, j);
      if (value != 0) {
        add_edge(i, j, value);
      }
    }
  }
}

uint64_t OFM::get_n() {
  node_lock.lock_shared();
  uint64_t size = nodes.size();
  node_lock.unlock_shared();
  return size;
}

uint64_t OFM::get_size() {
  node_lock.lock_shared();
  uint64_t size = nodes.capacity() * sizeof(node_t);
  size += sizeof(OFM);
  size += edges.N * sizeof(edge_t);
  node_lock.unlock_shared();
  return size;
}


void print_array(edge_list_t *edges) {
  printf("N = %lu, logN = %d\n", edges->N, edges->logN);
  for (uint64_t i = 0; i < edges->N; i++) {
    if (edges->dests[i]==NULL_VAL) {
      printf("%lu-x ", i);
    } else if ((edges->dests[i]==SENT_VAL) || i == 0) {
      uint32_t value = edges->vals[i];
      if (value == NULL_VAL) {
        value = 0;
      }
      printf("\n%lu-s(%u):(?, ?) ", i, value);
    } else {
      printf("%lu-(%d, %u) ", i, edges->dests[i], edges->vals[i]);
    }
  }
  printf("\n\n");
}

void OFM::print_array(uint64_t worker_num) {
  printf("worker num: %lu, N = %lu, logN = %d, density_limit = %f\n", worker_num, edges.N, edges.logN, edges.density_limit);
  for (uint64_t i = 0; i < edges.N; i++) {
    if (edges.dests[i]==NULL_VAL) {
      printf("%lu-x ", i);
    } else if ((edges.dests[i] == SENT_VAL) || i == 0) {
      uint32_t value = edges.vals[i];
      if (value == NULL_VAL) {
        value = 0;
      }
      printf("\n worker num: %lu, %lu-s(%u):(%d, %d)(%d) (%s, %d ", worker_num, i, value, nodes[value].beginning,
             nodes[value].end, nodes[value].num_neighbors, nodes[value].lock.check_unlocked() ? "Unlocked": "Locked", nodes[value].lock.reason);
      #ifndef NDEBUG
        printf(" by %u)", nodes[value].lock.owner);
      #else
        printf(")");
      #endif
    } else {
      printf("%lu-(%d, %u) ", i, edges.dests[i], edges.vals[i]);
    }
  }
  printf("\n\n");
}

// get density of a node
// should already be locked if you are calling get density
double get_density(edge_list_t *list, uint64_t index, uint64_t len) {
  uint64_t full = 0;
  uint32_t volatile  * volatile dests = list->dests;
  for (uint64_t i = index; i < index+len; i+=4) {
      uint64_t add = (dests[i]!=NULL_VAL) + (dests[i+1]!=NULL_VAL) + (dests[i+2]!=NULL_VAL) + (dests[i+3]!=NULL_VAL);
      full += add;
  }
  double full_d = (double)full;
  return full_d / len;
}

uint64_t get_density_count(edge_list_t *list, uint64_t index, uint64_t len) {
 /* 
  uint32_t full = 0;
  uint32_t i = index;
  while (i < index + len) {
    if (!is_null(list->items[i].e)) {
      full++;
      i++;
    } else {
      i = next_leaf(i, list->logN);
    }
  }
  return full;
*/
  /*
  cilk::reducer< cilk::op_add<uint32_t> > full;
  cilk_for (uint32_t i = index; i < index+len; i++) {
    if (!is_null(list->items[i].e)) {
      (*full)++;
    }
  }
  return full.get_value();
  */
  // fater without paralleliszation since it gets properly vectorized
  uint32_t volatile  * volatile dests = list->dests;
  uint64_t full = 0;
  for (uint64_t i = index; i < index+len; i+=4) {
      uint64_t add = (dests[i]!=NULL_VAL) + (dests[i+1]!=NULL_VAL) + (dests[i+2]!=NULL_VAL) + (dests[i+3]!=NULL_VAL);
      //__sync_fetch_and_add(&full, add);
      full+=add;
  }
  return full;
  /*
  cilk::reducer< cilk::op_add<uint32_t> > full;
  cilk_for (uint32_t i = index; i < index+len; i+=4) {
      uint32_t add = !is_null(list->items[i].e) + !is_null(list->items[i+1].e) + !is_null(list->items[i+2].e) + !is_null(list->items[i+3].e);
      (*full)+=add;
  }
  return full.get_value();
  */
  
}

uint64_t get_density_count_par(edge_list_t *list, uint64_t index, uint64_t len, std::vector<uint64_t> &sub_counts) {
  cilk::reducer< cilk::op_add<uint64_t> > total;
  uint32_t volatile  * volatile dests = list->dests;
  cilk_for(uint64_t j = index; j < index+len; j+= REDISTRIBUTE_PAR_SIZE) {
    uint64_t full = 0;
    for (uint64_t i = j; i < j+REDISTRIBUTE_PAR_SIZE; i+=4) {
        uint64_t add = (dests[i]!=NULL_VAL) + (dests[i+1]!=NULL_VAL) + (dests[i+2]!=NULL_VAL) + (dests[i+3]!=NULL_VAL);
        full+=add;
    }
    (*total)+=full;
    sub_counts[(j-index)/REDISTRIBUTE_PAR_SIZE] = full;
  }
  /*
  for (int i = 0; i < sub_counts.size(); i++) {
    printf("subcount[%d] = %u\n", i, sub_counts[i]);
  }
  print_array(list);
  */
  return total.get_value();  
}

bool check_no_full_leaves(edge_list_t *list, uint64_t index, uint64_t len) {
  for (uint64_t i = index; i < index + len; i+= list->logN) {
    bool full = true;
    for (uint64_t j = i; j < i + list->logN; j++) {
       if (list->dests[j]==NULL_VAL) {
        full = false;
      }
    }
    if (full) {
      return false;
    }
  }
  return true;
}

// height of this node in the tree
uint32_t get_depth(edge_list_t *list, uint64_t len) { return bsr_word(list->N / len); }

// when adjusting the list size, make sure you're still in the
// density bound
pair_double density_bound(edge_list_t *list, uint64_t depth) {
  pair_double pair;

  // between 1/4 and 1/2
  // pair.x = 1.0/2.0 - (( .25*depth)/list->H);
  // between 1/8 and 1/4
  pair.x = 1.0 / 4.0 - ((.125 * depth) / list->H);
  pair.y = 3.0 / 4.0 + ((.25 * depth) / list->H);
  if (pair.y > list->density_limit) {
    pair.y = list->density_limit+.001;
  }
  return pair;
}

//TODO make it so the first element is known to always be the first element and don't special case it so much
//assumes the node_lock is held
void OFM::fix_sentinel(uint32_t node_index, uint64_t in) {
  // we know the first sentinal will never move so we just ignore it
  assert(node_index > 0);
  //node_lock.lock_shared();
  if (in >= 1UL<<32) {
	  printf("fix_sentinel called with something too big: %lu, the graph exceeds the capabilities\n",in);
	  exit(-1);
  }
	if (node_index == 0) {
	    printf("got 0 for node index, should never happen\n");
	    while(1){}
	}
  nodes[node_index - 1].end = in;

  nodes[node_index].beginning = in;
  if (node_index == nodes.size() - 1) {
    nodes[node_index].end = edges.N - 1;
  }
  //node_lock.unlock_shared();
}


// Evenly redistribute elements in the ofm, given a range to look into
// index: starting position in ofm structure
// len: area to redistribute
// should already be locked
void OFM::redistribute(uint64_t index, uint64_t len) {
  //printf("len = %u\n", len);
  assert(find_leaf(&edges, index) == index);
  
  //printf("REDISTRIBUTE START: index:%u, len %u, worker %lu\n", index, len, get_worker_num());
  //print_array(get_worker_num());
  // std::vector<edge_t> space(len); //
  //TODO if its small use the stack
  // for small cases put on the stack
  uint32_t *space_vals;
  uint32_t *space_dests;
  uint32_t volatile  * volatile vals = edges.vals;
  uint32_t volatile  * volatile dests = edges.dests;
  uint64_t j = 0;
  if (len == edges.logN) {
    return;
    j = index;
    //printf("index = %u\n",index);
    //print_array(0);
    for (uint64_t i = index; i < index + len; i++) {
      vals[j] = vals[i];
      dests[j] = dests[i];
      if (dests[j]==SENT_VAL) {
        // fixing pointer of node that goes to this sentinel
        uint64_t node_index = vals[j];
        fix_sentinel(node_index, j);
      }
      j += (dests[j]!=NULL_VAL);
    }
    for (uint64_t i = j; i < index+len; i++) {
      vals[i] = 0;
      dests[i] = NULL_VAL;
    }
    //print_array(0);
    return;
  } else {
    space_vals = (uint32_t *)malloc(len * sizeof(*(edges.vals)));
	if (space_vals == NULL) {
		printf("bad malloc in redistribute 1\n");
		exit(-1);
	}
    space_dests = (uint32_t *)malloc(len * sizeof(*(edges.dests)));
	if (space_dests == NULL) {
		printf("bad malloc in redistribute 2\n");
		exit(-1);
	}

  }

  // move all items in ofm in the range into
  // a temp array
  /*
  int i = index;
  while (i < index + len) {
    if (!is_null(edges.items[i])) {
      space[j] = edges.items[i];
      edges.items[i].value = 0;
      edges.items[i].dest = 0;
      i++;
      j++;
    } else {
      i = next_leaf(i, edges.logN);
    }
  }
  */
  //TODO could parralize if get_density_count gave us more data, but doesn't seem to be a bottle neck
  // could get better cache behavior if I go back and forth with reading and writing
  for (uint64_t i = index; i < index + len; i++) {
    space_vals[j] = vals[i];
    space_dests[j] = dests[i];
    // counting non-null edges
    j += (space_dests[j]!=NULL_VAL);
    // setting section to null
    vals[i] = 0;
    dests[i] = NULL_VAL;
  }
  /*
  if (((double)j)/len > ((double)(edges.logN-1)/edges.logN)) {
    printf("too dense in redistribute, j = %u, len = %u, index = %u for worker %lu\n",j, len, index, get_worker_num() );
    print_array(get_worker_num());
  }*/
  assert( ((double)j)/len <= ((double)(edges.logN-1)/edges.logN));

  uint64_t num_leaves = len >> edges.loglogN;
  uint64_t count_per_leaf = j / num_leaves;
  uint64_t extra = j % num_leaves;

  // parallizing does not make it faster
  for (uint64_t i = 0; i < num_leaves; i++) {
    uint64_t count_for_leaf = count_per_leaf + (i < extra);
    uint64_t in = index + (edges.logN * (i));
    uint64_t j2 = count_per_leaf*i + min(i,extra);
    //TODO could be parallized, but normally only up to size 32
    uint64_t j3 = j2;
    for (uint64_t k = in; k < count_for_leaf+in; k++) {
      vals[k] = space_vals[j2];
      j2++;
    }
    for (uint64_t k = in; k < count_for_leaf+in; k++) {
      dests[k] = space_dests[j3];
      if (dests[k]==SENT_VAL) {
        // fixing pointer of node that goes to this sentinel
        uint32_t node_index = vals[k];
        fix_sentinel(node_index, k);
      }
      j3++;
    }

  }
  free(space_dests);
  free(space_vals);
  /*
  if (!check_no_full_leaves(&edges, index, len)) {
    printf("some leaves are full, current density is %f, index = %u, len = %u\n", get_density(&edges, index, len), index, len);
    print_array(get_worker_num());
  }*/
  assert(check_no_full_leaves(&edges, index, len));
  //printf("REDISTRIBUTE END: index:%u, len %u, worker %lu\n", index, len, get_worker_num());
}

void OFM::redistribute_par(uint64_t index, uint64_t len, std::vector<uint64_t> &sub_counts, uint64_t num_elements, bool for_double) {
  assert(find_leaf(&edges, index) == index);
  //printf("par len = %u\n", len);

  uint32_t *space_vals = (uint32_t *)aligned_alloc(64, len * sizeof(*(edges.vals)));
  uint32_t *space_dests = (uint32_t *)aligned_alloc(64, len * sizeof(*(edges.dests)));
  if (space_vals == NULL) {
    printf("bad malloc for space_vals in redistribute_par\n");
  }
  if (space_dests == NULL) {
    printf("bad malloc for space_dests in redistribute_par\n");
  }
  uint64_t j = 0;

  // move all items in ofm in the range into
  // a temp array


  //TODO parallel prefix sum
  for (uint64_t i = 1; i < sub_counts.size(); i++) {
    sub_counts[i] += sub_counts[i-1];
    //printf("sub_counts[%d] = %u\n", i, sub_counts[i]);
  }

  // could get better cache behavior if I go back and forth with reading and writing
  uint32_t volatile  * volatile vals = edges.vals;
  uint32_t volatile  * volatile dests = edges.dests;
 /* 
  for (uint32_t i = index; i < index + len; i++) {
    j += (dests[i]!=NULL_VAL);
    printf("start = %u, i = %u\n", j, i);
  }
  */
  j = num_elements;
  uint64_t end = index+len;
  if (for_double) {
    end = end/2;
  }
  cilk_for(uint64_t k = index; k < end; k+=REDISTRIBUTE_PAR_SIZE) {
    //printf("index = %u, len = %u, k = %u REDISTRIBUTE_PAR_SIZE = %u\n", index, len, k , REDISTRIBUTE_PAR_SIZE);
      uint64_t start;
      if (k == index) {
        start = 0;
      } else {
        start = sub_counts[(k-index)/REDISTRIBUTE_PAR_SIZE-1];
      }

      for (uint64_t i = k; i < k+REDISTRIBUTE_PAR_SIZE; i+=8) {
        //printf("start = %u, i = %u\n", start, i);
        for (uint64_t j = i; j < i+8; j++) {
          if (dests[j] != NULL_VAL) {
            space_vals[start] = vals[j];
            space_dests[start] = dests[j];
            start += 1;
          }
#ifndef __SSE4_2__
	    vals[j] = 0;
	    dests[j] = NULL_VAL;
#endif
        }
        // setting section to null
#ifdef __SSE4_2__
        _mm256_storeu_si256((__m256i *) (&edges.vals[i]), _mm256_setzero_si256());
        _mm256_storeu_si256((__m256i *) (&edges.dests[i]), _mm256_set_epi32(NULL_VAL, NULL_VAL, NULL_VAL, NULL_VAL, NULL_VAL, NULL_VAL, NULL_VAL, NULL_VAL));  
#endif
      }

  }
/*  
  j = 0;
  for(uint32_t k = index; k < index+len; k+=REDISTRIBUTE_PAR_SIZE) {
    for (uint32_t i = k; i < k + REDISTRIBUTE_PAR_SIZE; i++) {
      space_vals[j] = vals[i];
      space_dests[j] = dests[i];
      // counting non-null edges
      j += (space_dests[j]!=NULL_VAL);
      // setting section to null
      vals[i] = 0;
      dests[i] = NULL_VAL;
    }
  }
*/
  assert( ((double)j)/len <= ((double)(edges.logN-1)/edges.logN));

  uint64_t num_leaves = len >> edges.loglogN;
  uint64_t count_per_leaf = j / num_leaves;
  uint64_t extra = j % num_leaves;

  // parallizing does not make it faster
  cilk_for (uint64_t i = 0; i < num_leaves; i++) {
    uint64_t count_for_leaf = count_per_leaf + (i < extra);
    uint64_t in = index + ((i) << edges.loglogN);
    uint64_t j2 = count_per_leaf*i + min(i,extra);
    //TODO could be parallized, but normally only up to size 32
    uint64_t j3 = j2;
    assert(j3 < len);
    
    for (uint64_t k = in; k < count_for_leaf+in; k++) {
      vals[k] = space_vals[j2];
      j2++;
    }
    for (uint64_t k = in; k < count_for_leaf+in; k++) {
      dests[k] = space_dests[j3];
      if (dests[k]==SENT_VAL) {
        // fixing pointer of node that goes to this sentinel
        uint64_t node_index = vals[k];
         
        fix_sentinel(node_index, k);
      }
      j3++;
      assert(j3 < len);
    }

  }
  free(space_dests);
  free(space_vals);

  assert(check_no_full_leaves(&edges, index, len));
}


//TODO pass in subcounts and do redistibute_par when big
void OFM::double_list(uint64_t task_id, std::vector<uint64_t> &sub_counts, uint64_t num_elements) {
  //printf("doubling list by worker %lu\n", get_worker_num());
  grab_all_locks(task_id, true, DOUBLE);
  uint64_t new_N = 2UL * edges.N;
  edges.loglogN = bsr_word(bsr_word(edges.N) + 1);
  edges.logN = (1 << edges.loglogN);
  edges.mask_for_leaf = ~(edges.logN - 1);
  assert(edges.logN > 0);
  edges.density_limit = ((double) edges.logN - 1)/edges.logN;
  edges.H = bsr_word(new_N / edges.logN);
  for (uint32_t i = 0; i <= edges.H; i++) {
    upper_density_bound[i] = density_bound(&edges, i).y;
    lower_density_bound[i] = density_bound(&edges, i).x;
  }

  edges.N = new_N;


  /* 
  uint32_t *new_dests = (uint32_t *)malloc(new_N * sizeof(*(edges.dests)));
  uint32_t *new_vals = (uint32_t *)malloc(new_N * sizeof(*(edges.vals)));
  cilk_for (uint32_t i = 0; i < edges.N / 2; i++) {
    new_vals[i] = vals[i]; // setting second half to null
    new_dests[i] = dests[i]; // setting second half to null
  }
  edges.dests = new_dests;
  edges.vals = new_vals;
  vals = edges.vals;
  dests = edges.dests;
  */
  uint32_t *new_dests = (uint32_t *)aligned_alloc(64, new_N * sizeof(*(edges.dests)));
	if (new_dests == NULL) {
		printf("bad ralloc in double 1\n");
		exit(-1);
	}
  if (sub_counts.size() == 0) {
    for (uint64_t i = 0; i < new_N / 2; i++) {
      new_dests[i] = edges.dests[i];
    }
  } else {
    cilk_for (uint64_t i = 0; i < new_N / 2; i+= REDISTRIBUTE_PAR_SIZE) {
      for (uint64_t j = i; j < i+REDISTRIBUTE_PAR_SIZE; j++) {
        new_dests[j] = edges.dests[j];
      }
      //memcpy((void*)&new_dests[i], (void*)&edges.dests[i], REDISTRIBUTE_PAR_SIZE * sizeof(uint32_t));
    }
  }

  free((void*)edges.dests);
  edges.dests = new_dests;
  //edges.dests = (uint32_t *)realloc((void*)edges.dests, new_N * sizeof(*(edges.dests)));
  uint32_t *new_vals = (uint32_t *)aligned_alloc(64, new_N * sizeof(*(edges.vals)));
	if (new_vals == NULL) {
		printf("bad ralloc in double 1\n");
		exit(-1);
	}
  if (sub_counts.size() == 0) {
    for (uint64_t i = 0; i < new_N / 2; i++) {
      new_vals[i] = edges.vals[i];
    }
  } else {
    cilk_for (uint64_t i = 0; i < new_N / 2; i+= REDISTRIBUTE_PAR_SIZE) {
      for (uint64_t j = i; j < i+REDISTRIBUTE_PAR_SIZE; j++) {
        new_vals[j] = edges.vals[j];
      }
      //memcpy((void*)&new_vals[i], (void*)&edges.vals[i], REDISTRIBUTE_PAR_SIZE * sizeof(*(edges.vals)));
    }
  }
  free((void*)edges.vals);
  edges.vals = new_vals;
  //edges.vals = (uint32_t *)realloc((void*)edges.vals, new_N * sizeof(*(edges.vals)));
  uint32_t volatile * volatile vals = edges.vals;
  uint32_t volatile * volatile dests = edges.dests;

  assert(edges.dests != NULL && edges.vals != NULL);
  //memset((void*)(edges.items+(new_N / 2)), 0, sizeof(edge_t) * (new_N / 2) );
  
  cilk_for (uint64_t i = new_N / 2; i < edges.N; i++) {
    vals[i] = 0; // setting second half to null
    dests[i] = NULL_VAL; // setting second half to null
  }
  
  //printf("List doubled: N - %u, logN = %u, H = %u\n", edges.N, edges.logN, edges.H);
  if (sub_counts.size() == 0) {
    redistribute(0, new_N);
  } else {
    redistribute_par(0, new_N, sub_counts, num_elements, true);
  }
  release_all_locks(task_id, true, GENERAL);
}

void OFM::half_list(uint64_t task_id, std::vector<uint64_t> &sub_counts, uint64_t num_elements) {
  //printf("doubling list by worker %lu\n", get_worker_num());
  grab_all_locks(task_id, true, DOUBLE);
  uint64_t new_N = edges.N / 2;
  edges.loglogN = bsr_word(bsr_word(edges.N) - 1);
  edges.logN = (1 << edges.loglogN);
  edges.mask_for_leaf = ~(edges.logN - 1);
  assert(edges.logN > 0);
  edges.density_limit = ((double) edges.logN - 1)/edges.logN;
  edges.H = bsr_word(new_N / edges.logN);
  for (uint32_t i = 0; i <= edges.H; i++) {
    upper_density_bound[i] = density_bound(&edges, i).y;
    lower_density_bound[i] = density_bound(&edges, i).x;
  }


  /* 
  uint32_t *new_dests = (uint32_t *)malloc(new_N * sizeof(*(edges.dests)));
  uint32_t *new_vals = (uint32_t *)malloc(new_N * sizeof(*(edges.vals)));
  cilk_for (uint32_t i = 0; i < edges.N / 2; i++) {
    new_vals[i] = vals[i]; // setting second half to null
    new_dests[i] = dests[i]; // setting second half to null
  }
  edges.dests = new_dests;
  edges.vals = new_vals;
  vals = edges.vals;
  dests = edges.dests;
  */
  uint32_t *new_dests = (uint32_t *)malloc(new_N * sizeof(*(edges.dests)));
	if (new_dests == NULL) {
		printf("bad malloc in half 1\n");
		exit(-1);
	}
  uint32_t *new_vals = (uint32_t *)malloc(new_N * sizeof(*(edges.vals)));
	if (new_vals == NULL) {
		printf("bad malloc in half 2\n");
		exit(-1);
	}
  uint32_t * vals = (uint32_t *)edges.vals;
  uint32_t * dests = (uint32_t *)edges.dests;

  assert(edges.dests != NULL && edges.vals != NULL);
  
  
  //printf("List doubled: N - %u, logN = %u, H = %u\n", edges.N, edges.logN, edges.H);
  if (sub_counts.size() == 0) {
    uint64_t start = 0;
    for (uint64_t i = 0; i < edges.N; i++) {
      if (dests[i] != NULL_VAL) {
        new_vals[start] = vals[i];
        new_dests[start] = dests[i];
        start += 1;
      }
    }
    free(vals);
    free(dests);
    edges.vals = new_vals;
    edges.dests = new_dests;
    edges.N = new_N;
    redistribute(0, new_N);
  } else {
    cilk_for(uint32_t k = 0; k < sub_counts.size(); k+=1) {
      //printf("index = %u, len = %u, k = %u REDISTRIBUTE_PAR_SIZE = %u\n", index, len, k , REDISTRIBUTE_PAR_SIZE);
      uint64_t start;
      if (k == 0) {
        start = 0;
      } else {
        start = sub_counts[k-1];
      }

      for (uint64_t i = REDISTRIBUTE_PAR_SIZE*k; i < REDISTRIBUTE_PAR_SIZE*(k+1); i+=8) {
        //printf("start = %u, i = %u\n", start, i);
        for (uint64_t j = i; j < i+8; j++) {
          if (dests[j] != NULL_VAL) {
            new_vals[start] = vals[j];
            new_dests[start] = dests[j];
            start += 1;
          }
        }
      }
    }
    free(vals);
    free(dests);
    edges.vals = new_vals;
    edges.dests = new_dests;
    edges.N = new_N;
    redistribute(0, new_N);
  }
  release_all_locks(task_id, true, GENERAL);
}

// index is the beginning of the sequence that you want to slide right.
// we wil always hold locks to the end of the leaf so we don't need to lock here
void OFM::slide_right(uint64_t index, uint32_t *vals, uint32_t *dests) {
  // edges.list_lock.lock_shared();
  uint32_t el_val = vals[index];
  uint32_t el_dest = dests[index];
  dests[index] = NULL_VAL;
  vals[index] = 0;
  //uint64_t original_index = index;
  //printf("start of slide right, original_index: %d, worker number: %lu, \n", original_index, get_worker_num());
  // edges.lock_array[current_lock].print();
  index++;
  //uint32_t leaf = find_leaf(&edges, index);
  while (index < edges.N && (dests[index]!=NULL_VAL)) {
    assert(find_leaf(&edges, index) == leaf);
    uint32_t temp_val = vals[index];
    uint32_t temp_dest = dests[index];
    vals[index] = el_val;
    dests[index] = el_dest;
    if (el_dest == SENT_VAL) {
      // fixing pointer of node that goes to this sentinel
      uint32_t node_index = el_val;
      fix_sentinel(node_index, index);
    }
    el_val = temp_val;
    el_dest = temp_dest;
    index++;
  }
  if (el_dest == SENT_VAL) {
    // fixing pointer of node that goes to this sentinel
    uint32_t node_index = el_val;
    if (node_index == 0) {
	    printf("got 0 for node index, should never happen\n");
	    while(1){}
    }
    fix_sentinel(node_index, index);
  }
  //printf("middle of slide right, original_index: %d, worker number: %lu, current lock = %d\n", original_index, get_worker_num(), current_lock);

  // TODO There might be an issue with this going of the end sometimes
  assert(index != edges.N);

  vals[index] = el_val;
  dests[index] = el_dest;
  assert(find_leaf(&edges, index) == leaf);
  assert(check_no_full_leaves(&edges, original_index, edges.logN));
}

// index is the beginning of the sequence that you want to slide left.
// the element we start at will be deleted
// we wil always hold locks to the end of the leaf so we don't need to lock here
void OFM::slide_left(uint64_t index, uint32_t *vals, uint32_t *dests) {
  // edges.list_lock.lock_shared();
  //uint64_t original_index = index;
  //printf("start of slide right, original_index: %d, worker number: %lu, \n", original_index, get_worker_num());
  // edges.lock_array[current_lock].print();
  //uint64_t leaf = find_leaf(&edges, index);
  while (index+1 < edges.N) {
    assert(find_leaf(&edges, index) == leaf);
    uint32_t temp_val = vals[index+1];
    uint32_t temp_dest = dests[index+1];
    vals[index] = temp_val;
    dests[index] = temp_dest;
    if (temp_dest == SENT_VAL) {
      // fixing pointer of node that goes to this sentinel
      uint32_t node_index = temp_val;
      fix_sentinel(node_index, index);
    }
    if (dests[index]==NULL_VAL) {
      break;
    }
    index++;
  }

  assert(index != edges.N);

  assert(find_leaf(&edges, index) == leaf);
  assert(check_no_full_leaves(&edges, original_index, edges.logN));
}

// important: make sure start, end don't include sentinels
// returns the index of the smallest element bigger than you in the range
// [start, end)
// if no such element is found, returns end (because insert shifts everything to
// the right)
// assumes we already hold the list_lock and the relevant lock_array locks
uint64_t binary_search(edge_list_t *list, uint32_t elem_dest, uint32_t elem_val, uint64_t start,
                       uint64_t end) {
  uint32_t volatile  * volatile dests = list->dests;

  uint64_t mid = (start + end) / 2;
  // print_array(list);
  while (start + 1 < end) {
    
    //__builtin_prefetch ((void *)&dests[(mid+end)/2], 0, 3);
    //__builtin_prefetch ((void *)&dests[(start + mid)/2], 0, 3);
    // printf("start = %d, end = %d, dest = %d, mid = %d, val =%u\n", start, end, elem_dest, mid, elem_val);
    /*
    if (mid % list->logN > list->logN / 2) {
      // if beginning of next leaf is before end of binary search
      uint32_t temp = next_leaf(mid, list->loglogN);
      if(temp < end) {
        mid = temp;
      }
    }
    */
    uint32_t item_dest = dests[mid];

    //if is_null
    if (item_dest==NULL_VAL) {
      // first check the next leaf
      uint64_t check = next_leaf(mid, list->loglogN);
      //TODO deal with check is null
      if (check > end) {
        end = mid;
        mid = (start + end) / 2;
        continue;
      }
      // if is_null
      if (dests[check]==NULL_VAL) {
        uint64_t early_check = find_prev_valid(dests, mid);
        // if we found the sentinel, go right after it
        if (dests[early_check] == SENT_VAL) {
          return early_check + 1;
        }
        if (early_check < start) {
          start = mid;
          mid = (start + end) / 2;
          //__builtin_prefetch ((void *)&dests[mid], 0, 3);
          continue;
        }  
        check = early_check;
      } 
      // printf("check = %d\n", check);
      uint32_t dest = dests[check];
      if (elem_dest == dest) {
        // cleanup before return
        return check;
      } else if (elem_dest < dest) {
        end = find_prev_valid(dests, mid) + 1;

      } else {
        if (check == start) {
          start = check + 1;
        } else {
          start = check;
        }
        // otherwise, searched for item is more than current and we set start
      }
      mid = (start + end) / 2;
      //__builtin_prefetch ((void *)&dests[mid], 0, 3);
      continue;
    }

    if (elem_dest < item_dest) {
      end = mid; // if the searched for item is less than current item, set end
      mid = (start + end) / 2;
    } else if (elem_dest > item_dest) {
      start = mid;
      mid = (start + end) / 2;
      // otherwise, sesarched for item is more than current and we set start
    } else if (elem_dest == item_dest) {  // if we found it, return
      // cleanup before return
      return mid;
    }
  }
  if (end < start) {
    start = end;
  }
  assert(start >= 0);
  tbassert(end < list->N, "end: %u, list->N: %u\n", end, list->N);

  //trying to encourage the packed left property so if they are both null go to the left
  if ((dests[start]==NULL_VAL) && (dests[end]==NULL_VAL)) {
    end = start;
  }

  // handling the case where there is one element left
  // if you are leq, return start (index where elt is)
  // otherwise, return end (no element greater than you in the range)
  // printf("start = %d, end = %d, n = %d\n", start,end, list->N);
  if (elem_dest <= dests[start] && (dests[start]!=NULL_VAL)) {
    end = start;
  }
  // cleanup before return

  return end;
}

uint32_t OFM::find_value(uint32_t src, uint32_t dest) {
  uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2);

  node_lock.lock_shared(task_id);
  uint32_t e_value = 0;
  uint32_t e_dest = dest;
  //printf("src: %d, dest: %d\n", src, dest);
  //printf("beginning: %d, end:%d\n", nodes[src].beginning+1, nodes[src].end);

  //the lock has been deleted
  // not sure why this has to be exclusive
  if (!nodes[src].lock.lock(task_id, GENERAL)) {
    node_lock.unlock_shared(task_id);
    return find_value(src,dest);
  }
  uint64_t loc =
      binary_search(&edges, e_dest, e_value, nodes[src].beginning + 1, nodes[src].end);
  //printf("loc = %d, looking for %u, %u, found %u, %u\n",loc, src, dest, edges.dests[loc], edges.vals[loc]);
  e_dest = edges.dests[loc];
  e_value = edges.vals[loc];
  nodes[src].lock.unlock(task_id, GENERAL);
  // edges.list_lock.unlock_shared();
  node_lock.unlock_shared(task_id);
  // print_array();
  // printf("loc: %d\n", loc);
  //TODO probably don't need the first check since we will never look for null
  if ((e_dest != NULL_VAL) && e_dest == dest) {
    return e_value;
  } else {
    return 0;
  }
}

//assumes node_lock is held
uint32_t OFM::find_contaning_node(uint64_t index) {
  uint64_t start = 0; 
  uint64_t end = nodes.size()-1;
  while (end - start > 1) {
    uint64_t middle = (end + start) / 2;
    uint64_t node_start = nodes[middle].beginning;
    uint64_t node_end = nodes[middle].end;
    if ( index >= node_start && index < node_end){
      return middle;
    } else if (index < node_start) {
      end = middle;
    } else if (index >= node_end) {
      start = middle;
    } else {
      printf("should not happen\n");
      assert(false);
    }
  }
  if ( index >= nodes[start].beginning && index < nodes[start].end){
    return start;
  } else if ( index >= nodes[end].beginning && index < nodes[end].end) {
    return end;
  } else if (index >= nodes[nodes.size() - 1].end) {
      return nodes.size() - 1;
  } else {
    //printf("no containing node trying again\n");
    return find_contaning_node(index);
  }
}

pair_int OFM::which_locks_in_range(uint64_t index, uint64_t len, uint64_t guess) {
  uint32_t start;
which_locks_in_range_start:
  
  if (index >= nodes[guess].beginning && index < nodes[guess].end) {
      start = guess;
  } else {
    bool found = false;
    uint64_t range = 50;
    uint64_t st = (guess > range) ? guess - range : 0;
    uint64_t end = (guess + range < nodes.size()) ? guess + range : nodes.size() - 1;
    if (index >= nodes[st].beginning && index < nodes[end].end) {
      for (uint64_t est = st; est < end; est++) {
        if (index >= nodes[est].beginning && index < nodes[est].end) {
          start = est;
          found = true;
          break;
        }
      }
    }
    if (!found) {
      start = find_contaning_node(index);
    }
  }

  if (edges.dests[index]==SENT_VAL) {
    // we are unlocked so this can cuase a rare segfault so just start over
    if (start == 0) {
      goto which_locks_in_range_start;
    }
    start--;
  }
  uint64_t end_index = next_leaf(index + len, edges.loglogN) - 1;
  uint32_t end;
  if (end_index < nodes[start].end) {
    end = start;
  } else {
    end = find_contaning_node(end_index);
  }
  //printf("grabbing locks %d through %d, by %lu with reason %d\n", start, end, get_worker_num(), reason);
  if (end < start) {
    goto which_locks_in_range_start;
  }
  return {start, end};
}


pair_int OFM::grab_locks_in_range_with_resets(uint64_t task_id, uint64_t index, uint64_t len, REASONS reason, uint64_t guess) {
  uint64_t orig_n = edges.N;
  start_grab_locks_in_range_with_resets:
  if (orig_n != edges.N) {
    return {0xFFFFFFFF,0};
  }
  pair_int ends = which_locks_in_range(index, len, guess);
  //printf("grabbing locks %d through %d, by %lu with reason %d\n", ends.x, ends.y, get_worker_num(), reason);
  for (uint32_t i = ends.x; i <= ends.y; i++) {
    if (!nodes[i].lock.try_lock(task_id, reason)) {
      if ( i > 0) {
        release_locks_in_range(task_id, {ends.x, i - 1}, SAME);
      }
      //usleep(1);
      //return grab_locks_in_range_with_resets(task_id, index, len, reason);
      goto start_grab_locks_in_range_with_resets;
    }
  }
  //printf("grabed locks %d through %d, by %lu with reason %d\n", ends.x, ends.y, get_worker_num(), reason);
  // make sure nothing changed before we managed to grab the locks
  if (!((nodes[ends.x].beginning <= index) && ((index+len <= nodes[ends.y].end) || (ends.y == nodes.size() -1)))) {
    release_locks_in_range(task_id, ends, SAME);
    goto start_grab_locks_in_range_with_resets;
  }
  return ends;
}

void OFM::release_locks_in_range(uint64_t task_id, pair_int locks, REASONS reason) {
  uint64_t start = locks.x;
  uint64_t end = locks.y;
  //printf("releasing locks %d through %d, by worker %lu with reason %d\n", start, end, get_worker_num(), reason);
  for (uint32_t i = start; i <= end; i++) {
    nodes[i].lock.unlock(task_id, reason);
  }
}

pair_int OFM::which_locks_for_leaf(uint32_t src) {
  uint64_t start_index = find_leaf(&edges, nodes[src].beginning);
  uint64_t end_index = next_leaf(nodes[src].end, edges.loglogN);
  uint32_t first_node = src;
  while (nodes[first_node].beginning > start_index) {
    first_node--;
  }
  uint32_t last_node = src;
  while (nodes[last_node].end < end_index && last_node < nodes.size() -1) {
    last_node++;
  }
  return {first_node, last_node};
}

pair_int OFM::grab_locks_for_leaf_with_resets(uint64_t task_id, uint32_t src, REASONS reason) {
start_grab_locks_for_leaf_with_resets:
  pair_int ends = which_locks_for_leaf(src);
  //printf("grabbing locks %d through %d, by %lu with reason %d\n", ends.x, ends.y, get_worker_num(), reason);
  for (uint32_t i = ends.x; i <= ends.y; i++) {
    if (!nodes[i].lock.try_lock(task_id, reason)) {
      if ( i > 0) {
        release_locks_in_range(task_id, {ends.x, i - 1}, SAME);
      }
      usleep(1);
     //sched_yield();
      goto start_grab_locks_for_leaf_with_resets;
    }
  }
  //printf("grabed locks %d through %d, by %lu with reason %d\n", ends.x, ends.y, get_worker_num(), reason);
  pair_int ends_check = which_locks_for_leaf(src);
  if (ends.x != ends_check.x || ends.y != ends_check.y) {
    release_locks_in_range(task_id, ends);
    goto start_grab_locks_for_leaf_with_resets;
  }
  return ends;
}

bool OFM::check_every_lock_in_leaf(uint64_t task_id, uint64_t index) {
  uint64_t start = find_leaf(&edges, index);
  uint64_t end = next_leaf(index, edges.loglogN);
  for (uint64_t i = start; i < end; i++) {
    uint32_t node = find_contaning_node(i);
    assert(nodes[node].lock.i_own_lock(task_id));
    if (!nodes[node].lock.i_own_lock(task_id)) {
      return false;
    }
  }
  return true;
}

bool OFM::check_every_lock_in_node(uint64_t task_id, uint64_t index, uint64_t len) {
  uint64_t start = find_leaf(&edges, index);
  for (uint64_t i = start; i < start + len; i++) {
    uint32_t node = find_contaning_node(i);
    if (!nodes[node].lock.i_own_lock(task_id)) {
      return false;
    }
    assert(nodes[node].lock.i_own_lock(task_id));
  }
  return true;
}


// insert elem at index
// assumes the lock on index is held by the parent
// and releases it when it is done with it
void OFM::insert(uint64_t task_id, uint64_t index, uint32_t elem_dest, uint32_t elem_val, uint32_t src, pair_int held_locks) {
  assert(check_every_lock_in_leaf(task_id, index));
  assert(check_no_full_leaves(&edges, find_leaf(&edges, index), edges.logN));
  
  /*
  if (index > 0 && edges.dests[index - 1] == NULL_VAL && find_leaf(&edges, index) != index) {
    printf("insert conditional || index = %d src = %u elem dest = %u, elem value = %u, next leaf: %d\n", index, src, elem_dest, elem_val, next_leaf(index, edges.loglogN));
    print_array(get_worker_num());
  }*/


  //printf("index = %lu elem dest = %u, elem value = %u, worker %lu\n", index, elem_dest, elem_val, get_worker_num());
  //print_array(get_worker_num());
  uint64_t orig_n = edges.N;
  
  if ((elem_dest != SENT_VAL) && (elem_val != 0)) {
    // -1 for the case that the location is the spot of the next sentinal
    assert(src == find_contaning_node(index-1));
  }
  // edges.list_lock.lock_shared();
  assert(nodes[src].lock.i_own_lock(task_id));
  uint32_t level = edges.H;
  uint64_t len = edges.logN;

  uint32_t * vals = (uint32_t *) edges.vals;
  uint32_t * dests = (uint32_t *) edges.dests;
  // always deposit on the left
  if (dests[index]==NULL_VAL) {
    // printf("added to empty\n");
    vals[index] = elem_val;
    dests[index] = elem_dest;
    assert(check_every_lock_in_leaf(task_id, index));

  } else {
    assert(index < edges.N - 1);
    // if the edge already exists in the graph, update its value
    // do not make another edge
    if ((elem_dest != SENT_VAL) && dests[index] == elem_dest) {
      vals[index] = elem_val;
      assert(check_every_lock_in_leaf(task_id, index));
      // edges.list_lock.unlock_shared();
      assert(check_no_full_leaves(&edges, find_leaf(&edges, index), edges.logN));
      release_locks_in_range(task_id, held_locks);
      assert(check_no_node_locks_for_me(task_id));
      return;
    } else {
      // slide right assumes we hold the lock till the end of the leaf
      assert(check_every_lock_in_leaf(task_id, index));
      slide_right(index, vals, dests);
      // printf("after sliding, index = %d\n", index);
      // print_array();
      vals[index] = elem_val;
      dests[index] = elem_dest;
      assert(check_every_lock_in_leaf(task_id, index));
      // print_array();
    }
  }
  
  //assert(vals[index]!=0);
  //print_array();
  uint32_t node_index = find_leaf(&edges, index);
  double density = get_density(&edges, node_index, len);
  //printf("density = %f, %d\n", density, density == 1);

  // spill over into next level up, node is completely full.
  if (density == 1) {
    //printf("first rebalence\n");
    len*=2;
    level--;
    node_index = find_node(node_index, len);
  } else {
    redistribute(node_index, len);
  }
  assert(edges.N == orig_n);
  assert(check_every_lock_in_leaf(task_id, index));
  // being fancy here since some of the locks are done with and others might be used for the rebalence
  pair_int locks_for_rebalance = which_locks_in_range(node_index, len, src);
  //printf("relasing locks %d through %d with reason %d\n", locks_for_leaf.x, locks_for_rebalance.x -1, GENERAL);
  //printf("before with worker %lu\n", get_worker_num());
  //printf("worker %lu unlocking %d through %d with GENERAL \n", get_worker_num(), held_locks.x, locks_for_rebalance.x);
  cilk_for (uint32_t i = held_locks.x; i < locks_for_rebalance.x; i++) {
    nodes[i].lock.unlock(task_id);
  }
  //printf("worker %lu relasing locks %d through %d with reason REBALANCE\n",get_worker_num(), max(locks_for_rebalance.x, locks_for_leaf.x), locks_for_leaf.y);
  cilk_for (uint32_t i = max(locks_for_rebalance.x, held_locks.x) ; i <= min(held_locks.y, locks_for_rebalance.y); i++) {
    nodes[i].lock.unlock(task_id, REBALANCE);
  }
  //printf("worker %lu unlocking %d through %d with GENERAL \n", get_worker_num(), locks_for_rebalance.y +1,locks_for_leaf.y);

  cilk_for (uint32_t i = locks_for_rebalance.y +1; i <= held_locks.y; i++) {
    nodes[i].lock.unlock(task_id);
  }

  //printf("unlocking %d by worker %lu\n",src, get_worker_num());
  //printf("grabbing %d, %d\n", node_index, len);
  //printf("about to grab\n");
  pair_int lock_span = {max(locks_for_rebalance.x, held_locks.x),  max(locks_for_rebalance.y, held_locks.y)};
  assert(check_no_node_locks_held_by_me(task_id));

  held_locks = grab_locks_in_range_with_resets(task_id, node_index, len, REBALANCE, src);
  lock_span.x = min(held_locks.x, lock_span.x);
  lock_span.y = max(held_locks.y, lock_span.y);
  // if somebody doubled and we let them with a reset
  if (edges.N != orig_n) {
    //printf("worker %lu: somebody else doubled on us so we assume it has been distributed properly 1\n", get_worker_num());
    //printf("worker %lu: edges.N = %u, orig_n = %u\n",get_worker_num(), edges.N, orig_n);
    release_locks_in_range(task_id, held_locks, GENERAL);
    for (uint32_t i = lock_span.x; i <= lock_span.y; i++) {
      if (i < held_locks.x || i > held_locks.y) {
        nodes[i].lock.lock(task_id, REBALANCE);
        nodes[i].lock.unlock(task_id);
      }
    }
    assert(check_no_node_locks_for_me(task_id));
    return;
  }
  // printf("node_index3 = %d\n", node_index);
  // print_array();

  // get density of the leaf you are in
  double density_b = upper_density_bound[level];
  uint64_t density_count = get_density_count(&edges, node_index, len);
  density = ((double)density_count)/len;
  // printf("density %f, upperbound %f, len = %d, N = %d, logN = %d\n", density,
  // density_b.y, len, list->N, list->logN);

  // while density too high, go up the implicit tree
  // go up to the biggest node above the density bound
  //printf("node_index = %d, desnsity = %f, density bound = %f\n", node_index, density, density_b.y);
  std::vector<uint64_t> sub_counts(0);
  while (density >= density_b) {
    //printf("node_index = %d, desnsity = %f, density bound = %f, len = %d, worker = %lu\n", node_index, density, density_b.y, len, get_worker_num());
    len *= 2;
    if (len <= edges.N) {
      release_locks_in_range(task_id, held_locks, REBALANCE);
      level--;
      uint32_t new_node_index = find_node(node_index, len);
      held_locks = grab_locks_in_range_with_resets(task_id, new_node_index, len, REBALANCE, src);
      if (held_locks.x < lock_span.x) {
        lock_span.x = held_locks.x;
      }
      if (held_locks.y > lock_span.y) {
        lock_span.y = held_locks.y;
      }
      if (edges.N != orig_n) {
        //printf("worker %lu: somebody else doubled on us so we assume it has been distributed properly 2\n", get_worker_num());
        //printf("worker %lu: edges.N = %u, orig_n = %u\n",get_worker_num(), edges.N, orig_n);
        release_locks_in_range(task_id, held_locks, GENERAL);
        for (uint32_t i = lock_span.x; i <= lock_span.y; i++) {
          if (i < held_locks.x || i > held_locks.y) {
            nodes[i].lock.lock(task_id, REBALANCE);
            nodes[i].lock.unlock(task_id);
          }
        }
        assert(check_no_node_locks_for_me(task_id));
        return;
      }
      if (len <= REDISTRIBUTE_PAR_SIZE) {
        density_count = get_density_count(&edges, new_node_index, len);
      } else {
        sub_counts.resize(len/REDISTRIBUTE_PAR_SIZE);
        density_count = get_density_count_par(&edges, new_node_index, len, sub_counts);
      }
      // to help prevent double doubling by knowing how big it was on the last count
      orig_n = edges.N;
      node_index = new_node_index;
      density_b = upper_density_bound[level];
      density = ((double) density_count )/len;
    } else {
      // if you reach the root, double the list
      release_locks_in_range(task_id, held_locks, DOUBLE);
      //if (edges.N == orig_n) {
        // for (int i = 0; i < lock_count; i++) {
        //  edges.list_lock.unlock_shared();
        //}
        //TODO don't double double
        //printf("second double due to worker %lu\n", get_worker_num());
        //print_array(get_worker_num());
        assert(check_no_node_locks_held_by_me(task_id));
        double_list(task_id, sub_counts, density_count);
        // -1 for the lock at the begining of the function
        // we want to leave with the same number as we entered with
        //for (int i = 0; i < lock_count-1; i++) {
        //  edges.list_lock.lock_shared();
        //}
     // }
        assert(check_no_node_locks_for_me(task_id));
      return;
    }
    // printf("density %f, upperbound %f, len = %d, N = %d, logN = %d\n",
    // density, density_b.y, len, list->N, list->logN);
  }
  assert(((double)get_density_count(&edges, node_index, len))/ len == density);
  assert(density < density_b);
  assert(density <= (((double) edges.logN - 1) / edges.logN));
  //print_array(get_worker_num());
  if(len > edges.logN) {
     if (len <= REDISTRIBUTE_PAR_SIZE) {
      redistribute(node_index, len);
     } else {
      redistribute_par(node_index, len, sub_counts, density_count);
     }
    
  }
  //printf("this relase? %d, %d\n", node_index, len);
  assert(check_no_full_leaves(&edges, node_index, len));
  release_locks_in_range(task_id, held_locks, GENERAL);

  for (uint32_t i = lock_span.x; i < held_locks.x; i++) {
    nodes[i].lock.lock(task_id, REBALANCE);
    nodes[i].lock.unlock(task_id);
  }
  for (uint32_t i = held_locks.y+1; i <= lock_span.y; i++) {
    nodes[i].lock.lock(task_id, REBALANCE);
    nodes[i].lock.unlock(task_id);
  }
  // printf("node_index5 = %d\n", node_index);
  // print_array();
  // edges.list_lock.unlock_shared();
    assert(check_no_node_locks_for_me(task_id));
  return;
}

// remove elem at index
// assumes the lock on index is held by the parent
// and releases it when it is done with it
void OFM::remove(uint64_t task_id, uint64_t index, uint32_t elem_dest, uint32_t src, pair_int held_locks) {
  assert(check_every_lock_in_leaf(task_id, index));
  assert(check_no_full_leaves(&edges, find_leaf(&edges, index), edges.logN));
  
  /*
  if (index > 0 && edges.dests[index - 1] == NULL_VAL && find_leaf(&edges, index) != index) {
    printf("insert conditional || index = %d src = %u elem dest = %u, elem value = %u, next leaf: %d\n", index, src, elem_dest, elem_val, next_leaf(index, edges.loglogN));
    print_array(get_worker_num());
  }*/


  //printf("index = %d src = %u elem dest = %u, elem value = %u, worker %lu\n", index, src, elem.dest, elem.value, get_worker_num());
  //print_array(get_worker_num());
  uint64_t orig_n = edges.N;
  
  if ((elem_dest != SENT_VAL) ) {
    // -1 for the case that the location is the spot of the next sentinal
    assert(src == find_contaning_node(index-1));
  }
  // edges.list_lock.lock_shared();
  assert(nodes[src].lock.i_own_lock(task_id));
  uint32_t level = edges.H;
  uint64_t len = edges.logN;

  uint32_t * vals = (uint32_t *) edges.vals;
  uint32_t * dests = (uint32_t *) edges.dests;
  // always deposit on the left
  assert(index < edges.N - 1);
  // if the edge already exists in the graph, update its value
  // do not make another edge
  // slide right assumes we hold the lock till the end of the leaf
  assert(check_every_lock_in_leaf(task_id, index));
  slide_left(index, vals, dests);
  // printf("after sliding, index = %d\n", index);
  // print_array();
  assert(check_every_lock_in_leaf(task_id, index));
  // print_array();

  
  //assert(vals[index]!=0);
  //print_array();
  uint32_t node_index = find_leaf(&edges, index);
  double density = get_density(&edges, node_index, len);
  //printf("density = %f, %d\n", density, density == 1);

  // spill over into next level up, node is completely full.
  if (density == 0) {
    //printf("first rebalence\n");
    len*=2;
    level--;
    node_index = find_node(node_index, len);
  } else {
    redistribute(node_index, len);
  }
  assert(edges.N == orig_n);
  assert(check_every_lock_in_leaf(task_id, index));
  // being fancy here since some of the locks are done with and others might be used for the rebalence
  pair_int locks_for_rebalance = which_locks_in_range(node_index, len, src);
  //printf("relasing locks %d through %d with reason %d\n", locks_for_leaf.x, locks_for_rebalance.x -1, GENERAL);
  //printf("before with worker %lu\n", get_worker_num());
  //printf("worker %lu unlocking %d through %d with GENERAL \n", get_worker_num(), locks_for_leaf.x, locks_for_rebalance.x);
  cilk_for (uint32_t i = held_locks.x; i < locks_for_rebalance.x; i++) {
    nodes[i].lock.unlock(task_id);
  }
  //printf("worker %lu relasing locks %d through %d with reason REBALANCE\n",get_worker_num(), max(locks_for_rebalance.x, locks_for_leaf.x), locks_for_leaf.y);
  for (uint32_t i = max(locks_for_rebalance.x, held_locks.x) ; i <= min(held_locks.y, locks_for_rebalance.y); i++) {
    nodes[i].lock.unlock(task_id, REBALANCE);
  }
  //printf("worker %lu unlocking %d through %d with GENERAL \n", get_worker_num(), locks_for_rebalance.y +1,locks_for_leaf.y);

  for (uint32_t i = locks_for_rebalance.y +1; i <= held_locks.y; i++) {
    nodes[i].lock.unlock(task_id);
  }

  //printf("unlocking %d by worker %lu\n",src, get_worker_num());
  //printf("grabbing %d, %d\n", node_index, len);
  //printf("about to grab\n");
  pair_int lock_span = {max(locks_for_rebalance.x, held_locks.x),  max(locks_for_rebalance.y, held_locks.y)};
  assert(check_no_node_locks_held_by_me(task_id));

  held_locks = grab_locks_in_range_with_resets(task_id, node_index, len, REBALANCE, src);
  lock_span.x = min(held_locks.x, lock_span.x);
  lock_span.y = max(held_locks.y, lock_span.y);
  // if somebody doubled and we let them with a reset
  if (edges.N != orig_n) {
    //printf("worker %lu: somebody else doubled on us so we assume it has been distributed properly 1\n", get_worker_num());
    //printf("worker %lu: edges.N = %u, orig_n = %u\n",get_worker_num(), edges.N, orig_n);
    release_locks_in_range(task_id, held_locks, GENERAL);
    for (uint32_t i = lock_span.x; i <= lock_span.y; i++) {
      if (i < held_locks.x || i > held_locks.y) {
        nodes[i].lock.lock(task_id, REBALANCE);
        nodes[i].lock.unlock(task_id);
      }
    }
    assert(check_no_node_locks_for_me(task_id));
    return;
  }
  // printf("node_index3 = %d\n", node_index);
  // print_array();

  // get density of the leaf you are in
  double density_b = lower_density_bound[level];
  uint64_t density_count = get_density_count(&edges, node_index, len);
  density = ((double)density_count)/len;
  // printf("density %f, upperbound %f, len = %d, N = %d, logN = %d\n", density,
  // density_b.y, len, list->N, list->logN);

  // while density too high, go up the implicit tree
  // go up to the biggest node above the density bound
  //printf("node_index = %d, desnsity = %f, density bound = %f\n", node_index, density, density_b.y);
  std::vector<uint64_t> sub_counts(0);
  while (density <= density_b) {
    //printf("node_index = %d, desnsity = %f, density bound = %f, len = %d, worker = %lu\n", node_index, density, density_b.y, len, get_worker_num());
    len *= 2;
    if (len <= edges.N) {
      release_locks_in_range(task_id, held_locks, REBALANCE);
      level--;
      uint32_t new_node_index = find_node(node_index, len);
      held_locks = grab_locks_in_range_with_resets(task_id, new_node_index, len, REBALANCE, src);
      if (held_locks.x < lock_span.x) {
        lock_span.x = held_locks.x;
      }
      if (held_locks.y > lock_span.y) {
        lock_span.y = held_locks.y;
      }
      if (edges.N != orig_n) {
        //printf("worker %lu: somebody else doubled on us so we assume it has been distributed properly 2\n", get_worker_num());
        //printf("worker %lu: edges.N = %u, orig_n = %u\n",get_worker_num(), edges.N, orig_n);
        release_locks_in_range(task_id, held_locks, GENERAL);
        for (uint32_t i = lock_span.x; i <= lock_span.y; i++) {
          if (i < held_locks.x || i > held_locks.y) {
            nodes[i].lock.lock(task_id, REBALANCE);
            nodes[i].lock.unlock(task_id);
          }
        }
        assert(check_no_node_locks_for_me(task_id));
        return;
      }
      if (len <= REDISTRIBUTE_PAR_SIZE) {
        density_count = get_density_count(&edges, new_node_index, len);
      } else {
        sub_counts.resize(len/REDISTRIBUTE_PAR_SIZE);
        density_count = get_density_count_par(&edges, new_node_index, len, sub_counts);
      }
      // to help prevent double doubling by knowing how big it was on the last count
      orig_n = edges.N;
      node_index = new_node_index;
      density_b = lower_density_bound[level];
      density = ((double) density_count )/len;
    } else {
      // if you reach the root, double the list
      release_locks_in_range(task_id, held_locks, DOUBLE);
      //if (edges.N == orig_n) {
        // for (int i = 0; i < lock_count; i++) {
        //  edges.list_lock.unlock_shared();
        //}
        //TODO don't double double
        //printf("second double due to worker %lu\n", get_worker_num());
        //print_array(get_worker_num());
        assert(check_no_node_locks_held_by_me(task_id));
        half_list(task_id, sub_counts, density_count);
        // -1 for the lock at the begining of the function
        // we want to leave with the same number as we entered with
        //for (int i = 0; i < lock_count-1; i++) {
        //  edges.list_lock.lock_shared();
        //}
     // }
        assert(check_no_node_locks_for_me(task_id));
      return;
    }
    // printf("density %f, upperbound %f, len = %d, N = %d, logN = %d\n",
    // density, density_b.y, len, list->N, list->logN);
  }
  assert(((double)get_density_count(&edges, node_index, len))/ len == density);
  assert(density > density_b);
  //print_array(get_worker_num());
  if(len > edges.logN) {
     if (len <= REDISTRIBUTE_PAR_SIZE) {
      redistribute(node_index, len);
     } else {
      redistribute_par(node_index, len, sub_counts, density_count);
     }
    
  }
  //printf("this relase? %d, %d\n", node_index, len);
  assert(check_no_full_leaves(&edges, node_index, len));
  release_locks_in_range(task_id, held_locks, GENERAL);

  for (uint32_t i = lock_span.x; i < held_locks.x; i++) {
    nodes[i].lock.lock(task_id, REBALANCE);
    nodes[i].lock.unlock(task_id);
  }
  for (uint32_t i = held_locks.y+1; i <= lock_span.y; i++) {
    nodes[i].lock.lock(task_id, REBALANCE);
    nodes[i].lock.unlock(task_id);
  }
  // printf("node_index5 = %d\n", node_index);
  // print_array();
  // edges.list_lock.unlock_shared();
    assert(check_no_node_locks_for_me(task_id));
  return;
}




#include <immintrin.h>
#include <iostream>
#include <iomanip>    

template<class T> inline void Log(const __m256i & value) {
    const size_t n = sizeof(__m256i) / sizeof(T);
    T buffer[n];
    _mm256_storeu_si256((__m256i*)buffer, value);
    for (size_t i = 0; i < n; i++)
        std::cout << buffer[i] << " ";
    printf("\n");
}


/*
//I don't think this is quite correct, but it did work at some point
std::vector<uint32_t>
OFM::sparse_matrix_vector_multiplication(std::vector<uint32_t> const &v) {
  std::vector<uint32_t> result(nodes.size(), 0);

  uint32_t num_vertices = nodes.size();
  
  bool vector = true;

  cilk_for (uint32_t i = 0; i < num_vertices; i++) {
    nodes[i].lock.lock_shared();
    // printf("i: %d\n", i);
    // +1 to avoid sentinel
    uint32_t j = nodes[i].beginning + 1;
    __m256i temp_sum = _mm256_setzero_si256();
    while (j < nodes[i].end) {
       if (!vector || (nodes[i].end - j < 4 || j >= next_leaf(j, edges.logN)-4)) { 
        if (!is_null(edges.items[j].e)) {
          result[i] += edges.items[j].e.value * v[edges.items[j].e.dest];
          j++;
        } else {
          j = next_leaf(j, edges.logN);
        }
      
      } else {
        __m256i edgegroup = _mm256_loadu_si256((__m256i*) &edges.items[j]);
        __m256i values_group = _mm256_srli_epi64(edgegroup, 32);
        __m256i mask = _mm256_setr_epi64x(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);
        __m256i dest_group_targets = _mm256_and_si256(edgegroup, mask);
        __m256i dest_group = _mm256_i32gather_epi32((int*)v.data(), dest_group_targets, sizeof(v[0]));
        __m256i partial_sums = _mm256_mul_epi32(dest_group, values_group);
        temp_sum = _mm256_add_epi64(temp_sum, partial_sums);
        //temp_sums = _mm256_and_si256(temp_sums, mask);    


        __m256i mask2 = _mm256_setr_epi64x(0xFFFFFFFFFFFFFFFFUL, 0xFFFFFFFFFFFFFFFFUL, 0xFFFFFFFFFFFFFFFFUL, 0xFFFFFFFF00000000UL);
        if (_mm256_testc_si256(mask2, edgegroup)) {
          j = next_leaf(j, edges.logN);
        } else {
          j+=4;
        }
      }
    }
    if (vector) {
      result[i] += _mm256_extract_epi32(temp_sum, 0);
      result[i] += _mm256_extract_epi32(temp_sum, 2);
      result[i] += _mm256_extract_epi32(temp_sum, 4);
      result[i] += _mm256_extract_epi32(temp_sum, 6);
    } 
  } 
  return result;
}
*/

void OFM::print_graph() {
  uint32_t num_vertices = nodes.size();
  for (uint32_t i = 0; i < num_vertices; i++) {
    // printf("i: %d\n", i);
    // +1 to avoid sentinel
    uint32_t matrix_index = 0;

    for (uint64_t j = nodes[i].beginning + 1; j < nodes[i].end; j++) {
      if (edges.dests[j]!=NULL_VAL) {
        while (matrix_index < edges.dests[j]) {
          printf("000 ");
          matrix_index++;
        }
        printf("%03d ", edges.vals[j]);
        matrix_index++;
      }
    }
    for (uint64_t j = matrix_index; j < num_vertices; j++) {
      printf("000 ");
    }
    printf("\n");
  }
}

void OFM::add_node() {
  uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2);
  node_lock.lock(task_id);
  node_t node;
  uint32_t len = nodes.size();
  edge_t sentinel;
  sentinel.dest = SENT_VAL; // placeholder
  sentinel.value = len;       // back pointer

  if (len > 0) {
    node.beginning = nodes[len - 1].end;
    if(node.beginning == edges.N - 1) {
      uint32_t volatile  * volatile dests = edges.dests;
      uint32_t leaf = find_leaf(&edges, node.beginning);
      // moving the beginning of the node you are inserting left
      //TODO jump back to the last leaf and look forward from there
      if (nodes[len-1].num_neighbors == 0) {
        node.beginning = nodes[len - 1].beginning + 1;
      } else {
        while(dests[node.beginning - 1] == NULL_VAL && node.beginning != leaf /* && next_leaf(node.beginning, edges.loglogN) != node.beginning*/) {
          node.beginning -= 1;
        }
      }
    }
    node.end = node.beginning + 1;
    // fix previous node to set its end to your beginning
    nodes[len - 1].end = node.beginning;
  } else {
    node.beginning = 0;
    node.end = 1;
    // be special to the first one since we know it never moves do it doesn't need to look like a sentinal since it doesn't have to be fixed ever
    sentinel.value = NULL_VAL;
    sentinel.dest = 0;
  }
  // printf("sentinel dest: %u, value: %u\n", sentinel.dest, sentinel.value);
  node.num_neighbors = 0;
  node.lock.name("nodeLock");
  node.lock.number(nodes.size());

  nodes.push_back(node);
  pair_int held_locks = grab_locks_for_leaf_with_resets(task_id, nodes.size() - 1);
  // edges.list_lock.lock_shared();
  node_lock.make_shared(task_id);
  uint32_t loc = node.beginning;
  tbassert(loc < edges.N, "loc: %d, edges.N: %d\n", loc, edges.N);
  insert(task_id, loc, sentinel.dest, sentinel.value, nodes.size() - 1, held_locks);
  // printf("end of insert reason: %d\n", edges.lock_array[find_lock_index(&edges, loc)].reason);
  // edges.list_lock.unlock_shared();
  node_lock.unlock_shared(task_id);
  // print_array(len);
}

// add a node to the graph
void OFM::add_node_fast(uint64_t task_id) {
  node_lock.lock(task_id);
  node_t node;
  uint32_t len = nodes.size();
  edge_t sentinel;
  sentinel.dest = SENT_VAL; // placeholder
  sentinel.value = len;       // back pointer

  if (len > 0) {
    node.beginning = nodes[len - 1].end;
    if(node.beginning == edges.N - 1) {
      uint32_t volatile  * volatile dests = edges.dests;
      uint32_t leaf = find_leaf(&edges, node.beginning);
      // moving the beginning of the node you are inserting left
      //TODO jump back to the last leaf and look forward from there
      if (nodes[len-1].num_neighbors == 0) {
        node.beginning = nodes[len - 1].beginning + 1;
      } else {
        while(dests[node.beginning - 1] == NULL_VAL && node.beginning != leaf /* && next_leaf(node.beginning, edges.loglogN) != node.beginning*/) {
          node.beginning -= 1;
        }
      }
    }
    node.end = node.beginning + 1;
    // fix previous node to set its end to your beginning
    nodes[len - 1].end = node.beginning;
  } else {
    node.beginning = 0;
    node.end = 1;
    // be special to the first one since we know it never moves do it doesn't need to look like a sentinal since it doesn't have to be fixed ever
    sentinel.value = NULL_VAL;
    sentinel.dest = 0;
  }
  // printf("sentinel dest: %u, value: %u\n", sentinel.dest, sentinel.value);
  node.num_neighbors = 0;
  node.lock.name("nodeLock");
  node.lock.number(nodes.size());

  nodes.push_back(node);
  pair_int held_locks = grab_locks_for_leaf_with_resets(task_id, nodes.size() - 1);
  // edges.list_lock.lock_shared();
  node_lock.make_shared(task_id);
  uint32_t loc = node.beginning;
  tbassert(loc < edges.N, "loc: %d, edges.N: %d\n", loc, edges.N);
  insert(task_id, loc, sentinel.dest, sentinel.value, nodes.size() - 1, held_locks);
  // printf("end of insert reason: %d\n", edges.lock_array[find_lock_index(&edges, loc)].reason);
  // edges.list_lock.unlock_shared();
  node_lock.unlock_shared(task_id);
  // print_array(len);
}

//TODO deal with the case that multiple threads do the binary search and try and make the same region exclusive
void OFM::add_edge(uint32_t src, uint32_t dest, uint32_t value) {
  // printf("src = %d, dest = %d, val = %d\n", src,dest,value);
  if (value != 0) {
    uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2);
    assert(check_no_locks_for_me(task_id));
    node_lock.lock_shared(task_id);
    pair_int held_locks = grab_locks_for_leaf_with_resets(task_id, src);
    assert(nodes[src].lock.i_own_lock(task_id));
    node_t node = nodes[src];


    edge_t e;
    e.dest = dest;
    e.value = value;
    
    // edges.list_lock.lock_shared();
    //uint32_t node_begin = node.beginning;
    //uint32_t node_end = node.end;
    __sync_fetch_and_add(&nodes[src].num_neighbors,1);
    // printf("looking in region %u, %u by worker %lu\n", node.beginning + 1, node.end, get_worker_num());
    uint64_t loc_to_add =
        binary_search(&edges, e.dest, e.value, node.beginning + 1, node.end);

    pair_int needed_locks = which_locks_in_range(find_leaf(&edges, loc_to_add), edges.logN, src);
    for (uint32_t i = held_locks.x; i < needed_locks.x; i++) {
      nodes[i].lock.unlock(task_id, GENERAL);
      held_locks.x = i+1;
    }
    for (uint32_t i = needed_locks.y+1; i <= held_locks.y; i++) {
      nodes[i].lock.unlock(task_id, GENERAL);
      
    }
    if (needed_locks.y+1 <= held_locks.y) {
      held_locks.y = needed_locks.y;
    }
    // print_array();
    //printf("loc_to_add: %u by worker %lu\n", loc_to_add, get_worker_num());
    // printf("src: %d, dest: %d by worker %lu\n", src, dest, get_worker_num());
    // print_array();

    
    assert(nodes[src].lock.i_own_lock(task_id));
    assert(check_every_lock_in_leaf(task_id, loc_to_add));
    /*
    if (!check_no_full_leaves(&edges, find_leaf(&edges, loc_to_add), edges.logN)) {
      printf("worker %lu is trying to insert into a full leaf\n", get_worker_num());
      print_array(get_worker_num());
    }
    */
    assert(check_no_full_leaves(&edges, find_leaf(&edges, loc_to_add), edges.logN));
    insert(task_id, loc_to_add, e.dest, e.value, src, held_locks);
    //printf("worker %lu done\n", get_worker_num());
    // edges.list_lock.unlock_shared();
    //print_array();
    node_lock.unlock_shared(task_id);
    assert(check_no_locks_for_me(task_id));
  }
}

void OFM::remove_edge(uint32_t src, uint32_t dest) {
  //printf(" trying to remove src = %d, dest = %d\n", src,dest);
  uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2);
  assert(check_no_locks_for_me(task_id));
  node_lock.lock_shared(task_id);
  pair_int held_locks = grab_locks_for_leaf_with_resets(task_id, src);
  assert(nodes[src].lock.i_own_lock(task_id));
  node_t node = nodes[src];

  
  // edges.list_lock.lock_shared();
  //uint32_t node_begin = node.beginning;
  //uint32_t node_end = node.end;
  // printf("looking in region %u, %u by worker %lu\n", node.beginning + 1, node.end, get_worker_num());
  uint64_t loc_to_remove =
      binary_search(&edges, dest, 0, node.beginning + 1, node.end);
  if (edges.dests[loc_to_remove] != dest) {
    //if the edge isn't there
    //printf("not removed\n");
    release_locks_in_range(task_id, held_locks);
    node_lock.unlock_shared(task_id);
    return;
  }
  //printf("removed\n");
  __sync_fetch_and_add(&nodes[src].num_neighbors,-1);
  pair_int needed_locks = which_locks_in_range(find_leaf(&edges, loc_to_remove), edges.logN, src);
  for (uint32_t i = held_locks.x; i < needed_locks.x; i++) {
    nodes[i].lock.unlock(task_id, GENERAL);
    held_locks.x = i+1;
  }
  for (uint32_t i = needed_locks.y+1; i <= held_locks.y; i++) {
    nodes[i].lock.unlock(task_id, GENERAL);
    
  }
  if (needed_locks.y+1 <= held_locks.y) {
    held_locks.y = needed_locks.y;
  }
  // print_array();
  //printf("loc_to_add: %u by worker %lu\n", loc_to_add, get_worker_num());
  // printf("src: %d, dest: %d by worker %lu\n", src, dest, get_worker_num());
  // print_array();

  
  assert(nodes[src].lock.i_own_lock(task_id));
  assert(check_every_lock_in_leaf(task_id, loc_to_add));
  /*
  if (!check_no_full_leaves(&edges, find_leaf(&edges, loc_to_add), edges.logN)) {
    printf("worker %lu is trying to insert into a full leaf\n", get_worker_num());
    print_array(get_worker_num());
  }
  */
  assert(check_no_full_leaves(&edges, find_leaf(&edges, loc_to_add), edges.logN));
  remove(task_id, loc_to_remove, dest, src, held_locks);
  //printf("worker %lu done\n", get_worker_num());
  // edges.list_lock.unlock_shared();
  //print_array();
  node_lock.unlock_shared(task_id);
  assert(check_no_locks_for_me(task_id));
}

void OFM::remove_edge_fast(uint32_t src, uint32_t dest, uint64_t task_id) {
  //printf(" trying to remove src = %d, dest = %d\n", src,dest);
  pair_int held_locks = grab_locks_for_leaf_with_resets(task_id, src);
  node_t node = nodes[src];

  
  // edges.list_lock.lock_shared();
  //uint32_t node_begin = node.beginning;
  //uint32_t node_end = node.end;
  // printf("looking in region %u, %u by worker %lu\n", node.beginning + 1, node.end, get_worker_num());
  uint64_t loc_to_remove =
      binary_search(&edges, dest, 0, node.beginning + 1, node.end);
  if (edges.dests[loc_to_remove] != dest) {
    //if the edge isn't there
    //printf("not removed\n");
    release_locks_in_range(task_id, held_locks);
    return;
  }
  //printf("removed\n");
  __sync_fetch_and_add(&nodes[src].num_neighbors,-1);
  pair_int needed_locks = which_locks_in_range(find_leaf(&edges, loc_to_remove), edges.logN, src);
  for (uint32_t i = held_locks.x; i < needed_locks.x; i++) {
    nodes[i].lock.unlock(task_id, GENERAL);
    held_locks.x = i+1;
  }
  for (uint32_t i = needed_locks.y+1; i <= held_locks.y; i++) {
    nodes[i].lock.unlock(task_id, GENERAL);
    
  }
  if (needed_locks.y+1 <= held_locks.y) {
    held_locks.y = needed_locks.y;
  }
  // print_array();
  //printf("loc_to_add: %u by worker %lu\n", loc_to_add, get_worker_num());
  // printf("src: %d, dest: %d by worker %lu\n", src, dest, get_worker_num());
  // print_array();

  
  assert(nodes[src].lock.i_own_lock(task_id));
  assert(check_every_lock_in_leaf(task_id, loc_to_add));
  /*
  if (!check_no_full_leaves(&edges, find_leaf(&edges, loc_to_add), edges.logN)) {
    printf("worker %lu is trying to insert into a full leaf\n", get_worker_num());
    print_array(get_worker_num());
  }
  */
  assert(check_no_full_leaves(&edges, find_leaf(&edges, loc_to_add), edges.logN));
  remove(task_id, loc_to_remove, dest, src, held_locks);
  //printf("worker %lu done\n", get_worker_num());
  // edges.list_lock.unlock_shared();
  //print_array();
}


//TODO can't set edges to zero
void OFM::add_edge_update(uint32_t src, uint32_t dest, uint32_t value) {
  if (value != 0) {
    uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2);
    assert(check_no_locks_for_me(task_id));
    node_lock.lock_shared(task_id);
    pair_int held_locks = grab_locks_for_leaf_with_resets(task_id, src);
    assert(nodes[src].lock.i_own_lock(task_id));
    node_t node = nodes[src];


    edge_t e;
    e.dest = dest;
    e.value = value;
    
    // edges.list_lock.lock_shared();
    //uint32_t node_begin = node.beginning;
    //uint32_t node_end = node.end;
    // printf("looking in region %u, %u by worker %lu\n", node.beginning + 1, node.end, get_worker_num());
    uint64_t loc_to_add =
        binary_search(&edges, e.dest, e.value, node.beginning + 1, node.end);


    if (edges.dests[loc_to_add] == dest) {
      edges.vals[loc_to_add] = value;
      release_locks_in_range(task_id, held_locks);
      node_lock.unlock_shared(task_id);
      return;
    }
    pair_int needed_locks = which_locks_in_range(find_leaf(&edges, loc_to_add), edges.logN, src);
    for (uint32_t i = held_locks.x; i < needed_locks.x; i++) {
      nodes[i].lock.unlock(task_id, GENERAL);
      held_locks.x = i+1;
    }
    for (uint32_t i = needed_locks.y+1; i <= held_locks.y; i++) {
      nodes[i].lock.unlock(task_id, GENERAL);
      
    }
    if (needed_locks.y+1 <= held_locks.y) {
      held_locks.y = needed_locks.y;
    }
    // printf("loc_to_add: %d\n", loc_to_add);
    // printf("src: %d, dest: %d\n", src, dest);
    __sync_fetch_and_add(&nodes[src].num_neighbors,1);
    // print_array();
    insert(task_id, loc_to_add, e.dest, e.value, src, held_locks);
    // edges.list_lock.unlock_shared();
    // print_array();
    node_lock.unlock_shared(task_id);
    assert(check_no_locks_for_me(task_id));
  }
}
void __attribute__ ((noinline)) OFM::add_edge_update_fast(uint32_t src, uint32_t dest, uint32_t value, uint64_t task_id) {
    //printf(" trying to add src = %d, dest = %d\n", src,dest);
    __builtin_prefetch(&nodes[src]);
    __builtin_prefetch(&edges);
    if (value != 0) {
      assert(check_no_locks_for_me(task_id));
      pair_int held_locks = grab_locks_for_leaf_with_resets(task_id, src);
      assert(nodes[src].lock.i_own_lock(task_id));
      node_t node = nodes[src];


      edge_t e;
      e.dest = dest;
      e.value = value;
      
      // edges.list_lock.lock_shared();
      //uint32_t node_begin = node.beginning;
      //uint32_t node_end = node.end;
      // printf("looking in region %u, %u by worker %lu\n", node.beginning + 1, node.end, get_worker_num());
      uint64_t loc_to_add =
          binary_search(&edges, e.dest, e.value, node.beginning + 1, node.end);


      if (edges.dests[loc_to_add] == dest) {
        //printf(" just updating value\n");
        edges.vals[loc_to_add] = value;
        release_locks_in_range(task_id, held_locks);
        return;
      }

      //printf("adding the edge\n");
      pair_int needed_locks = which_locks_in_range(find_leaf(&edges, loc_to_add), edges.logN, src);
      for (uint32_t i = held_locks.x; i < needed_locks.x; i++) {
        nodes[i].lock.unlock(task_id, GENERAL);
        held_locks.x = i+1;
      }
      for (uint32_t i = needed_locks.y+1; i <= held_locks.y; i++) {
        nodes[i].lock.unlock(task_id, GENERAL);
        
      }
      if (needed_locks.y+1 <= held_locks.y) {
        held_locks.y = needed_locks.y;
      }
      // printf("loc_to_add: %d\n", loc_to_add);
      // printf("src: %d, dest: %d\n", src, dest);
      __sync_fetch_and_add(&nodes[src].num_neighbors,1);
      // print_array();
      insert(task_id, loc_to_add, e.dest, e.value, src, held_locks);
      // edges.list_lock.unlock_shared();
      // print_array();
      assert(check_no_locks_for_me(task_id));
    }
}
void OFM::add_edge_batch_update_no_val(uint32_t *srcs, uint32_t *dests, uint32_t edge_count) {
  uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2*edge_count);
  node_lock.lock_shared(0);
  cilk_for (uint32_t i = 0; i < edge_count; i++) {
    uint32_t src = srcs[i];
    uint32_t dest = dests[i];
    add_edge_update_fast(src,dest,1,task_id+2*i);
  }
  node_lock.unlock_shared(0);
}
void OFM::remove_edge_batch(uint32_t *srcs, uint32_t *dests, uint32_t edge_count) {
  uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2*edge_count);
  node_lock.lock_shared(0);
  cilk_for (uint32_t i = 0; i < edge_count; i++) {
    uint32_t src = srcs[i];
    uint32_t dest = dests[i];
    remove_edge_fast(src,dest,task_id+2*i);
  }
  node_lock.unlock_shared(0);
}

void OFM::add_edge_batch_update(uint32_t *srcs, uint32_t *dests, uint32_t *values, uint32_t edge_count) {
  uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2*edge_count);
  node_lock.lock_shared(0);
  cilk_for (uint32_t i = 0; i < edge_count ; i++) {
    uint32_t src = srcs[i];
    uint32_t dest = dests[i];
    uint32_t value = values[i];
    add_edge_update_fast(src,dest,value,task_id+2*i);
  }
  node_lock.unlock_shared(0);
}


OFM::OFM(uint32_t init_n) {
  next_task_id = 1;
  //making sure logN is at least 4
  edges.N = max(2 << bsr_word(init_n*2), 16);
  // printf("%d\n", bsf_word(list->N));
  edges.loglogN = bsr_word(bsr_word(edges.N) + 1);
  edges.logN = (1 << edges.loglogN);
  edges.mask_for_leaf = ~(edges.logN - 1);
  assert(edges.logN > 0);
  edges.density_limit = ((double) edges.logN - 1)/edges.logN;
  edges.H = bsr_word(edges.N / edges.logN);
  for (uint32_t i = 0; i <= edges.H; i++) {
    upper_density_bound[i] = density_bound(&edges, i).y;
    lower_density_bound[i] = density_bound(&edges, i).x;
  }

  // printf("N = %d, logN = %d, loglogN = %d, H = %d\n", list->N, list->logN,
  // list->loglogN, list->H);
  
  // edges.list_lock.name("listLock");
  node_lock.name("nodeLock");  
  edges.vals = (uint32_t *)malloc(edges.N * sizeof(*(edges.vals)));
	if (edges.vals == NULL) {
		printf("bad malloc in create 1\n");
		exit(-1);
	}
  edges.dests = (uint32_t *)malloc(edges.N * sizeof(*(edges.dests)));
	if (edges.dests == NULL) {
		printf("bad malloc in create 2\n");
		exit(-1);
	}
  for (uint32_t i = 0; i < edges.N; i++) {
    edges.vals[i] = 0;
    edges.dests[i] = NULL_VAL;
  }
  //TODO might be an issue if we grow it one at a time and let nodes be moved during operation
  nodes.reserve(init_n);
  uint64_t task_id = __sync_fetch_and_add(&next_task_id, 2*init_n);
  for (uint32_t i = 0; i < init_n; i++) {
    add_node_fast(task_id + (i * 2));
  }
}

OFM::OFM(OFM &other) {
  nodes = other.nodes;
  next_task_id = 1;
  edges.N = other.edges.N;
  edges.loglogN = other.edges.loglogN;
  edges.logN = other.edges.logN;
  edges.mask_for_leaf = other.edges.mask_for_leaf;
  edges.density_limit = other.edges.density_limit;
  edges.H = other.edges.H;
  for (uint32_t i = 0; i <= edges.H; i++) {
    upper_density_bound[i] = other.upper_density_bound[i];
    lower_density_bound[i] = other.lower_density_bound[i];
  }
  node_lock.name("nodeLock");
  edges.vals = (uint32_t *)malloc(edges.N * sizeof(*(edges.vals)));
	if (edges.vals == NULL) {
		printf("bad malloc in create copy 1\n");
		exit(-1);
	}
  memcpy((void *)edges.vals, (void *)other.edges.vals, edges.N * sizeof(*(edges.vals)));
  edges.dests = (uint32_t *)malloc(edges.N * sizeof(*(edges.dests)));
	if (edges.dests == NULL) {
		printf("bad malloc in create copy 2\n");
		exit(-1);
	}
  memcpy((void *)edges.dests, (void *)other.edges.dests, edges.N * sizeof(*(edges.dests)));
}

//TODO free lock array
OFM::~OFM() { 
  free((void*)edges.vals);
  free((void*)edges.dests);
}
