#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <stdint.h>
#include <queue>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cilk/cilk.h>

#include "helpers.h"
#include "bitmap.h"
#include "platform_atomics.h"
#include "pvector.h"
#include "sliding_queue.h"
#include <cilk/reducer_opadd.h>
using namespace std;
class Graph
{
    public:

    virtual ~Graph() = 0;
    // no data
    virtual uint64_t get_n() = 0;
    virtual uint32_t find_value(uint32_t src, uint32_t dest) = 0;
    virtual vector<uint32_t> sparse_matrix_vector_multiplication(std::vector<uint32_t> const &v) = 0;
    virtual void print_graph() = 0;
    virtual void add_node() = 0;
    virtual void add_edge(uint32_t src, uint32_t dest, uint32_t value) = 0;
    virtual void add_edge_update(uint32_t src, uint32_t dest, uint32_t value) = 0;
    virtual uint64_t get_size() = 0;
    virtual uint64_t triangle_count() = 0;
    virtual vector<double> pagerank(std::vector<double> const &node_values) = 0;
    virtual vector<double> pagerank_pull(std::vector<double> &node_values) = 0;
    virtual vector<tuple<uint32_t, uint32_t, uint32_t> > get_edges() = 0;
    virtual void add_edge_batch_update(uint32_t *srcs, uint32_t *dests, uint32_t *values, uint32_t edge_count) = 0;
    virtual uint32_t num_neighbors(uint32_t node){return 0;};

    // takes in start node, returns vector of distances (length n)
    // max_int if not connected
    virtual vector<int32_t> bfs(int32_t start_node) = 0;
    virtual pvector<int32_t> parallel_bfs(int32_t source,int32_t total_edges, int alpha = 5, int beta = 18) = 0;
    virtual void convert(Graph* g) = 0;
    // virtual vector<uint32_t> find_neighbors(uint32_t v) = 0;

    void make_symetric() {
      for (uint32_t i = 0; i < get_n(); i++) {
        for (uint32_t j = 0; j < get_n(); j++) {
          uint32_t value = find_value(i,j);
          if (value > 0) {
            add_edge_update(i,j,value);      
            add_edge_update(j,i,value);
          }
        }
      }
    }
    
};

// needed so the correct destructor is called when using base pointer
Graph::~Graph() {

}
//TODO might be a bug in how we initialize the parent vector, gap has it each to its own negative degree, not sure it matters
#define PARALLEL_BFS \
pvector<int32_t> parallel_bfs(int32_t source, int32_t total_edges, int alpha, int beta) {\
  int n_workers = __cilkrts_get_nworkers();\
  pvector<int32_t> parent(get_n());\
  for (int32_t n = 0; n < get_n(); n++) {\
    parent[n] = num_neighbors(n) != 0 ? - num_neighbors(n) : -1;\
  }\
  parent[source] = source;\
  SlidingQueue<int32_t> queue(get_n());\
  queue.push_back(source);\
  queue.slide_window();\
  Bitmap curr(get_n());\
  curr.reset();\
  Bitmap front(get_n());\
  front.reset();\
  int64_t edges_to_check = total_edges;\
  int64_t scout_count = num_neighbors(source);\
  while (!queue.empty()) {\
    if (scout_count > edges_to_check / alpha) {\
      int64_t awake_count, old_awake_count;\
      cilk_for (int32_t i = queue.shared_out_start; i < queue.shared_out_end; i++) {\
        int32_t u = queue.shared[i];\
        front.set_bit_atomic(u);\
      }\
      awake_count = queue.size();\
      queue.slide_window();\
      do {\
        old_awake_count = awake_count;\
        awake_count = 0;\
        curr.reset();\
        vector<int64_t> awake_count_vector(n_workers, 0);\
        cilk_for (int32_t u=0; u < get_n(); u++) {\
          uint32_t worker_num = __cilkrts_get_worker_number();\
          if (parent[u] < 0) {\
            for (iterator it = begin(u); it != end(u); ++it) {\
              edge_t edge = *it;\
              int32_t v = edge.dest;\
              if (front.get_bit(v)) {\
                parent[u] = v;\
                awake_count_vector[worker_num]+=1;\
                curr.set_bit(u);\
                break;\
              }\
            }\
          }\
        }\
        for (auto &item : awake_count_vector) {\
          awake_count+=item;\
        }\
        front.swap(curr);\
      } while ((awake_count >= old_awake_count) || (awake_count > get_n() / beta));\
      QueueBuffer<int32_t> *queue_array = (QueueBuffer<int32_t> *)malloc(4*sizeof(QueueBuffer<int32_t>) * n_workers);\
      for (int i = 0; i < n_workers; i++) {\
        new(&queue_array[i*4]) QueueBuffer(queue);\
      }\
      cilk_for (int32_t n=0; n < get_n(); n++) {\
        if (front.get_bit(n)) {\
          queue_array[__cilkrts_get_worker_number()*4].push_back(n);\
        }\
      }\
      cilk_for (int i = 0; i < n_workers; i++) {\
        queue_array[i*4].flush();\
      }\
      queue.slide_window();\
      scout_count = 1;\
    } else {\
      edges_to_check -= scout_count;\
      scout_count = 0;\
      vector<int64_t> scout_count_vector(n_workers, 0);\
      QueueBuffer<int32_t> *queue_array = (QueueBuffer<int32_t> *)malloc(4*sizeof(QueueBuffer<int32_t>) * n_workers);\
      for (int i = 0; i < n_workers; i++) {\
        new(&queue_array[i*4]) QueueBuffer(queue);\
      }\
      cilk_for (int32_t i = queue.shared_out_start; i < queue.shared_out_end; i++) {\
        int32_t u = queue.shared[i];\
        uint32_t worker_num = __cilkrts_get_worker_number();\
        for (iterator it = begin(u); it != end(u); ++it) {\
          edge_t edge = *it;\
          int32_t v = edge.dest;\
          int32_t curr_val = parent[v];\
          if (curr_val < 0) {\
            if (compare_and_swap(parent[v], curr_val, u)) {\
              queue_array[4*worker_num].push_back(v);\
              scout_count_vector[worker_num] += -curr_val;\
            }\
          }\
        }\
      }\
      cilk_for (int i = 0; i < n_workers; i++) {\
        queue_array[i*4].flush();\
      }\
      for (auto &item : scout_count_vector) {\
        scout_count+=item;\
      }\
      queue.slide_window();\
    }\
  }\
  for (int32_t n = 0; n < get_n(); n++) {\
    if (parent[n] < -1) {\
      parent[n] = -1;\
    }\
  }\
  return parent;\
}
  
//should technially lock in here if we want to parallel read and write
//not needed if we only want parallel read or parallel write, but not both at the same time
#define BFS \
vector<int32_t> bfs(int32_t start_node) {\
  uint64_t n = get_n();\
  vector<int32_t> out (n, -1);\
  queue<uint32_t> next;\
  next.push(start_node);\
  out[start_node] = start_node;\
  while(next.size() > 0) {\
    uint32_t active = next.front();\
    next.pop();\
    for (iterator it = begin(active); it != end(active); ++it) {\
      edge_t edge = *it;\
      if (out[edge.dest] == -1) {\
        next.push(edge.dest);\
        out[edge.dest] = active;\
      }\
    }\
  }\
  return out;\
}

#define PAGERANK \
vector<double> pagerank(std::vector<double> const &node_values) {\
  uint64_t n = get_n();\
  vector<double> output(n, 0);\
  int workers = __cilkrts_get_nworkers();\
  vector<double> output_par_(n*workers, 0);\
  double *output_par = output_par_.data();\
  cilk_for(int i = 0; i < n; i++) {\
    int worker_num = __cilkrts_get_worker_number();\
    if (num_neighbors(i) > 0) {\
      double contrib = (node_values[i] / num_neighbors(i));\
      for (iterator it = begin(i); it != end(i); ++it) {\
        edge_t edge = *it;\
        output_par[worker_num*n + edge.dest] += contrib;\
      }\
    }\
  }\
  cilk_for (uint32_t i = 0; i < n; i++) {\
    for (int j = 0; j < workers; j++) {\
      output[i] += output_par[j*n + i];\
    }\
  }\
  return output;\
}\
vector<double> pagerank_pull(std::vector<double> &node_values) {\
  uint64_t n = get_n();\
  cilk_for(int i = 0; i < n; i++) {\
    node_values[i] = node_values[i]/num_neighbors(i);\
  }\
  vector<double> output(n, 0);\
  cilk_for(int i = 0; i < n; i++) {\
    for (iterator it = begin(i); it != end(i); ++it) {\
      edge_t edge = *it;\
      output[i] += node_values[edge.dest];\
    }\
  }\
  return output;\
}\
vector<double> pagerank_simple(std::vector<double> const &node_values) {\
  uint64_t n = get_n();\
  vector<double> output(n, 0);\
  for(int i = 0; i < n; i++) {\
    if (num_neighbors(i) > 0) {\
      double contrib = (node_values[i] / num_neighbors(i));\
      for (iterator it = begin(i); it != end(i); ++it) {\
        edge_t edge = *it;\
        output[edge.dest] += contrib;\
      }\
    }\
  }\
  return output;\
}\





#define SPMV \
std::vector<uint32_t> sparse_matrix_vector_multiplication(std::vector<uint32_t> const &v) {\
    std::vector<uint32_t> result(get_n(), 0);\
    cilk_for(int i = 0; i < get_n(); i++) {\
      uint32_t temp = 0;\
      for (iterator it = begin(i); it != end(i); ++it) {\
        edge_t edge = *it;\
          temp += edge.value * v[edge.dest];\
      }\
      result[i] = temp;\
    }\
    return result;\
}


//really bad algorythm, but I am not sure what to do or the unsorted data structures
#define TRIANGLE_COUNT \
uint64_t triangle_count() {\
  uint64_t total = 0;\
  for(uint32_t i = 0; i < get_n(); i++) {\
    for(uint32_t j = 0; j < get_n(); j++) {\
      if (j > i) break;\
      for(uint32_t k = 0; k < get_n(); k++) {\
        if (k > j) break;\
        if (find_value(i,j) > 0 && find_value(j,k) > 0 && find_value(k,i) > 0) {\
          total++;\
        }\
      }\
    }\
  }\
  return total;\
}
//ONLY WORKS FOR SORTED NEIGHBORS
#define TRIANGLE_COUNT_SORTED \
uint64_t triangle_count() {\
  cilk::reducer< cilk::op_add<uint64_t> > total(0);\
  cilk_for(uint32_t u = 0; u < get_n(); u++) {\
    uint64_t temp = 0;\
    for (iterator it_u = begin(u); it_u != end(u); ++it_u) {\
      uint32_t v = (*it_u).dest;\
      if (v > u) break;\
      iterator it = begin(u);\
      for (iterator it_v = begin(v); it_v != end(v); ++it_v) {\
        uint32_t w = (*it_v).dest;\
        if (w > v) break;\
        while ((*it).dest < w) ++it;\
        if (w == (*it).dest) {\
          temp++;\
        }\
      }\
    }\
    *total+=temp;\
  }\
  return total.get_value();\
}

vector<string> split(string &s, char delim) {
    std::stringstream ss(s);
    std::string item;
    vector<string> result;
    while (std::getline(ss, item, delim)) {
	result.push_back(item);
    }
    return result;
}
#endif
