/*
 * adjacency matrix
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <vector>
#include <cstring>
#include "helpers.h"

#include "Graph.hpp"

class CSR : public Graph {
    public:
    // data members
    uint32_t *nodes;
    uint32_t *edges;
    uint32_t *values;
    SpinLock *locks;
    uint32_t N;
    uint32_t M;
 
    // function headings
    CSR(uint32_t init_n);
    ~CSR();

    uint64_t get_size();    
    uint64_t get_size_ideal();    
    uint64_t get_n();
    uint32_t find_value(uint32_t src, uint32_t dest);
    uint32_t find_value_index(uint32_t src, uint32_t dest);
    void print_graph();
    void print_arrays();
    void add_node();
    void add_edge(uint32_t src, uint32_t dest, uint32_t value);
    void add_edge_update(uint32_t src, uint32_t dest, uint32_t value);
    void convert(Graph* g);
    void add_edge_batch_update(uint32_t *srcs, uint32_t *dests, uint32_t *vals, uint32_t edge_count);

    void add_file3(string filename);
    vector<tuple<uint32_t, uint32_t, uint32_t> > get_edges() {
      printf("get_edges not implemented for csr\n");
      exit(1);
    }

    uint32_t num_neighbors(uint32_t node) {
      uint32_t start = nodes[node];
      uint32_t end;
      if (node < N - 1) { end = nodes[node+1]; }
      else { end = M; }
      return end - start;
    }

    class iterator {
    public:
      uint32_t index;
      uint32_t *edges;
      uint32_t *values;
      iterator(CSR *G, uint32_t node, bool start) {
        if (!start) {
          uint32_t end;
          if (node < G->N - 1) { end = G->nodes[node+1]; }
          else { end = G->M; }
          index = end;
          return;
        }
        index = G->nodes[node];
        edges = G->edges;
        values = G->values;
        return;
      }
      bool operator==(const iterator& other) const {
        return index == other.index;
      }
      bool operator!=(const iterator& other) const {
        return index != other.index;
      }
      iterator& operator++() {
        index+=1;
        return *this;
      }
      edge_t operator*() const {
        return {values[index], edges[index]};
      }
    };
    iterator begin(uint32_t node) {
      return iterator(this, node, true);
    }   
    iterator end(uint32_t node) {
      return iterator(this, node, false);
    }   

 
    BFS
    PAGERANK
    SPMV
    TRIANGLE_COUNT_SORTED
    PARALLEL_BFS
};  

// for soc-
// starting at 1
void CSR::add_file3(string filename) {
  vector<tuple<uint32_t, uint32_t, uint32_t>> edges_to_add;
  ifstream myfile(filename.c_str()); 
  string line;
  if (myfile.is_open()) {
    while ( getline (myfile,line) ) {
      vector<string> elems = split(line, '\t');
        int src = atoi(elems[0].c_str())-1;
        int dest = atoi(elems[1].c_str())-1;
        edges_to_add.push_back(make_tuple( src, dest, 1 ));
        // if (line_num++ > 400000000) {
        //  break;
        // }
    }
  myfile.close();
  // return 0;
  } else {
    printf("file was not opened\n");
  }
    // populate edges
  sort(edges_to_add.begin(), edges_to_add.end());

  uint32_t current_node = 0;
  for(uint32_t i = 0; i < edges_to_add.size(); i++) {
    uint32_t src = get<0>(edges_to_add[i]);
    uint32_t dest = get<1>(edges_to_add[i]);
    uint32_t val = get<2>(edges_to_add[i]);
    while (src >= current_node) {
      add_node();
      current_node++;
    }
    edges = (uint32_t *) realloc(edges, (M+1)*sizeof(uint32_t));
    values = (uint32_t *) realloc(values, (M+1)*sizeof(uint32_t));
    edges[M] = dest;
    values[M] = val;
    M += 1;
  }
}

void CSR::convert(Graph* g) {
  free(nodes);
  free(edges);
  free(values);
  nodes = nullptr;
  edges = nullptr;
  values = nullptr;
  vector<tuple<uint32_t, uint32_t, uint32_t> > edges_to_add = g->get_edges();
  N = 0;
  M = 0;
  // populate edges
  sort(edges_to_add.begin(), edges_to_add.end());

  uint32_t current_node = 0;
  for(uint32_t i = 0; i < edges_to_add.size(); i++) {
    uint32_t src = get<0>(edges_to_add[i]);
    uint32_t dest = get<1>(edges_to_add[i]);
    uint32_t val = get<2>(edges_to_add[i]);
    while (src >= current_node) {
      add_node();
      current_node++;
    }
    edges = (uint32_t *) realloc(edges, (M+1)*sizeof(uint32_t));
    values = (uint32_t *) realloc(values, (M+1)*sizeof(uint32_t));
    edges[M] = dest;
    values[M] = val;
    M += 1;
  }
}

uint64_t CSR::get_n() {
  return N;
}

uint64_t CSR::get_size() {
    return (N + M + M * sizeof(uint32_t));
}

void CSR::print_arrays() {
  printf("NODES\n");
  for(uint32_t i = 0; i < N; i++) {
    printf("%d ", nodes[i]);
  }
  printf("\nEDGES\n");
  for(uint32_t i = 0; i < M; i++) {
    printf("%d ", edges[i]);
  }
  printf("\nVALUES\n");
  for(uint32_t i = 0; i < M; i++) {
    printf("%d ", values[i]);
  }
  printf("\n");
}

// find value of (src, dest)
uint32_t CSR::find_value(uint32_t src, uint32_t dest) {
    uint32_t end = 0;
    if (src ==  N - 1) { end = M; }
    else { end = nodes[src + 1]; }
    uint32_t start = nodes[src];
    if(nodes[src] == end) {
      return 0;
    }

    if(edges[start] == dest) {
      return values[start];
    }

    while (start + 1 < end) { 
	uint32_t mid = (start + end) / 2;
	if (edges[mid] == dest) { return values[mid]; }
	else if (edges[mid] > dest) {
	  end = mid;
	}
	else {
	  start = mid;
        }
    } 
    // not found
    return 0;
}

uint32_t CSR::find_value_index(uint32_t src, uint32_t dest) {
    uint32_t end = 0;
    if (src ==  N - 1) { end = M; }
    else { end = nodes[src + 1]; }
    uint32_t start = nodes[src];
    if(nodes[src] == end) {
      return UINT32_MAX;
    }

    if(edges[start] == dest) {
      return start;
    }

    while (start + 1 < end) { 
      uint32_t mid = (start + end) / 2;
      if (edges[mid] == dest) { 
        return mid; 
      }
      else if (edges[mid] > dest) {
        end = mid;
      }
      else {
        start = mid;
      }
    } 
    // not found
    return UINT32_MAX;
}

void CSR::add_node() {
    SpinLock lock;
    locks = (SpinLock *) realloc(locks, (N+1)*sizeof(SpinLock));
    nodes = (uint32_t *) realloc(nodes, (N+1)*sizeof(uint32_t));
    locks[N] = lock;
    nodes[N] = M;
    N += 1;
}

// src, dest < N
void CSR::add_edge(uint32_t src, uint32_t dest, uint32_t value) {
    uint32_t end = 0;
    if (src ==  N - 1) { end = M; }
    else { end = nodes[src + 1]; }
    uint32_t start = nodes[src];
    //printf("adding edge (%u, %u, %u): start = %d, end = %d\n", src, dest, value, start, end);
    if(M == 0) {
      edges = (uint32_t *) realloc(edges, (M+1)*sizeof(uint32_t));
      values = (uint32_t *) realloc(values, (M+1)*sizeof(uint32_t));
      edges[M] = dest;
      values[M] = value;
      M += 1;
      for(uint32_t i = src + 1; i < N; i++) {
        nodes[i]++;
      }
      return;
    }
    if (dest < edges[start] || start == end) { 
            edges = (uint32_t *) realloc(edges, (M+1)*sizeof(uint32_t));
            values = (uint32_t *) realloc(values, (M+1)*sizeof(uint32_t));
            memmove(edges + start + 1, edges + start, (M - start)*sizeof(uint32_t));
            edges[start] = dest;
            memmove(values + start + 1, values + start, (M - start)*sizeof(uint32_t));
            values[start] = value;
            M += 1;
        for(uint32_t i = src + 1; i < N; i++) {
          nodes[i]++;
        }
	      return;
    } else if (end == 0 || dest > edges[end - 1]) {
            edges = (uint32_t *) realloc(edges, (M+1)*sizeof(uint32_t));
            values = (uint32_t *) realloc(values, (M+1)*sizeof(uint32_t));
            memmove(edges + end + 1, edges + end, (M - end)*sizeof(uint32_t));
            edges[end] = dest;
            memmove(values + end + 1, values + end, (M - end)*sizeof(uint32_t));
            values[end] = value;
            M += 1;

        for(uint32_t i = src + 1; i < N; i++) {
          nodes[i]++;
        }
	      return;
    }
    // could be optimized by making it a binary search
    // but it doesn't change the complexity
    else {
      for (uint32_t i = start; i < end; i++) {
        if (edges[i] > dest) {
          edges = (uint32_t *) realloc(edges, (M+1)*sizeof(uint32_t));
          values = (uint32_t *) realloc(values, (M+1)*sizeof(uint32_t));
          memmove(edges + i + 1, edges + i, (M - i)*sizeof(uint32_t));
          edges[i] = dest;
          memmove(values + i + 1, values + i, (M - i)*sizeof(uint32_t));
          values[i] = value;
          M += 1;
          break;
        } else if (edges[i] == dest) {
          values[i] = value;
          return;
        }
      }
      for(uint32_t i = src + 1; i < N; i++) {
        nodes[i]++;
      }
  }
}

void CSR::add_edge_update(uint32_t src, uint32_t dest, uint32_t value) {
    uint32_t index = find_value_index(src, dest);
    if (index < UINT32_MAX) {
      values[index] = value;
      return;
    }
    uint32_t end = 0;
    if (src ==  N - 1) { end = M; }
    else { end = nodes[src + 1]; }
    uint32_t start = nodes[src];
    //printf("start = %d, end = %d\n", start, end);
    if(M == 0) {
      edges = (uint32_t *) realloc(edges, (M+1)*sizeof(uint32_t));
      values = (uint32_t *) realloc(values, (M+1)*sizeof(uint32_t));
      edges[M] = dest;
      values[M] = value;
      M += 1;
      for(uint32_t i = src + 1; i < N; i++) {
        nodes[i]++;
      }
      return;
    }
    if (dest < edges[start] || start == end) { 
      edges = (uint32_t *) realloc(edges, (M+1)*sizeof(uint32_t));
      values = (uint32_t *) realloc(values, (M+1)*sizeof(uint32_t));
      memmove(edges + start + 1, edges + start, (M - start)*sizeof(uint32_t));
      edges[start] = dest;
      memmove(values + start + 1, values + start, (M - start)*sizeof(uint32_t));
      values[start] = value;
      M += 1;
        for(uint32_t i = src + 1; i < N; i++) {
          nodes[i]++;
        }
        return;
    } else if (end == 0 || dest > edges[end - 1]) {
            edges = (uint32_t *) realloc(edges, (M+1)*sizeof(uint32_t));
            values = (uint32_t *) realloc(values, (M+1)*sizeof(uint32_t));
            memmove(edges + end + 1, edges + end, (M - end)*sizeof(uint32_t));
            edges[end] = dest;
            memmove(values + end + 1, values + end, (M - end)*sizeof(uint32_t));
            values[end] = value;
            M += 1;

        for(uint32_t i = src + 1; i < N; i++) {
          nodes[i]++;
        }
        return;
    }
    // could be optimized by making it a binary search
    // but it doesn't change the complexity
    else {
      for (uint32_t i = start; i < end; i++) {
        if (edges[i] > dest) {
          edges = (uint32_t *) realloc(edges, (M+1)*sizeof(uint32_t));
          values = (uint32_t *) realloc(values, (M+1)*sizeof(uint32_t));
          memmove(edges + i + 1, edges + i, (M - i)*sizeof(uint32_t));
          edges[i] = dest;
          memmove(values + i + 1, values + i, (M - i)*sizeof(uint32_t));
          values[i] = value;
          M += 1;
          break;
        } else if (edges[i] == dest) {
          values[i] = value;
          return;
        }
      }
      for(uint32_t i = src + 1; i < N; i++) {
        nodes[i]++;
      }
  }
}

void CSR::add_edge_batch_update(uint32_t *srcs, uint32_t *dests, uint32_t *vals, uint32_t edge_count) {
  for(int i = 0; i < edge_count; i++) {
    add_edge_update(srcs[i], dests[i], vals[i]);
  }
}


void CSR::print_graph() {
  for(uint32_t i = 0; i < N; i++) {
    vector<uint32_t> edgelist(N, 0);
    uint32_t end = 0;
    if (i == N - 1) { end = M; }
    else { end = nodes[i + 1]; }
    uint32_t start = nodes[i];
           
    for(uint32_t j = start; j < end; j++) {
      edgelist[edges[j]] = values[j];     
    }
    for(uint32_t j = 0; j < N; j++) {
          printf("%03d ", edgelist[j]);
    }
    printf("\n");
  }
}

// constructor
CSR::CSR(uint32_t init_n) {
  edges = nullptr;
  values = nullptr;
  nodes = nullptr;
  locks = nullptr;
  N = 0;
  M = 0;
  for(uint32_t i = 0; i < init_n; i++) {
    add_node();
  }
}

CSR::~CSR() {
  free(edges);
  free(values);
  free(nodes);
  free(locks);
}

/*
int main() {
    // graph_t g;
    // setup(&g);
    CSR g = CSR(5);
    while (1) {
        uint32_t src, dest, value;
        scanf("%d %d %d", &src, &dest, &value);
        // printf("src:%d, dst: %d, val:%d\n", src, dest, value);
        g.add_edge(src, dest, value);
        g.print_graph();

        
        std::vector<uint32_t> temp;
        temp.push_back(10);
        temp.push_back(1);
        temp.push_back(0);
        temp.push_back(0);
        temp.push_back(10);
        
        std::vector<uint32_t> res = g.sparse_matrix_vector_multiplication(temp);
        for(int i = 0; i < 5; i++) { printf("%d ", res[i]); }
        printf("\n");
    }
}
*/
