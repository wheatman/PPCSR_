/*
 * adjacency matrix
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <vector>
#include <cstring>
#include <assert.h>


#include "Graph.hpp"

class AdjacencyMatrix : public Graph {
    public:
    // data members
    uint64_t n; // num vertices
    uint32_t* values;
    vector<uint32_t> num_neighbors;
    
    // function headings
    AdjacencyMatrix(uint32_t init_n);
    ~AdjacencyMatrix();

    uint64_t get_size();    
    uint64_t get_n();
    uint32_t find_value(uint32_t src, uint32_t dest);
    vector<uint32_t> sparse_matrix_vector_multiplication(std::vector<uint32_t> const &v);
    void print_graph();
    void add_node();
    void add_edge(uint32_t src, uint32_t dest, uint32_t value);
    void add_edge_update(uint32_t src, uint32_t dest, uint32_t value);
    // vector<uint32_t> find_neighbors(uint32_t v);
    void convert(Graph* g); // convert a graph representation to this graph reprsentation
    vector<double> pagerank(std::vector<double> const &node_values);
    vector<double> pagerank_pull(std::vector<double> &node_values) {
      return pagerank(node_values);
    }
    vector<int32_t> bfs(int32_t start_node);
    pvector<int32_t> parallel_bfs(int32_t start_node, int32_t total_edges, int alpha, int beta) {
      printf("not implemented\n");
      exit(-1);
    }

    void add_edge_batch_update(uint32_t *srcs, uint32_t *dests, uint32_t *values, uint32_t edge_count);
    vector<tuple<uint32_t, uint32_t, uint32_t> > get_edges() {
      uint64_t n = get_n();
      vector<tuple<uint32_t, uint32_t, uint32_t>> output;
      for(int i = 0; i < n; i++) { 
        for(int j = 0; j < n; j++) {
          if (find_value(i,j) > 0) {
            output.push_back(make_tuple(i, j, find_value(i,j)));
          }
        }
      } 
      return output;
    }
    uint64_t triangle_count() {
      uint64_t total = 0;
      for (uint64_t i = 0; i < n; i++) {
        for (uint64_t j = 0; j < n; j++) {
          if (j > i) break;
          for (uint64_t k = 0; k < n; k++) {
            if (k > j) break;
            if (find_value(i,j) > 0 && find_value(j,k) > 0 && find_value(k,i) > 0) {
              total++;
            }
          }
        }   
      }
      return total;
    }

    
};

vector<int32_t> AdjacencyMatrix::bfs(int32_t start_node) {
  uint64_t n = get_n();
  vector<int32_t> out (n, -1);
  queue<uint32_t> next; // queue
  next.push(start_node);
  out[start_node] = start_node;
  
  while(next.size() > 0) {
    uint32_t active = next.front();
    next.pop();
    assert(out[active] != UINT32_MAX);
    for(uint32_t i = 0; i < n; i++) {
      if(find_value(active, i) != 0 && out[i] == -1) {
        next.push(i);
        out[i] = active;
      }
    }
  }
  return out;
}

vector<double> AdjacencyMatrix::pagerank(std::vector<double> const &node_values) {
  vector<double> output(n, 0);
  for(uint32_t i = 0; i < n; i++) {
    double contrib = 0;
    if (num_neighbors[i] > 0) {
      contrib = (node_values[i] / num_neighbors[i]);
    }

    for(uint32_t j = 0; j < n; j++) {
      if(find_value(i, j) != 0) {
        output[j] += contrib;
      }
    }
  }
  return output;
}

uint64_t AdjacencyMatrix::get_n() {
  return n;
}

void AdjacencyMatrix::convert(Graph* g) {
  n = g->get_n();
  values = (uint32_t*) realloc(values,  n * n * sizeof(uint32_t));
  for(uint32_t i = 0; i < n; i ++) {
    for(uint32_t j = 0; j < n; j++) {
      // find_value returns 0 if not found.
      add_edge(i, j, g->find_value(i, j));
    }
  }
}

uint64_t AdjacencyMatrix::get_size() {
    return n * n * sizeof(*values) + num_neighbors.capacity()*sizeof(uint32_t);
}

// find value of (src, dest)
uint32_t AdjacencyMatrix::find_value(uint32_t src, uint32_t dest) {
    return values[src*n + dest];
}

void AdjacencyMatrix::add_node() {
    uint64_t new_size = (n+1) * (n+1);
    uint32_t* new_graph = (uint32_t*)malloc(sizeof(uint32_t) * new_size);
    /// printf("n: %d\n", n);
    // clear the new graph
    memset(new_graph, 0, sizeof(uint32_t) * new_size);

    // copy old graph into new graph
    for(uint32_t i = 0; i < n; i++) {
        for(uint32_t j = 0; j < n; j++) {
            new_graph[i * ( n + 1 ) + j] = values[i * n + j];
        }
    }

    // update graph with new rep
    if(values != NULL) {
        free(values);
    }

    values = new_graph;
    n++;
    num_neighbors.push_back(0);
}

// src, dest < N
void AdjacencyMatrix::add_edge(uint32_t src, uint32_t dest, uint32_t value) {
    values[src*n + dest] = value;
    num_neighbors[src]++;
}
void AdjacencyMatrix::add_edge_update(uint32_t src, uint32_t dest, uint32_t value) {
  if (values[src*n + dest] == 0) {
    num_neighbors[src]++;
  }
  values[src*n + dest] = value;
}
void AdjacencyMatrix::add_edge_batch_update(uint32_t *srcs, uint32_t *dests, uint32_t *values, uint32_t edge_count) {
  cilk_for(int i = 0; i < edge_count; i++) {
    add_edge_update(srcs[i], dests[i], values[i]);
  }
}

void AdjacencyMatrix::print_graph() {
    for(uint32_t i = 0; i < n; i++) {
        for(uint32_t j = 0; j < n; j++) {
            printf("%03d ", values[i*n + j]);
        }
        printf("\n");
    }
}

// constructor
AdjacencyMatrix::AdjacencyMatrix(uint32_t init_n) {
  n = (uint64_t) init_n;
  uint64_t new_size = n * n;
  //printf("new_size: %lu\n", new_size);
  values = (uint32_t*)malloc(sizeof(uint32_t) * new_size);
  assert(values != NULL);
  memset(values, 0, sizeof(uint32_t) * new_size);
  for (uint32_t i = 0; i < init_n; i++) {
    num_neighbors.push_back(0);
  }
}

AdjacencyMatrix::~AdjacencyMatrix() {
  free(values);
}

/*
// find all neighbors of a vertex v
std::vector<uint32_t> AdjacencyMatrix::find_neighbors(uint32_t v) {
  std::vector<uint32_t> result;
  for(int i = 0; i < n; i++) {
    if(values[v*n + i] != 0) {
      result.push_back(i);
    }
  }
  return result;
}
*/

std::vector<uint32_t> AdjacencyMatrix::sparse_matrix_vector_multiplication(std::vector<uint32_t> const &v) {
    std::vector<uint32_t> result;
    for(uint32_t i = 0; i < n; i++) {
        uint32_t temp = 0;
        for(uint32_t j = 0; j < n; j++) {
            temp += values[i*n + j] * v[j];
        }
        result.push_back(temp);
    }
    return result;
}

/*
int main() {
    // graph_t g;
    // setup(&g);
    AdjacencyMatrix g = AdjacencyMatrix(5);
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
