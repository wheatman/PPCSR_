/*
 * adjacency vector
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <vector>
#include <assert.h>
#include "Graph.hpp"
#include <string>
#include <sstream>
#include <iterator>
#include <iostream>
#include <fstream>
#include <tuple>
#include "helpers.h"


typedef struct _adj_edge_vec {
    uint32_t dest;
    uint32_t val;
    // struct _adj_edge_vec* next;
} adj_edge_vec_t;

typedef struct _adj_node_vec {
    std::vector<adj_edge_vec_t> edges;
    SpinLock lock;
    // adj_edge_vec_t* head;
    // uint32_t num_neighbors;
} adj_node_vec_t;

class AdjacencyVector : public Graph {
    public:
    // data members
    std::vector<adj_node_vec_t> nodes;
    
    // function headings
    AdjacencyVector(uint32_t init_n);
    // ~AdjacencyVector(); // destructor

    uint64_t get_size();
    uint32_t find_value(uint32_t src, uint32_t dest);
    void print_graph();
    void add_node();
    void add_edge(uint32_t src, uint32_t dest, uint32_t value);
    void add_edge_update(uint32_t src, uint32_t dest, uint32_t value);
    void add_file(string filename);
    void add_file2(string filename);
    void add_file3(string filename);
    uint64_t get_n();
    void convert(Graph* g);
    vector<tuple<uint32_t, uint32_t, uint32_t> > get_edges();
    void add_edge_batch_update(uint32_t *srcs, uint32_t *dests, uint32_t *values, uint32_t edge_count);
    uint32_t num_neighbors(uint32_t node) {
      return nodes[node].edges.size();
    }
    
    class iterator {
    public:
      uint32_t index;
      std::vector<adj_edge_vec_t> edges;
      iterator(AdjacencyVector *G, uint32_t node, bool start) {
        if (!start) {
          index = G->num_neighbors(node);
          return;
        }
        index = 0;
        edges = G->nodes[node].edges;
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
        return {edges[index].val, edges[index].dest};
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
    TRIANGLE_COUNT  
    PARALLEL_BFS  
};

void AdjacencyVector::convert(Graph* g) {
  printf("convert not impemented for AdjacencyVector exiting now");
  exit(1);
}


uint64_t AdjacencyVector::get_n() {
  return nodes.size();
}

vector<tuple<uint32_t, uint32_t, uint32_t> > AdjacencyVector::get_edges() {
  uint64_t n = get_n();
  vector<tuple<uint32_t, uint32_t, uint32_t>> output;

  for(int i = 0; i < n; i++) {
    for(auto &temp : nodes[i].edges) {
      output.push_back(make_tuple(i, temp.dest,temp.val));
    }
  }
  return output;
}

void AdjacencyVector::add_file(string filename) {
  ifstream myfile(filename.c_str()); 
  bool first = false;
  int n, m;
  string line;
  if (myfile.is_open()) {
    while ( getline (myfile,line) ) {
      if(line[0] == '%') {  continue;  }
      vector<string> elems = split(line, ' ');
      if (!first) {
        n = atoi(elems[0].c_str() );
        for(int i = 0; i < n; i++) {
	  add_node();
        }
        m = atoi(elems[2].c_str());
        first = true;
      } else {
        int src = atoi(elems[0].c_str());
        int dest = atoi(elems[1].c_str());
        
        add_edge( src, dest, 1 );
      }
    }
  myfile.close();
  // return 0;
  }
}

//for twitter
// I don't think this is correct
void AdjacencyVector::add_file2(string filename) {
  ifstream myfile(filename.c_str()); 
  string line;
  if (myfile.is_open()) {
    int line_num = 0;
    while ( getline (myfile,line) ) {
      if(line[0] == '%') {  continue;  }
      vector<string> elems = split(line, ' ');
        int src = atoi(elems[0].c_str())-1;

        while (src >= get_n()+1) {
          add_node();
        }
        int dest = atoi(elems[1].c_str())-1;

        while (dest >= get_n()+1) {
          add_node();
        }

        add_edge( src, dest, 1 );
        if (line_num++ > 400000000) {
          break;
        }
    }
  myfile.close();
  // return 0;
  }
}

// for soc-
// starting at 1
void AdjacencyVector::add_file3(string filename) {
  ifstream myfile(filename.c_str()); 
  string line;
  if (myfile.is_open()) {
    while ( getline (myfile,line) ) {
      vector<string> elems = split(line, '\t');
        int src = atoi(elems[0].c_str())-1;

        while (src >= get_n()) {
          add_node();
        }
        int dest = atoi(elems[1].c_str())-1;

        while (dest >= get_n()) {
          add_node();
        }
        add_edge( src, dest, 1 );
    }
  myfile.close();
  // return 0;
  } else {
    printf("file was not opened\n");
  }
}

uint64_t AdjacencyVector::get_size() {
  uint64_t size = nodes.capacity() * sizeof(adj_node_vec_t);
  
  for(int i = 0; i < nodes.size(); i++) {
    size += nodes[i].edges.capacity() * sizeof(adj_edge_vec_t);
  }
  return size;
}

uint32_t AdjacencyVector::find_value(uint32_t src, uint32_t dest) {
  for(auto &e : nodes[src].edges) {
    if (e.dest == dest) { return e.val; }
  }
  return 0;
}

// add a disconnected new node
void AdjacencyVector::add_node() {
    adj_node_vec_t node;
    nodes.push_back(node);
}

// src, dest < N
void AdjacencyVector::add_edge(uint32_t src, uint32_t dest, uint32_t value) {
    nodes[src].lock.lock();
    adj_edge_vec_t e;
    // printf("add edge: %d, %d, %d\n", src, dest, value);
    e.val = value;
    e.dest = dest;
    nodes[src].edges.push_back(e);
    nodes[src].lock.unlock();
}

void AdjacencyVector::add_edge_update(uint32_t src, uint32_t dest, uint32_t value) {
  nodes[src].lock.lock();
  for(auto &e : nodes[src].edges) {
    if (e.dest == dest) { 
      e.val = value; 
      nodes[src].lock.unlock(); 
      return; 
    }
  }
  adj_edge_vec_t e;
  // printf("add edge: %d, %d, %d\n", src, dest, value);
  e.val = value;
  e.dest = dest;
  nodes[src].edges.push_back(e);
  nodes[src].lock.unlock();
}

void AdjacencyVector::add_edge_batch_update(uint32_t *srcs, uint32_t *dests, uint32_t *values, uint32_t edge_count) {
  cilk_for(int i = 0; i < edge_count; i++) {
    add_edge_update(srcs[i], dests[i], values[i]);
  }
}

void AdjacencyVector::print_graph() {
    for(int i = 0; i < nodes.size(); i++) { // iterate over the nodes
        vector<uint32_t> edgelist (nodes.size(), 0);
        for(auto &e : nodes[i].edges) {
          edgelist[e.dest] = e.val;
        }

        for(int j = 0; j < nodes.size(); j++) {
            printf("%03d ", edgelist[j]);
        }
        printf("\n");
    }
}

AdjacencyVector::AdjacencyVector(uint32_t init_n) {
  for (int i = 0; i < init_n; i++) {
      add_node();
      // print_graph();
  }
}
