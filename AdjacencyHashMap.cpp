/*
 * adjacency list
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
#include <unordered_map>
#include "helpers.h"
/*
typedef struct _adj_edge {
    uint32_t dest;
    uint32_t val;
} adj_edge_t;
*/


typedef struct _ahm_node {
  std::unordered_map<uint32_t, uint32_t> edges;
  SpinLock lock;
} ahm_node_t;

class AdjacencyHashMap : public Graph {
    public:
    // data members
    std::vector<ahm_node_t> nodes;
    // std::vector<std::unordered_map<uint32_t, uint32_t>> nodes;
    // std::vector<SpinLock> locks;
    
    // function headings
    AdjacencyHashMap(uint32_t init_n);
    // ~AdjacencyHashMap(); // destructor

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
      std::unordered_map<uint32_t, uint32_t>::iterator it;
      iterator(AdjacencyHashMap *G, uint32_t node, bool start) {
        if (!start) {
          it = G->nodes[node].edges.end();
          return;
        }
        it = G->nodes[node].edges.begin();
        return;
      }
    bool operator==(const iterator& other) const {
      return it == other.it;
    }
    bool operator!=(const iterator& other) const {
      return it != other.it;
    }
    iterator& operator++() {
      ++it;
      return *this;
    }
    edge_t operator*() const {
      return {it->second, it->first};
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

void AdjacencyHashMap::convert(Graph* g) {
  printf("convert not impemented for AdjacencyHashMap exiting now");
  exit(1);
}


uint64_t AdjacencyHashMap::get_n() {
  return nodes.size();
}

vector<tuple<uint32_t, uint32_t, uint32_t> > AdjacencyHashMap::get_edges() {
  uint64_t n = get_n();
  vector<tuple<uint32_t, uint32_t, uint32_t>> output;

  for(int i = 0; i < n; i++) {
    for(auto &temp : nodes[i].edges) {
      output.push_back(make_tuple(i, temp.first,temp.second));
    }
  }
  return output;
}

void AdjacencyHashMap::add_file(string filename) {
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
void AdjacencyHashMap::add_file2(string filename) {
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
void AdjacencyHashMap::add_file3(string filename) {
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

uint64_t AdjacencyHashMap::get_size() {
  printf("can't get capacity of a hash map in c++\n");
  printf("here's a lower bound\n");

  // node list
  uint64_t size = nodes.capacity() * sizeof(std::unordered_map<uint32_t, uint32_t>);

  uint64_t map_size = 0;
  // edges
  for(int i = 0; i < nodes.size(); i++) {
    for(uint32_t j = 0; j < nodes[i].edges.bucket_count(); j++) {
      uint32_t bucket_size = nodes[i].edges.bucket_size(j);
      if (bucket_size == 0) {
        map_size++;
      } else {
        map_size += bucket_size;
      }
    }
  }
  map_size *= sizeof(uint32_t);
  return size + map_size;
}

uint32_t AdjacencyHashMap::find_value(uint32_t src, uint32_t dest) {
  auto search = nodes[src].edges.find(dest);
  if (search != nodes[src].edges.end()) {
    return search->second;
  }
  return 0;
}

// add a disconnected new node
void AdjacencyHashMap::add_node() {
  ahm_node_t node;
  // node.edges = std::unordered_map<uint32_t, uint32_t>();
  // node.lock = SpinLock();
  nodes.push_back(node);
}

void AdjacencyHashMap::add_edge_batch_update(uint32_t *srcs, uint32_t *dests, uint32_t *values, uint32_t edge_count)
{
  cilk_for(int i = 0; i < edge_count; i++) {
      add_edge_update(srcs[i], dests[i], values[i]);
  }
}

// src, dest < N
void AdjacencyHashMap::add_edge(uint32_t src, uint32_t dest, uint32_t value) {
  nodes[src].lock.lock();
  nodes[src].edges[dest] = value;
  nodes[src].lock.unlock();
}

void AdjacencyHashMap::add_edge_update(uint32_t src, uint32_t dest, uint32_t value) {
  nodes[src].lock.lock();
  nodes[src].edges[dest] = value;
  nodes[src].lock.unlock();
}

void AdjacencyHashMap::print_graph() {
    for(int i = 0; i < nodes.size(); i++) { // iterate over the nodes
        vector<uint32_t> edgelist (nodes.size(), 0);
        for(auto &e : nodes[i].edges) {
          edgelist[e.first] = e.second;
        }

        for(int j = 0; j < nodes.size(); j++) {
            printf("%03d ", edgelist[j]);
        }
        printf("\n");
    }
}

AdjacencyHashMap::AdjacencyHashMap(uint32_t init_n) {
  for (int i = 0; i < init_n; i++) {
      add_node();
      // print_graph();
  }
}
