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
#include "SpinLock.cpp"
#include <atomic>
#include <cilk/cilk.h>
#include "helpers.h"

typedef struct _adj_edge {
    uint32_t dest;
    uint32_t val;
    struct _adj_edge* next;
} adj_edge_t;

typedef struct _adj_node {
    adj_edge_t* head;
    uint32_t num_neighbors;
    SpinLock lock;
} adj_node_t;

class AdjacencyList : public Graph {
    public:
    // data members
    std::vector<adj_node_t> nodes;
    
    // function headings
    AdjacencyList(uint32_t init_n);
    ~AdjacencyList(); // destructor

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
      return nodes[node].num_neighbors;
    }    
 
    class iterator {
    public:
      adj_edge_t *place;
      iterator(AdjacencyList *G, uint32_t node, bool start) {
        if (!start) {
          place = nullptr;
          return;
        }
        place = G->nodes[node].head;
        return;
      }
      bool operator==(const iterator& other) const {
        return place == other.place;
      }
      bool operator!=(const iterator& other) const {
        return place != other.place;
      }
      iterator& operator++() {
        place = place->next;
        return *this;
      }
      edge_t operator*() const {
        return {place->val, place->dest};
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
void AdjacencyList::convert(Graph* g) {
  printf("convert not impemented for AdjencyList exiting now");
  exit(1);
}

uint64_t AdjacencyList::get_n() {
  return nodes.size();
}

vector<tuple<uint32_t, uint32_t, uint32_t> > AdjacencyList::get_edges() {
  uint64_t n = get_n();
  vector<tuple<uint32_t, uint32_t, uint32_t>> output;

  for(int i = 0; i < n; i++) {
    adj_edge_t* temp = nodes[i].head;
    while (temp != NULL) {
      output.push_back(make_tuple(i, temp->dest,temp->val));
      temp = temp->next;
    }
  }
  return output;
}

void AdjacencyList::add_file(string filename) {
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
void AdjacencyList::add_file2(string filename) {
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
void AdjacencyList::add_file3(string filename) {
  ifstream myfile(filename.c_str()); 
  string line;
  if (myfile.is_open()) {
    //int line_num = 0;
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

uint64_t AdjacencyList::get_size() {
  uint64_t size = nodes.capacity() * sizeof(adj_node_t);
  
  for(int i = 0; i < nodes.size(); i++) {
      adj_edge_t* temp = nodes[i].head;
      while (temp != NULL) {
        size += sizeof(adj_edge_t);
        temp = temp->next;
      }
  }
  return size;
}

uint32_t AdjacencyList::find_value(uint32_t src, uint32_t dest) {
    nodes[src].lock.lock();
    adj_edge_t* e = nodes[src].head;
    while(e != NULL) {
        if(e->dest == dest) { 
          uint32_t val = e->val;
          nodes[src].lock.unlock();
          return val; 
        }
        e = e->next;
    }
    nodes[src].lock.unlock();
    return 0;
}

// add a disconnected new node
void AdjacencyList::add_node() {
    adj_node_t node;
    node.head = NULL;
    node.num_neighbors = 0;
    nodes.push_back(node);
}

// src, dest < N
void AdjacencyList::add_edge(uint32_t src, uint32_t dest, uint32_t value) {
    nodes[src].lock.lock();
    adj_edge_t* e = (adj_edge_t *)malloc(sizeof(adj_edge_t));
    // printf("add edge: %d, %d, %d\n", src, dest, value);
    e->val = value;
    e->dest = dest;
    e->next = nodes[src].head;
    nodes[src].head = e;
    nodes[src].num_neighbors++;
    assert(nodes[src].head != NULL);
    nodes[src].lock.unlock();
}

void AdjacencyList::add_edge_update(uint32_t src, uint32_t dest, uint32_t value) {
    nodes[src].lock.lock();
    adj_edge_t* e = nodes[src].head;
    while(e != NULL) {
        if(e->dest == dest) {
          e->val = value;
          nodes[src].lock.unlock();
          return; 
        }
        e = e->next;
    }
    e = (adj_edge_t *)malloc(sizeof(adj_edge_t));
    // printf("add edge: %d, %d, %d\n", src, dest, value);
    e->val = value;
    e->dest = dest;
    e->next = nodes[src].head;
    nodes[src].head = e;
    nodes[src].num_neighbors++;
    assert(nodes[src].head != NULL);
    nodes[src].lock.unlock();
}


void AdjacencyList::print_graph() {
    for(int i = 0; i < nodes.size(); i++) { // iterate over the nodes
        vector<uint32_t> edgelist (nodes.size(), 0);
        adj_edge_t * temp = nodes[i].head;
        while(temp != NULL) {
            edgelist[temp->dest] = temp->val;
            temp = temp->next;
        }
        for(int j = 0; j < nodes.size(); j++) {
            printf("%03d ", edgelist[j]);
        }
        printf("\n");
    }
}

AdjacencyList::AdjacencyList(uint32_t init_n) {
  for (int i = 0; i < init_n; i++) {
      add_node();
      // print_graph();
  }
}

AdjacencyList::~AdjacencyList() {
  for (int i = 0; i < nodes.size(); i++) {
      adj_edge_t* e = nodes[i].head;
      while (e) {
        adj_edge_t* e_old = e;
        e = e->next;
        free(e_old);
      }
  }
}

void AdjacencyList::add_edge_batch_update(uint32_t *srcs, uint32_t *dests, uint32_t *values, uint32_t edge_count) {
  cilk_for(int i = 0; i < edge_count; i++) {
    add_edge_update(srcs[i], dests[i], values[i]);
  }
}
