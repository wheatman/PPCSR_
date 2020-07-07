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

#include "SpinLock.cpp"
#include <atomic>
#include <cilk/cilk.h>
#include "helpers.h"


#define BLOCK_SIZE 8 // cache line size / words per edge
	    // 64 / 8
	    // potentially could do better with a different block size

typedef struct _blocked_adj_edge {
  uint32_t dest;
  uint32_t val;
} blocked_adj_edge_t;

typedef struct _block_t {
    uint8_t count;
    blocked_adj_edge_t edges[BLOCK_SIZE]; 
    struct _block_t* next;
} block_t;

typedef struct _blocked_adj_node {
  block_t* head;
  uint32_t num_neighbors;
  SpinLock lock;
} blocked_adj_node_t;

class BlockedAdjacencyList : public Graph {
    public:
    // data members
    std::vector<blocked_adj_node_t> nodes;
    
    // function headings
    BlockedAdjacencyList(uint32_t init_n);
    ~BlockedAdjacencyList(); // destructor

    uint64_t get_size();
    uint32_t find_value(uint32_t src, uint32_t dest);
    void print_graph();
    void print_lists();
    void add_node();
    void add_edge(uint32_t src, uint32_t dest, uint32_t value);
    void add_edge_update(uint32_t src, uint32_t dest, uint32_t value);
    uint64_t get_n();
    void convert(Graph* g);
    vector<tuple<uint32_t, uint32_t, uint32_t> > get_edges();
    void add_file3(string filename);
    void add_edge_batch_update(uint32_t *srcs, uint32_t *dests, uint32_t *values, uint32_t edge_count);
    uint32_t num_neighbors(uint32_t node) {
      return nodes[node].num_neighbors;
    }

    class iterator {
    public:
      block_t *block;
      uint32_t index;
      iterator(BlockedAdjacencyList *G, uint32_t node, bool start) {
      if (!start) {
        block = nullptr;
        return;
      }
      block = G->nodes[node].head;
      index = 0;
      return;
    }
    bool operator==(const iterator& other) const {
      return (block == other.block) && (index == other.index);
    }
    bool operator!=(const iterator& other) const {
      // fast but not correct in general, but we only compare not equal to find the end, and end is nullptr
      return block != other.block;
    }
    iterator& operator++() {
      if (index + 1 < block->count) {
        index+=1;
        return *this;
      }
      block = block->next;
      index = 0;
      return *this;
    }
    edge_t operator*() const {
      return {block->edges[index].val, block->edges[index].dest};
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

void BlockedAdjacencyList::convert(Graph* g) {
  // destruct current data structure
  for (int i = 0; i < nodes.size(); i++) {
      block_t* e = nodes[i].head;
      while (e) {
        block_t* e_old = e;
        e = e->next;
        free(e_old);
      }
  }

  int n = g->get_n();
  nodes.clear();

  // populate nodes
  for(int i = 0; i < n; i++) { add_node(); }

  vector<tuple<uint32_t, uint32_t, uint32_t> > edges = g->get_edges();
  std::random_shuffle ( edges.begin(), edges.end() );

  // populate edges
  for(int i = 0; i < edges.size(); i++) {
    add_edge(get<0>(edges[i]), get<1>(edges[i]), get<2>(edges[i]));
  }
}

// for soc-
// starting at 1
void BlockedAdjacencyList::add_file3(string filename) {
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

vector<tuple<uint32_t, uint32_t, uint32_t> > BlockedAdjacencyList::get_edges() {
  uint64_t n = get_n();
  vector<tuple<uint32_t, uint32_t, uint32_t>> output;

  for(int i = 0; i < n; i++) {
    block_t* temp = nodes[i].head;
    while (temp != NULL) {
      for(int j = 0; j < temp->count; j++) {
        output.push_back(make_tuple(i, temp->edges[j].dest,temp->edges[j].val));
      }
      
      temp = temp->next;
    }
  }
  return output;
}


uint64_t BlockedAdjacencyList::get_n() {
  return nodes.size();
}

uint64_t BlockedAdjacencyList::get_size() {
  // size of nodes
  uint64_t size = nodes.capacity() * sizeof(blocked_adj_node_t);
 
  for(int i = 0; i < nodes.size(); i++) {
      block_t* temp = nodes[i].head;
      while (temp != NULL) {
	  size += sizeof(block_t);
	  temp = temp->next;
      }
  }
  return size;
}

uint32_t BlockedAdjacencyList::find_value(uint32_t src, uint32_t dest) {
    nodes[src].lock.lock();
    block_t* e = nodes[src].head;
    //printf("SRC: %d, DEST: %d\n", src, dest);
    while(e != NULL) {
      for (int i = 0; i < e->count; i++) {
        //printf("dest: %d, value: %d\n", e->edges[i].dest, e->edges[i].val);
        if(e->edges[i].dest == dest) {
          uint32_t val = e->edges[i].val;
          nodes[src].lock.unlock();
          return val;
	}
      }
      e = e->next;
    }
    nodes[src].lock.unlock();
    return 0;
}

// add a disconnected new node
void BlockedAdjacencyList::add_node() {
  blocked_adj_node_t node;
  node.head = NULL;
  node.num_neighbors = 0;
  nodes.push_back(node);
  node.lock.x = 0;
}

// src, dest < N
void BlockedAdjacencyList::add_edge(uint32_t src, uint32_t dest, uint32_t value) {
  if (value != 0) {
    nodes[src].lock.lock();
    block_t* blk = nodes[src].head;
    nodes[src].num_neighbors++; // add to neighbors

    if(blk != NULL && blk->count < BLOCK_SIZE) { // insert
      blk->edges[blk->count] = { dest, value };
      //printf("blk dest %d, val %d\n", blk->edges[blk->count].dest, blk->edges[blk->count].val);
      blk->count++;
      //assert(find_value(src, dest) == value);
    } else { // block is full
       block_t* new_block = (block_t *)malloc(sizeof(block_t));
       new_block->edges[0] = {dest, value};
       new_block->next = blk;
       new_block->count = 1;
       nodes[src].head = new_block;
       //assert(nodes[src].head != NULL);
       //assert(find_value(src, dest) == value);

   }
  nodes[src].lock.unlock();
  }
}

void BlockedAdjacencyList::add_edge_update(uint32_t src, uint32_t dest, uint32_t value) {
  if (value != 0) {
    nodes[src].lock.lock();
    block_t* blk = nodes[src].head;
    //printf("SRC: %d, DEST: %d\n", src, dest);
    while(blk != NULL) {
      for (int i = 0; i < blk->count; i++) {
        //printf("dest: %d, value: %d\n", e->edges[i].dest, e->edges[i].val);
        if(blk->edges[i].dest == dest) {
          blk->edges[i].val = value;
          nodes[src].lock.unlock();
          return;
        }
      }
      blk = blk->next;
    }
    blk = nodes[src].head;
    nodes[src].num_neighbors++; // add to neighbors

    if(blk != NULL && blk->count < BLOCK_SIZE) { // insert
      blk->edges[blk->count] = { dest, value };
      //printf("blk dest %d, val %d\n", blk->edges[blk->count].dest, blk->edges[blk->count].val);
      blk->count++;
      //assert(find_value(src, dest) == value);
    } else { // block is full
       block_t* new_block = (block_t *)malloc(sizeof(block_t));
       new_block->edges[0] = {dest, value};
       new_block->next = blk;
       new_block->count = 1;
       nodes[src].head = new_block;
       //assert(nodes[src].head != NULL);
       //assert(find_value(src, dest) == value);

   }
  nodes[src].lock.unlock();
  }
}

void BlockedAdjacencyList::print_graph() {
    for(int i = 0; i < nodes.size(); i++) { // iterate over the nodes
        vector<uint32_t> edgelist (nodes.size(), 0);
        block_t * temp = nodes[i].head;
        while(temp != NULL) {
          for (int i = 0; i < temp->count; i++) {
	    edgelist[temp->edges[i].dest] = temp->edges[i].val;
	  }
	  temp = temp->next;
	}

        for(int j = 0; j < nodes.size(); j++) {
            printf("%03d ", edgelist[j]);
        }
        printf("\n");
    }
}

void BlockedAdjacencyList::print_lists() {
    for(int i = 0; i < nodes.size(); i++) { // iterate over the nodes
        block_t * temp = nodes[i].head;
        while(temp != NULL) {
          for (int j = 0; j < temp->count; j++) {
            printf("(%d, %d, %d), ", i,  temp->edges[j].dest, temp->edges[j].val);
          }
          temp = temp->next;
        }
        printf("\n");
    }
}

BlockedAdjacencyList::BlockedAdjacencyList(uint32_t init_n) {
  for (int i = 0; i < init_n; i++) {
      add_node();
      // print_graph();
  }
}

BlockedAdjacencyList::~BlockedAdjacencyList() {
  for (int i = 0; i < nodes.size(); i++) {
      block_t* e = nodes[i].head;
      while (e) {
        block_t* e_old = e;
        e = e->next;
        free(e_old);
      }
  }
}

void BlockedAdjacencyList::add_edge_batch_update(uint32_t *srcs, uint32_t *dests, uint32_t *values, uint32_t edge_count) {
  cilk_for(int i = 0; i < edge_count; i++) {
    add_edge_update(srcs[i], dests[i], values[i]);
  }
}
