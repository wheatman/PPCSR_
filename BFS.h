#include "Map.cpp"
// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of 
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


struct BFS_F {
  int32_t* Parents;
  BFS_F(int32_t* _Parents) : Parents(_Parents) {}
  inline bool update (uint32_t s, uint32_t d) { //Update
    if(Parents[d] == -1) { Parents[d] = s; return 1; }
    else return 0;
  }
  inline bool updateAtomic (uint32_t s, uint32_t d){ //atomic version of Update
    return __sync_bool_compare_and_swap(&Parents[d],-1,s);
  }
  //cond function checks if vertex has been visited yet
  inline bool cond (uint32_t d) { return (Parents[d] == -1); } 
};

int32_t* BFS_with_edge_map(OFM &G, uint32_t src) {
  long start = src;
  long n = G.get_n();
  //creates Parents array, initialized to all -1, except for start
  int32_t* Parents = (int32_t *) malloc(n * sizeof(uint32_t));
  parallel_for(long i=0;i<n;i++) Parents[i] = -1;
  Parents[start] = start;
  VertexSubset frontier = VertexSubset(start, n); //creates initial frontier
  while(frontier.get_n() > 0){ //loop until frontier is empty
    //printf("frontier size  %lu\n", frontier.get_n());
    edgeMap(G, frontier, BFS_F(Parents));    
    frontier.move_next_to_current();
  } 
  return Parents;
}
