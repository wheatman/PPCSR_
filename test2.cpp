// test framework
#define NDEBUG
//#define cilk_for for
//#define usleep(a)
#define batch_size 10000
#include <set>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
// #include <pair>
#include "AdjacencyList.cpp"
#include "BlockedAdjacencyList.cpp"
#include "AdjacencyHashMap.cpp"
#include "AdjacencyMatrix.cpp"
#include "AdjacencyVector.cpp"
#include "CSR.cpp"
#include "Graph.hpp"
#include "OFM.cpp"
#include "BFS.h"
#include "BC.h"
#include "TC.h"
#include "Pagerank.h"
#include "Components.h"
#include "rmat_util.h"
#include <ctime>
#include <iostream>
#include <random>
#include <sys/time.h>
#define MAXVAL 100

#define batch_size_for_test 10

/*
static uint32_t x = 123456789;
static uint32_t y = 362436069;
static uint32_t z = 521288629;
static uint32_t w = 88675123;

uint32_t xor128(void) {
  uint32_t t;
  t = x ^ (x << 11);
  x = y; y = z; z = w;
  return w = w ^ (w >> 19) ^ (t ^ (t >> 8));
}
void srand() {
  x = 123456789;
  y = 362436069;
  z = 521288629;
  w = 88675123;
}
*/
uint32_t rand_in_range(uint32_t max) { return rand() % max; }

bool compare_matrices(Graph *g1, Graph *g2, int num_nodes) {
  for (int i = 0; i < num_nodes; i++) {
    for (int j = 0; j < num_nodes; j++) {
      if (g1->find_value(i, j) != g2->find_value(i, j)) {
        printf("i: %d, j: %d, g1_val: %d, g2_vl: %d\n", i, j,
               g1->find_value(i, j), g2->find_value(i, j));
        if (j + 1 < num_nodes) {
          printf("i: %d, j: %d, g1_val: %d, g2_vl: %d\n", i, j + 1,
                 g1->find_value(i, j + 1), g2->find_value(i, j + 1));
        }
        if (j - 1 >= 0) {
          printf("i: %d, j: %d, g1_val: %d, g2_vl: %d\n", i, j - 1,
                 g1->find_value(i, j - 1), g2->find_value(i, j - 1));
        }
        return false;
      }
    }
  }
  return true;
}

bool verify_structure(Graph *g, uint32_t num_nodes, uint32_t num_edges, int trials) {
      // printf("########## starting trial %d out of %d, nodes = %u, edges = %u "
      //        "##########\n",
      //        j, trials, num_nodes, num_edges);
      AdjacencyMatrix matrix = AdjacencyMatrix(num_nodes);
      // AdjacencyHashMap adjhash = AdjacencyHashMap(num_nodes);
      uint32_t srcs[batch_size_for_test] = {0};
      uint32_t dests[batch_size_for_test] = {0};
      uint32_t vals[batch_size_for_test] = {0};
      for (uint32_t i = 0; i <= num_edges / batch_size_for_test; i++) {
        for (uint32_t b = 0; b < batch_size_for_test; b++) {
          uint32_t src = rand_in_range(num_nodes);
          uint32_t dest = rand_in_range(num_nodes);
          while (matrix.find_value(src, dest) != 0) {
            src = rand_in_range(num_nodes);
            dest = rand_in_range(num_nodes);
          }
          uint32_t val = rand_in_range(MAXVAL) + 1;

          matrix.add_edge(src, dest, val);
          srcs[b] = src;
          dests[b] = dest;
          vals[b] = val;
        }
        for (uint32_t b = 0; b < batch_size_for_test; b++) {
          // printf("about to add edge %d, %d, %d by worker %lu\n",srcs[b],
          // dests[b], vals[b], get_worker_num());
          g->add_edge(srcs[b], dests[b], vals[b]);
        }
      }
      if (!compare_matrices(&matrix, g, num_nodes)) {
        printf("failed in add edges, nodes = %u, edges = %u\n", num_nodes,
               num_edges);
        matrix.print_graph();
        printf("\n");
        g->print_graph();
        return false;
      }
      for (uint32_t i = 0; i <= num_edges / batch_size_for_test; i++) {
        for (uint32_t b = 0; b < batch_size_for_test; b++) {
          uint32_t src = rand_in_range(num_nodes);
          uint32_t dest = rand_in_range(num_nodes);
          while (matrix.find_value(src, dest) != 0) {
            src = rand_in_range(num_nodes);
            dest = rand_in_range(num_nodes);
          }
          uint32_t val = rand_in_range(MAXVAL) + 1;

          matrix.add_edge_update(src, dest, val);
          srcs[b] = src;
          dests[b] = dest;
          vals[b] = val;
        }

        for (uint32_t b = 0; b < batch_size_for_test; b++) {
          g->add_edge_update(srcs[b], dests[b], vals[b]);
        }
      }
      if (!compare_matrices(&matrix, g, num_nodes)) {
        printf("failed in add edges update, nodes = %u, edges = %u\n",
               num_nodes, num_edges);
        matrix.print_graph();
        printf("\n");
        g->print_graph();
        return false;
      }

      vector<uint32_t> test_vector(num_nodes, 0);
      for (uint32_t i = 0; i < num_nodes; i++) {
        test_vector[i] = rand_in_range(MAXVAL);
      }
      vector<uint32_t> r1 =
          matrix.sparse_matrix_vector_multiplication(test_vector);
      vector<uint32_t> r2 =
          g->sparse_matrix_vector_multiplication(test_vector);
      if (r1 != r2) {
        printf("failed spvm, nodes = %u, edges = %u\n", num_nodes, num_edges);
        for (auto const &c : r1) {
          printf("%d ", c);
        }
        printf("\n");
        for (auto const &c : r2) {
          printf("%d ", c);
        }
        printf("\n");
        return false;
      }

      int start_node = rand_in_range(num_nodes);
      vector<int32_t> r3 = matrix.bfs(start_node);
      vector<int32_t> r4 = g->bfs(start_node);
      if (r3 != r4) {
        printf("failed bfs, nodes = %u, edges = %u\n", num_nodes, num_edges);
        matrix.print_graph();
        printf("\n");
        g->print_graph();
        for (auto const &c : r3) {
          printf("%d ", c);
        }
        printf("\n");
        for (auto const &c : r4) {
          printf("%d ", c);
        }
        printf("\n");
        return false;
      }

      vector<double> test_vector2(num_nodes, 0);
      for (uint32_t i = 0; i < num_nodes; i++) {
        test_vector2[i] = ((double)rand_in_range(MAXVAL)) / MAXVAL;
      }
      g->make_symetric(); 
      matrix.make_symetric(); 
      vector<double> r5 = matrix.pagerank(test_vector2);
      vector<double> r6 = g->pagerank_pull(test_vector2);
      bool different = false;
      for (uint32_t i = 0; i < r5.size(); i++) {
        if (fabs(r5[i] - r6[i]) > .0000000000001L) different = true;
      }
      if (different) {
        printf("failed pagerank, nodes = %u, edges = %u\n", num_nodes,
               num_edges);
        matrix.print_graph();
        printf("\n");
        g->print_graph();
        for (auto const &c : test_vector2) {
          printf("%f ", c);
        }
        printf("\n");
        for (auto const &c : r5) {
          printf("%f ", c);
        }
        printf("\n");
        for (auto const &c : r6) {
          printf("%f ", c);
        }
        printf("\n");
        return false;
      }
      uint64_t tc2 = g->triangle_count();
      uint64_t tc = matrix.triangle_count();
      if (tc != tc2) {
        printf("failed triangle count, got %lu, expected %lu\n",tc2, tc);
        return false; 
      }
      
  return true;
}

bool verify_adjacency_hashmap() {
  srand(0);
  uint32_t node_counts[4] = {5, 10, 30, 100};
  uint32_t edge_counts[4] = {5, 20, 100, 1000};
  int trials = 1000;
  for (int a = 0; a < 4; a++) {
    for (int j = 1; j < trials; j++) {
      uint32_t num_nodes = node_counts[a];
      uint32_t num_edges = edge_counts[a];
      /*
      printf("########## starting trial %d out of %d, nodes = %u, edges = %u "
             "##########\n",
             j, trials, num_nodes, num_edges);
       */
      // AdjacencyMatrix matrix = AdjacencyMatrix(num_nodes);
      AdjacencyHashMap adjhash = AdjacencyHashMap(num_nodes);

      if(!verify_structure(&adjhash, num_nodes, num_edges, trials)){
        return false;
      }
    }
  }
  return true;
}

bool verify_adjacency_list() {
  srand(0);
  uint32_t node_counts[4] = {5, 10, 30, 100};
  uint32_t edge_counts[4] = {5, 20, 100, 1000};
  int trials = 1000;
  for (int a = 0; a < 4; a++) {
    for (int j = 1; j < trials; j++) {
      uint32_t num_nodes = node_counts[a];
      uint32_t num_edges = edge_counts[a];
      AdjacencyList adjlist = AdjacencyList(num_nodes);
      if(!verify_structure(&adjlist, num_nodes, num_edges, trials)) {
        return false;
      }
    }
  }
  return true;
}
bool verify_blocked_adjacency_list() {
  srand(0);
  uint32_t node_counts[4] = {5, 10, 30, 100};
  uint32_t edge_counts[4] = {5, 20, 100, 1000};
  int trials = 1000;
  for (int a = 0; a < 4; a++) {
    for (int j = 1; j < trials; j++) {
      uint32_t num_nodes = node_counts[a];
      uint32_t num_edges = edge_counts[a];
      BlockedAdjacencyList adjlist = BlockedAdjacencyList(num_nodes);
      if(!verify_structure(&adjlist, num_nodes, num_edges, trials)) {
        return false;
      }
    }
  }
  return true;
}


bool verify_adjacency_vector() {
  srand(0);
  uint32_t node_counts[4] = {5, 10, 30, 100};
  uint32_t edge_counts[4] = {5, 20, 100, 1000};
  int trials = 1000;
  for (int a = 0; a < 4; a++) {
    for (int j = 1; j < trials; j++) {
      uint32_t num_nodes = node_counts[a];
      uint32_t num_edges = edge_counts[a];
      AdjacencyVector adjvec = AdjacencyVector(num_nodes);
      if(!verify_structure(&adjvec, num_nodes, num_edges, trials)) {
        return false;
      }
    }
  }
  return true;
}

bool verify_csr() {
  srand(0);
  uint32_t node_counts[5] = {5, 10, 30, 100, 1000};
  uint32_t edge_counts[5] = {5, 20, 100, 1000, 20000};
  int trials = 1000;
  for (int a = 0; a < 4; a++) {
    for (int j = 1; j < trials; j++) {
      uint32_t num_nodes = node_counts[a];
      uint32_t num_edges = edge_counts[a];
      CSR *csr = new CSR(num_nodes);
      if(!verify_structure(csr, num_nodes, num_edges, trials)) {
        return false;
      }
    }
  }
  return true;
}


bool verify_pcsr() {
  srand(0);
  uint32_t node_counts[4] = {5, 10, 30, 100};
  uint32_t edge_counts[4] = {5, 20, 100, 1000};
  int trials = 1000;
  for (int a = 0; a < 4; a++) {
    for (int j = 1; j < trials; j++) {
      uint32_t num_nodes = node_counts[a];
      uint32_t num_edges = edge_counts[a];
      OFM *ofm = new OFM(num_nodes);
      if(!verify_structure(ofm, num_nodes, num_edges, trials)) {
        return false;
      }
    }
  }
  return true;
}

/******* timing functions ********/
bool time_structure(Graph* g, uint32_t num_nodes, uint32_t num_edges) {
  srand(0);

  struct timeval start, end;
  uint32_t duration;

  
  // add edge update
  uint32_t *srcs = (uint32_t *) malloc(sizeof(uint32_t) * num_edges);
  uint32_t *dests = (uint32_t *) malloc(sizeof(uint32_t) * num_edges);
  uint32_t *vals = (uint32_t *) malloc(sizeof(uint32_t) * num_edges);

  for (uint32_t k = 0; k < num_edges; k++) {
    srcs[k] = rand_in_range(num_nodes);
    dests[k] = rand_in_range(num_nodes);
    vals[k] = rand_in_range(MAXVAL) + 1;
  }
  gettimeofday(&start, NULL);
  for(int i = 0; i < num_edges; i+=batch_size) {
    g->add_edge_batch_update(&srcs[i], &dests[i], &vals[i], batch_size);
  }

  gettimeofday(&end, NULL);
  free(srcs);
  free(dests);
  free(vals);
  duration =
            (end.tv_sec - start.tv_sec) * 1000000u + (end.tv_usec - start.tv_usec);
  if (num_edges > 0) {
    printf("add_edge_update :%f: \n", ((double)duration) / (1000000));
  } else {
    printf("add_edge_update :N/A: ");
  }


    vector<uint32_t> test_vector(num_nodes, 0);
    for (uint32_t i = 0; i < num_nodes; i++) {
      test_vector[i] = rand_in_range(MAXVAL);
    }

  // spvm
    
    gettimeofday(&start, NULL);
    // start = std::clock();
    for (int spvm_iters = 0; spvm_iters < 100; spvm_iters++) {
      test_vector = g->sparse_matrix_vector_multiplication(test_vector);
    }

    gettimeofday(&end, NULL);
    duration = (end.tv_sec-start.tv_sec)*1000000u+(end.tv_usec-start.tv_usec);
    printf("spvm: %f: ", ((double)duration)/1000000);
    printf("spvm[0]: %d: ", test_vector[0]);
    gettimeofday(&start, NULL);

  // bfs
  
    int start_node = rand_in_range(num_nodes);
    gettimeofday(&start, NULL);
    for (int bfs_iters = 0; bfs_iters < 100; bfs_iters++) {
      vector<int32_t> r4 = g->bfs(start_node);
      start_node = r4[0] % num_nodes;
    }
    gettimeofday(&end, NULL);
    duration = (end.tv_sec-start.tv_sec)*1000000u+(end.tv_usec-start.tv_usec);
    printf("bfs in %f,", ((double)duration)/1000000);
    printf("sn = %d,", start_node);
    

  // pagerank
  vector<double> test_vector2(num_nodes, 0);
  for (uint32_t i = 0; i < num_nodes; i++) {
    test_vector2[i] = ((double)rand_in_range(MAXVAL)) / MAXVAL;
  }
  gettimeofday(&start, NULL);
  for (int pr_iters = 0; pr_iters < 100; pr_iters++) {
    vector<double> test_vector2 = g->pagerank(test_vector2);
  }
  gettimeofday(&end, NULL);
  duration = (end.tv_sec-start.tv_sec)*1000000u+(end.tv_usec-start.tv_usec);
  printf("pagerank: %f:", ((double)duration)/1000000);
  printf("pr[0]: %f\n", test_vector2[0]);

  return true;

}



bool time_adjlist(uint32_t num_nodes, uint32_t num_edges, int trials) {
  srand(0);
  struct timeval start, end;
  uint32_t duration;
  for (int j = 0; j < trials; j++) {
    gettimeofday(&start, NULL);
    AdjacencyList *al = new AdjacencyList(num_nodes);
    gettimeofday(&end, NULL);
    duration =
        (end.tv_sec - start.tv_sec) * 1000000u + (end.tv_usec - start.tv_usec);
    printf("nodes: %d: edges : %d: ", num_nodes, num_edges);
    printf("created structure in :%f: ", ((double)duration) / 1000000);

    time_structure(al, num_nodes, num_edges);
    delete al;
  }
  return true;
}
bool time_blocked_adjlist(uint32_t num_nodes, uint32_t num_edges, int trials) {
  srand(0);
  struct timeval start, end;
  uint32_t duration;
  for (int j = 0; j < trials; j++) {
    gettimeofday(&start, NULL);
    BlockedAdjacencyList *bal = new BlockedAdjacencyList(num_nodes);
    gettimeofday(&end, NULL);
    duration =
        (end.tv_sec - start.tv_sec) * 1000000u + (end.tv_usec - start.tv_usec);
    printf("nodes: %d: edges : %d: ", num_nodes, num_edges);
    printf("created structure in :%f: ", ((double)duration) / 1000000);

    time_structure(bal, num_nodes, num_edges);
    delete bal;
  }
  return true;
}

bool time_pcsr(uint32_t num_nodes, uint32_t num_edges, int trials) {
  srand(0);
  struct timeval start, end;
  uint32_t duration;
  for (int j = 0; j < trials; j++) {
    gettimeofday(&start, NULL);
    OFM *ofm = new OFM(num_nodes);
    gettimeofday(&end, NULL);
    duration =
        (end.tv_sec - start.tv_sec) * 1000000u + (end.tv_usec - start.tv_usec);
    printf("nodes: %d: edges : %d: ", num_nodes, num_edges);
    printf("created structure in :%f: ", ((double)duration) / 1000000);
    gettimeofday(&start, NULL);

    time_structure(ofm, num_nodes, num_edges);
    delete ofm;
  }
  return true;
}

bool time_csr(uint32_t num_nodes, uint32_t num_edges, int trials) {
  srand(0);
  for (int j = 0; j < trials; j++) {
    OFM *ofm = new OFM(num_nodes);

    for (uint32_t k = 0; k < num_edges; k++) {
      uint32_t src = rand_in_range(num_nodes);
      uint32_t dest = rand_in_range(num_nodes);
      uint32_t val = rand_in_range(MAXVAL) + 1;

      ofm->add_edge_update(src, dest, val);
    }
    CSR csr = CSR(1);
    csr.convert(ofm);
    delete ofm;
    printf("nodes: %d: edges : %d: ", num_nodes, num_edges);
    printf("created structure in :N/A: ");
    time_structure(&csr, num_nodes, 0);

  }
  return true;
}

bool time_adjhashmap(uint32_t num_nodes, uint32_t num_edges, int trials) {
  srand(0);

  struct timeval start, end;
  uint32_t duration;
  for (int j = 0; j < trials; j++) {
    // start = std::clock();
    gettimeofday(&start, NULL);
    AdjacencyHashMap *ahm = new AdjacencyHashMap(num_nodes);
    gettimeofday(&end, NULL);
    // duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
    duration = (end.tv_sec - start.tv_sec) * 1000000u + (end.tv_usec - start.tv_usec);
    printf("nodes: %d: edges : %d: ", num_nodes, num_edges);
    printf("created structure in :%f: ", ((double)duration) / 1000000);
    gettimeofday(&start, NULL);

    time_structure(ahm, num_nodes, num_edges);
    delete ahm;
  }
  return true;
}
bool time_adjvector(uint32_t num_nodes, uint32_t num_edges, int trials) {
  srand(0);

  struct timeval start, end;
  uint32_t duration;
  for (int j = 0; j < trials; j++) {
    // start = std::clock();
    gettimeofday(&start, NULL);
    AdjacencyVector *av = new AdjacencyVector(num_nodes);
    gettimeofday(&end, NULL);
    duration = (end.tv_sec - start.tv_sec) * 1000000u + (end.tv_usec - start.tv_usec);
    printf("nodes: %d: edges : %d: ", num_nodes, num_edges);
    printf("created structure in :%f: ", ((double)duration) / 1000000);

    time_structure(av, num_nodes, num_edges);
    delete av;
  }
  return true;
}

vector<uint32_t> parent_to_depth(vector<int32_t> parents) {
  vector<uint32_t> depths(parents.size(), UINT32_MAX);
  for (int i = 0; i < parents.size(); i++) {
    for (int j = 0; j < parents.size(); j++) {
      if (parents[j] == j) {
        depths[j] = 0;
      } else if (parents[j] < 0) {
        continue;
      } else if (depths[parents[j]] < UINT32_MAX) {
        depths[j] = std::min(depths[parents[j]] + 1, depths[j]);
      }
    }
  }
  return depths;
}
vector<int32_t> pvector_to_vector(pvector<int32_t> &pvec, int len) {
  vector<int32_t> vec(len);
  for (int i = 0; i < len; i++) {
    vec[i] = pvec[i];
  }
  return vec;
}
bool time_structure_real_world_graphs(Graph* g, int graph_num) {
  srand(0);

  struct timeval start, end;
  uint32_t duration;

  
  // add edge update
  std::vector<uint32_t> srcs;
  std::vector<uint32_t> dests;
  std::vector<uint32_t> vals;
  string filename;
  if (graph_num == 0) {
    filename = "graphs/rand-out-soc-Slashdot0811.txt";
  }
  else if (graph_num == 1) {
    filename = "graphs/rand-soc-pokec-relationships.txt";
  }
  else if (graph_num == 2) {
    filename = "graphs/rand-out-soc-LiveJournal1.txt";
  } else if (graph_num == 3) {
    filename = "/efs/gapbs/benchmark/graphs/raw/twitter.el";
   } else if (graph_num == 4) {
    filename = "/efs/home/____/aspen/code/inputs/com-orkut.ungraph.txt";
  } else {
    printf("bad graph_num\n");
    return false;
  }
  ifstream myfile(filename.c_str());
  string line;
  if (myfile.is_open()) {
    while (getline(myfile, line)) {
      uint32_t elem_1;
      uint32_t elem_2;
      sscanf(line.c_str(), "%u   %u", &elem_1, &elem_2);
      
      int src = elem_1 - 1;

      while (src >= g->get_n()) {
        g->add_node();
      }
      int dest = elem_2 - 1;

      while (dest >= g->get_n()) {
        g->add_node();
      }
      // making the graph symetric
      srcs.push_back(src);
      srcs.push_back(dest);
      dests.push_back(dest);
      dests.push_back(src);
      uint32_t val = 1;//rand_in_range(MAXVAL);
      vals.push_back(val);
      vals.push_back(val);
    }
    myfile.close();
    // return 0;
  } else {
    printf("file was not opened\n");
  }
  uint32_t num_edges = srcs.size();
  gettimeofday(&start, NULL);
  int i = 0;
  for(; i < num_edges-batch_size; i+=batch_size) {
    g->add_edge_batch_update(&srcs[i], &dests[i], &vals[i], batch_size);
  }
  g->add_edge_batch_update(&srcs[i], &dests[i], &vals[i], num_edges % batch_size);

  uint32_t num_nodes = g->get_n();

  gettimeofday(&end, NULL);
  duration =
            (end.tv_sec - start.tv_sec) * 1000000u + (end.tv_usec - start.tv_usec);
  if (num_edges > 0) {
    printf("add_edge_update :%f:", ((double)duration) / (1000000));
  } else {
    printf("add_edge_update :N/A: ");
  }

    vector<uint32_t> test_vector(num_nodes, 0);
    for (uint32_t i = 0; i < num_nodes; i++) {
      test_vector[i] = rand_in_range(MAXVAL);
    }

  // spvm
   /* 
    gettimeofday(&start, NULL);
    // start = std::clock();
    for (int spvm_iters = 0; spvm_iters < 100; spvm_iters++) {
      test_vector = g->sparse_matrix_vector_multiplication(test_vector);
    }

    gettimeofday(&end, NULL);
    duration = (end.tv_sec-start.tv_sec)*1000000u+(end.tv_usec-start.tv_usec);
    printf("spvm: %f: ", ((double)duration)/1000000);
    printf("spvm[0]: %d: ", test_vector[0]);
    gettimeofday(&start, NULL);
*/
    
  // pagerank
/*
  vector<uint32_t> test_vector2(num_nodes, 0);
  for (uint32_t i = 0; i < num_nodes; i++) {
    test_vector2[i] = ((uint32_t)rand_in_range(MAXVAL)) / MAXVAL;
  }
  gettimeofday(&start, NULL);
  for (int pr_iters = 0; pr_iters < 100; pr_iters++) {
    vector<uint32_t> test_vector2 = g->pagerank(test_vector2);
  }
  gettimeofday(&end, NULL);
  duration = (end.tv_sec-start.tv_sec)*1000000u+(end.tv_usec-start.tv_usec);
  printf("pagerank: %f:", ((double)duration)/1000000);
  printf("pr[0]: %u\n", test_vector2[0]);
  */
/*
    printf("\nstarting make symetric\n");
    g->make_symetric();
    printf("done make symetric\n");
  */  

    //triangle count
    /*
    gettimeofday(&start, NULL);
    uint64_t total = g->triangle_count();
    gettimeofday(&end, NULL);
    duration = (end.tv_sec-start.tv_sec)*1000000u+(end.tv_usec-start.tv_usec);
    printf("triangle count in %f,", ((double)duration)/1000000);
    printf("count = %lu\n", total);
    */  
  
    // bfs
  
    int start_node = 10012;//rand_in_range(num_nodes);
    gettimeofday(&start, NULL);
    vector<int32_t> r4 = g->bfs(start_node);
    start_node = r4[0];
    gettimeofday(&end, NULL);
    duration = (end.tv_sec-start.tv_sec)*1000000u+(end.tv_usec-start.tv_usec);
    printf("bfs in %f,", ((double)duration)/1000000);
    printf("parent of 0 = %d,", start_node);
    printf("\n");

    start_node = 10012;//rand_in_range(num_nodes);
    //parallel_bfs
    gettimeofday(&start, NULL);
    for (int i = 0; i < 100; i++) {
      pvector<int32_t> r7 = g->parallel_bfs(start_node, num_edges);
      start_node = r7[0] % num_nodes;
    }
    gettimeofday(&end, NULL);
    duration = (end.tv_sec-start.tv_sec)*1000000u+(end.tv_usec-start.tv_usec);
    printf("parallel_bfs in %f,", ((double)duration)/1000000);
    printf("parent of 0 = %d\n", start_node);
    /* 
    vector<uint32_t> d1 = parent_to_depth(r4);
    vector<uint32_t> d2 = parent_to_depth(pvector_to_vector(r7, g->get_n()));
    if (d1 != d2) {
      printf("we have an issue\n");
      for (int i = 0; i < d1.size(); i++) {
        if (d1[i] != d2[i]) printf("DIFFER");
        printf("position %d, seq %u, par %u, parent seq %d, parent par %d\n", i, d1[i], d2[i], r4[i], r7[i]);
      }
    }
    */
    
  return true;
}

bool time_ofm_real_world_graphs(int graph_num, int iters = 10) {
  srand(0);

  struct timeval start, end;
  uint32_t duration;

  
  // add edge update
  std::vector<uint32_t> srcs;
  std::vector<uint32_t> dests;
  //std::vector<uint32_t> vals;
  string filename;
  bool fix = false;
  bool zero_index = false;
  int start_node = 0;//rand_in_range(num_nodes);
  if (graph_num == 0) {
    filename = "graphs/rand-out-soc-Slashdot0811.txt";
  }
  else if (graph_num == 1) {
    filename = "graphs/rand-soc-pokec-relationships.txt";
  }
  else if (graph_num == 2) {
    filename = "graphs/rand-out-soc-LiveJournal1.txt";
   } else if (graph_num == 3) {
    filename = "/efs/gapbs/benchmark/graphs/raw/twitter.el";
  } 
  else if (graph_num == 4) {
    filename = "/efs/home/____/aspen/code/inputs/com-orkut.ungraph.txt";
  } 
  else if (graph_num == 5) {
    filename = "rmat_sym_rand.el";
    fix = true;
    zero_index = true;
  } else if (graph_num == 6) {
    filename = "graphs/rand-out-soc-LiveJournal1.txt.shuf";
    zero_index = true;
    start_node = 9;
  } else if (graph_num == 7) {
    filename = "/efs/home/____/aspen/code/inputs/com-orkut.ungraph.txt.shuf";
    zero_index = true;
    start_node = 28;
  } else if (graph_num == 8) {
    filename = "../er_graph.el";
    start_node = 0;
    zero_index = true;
  } else if (graph_num == 9) {
    filename = "rmat_ligra.el.shuf";
    start_node = 19372;
    zero_index = true;
  } else if (graph_num == 10) {
	  printf("adding random edges, num_nodes = 50,000,000, num_edges = 2,000,000,000\n");
  } else {
    printf("bad graph_num, %d\n", graph_num);
    return false;
  }

  uint32_t num_nodes = 0;
  if (graph_num != 10) {
    ifstream myfile(filename.c_str());
    string line;
    if (myfile.is_open()) {
      while (getline(myfile, line)) {
        uint32_t elem_1 = 0;
        uint32_t elem_2 = 0;
        if (line[0] == '#') {
          continue;
        }
        sscanf(line.c_str(), "%u   %u", &elem_1, &elem_2);
        elem_1 += zero_index;
        elem_2 += zero_index;
        num_nodes = max(num_nodes, elem_1);
        num_nodes = max(num_nodes, elem_2);
        uint32_t src = elem_1 - 1;

        uint32_t dest = elem_2 - 1;
        if (src == dest) {
          continue;
        }

        // making the graph symetric
        srcs.push_back(src);
        dests.push_back(dest);
        //uint32_t val = 1;//rand_in_range(MAXVAL);
        //vals.push_back(val);
        if (!fix) {
          srcs.push_back(dest);
          dests.push_back(src);
          //vals.push_back(val);
        }
      }
      myfile.close();
      // return 0;
    } else {
      printf("file was not opened\n");
    }
  } else {
    srcs.resize(2000000000UL);
    dests.resize(2000000000UL);
    num_nodes=50000000;
    uint64_t degree=srcs.size() / num_nodes;
    for(uint32_t i = 0; i < num_nodes; i++) {
      std::mt19937 local_rng;
      local_rng.seed(i);
      for (uint64_t j = 0; j < degree; j++) {
        srcs[i+j*num_nodes] = i;
        dests[i+j*num_nodes] = local_rng() % num_nodes;
      }
    }
  }
  printf("done reading in the graph\n");
  uint64_t num_edges = srcs.size();
  int batch = 10000;
  if (num_edges > 1000000) {
    batch = 1000000;
  }
  
  if (num_edges > 100000000) {
    batch = 100000000;
  }
  OFM g = OFM(num_nodes);
  printf("done adding the nodes\n");
  gettimeofday(&start, NULL);
  uint64_t i = 0;
  for(; i < num_edges-batch; i+=batch) {
    g.add_edge_batch_update_no_val(&srcs[i], &dests[i], batch);
    printf("added %lu edges so far\n", i);
  }
  g.add_edge_batch_update_no_val(&srcs[i], &dests[i], num_edges % batch);

  num_nodes = g.get_n();

  gettimeofday(&end, NULL);
  duration =
            (end.tv_sec - start.tv_sec) * 1000000u + (end.tv_usec - start.tv_usec);
  if (num_edges > 0) {
    printf("add_edge_update :%f\n", ((double)duration) / (1000000));
  } else {
    printf("add_edge_update :N/A: ");
  }
  printf("number of edges = %lu, size of graph = %lu\n", g.num_edges(), g.get_size());

  gettimeofday(&start, NULL);
  for (int i = 0; i < iters; i++) {
    CSR_simple csr = g.convert_to_csr();
    free(csr.nodes);
    free(csr.dests);
    free(csr.values);
  }
  gettimeofday(&end, NULL);
  duration = (end.tv_sec - start.tv_sec) * 1000000u + (end.tv_usec - start.tv_usec);
  printf("convert :%f\n", ((double)duration) / (1000000*iters));

  vector<uint32_t> test_vector(num_nodes, 0);
  for (uint32_t i = 0; i < num_nodes; i++) {
    test_vector[i] = rand_in_range(MAXVAL);
  }
  gettimeofday(&start, NULL);
  uint32_t checker = 0;
  int32_t *bfs_values = NULL;
  for (int j = 0; j < iters; j++) {
    if (bfs_values != NULL) {
      free(bfs_values);
    }
    bfs_values = BFS_with_edge_map(g,start_node);
    checker += bfs_values[0];
  }
  gettimeofday(&end, NULL);
  duration = (end.tv_sec-start.tv_sec)*1000000u+(end.tv_usec-start.tv_usec);
  printf("parallel_bfs with edge_map in %f,", ((double)duration)/(1000000*iters));
  
  printf("parent of 0 = %d\n", checker / iters);
  free(bfs_values);
  vector<double> test_vector4(num_nodes);
  checker=0; 
  gettimeofday(&start, NULL);
  for (int j = 0; j < iters; j++) {
    test_vector4 = PR_S<double>(g, 5);
    checker += test_vector4[0];
  }
  gettimeofday(&end, NULL);
  duration = (end.tv_sec-start.tv_sec)*1000000u+(end.tv_usec-start.tv_usec);
  printf("pagerank with edge_map in %f, value of 0 = %u\n", ((double)duration)/(1000000*iters), checker/iters);

  checker = 0;
  gettimeofday(&start, NULL);
  uint32_t *cc_values = NULL;
  for (int j = 0; j < iters; j++) {
    if (cc_values != NULL) {
      free(cc_values);
    }
    cc_values = CC(g);
    checker+= cc_values[0];
  }
  gettimeofday(&end, NULL);
  duration = (end.tv_sec-start.tv_sec)*1000000u+(end.tv_usec-start.tv_usec);
  printf("CC with edge_map in %f, component of 0 is %u\n", ((double)duration)/(1000000*iters), checker/iters);
  if (cc_values != NULL) {
    free(cc_values);
  }
  
  std::vector<double> bc_values;
  checker=0;
  gettimeofday(&start, NULL);
  for (int j = 0; j < iters; j++) {
    bc_values = BC(g, start_node);
    checker += bc_values[0];
  }
  gettimeofday(&end, NULL);
  duration = (end.tv_sec-start.tv_sec)*1000000u+(end.tv_usec-start.tv_usec);
  printf("BC with edge_map in %f, value of 0 is %u\n", ((double)duration)/(1000000*iters), checker/iters);

  gettimeofday(&start, NULL);
  for (int j = 0; j < iters; j++) {
    TC(g);
  }
  gettimeofday(&end, NULL);
  duration = (end.tv_sec-start.tv_sec)*1000000u+(end.tv_usec-start.tv_usec);
  printf("TC with edge_map in %f\n", ((double)duration)/(1000000*iters));
  return true;
}

bool parallel_updates(int graph_num) {
  struct timeval start, end;
  uint32_t duration;

  
  // add edge update
  string filename;
  int fix = 0;
  if (graph_num == 0) {
    filename = "graphs/rand-out-soc-Slashdot0811.txt";
  }
  else if (graph_num == 1) {
    filename = "graphs/rand-soc-pokec-relationships.txt";
  }
  else if (graph_num == 2) {
    filename = "graphs/rand-out-soc-LiveJournal1.txt";
  } 
  else if (graph_num == 3) {
    filename = "/efs/gapbs/benchmark/graphs/raw/twitter.el";
  } else if (graph_num == 4) {
    filename = "/efs/home/____/aspen/code/inputs/com-orkut.ungraph.txt";
  } else if (graph_num == 5) {
    filename = "rmat.gr";
    fix = 1;
  } else if (graph_num == 6) {
    filename = "graphs/rand-out-soc-LiveJournal1.txt.shuf";
    fix = 1;
  } else if (graph_num == 7) {
    filename = "/efs/home/____/aspen/code/inputs/com-orkut.ungraph.txt.shuf";
    fix = 1;
  } else if (graph_num == 8) {
    filename = "/efs/home/____/aspen/code/inputs/rmat_sym_rand.el.shuf";
    fix = 1;
  } else if (graph_num == 9) {
    filename = "rmat_ligra.el.shuf";
    fix = 1;
  } else {
    printf("bad graph_num\n");
    return false;
  }
  printf("counting the lines in the file\n");
  ifstream line_counter(filename.c_str());
  uint64_t line_count = 0;
  string line_;
  if (line_counter.is_open()) {
    while (getline(line_counter, line_)) {
      if ('#' == line_[0]) {continue;}
      line_count++;
    }
  } else {
    printf("coult not open file %s\n", filename.c_str());
    exit(-1); 
  }
  if (line_count > UINT32_MAX) {
    printf("the file is too big, has %lu lines\n", line_count);
    exit(-1);
  }
  uint32_t *srcs = (uint32_t *) malloc(line_count*4);
  uint32_t *dests = (uint32_t *) malloc(line_count*4);
  ifstream myfile(filename.c_str());
  string line;
  uint32_t num_nodes = 0;
  printf("reading in the graph\n");
  uint64_t num_edges = 0;
  if (myfile.is_open()) {
    while (getline(myfile, line)) {
      if ('#' == line[0]) {continue;}
      uint32_t elem_1;
      uint32_t elem_2;
      sscanf(line.c_str(), "%u   %u", &elem_1, &elem_2);
      elem_1 += fix;
      elem_2 += fix;
      num_nodes = max(num_nodes, elem_1);
      num_nodes = max(num_nodes, elem_2);
      uint32_t src = elem_1 - 1;

      uint32_t dest = elem_2 - 1;
      if (src == dest) {
        continue;
      }
      if (src > 1000000000) {
        printf("bad read got src = %u\n", src);
      }
      if (dest > 1000000000) {
        printf("bad read got dest = %u\n", dest);
      }
      // making the graph symetric
      srcs[num_edges] = src;
      dests[num_edges] = dest;
      num_edges+=1;
    }
    myfile.close();
    // return 0;
  } else {
    printf("file was not opened\n");
  }
printf("done reading in the graph, N = %u, M = %lu (directed)\n", num_nodes, num_edges);
  
  printf("starting to make the graph\n");
  gettimeofday(&start, NULL);
  OFM G = OFM(num_nodes);
  int i = 0;
  uint32_t batch = 10000;
  if (num_edges > 1000000) {
    batch = 1000000;
  }
  for(; i < num_edges-batch; i+=batch) {
    G.add_edge_batch_update_no_val(&srcs[i], &dests[i], batch);
  }
  G.add_edge_batch_update_no_val(&srcs[i], &dests[i], num_edges % batch);
  i = 0;
  for(; i < num_edges-batch; i+=batch) {
    G.add_edge_batch_update_no_val(&dests[i], &srcs[i], batch);
  }
  G.add_edge_batch_update_no_val(&dests[i], &srcs[i], num_edges % batch);
  gettimeofday(&end, NULL);
  duration = (end.tv_sec - start.tv_sec) * 1000000u + (end.tv_usec - start.tv_usec);
  printf("done building the graph M = %lu, took %u\n", G.num_edges(), duration/1000000);
  free(srcs);
  free(dests);
  
  std::vector<uint32_t> update_sizes = {10, 100, 1000,10000,100000,1000000,10000000};//,100000000,1000000000};
  auto r = random_aspen();
  auto update_times = std::vector<double>();
  size_t n_trials = 1;
  for (size_t us=0; us<update_sizes.size(); us++) {
    uint64_t avg_insert = 0;
    uint64_t avg_delete = 0;
    cout << "Running bs: " << update_sizes[us] << endl;

    if (update_sizes[us] <= 10000000) {
      n_trials = 20;
    }
    else {
      n_trials = 5;
    }
    size_t updates_to_run = update_sizes[us];
    for (size_t ts=0; ts<n_trials; ts++) {
      //printf("the graph has size %lu the start, and %lu edges\n", G.get_size(), G.num_edges());
      //printf("starting copy\n");
      //OFM g = OFM(G);
      //printf("done copy\n");
      uint32_t num_nodes = G.get_n();

      std::vector<uint32_t> new_srcs(updates_to_run);
      std::vector<uint32_t> new_dests(updates_to_run);
      double a = 0.5;
      double b = 0.1;
      double c = 0.1;
      size_t nn = 1 << (log2_up(num_nodes) - 1);
      auto rmat = rMat<uint32_t>(nn, r.ith_rand(100+ts), a, b, c); 
      //std::vector<std::pair<uint32_t, uint32_t>> edges(updates_to_run);
      cilk_for(uint32_t i = 0; i < updates_to_run; i++) {
        std::pair<uint32_t, uint32_t> edge = rmat(i);
        //edges[i] = edge;
        new_srcs[i] = edge.first;
        new_dests[i] = edge.second;
        //new_srcs[i] = rand_in_range(num_nodes);
        //new_dests[i] = rand_in_range(num_nodes);
      }
      /*
      std::sort(edges.begin(), edges.end());
      cilk_for(uint32_t i = 0; i < updates_to_run; i++) {
        new_srcs[i] = edges[i].first;
        new_dests[i] = edges[i].second;
      }
      */
      gettimeofday(&start, NULL);
      G.add_edge_batch_update_no_val(new_srcs.data(), new_dests.data(), updates_to_run);
      gettimeofday(&end, NULL);
      duration =
                (end.tv_sec - start.tv_sec) * 1000000u + (end.tv_usec - start.tv_usec);
      avg_insert += duration;
      //printf("the graph has size %lu after the updates, and %lu edges\n", G.get_size(), G.num_edges());
      gettimeofday(&start, NULL);
      G.remove_edge_batch(new_srcs.data(), new_dests.data(), updates_to_run);
      gettimeofday(&end, NULL);
      duration =
                (end.tv_sec - start.tv_sec) * 1000000u + (end.tv_usec - start.tv_usec);
      avg_delete += duration;
    }
    double time_i = (double) avg_insert;
    double time_d = (double) avg_delete;
    time_i /= n_trials;
    time_i /= 1000000u;
    double insert_throughput = updates_to_run / time_i;
    time_d /= n_trials;
    time_d /= 1000000u;
    double delete_throughput = updates_to_run / time_d;
    printf("batch_size = %zu, average insert: %lu, throughput %e\n", updates_to_run, avg_insert/n_trials, insert_throughput);
    printf("batch_size = %zu, average delete: %lu, throughput %e\n", updates_to_run, avg_delete/n_trials, delete_throughput);
  }
  
  ofstream mfile;
  mfile.open("rmat_sym2.el");
  for(auto &edge : G.get_edges()) {
    mfile << std::get<0>(edge) << "\t" << std::get<1>(edge) << "\n";
  }
  mfile.close();
  
  return true;
}
void write_rmat_to_file(uint32_t num_nodes, uint32_t num_edges) {
      double a = 0.5;
      double b = 0.1;
      double c = 0.1;
      size_t nn = 1 << (log2_up(num_nodes) - 1);
      auto r = random_aspen();
      auto rmat = rMat<uint64_t>(nn, r.ith_rand(0), a, b, c); 
      printf("writting rmat graph to file, num_nodes = %zu. num_edges = %u\n", nn, num_edges);
      ofstream mfile;
      mfile.open("rmat.gr");
      for(uint32_t i = 0; i < num_edges; i++) {
        std::pair<uint64_t, uint64_t> edge = rmat(i);
        if (edge.first <= edge.second) continue;
        mfile << edge.first << "\t" << edge.second << "\n";
      }
      mfile.close();

}
bool worst_case(int graph_num) {
  struct timeval start, end;
  
  // add edge update
  string filename;
  bool fix = false;
  if (graph_num == 0) {
    filename = "graphs/rand-out-soc-Slashdot0811.txt";
  }
  else if (graph_num == 1) {
    filename = "graphs/rand-soc-pokec-relationships.txt";
  }
  else if (graph_num == 2) {
    filename = "graphs/rand-out-soc-LiveJournal1.txt";
  } 
  else if (graph_num == 3) {
    filename = "/efs/gapbs/benchmark/graphs/raw/twitter.el";
  } else if (graph_num == 4) {
    filename = "/efs/home/____/aspen/code/inputs/com-orkut.ungraph.txt";
  } else if (graph_num == 5) {
    filename = "rmat_sym_rand.el";
    fix = true;
  } else {
    printf("bad graph_num\n");
    return false;
  }
  printf("counting the lines in the file\n");
  ifstream line_counter(filename.c_str());
  uint64_t line_count = 0;
  string line_;
  if (line_counter.is_open()) {
    while (getline(line_counter, line_)) {
      if ('#' == line_[0]) {continue;}
      line_count++;
    }
  } else {
    printf("coult not open file %s\n", filename.c_str());
    exit(-1); 
  }
  if (line_count > UINT32_MAX) {
    printf("the file is too big, has %lu lines\n", line_count);
    exit(-1);
  }
  uint32_t *srcs = (uint32_t *) malloc(line_count*4);
  uint32_t *dests = (uint32_t *) malloc(line_count*4);
  ifstream myfile(filename.c_str());
  string line;
  uint32_t num_nodes = 0;
  printf("reading in the graph\n");
  uint64_t num_edges = 0;
  if (myfile.is_open()) {
    while (getline(myfile, line)) {
      if ('#' == line[0]) {continue;}
      uint32_t elem_1;
      uint32_t elem_2;
      sscanf(line.c_str(), "%u   %u", &elem_1, &elem_2);
      elem_1 += fix;
      elem_2 += fix;
      num_nodes = max(num_nodes, elem_1);
      num_nodes = max(num_nodes, elem_2);
      uint32_t src = elem_1 - 1;

      uint32_t dest = elem_2 - 1;
      if (src == dest) {
        continue;
      }
      if (src > 1000000000) {
        printf("bad read got src = %u\n", src);
      }
      if (dest > 1000000000) {
        printf("bad read got dest = %u\n", dest);
      }
      srcs[num_edges] = src;
      dests[num_edges] = dest;
      num_edges+=1;
    }
    myfile.close();
    // return 0;
  } else {
    printf("file was not opened\n");
  }
  printf("done reading in the graph, N = %u, M = %lu (directed)\n", num_nodes, num_edges);
  printf("starting to make the graph\n");
  for (int w = 32; w >= 1; w/=2) {
    uint32_t duration;
    std::vector<uint32_t> durations;
    std::vector<uint32_t> insert_num;
    char str[2];
    sprintf(str, "%d", w);
    __cilkrts_end_cilk(); 
    __cilkrts_set_param("nworkers", str);
    printf("workers = %d\n", w);
    OFM G = OFM(num_nodes);
    int cutoff = 1000;
    int less_than_cutoff = 0;
    int i = 0;
    uint32_t num = 0;
    for(; i < num_edges; i++) {
      gettimeofday(&start, NULL);
      G.add_edge_update(srcs[i], dests[i], 1);
      gettimeofday(&end, NULL);
      duration = (end.tv_sec - start.tv_sec) * 1000000u + (end.tv_usec - start.tv_usec);
      if (duration < cutoff) {
        less_than_cutoff++;
      } else {
        durations.push_back(duration);
        insert_num.push_back(num);
      }
      num+=1;
    }
    if (!fix) {
      i = 0;
      for(; i < num_edges; i++) {
        gettimeofday(&start, NULL);
        G.add_edge_update(dests[i], srcs[i], 1);
        gettimeofday(&end, NULL);
        duration = (end.tv_sec - start.tv_sec) * 1000000u + (end.tv_usec - start.tv_usec);
        if (duration < cutoff) {
          less_than_cutoff++;
        } else {
          durations.push_back(duration);
          insert_num.push_back(num);
        }
        num+=1;
      }
      less_than_cutoff /= 2;
    }
    printf("the graph has size %lu\n", G.get_size());
    printf("%f were less than the cutoff\n", ((double) less_than_cutoff)/num_edges);
    printf("the max time was %u\n", *std::max_element(durations.begin(), durations.end()));
  }
  free(srcs);
  free(dests);
  /*
  ofstream mfile;
  mfile.open("times.txt");
  for (uint32_t i = 0; i < durations.size(); i++) {
    mfile << durations[i] << "," << insert_num[i] << "\n";
  }
  */
  return true;
}

void rewrite_graph(std::string filename, int num_nodes, bool fix = false) {
  FILE *fr;
  FILE *fw;
  fr = fopen(filename.c_str(), "r");
  printf("starting to read file %s\n", filename.c_str());
  fw = fopen((filename+".shuf").c_str(), "w");
  size_t buf_size = 128;
  char *line = (char *) malloc(buf_size);
  std::vector<uint32_t> new_node_ids(num_nodes);
  for (int i = 0; i < num_nodes; i++) {
    new_node_ids[i] = i;
  }
  std::random_shuffle(new_node_ids.begin(), new_node_ids.end());
  if (fr) {
    while (getline(&line, &buf_size, fr) != -1) {
      if (line[0] == '#') continue;
      uint32_t elem_1;
      uint32_t elem_2;
      sscanf(line, "%u   %u", &elem_1, &elem_2);
      fprintf(fw, "%u   %u\n", new_node_ids[elem_1-(1-fix)], new_node_ids[elem_2-(1-fix)]);
    }
    // return 0;
  } else {
    printf("file was not opened\n");
  }
  fclose(fr);
  fclose(fw);
  printf("finished writing %s\n", (filename+".shuf").c_str());
}


int main(int argc, char *argv[]) {
  //write_rmat_to_file(1<<24, 1000000000);
  //rewrite_graph("graphs/rand-out-soc-LiveJournal1.txt", 4847571);
  //rewrite_graph("/efs/home/____/aspen/code/inputs/com-orkut.ungraph.txt", 3072627);
  //rewrite_graph("rmat_sym_rand.el", 8388608, true);
  //rewrite_graph("rmat_ligra.el", 8388608, true);
  //rewrite_graph("rmat_ligra.el", 8388608, true);
  //parallel_updates(atoi(argv[1]));
  time_ofm_real_world_graphs(atoi(argv[1]), atoi(argv[2]));
  //worst_case(atoi(argv[1]));
  return 0;
  if (argc == 2) {
    if (string(argv[1]) == "-v") {
      if (!verify_pcsr()) {
        printf("pcsr failed\n");
        return 0;
      } else {
        printf("pcsr verified\n");
      }
      if (!verify_adjacency_list()) {
        printf("adj list failed\n");
        return 0;
      } else {
        printf("adjlist verified\n");
      }
      if (!verify_blocked_adjacency_list()) {
        printf("blocked adj list failed\n");
        return 0;
      } else {
        printf("blocked adjlist verified\n");
      }
     if (!verify_adjacency_vector()) {
        printf("adj vector failed\n");
        return 0;
      } else {
        printf("adjvec verified\n");
      }
      if (!verify_adjacency_hashmap()) {
        printf("adj hashmap failed\n");
        return 0;
      } else {
        printf("ahm verified\n");
      }
      if (!verify_csr()) {
        printf("csr failed\n");
        return 0;
      } else {
        printf("csr verified\n");
      }
      return 0;
    }
  }
/*
  uint32_t starting_size = 8;
  assert(starting_size >= batch_size);
  if (argc==3) {
    if(string(argv[1]) == "-s") {
      starting_size = atoi(argv[2]);
    }
  }

  int trials = 1;
  uint32_t num_nodes =  100000;
  uint32_t starting_degree = starting_size;
  uint32_t max_degree = 0+1;
  uint32_t num_edges = num_nodes * starting_degree;



  while(num_edges < max_degree * 100000) {
    for (int i = 16; i >= 1; i/=2) {
      char str[2];
      sprintf(str, "%d", i);
      __cilkrts_end_cilk(); 
      __cilkrts_set_param("nworkers", str);
      int workers = __cilkrts_get_nworkers();
      printf("PCSR: workers: %d: ", workers);
      time_pcsr(num_nodes, num_edges, trials);
      continue;
      printf("Adjacency List: workers: %d: ", workers);
      time_adjlist(num_nodes, num_edges, trials);

      printf("Blocked Adjacency List: workers: %d: ", workers);
      time_blocked_adjlist(num_nodes, num_edges, trials);

      printf("Adjacency Hash Map: workers: %d: ", workers);
      time_adjhashmap(num_nodes, num_edges, trials);
 

      printf("Adjacency Vector: workers: %d: ", workers);
      time_adjvector(num_nodes, num_edges, trials);

   
      printf("CSR: workers: %d: ", workers);
      time_csr(num_nodes, num_edges, trials);



    }


    num_edges *= 2;
  }
*/
  //real_world_graphs
  OFM *ofm = new OFM(0);
      //time_structure_real_world_graphs(ofm, 2);
  //return 0;
  AdjacencyList *list = new AdjacencyList(0);
  BlockedAdjacencyList *blist = new BlockedAdjacencyList(0);
  AdjacencyHashMap *ahm = new AdjacencyHashMap(0);
  AdjacencyVector *av = new AdjacencyVector(0);
  CSR *csr = new CSR(0);
  printf("\n");
  for (int i = 2; i < 3; i++) {
    for (int j = 32; j >= 1; j/=2) {
      char str[2];
      sprintf(str, "%d", j);
      __cilkrts_end_cilk();
      __cilkrts_set_param("nworkers", str);
          
      int workers = __cilkrts_get_nworkers();
      printf("PCSR: workers: %d: ", workers);
      time_structure_real_world_graphs(ofm, i);
      continue;
      printf("Adjacency List: workers: %d: ", workers);
      time_structure_real_world_graphs(list, i);
      printf("Blocked Adjacency List: workers: %d: ", workers);
      time_structure_real_world_graphs(blist, i);
      printf("Adjacency Hash Map: workers: %d: ", workers);
      time_structure_real_world_graphs(ahm, i);
      printf("Adjacency Vector: workers: %d: ", workers);
      time_structure_real_world_graphs(av, i);
      printf("CSR: workers: %d: ", workers);
      time_structure_real_world_graphs(csr, i);
      printf("\n");
    }
  }

}
