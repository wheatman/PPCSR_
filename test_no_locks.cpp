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
#include "AdjacencyMatrix.cpp"
#include "Graph.hpp"
#include "OFM_no_locks.cpp"
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
        if (num_nodes < 20) {
          matrix.print_graph();
          printf("\n");
          g->print_graph();
        }
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
        if (num_nodes < 20) {
          matrix.print_graph();
          printf("\n");
          g->print_graph();
        }
        
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
      vector<double> r5 = matrix.pagerank(test_vector2);
      vector<double> r6 = g->pagerank(test_vector2);
      if (r5 != r6) {
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

static __inline__ uint64_t rdtsc(void)
{
    unsigned hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return ( (uint64_t)lo)|( ((uint64_t)hi)<<32 );
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
  std::vector<uint64_t> times(0);
  for(int i = 0; i < num_edges; i+=1) {
    uint64_t before = rdtsc();
    g->add_edge(srcs[i], dests[i], vals[i]);
    uint64_t after = rdtsc();
    times.push_back(after-before);
  }

  gettimeofday(&end, NULL);
  free(srcs);
  free(dests);
  free(vals);
  std::sort(times.begin(), times.end());
  duration =
            (end.tv_sec - start.tv_sec) * 1000000u + (end.tv_usec - start.tv_usec);
  if (num_edges > 0) {
    printf("add_edge_update :%f: %lu: %lu: %lu: \n", ((double)duration) / (1000000), times[times.size()*.90], times[times.size()*.99], times[times.size()-1]);
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
  // bfs
  /*
    int start_node = rand_in_range(num_nodes);
    gettimeofday(&start, NULL);
    for (int bfs_iters = 0; bfs_iters < 100; bfs_iters++) {
      vector<uint32_t> r4 = g->bfs(start_node);
      start_node = r4[0] % num_nodes;
    }
    gettimeofday(&end, NULL);
    duration = (end.tv_sec-start.tv_sec)*1000000u+(end.tv_usec-start.tv_usec);
    printf("bfs in %f,", ((double)duration)/1000000);
    printf("sn = %d,", start_node);
    */
/*
  // pagerank
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


int main(int argc, char *argv[]) {
  if (argc == 2) {
    if (string(argv[1]) == "-v") {
      if (!verify_pcsr()) {
        return 0;
      } else {
        printf("pcsr verified\n");
      }

    }
  }

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
  uint32_t max_degree = 1024+1;
  uint32_t num_edges = num_nodes * starting_degree;

  int workers = __cilkrts_get_nworkers();

  while(num_edges < max_degree * 100000) {

    printf("PCSR: workers: %d: ", workers);
    time_pcsr(num_nodes, num_edges, trials);

    num_edges *= 2;
  }
  return 0;
}
