#pragma once
#include "Map.cpp"
#include <vector>


//much of the sturcure taken from the aspen code
//https://github.com/ldhulipala/aspen/blob/master/code/algorithms/BC.h
// no license is given in this file and it is included in other files so just have the link
typedef uint32_t uintE;


template <typename E, typename EV>
inline E fetch_and_add(E *a, EV b) {
  volatile E newV, oldV;
  do {oldV = *a; newV = oldV + b;}
  while (!atomic_compare_and_swap(a, oldV, newV));
  return oldV;
}

struct BC_F {
  std::vector<double> &Scores;
  std::vector<bool> &Visited;
  BC_F(std::vector<double> &_Scores, std::vector<bool> &_Visited) : Scores(_Scores), Visited(_Visited) {}
  inline bool update(const uintE& s, const uintE& d) {
    double oldV = Scores[d];
    Scores[d] += Scores[s];
    return oldV == 0.0;
  }
  inline bool updateAtomic (const uintE& s, const uintE& d) {
    double to_add = Scores[s];
    double o_val = fetch_and_add(&Scores[d],to_add);
    return (o_val == 0);
  }
  inline bool cond (uintE d) { return Visited[d] == 0; }
};


// Vertex map function to mark visited vertex_subset
struct BC_Vertex_F {
  std::vector<bool>& Visited;
  BC_Vertex_F(std::vector<bool> &_Visited) : Visited(_Visited) {}
  inline bool operator() (uintE i) {
    Visited[i] = 1;
    return 1;
  }
};

// Vertex map function (used on backwards phase) to mark visited vertex_subset
// and add to Dependencies score
struct BC_Back_Vertex_F {
  std::vector<bool>& Visited;
  std::vector<double>& Dependencies;
  std::vector<double>& NumPaths;
  BC_Back_Vertex_F(std::vector<bool> &_Visited, std::vector<double> &_Dependencies, std::vector<double> &_NumPaths) :
    Visited(_Visited), Dependencies(_Dependencies), NumPaths(_NumPaths) {}
  inline bool operator() (uintE i) {
    Visited[i] = 1;
    Dependencies[i] += NumPaths[i];
    return 1;
  }
};

auto BC(OFM& G, const uintE& start, [[maybe_unused]] bool use_dense_forward=false) {
  size_t n = G.get_n();
  std::vector<double> NumPaths(n, 0);   
  std::vector<bool> Visited(n, false);
  Visited[start] = 1; NumPaths[start] = 1.0;
  VertexSubset Frontier = VertexSubset(start, n); //creates initial frontier

  std::vector<VertexSubset> Levels;
  long round = 0;
  while (Frontier.get_n() > 0) {
    //printf("%u, %u\n", round, Frontier.get_n());
    round++;
    Levels.push_back(Frontier);
    edgeMap(G, Frontier, BC_F(NumPaths,Visited));
    Frontier.move_next_to_current();
    vertexMap(Frontier, BC_Vertex_F(Visited), false); // mark visited
  }
  //printf("%f\n", NumPaths[0]); 
  Levels.push_back(Frontier);

  std::vector<double> Dependencies(n, 0);

  parallel_for(uint32_t i = 0; i < n; i++) {
    NumPaths[i] = 1/NumPaths[i];
  }

  parallel_for(uint32_t i = 0; i < n; i++) {
    Visited[i] = 0;
  }

  vertexMap(Levels[round-1], BC_Back_Vertex_F(Visited,Dependencies,NumPaths), false);

  for(long r=round-2;r>=0;r--) {
    edgeMap(G, Levels[r+1], BC_F(Dependencies,Visited), false);
    vertexMap(Levels[r], BC_Back_Vertex_F(Visited,Dependencies,NumPaths), false);
  }

  parallel_for(uint32_t i = 0; i < n; i++) {
    Dependencies[i] = (Dependencies[i]-NumPaths[i])/NumPaths[i];
  }
  return Dependencies;
}

