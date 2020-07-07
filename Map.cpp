#pragma once
#include "OFM.cpp"
#include "VertexSubset.cpp"

template<class F, bool output>
struct EDGE_MAP_SPARSE {
  OFM &G;
  VertexSubset &vs;
  F f;
  EDGE_MAP_SPARSE(OFM &G_, VertexSubset &vs_, F f_) : 
    G(G_), vs(vs_), f(f_){}
  inline bool update(uint32_t val) {
    G.map_sparse<F, VertexSubset, output>(f, vs, val);
    return false;
  }
};

template <class F>
void EdgeMapSparse(OFM &G, VertexSubset &vs, F f, bool output) {
  //printf("edge map sparse\n");
  vs.convert_to_sparse();
  if (output) {
    struct EDGE_MAP_SPARSE<F, true> v(G, vs, f);
    vs.map(v);
  } else {
    struct EDGE_MAP_SPARSE<F, false> v(G, vs, f);
    vs.map(v);
  }
}

template <class F>
void EdgeMapDense(OFM &G, VertexSubset &vs, F f, bool output) {
  //printf("edge map dense\n");
  vs.convert_to_dense();
  // needs a grainsize of at least 512 
  // so writes to the bitvector storing the next vertex set are going to different cache lines
  if (output) {
    if (vs.all) {
      parallel_for(uint64_t i_ = 0; i_ < G.get_n(); i_+=512) {
        uint64_t end = std::min(i_+512, (uint64_t) G.get_n());
        for (uint64_t i = i_; i < end; i++) {
          if (f.cond(i) == 1) {
            //printf("processing row %lu\n", i);
            G.map_dense<F, VertexSubset, true, true>(f, vs, i);
          }
        }
      }
    } else {
      parallel_for(uint64_t i_ = 0; i_ < G.get_n(); i_+=512) {
        uint64_t end = std::min(i_+512, (uint64_t) G.get_n());
        for (uint64_t i = i_; i < end; i++) {
          if (f.cond(i) == 1) {
            //printf("processing row %lu\n", i);
            G.map_dense<F, VertexSubset, true, false>(f, vs, i);
          }
        }
      }
    }
  } else {
    if (vs.all) {
      parallel_for(uint64_t i_ = 0; i_ < G.get_n(); i_+=512) {
        uint64_t end = std::min(i_+512, (uint64_t) G.get_n());
        for (uint64_t i = i_; i < end; i++) {
          if (f.cond(i) == 1) {
            //printf("processing row %lu\n", i);
            G.map_dense<F, VertexSubset, false, true>(f, vs, i);
          }
        }
      }
    } else {
      parallel_for(uint64_t i_ = 0; i_ < G.get_n(); i_+=512) {
        uint64_t end = std::min(i_+512, (uint64_t) G.get_n());
        for (uint64_t i = i_; i < end; i++) {
          if (f.cond(i) == 1) {
            //printf("processing row %lu\n", i);
            G.map_dense<F, VertexSubset, false, false>(f, vs, i);
          }
        }
      }
    }
  }
}

template <class F>
void edgeMap(OFM &G, VertexSubset &vs, F f, bool output = true, uint32_t threshold = 20) {
  //vs.print();
  //printf("%u, %u, %u\n", G.rows, threshold, vs.get_n());
  if (G.get_n()/threshold <= vs.get_n()) {
    return EdgeMapDense(G, vs, f, output);
  } else {
    return EdgeMapSparse(G, vs, f, output);
  }
}

template<class F>
struct VERTEX_MAP {
  VertexSubset &vs;
  F f;
  bool output;
  VERTEX_MAP(VertexSubset &vs_, F f_, bool output_) : vs(vs_), f(f_), output(output_) {}
  inline bool update(uint32_t val) {
    if (output) {
      if (f(val) == 1) {
        vs.insert(val);
      }
    } else {
      f(val);
    }
    return false;
  }
};

template <class F>
void vertexMap(VertexSubset &vs, F f, bool output = true) {
  struct VERTEX_MAP<F> v(vs, f, output);
  vs.map(v);
}
