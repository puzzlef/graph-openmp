#pragma once
#include "Graph.hxx"
#include "update.hxx"




#pragma region METHODS
#pragma region TRANSPOSE
/**
 * Transpose a graph.
 * @param a transposed graph (output)
 * @param x graph to transpose
 */
template <class H, class G>
inline void transposeW(H& a, const G& x) {
  a.reserve(x.span());
  x.forEachVertex([&](auto u, auto d) { a.addVertex(u, d); });
  x.forEachVertexKey([&](auto u) {
    x.forEachEdge(u, [&](auto v, auto w) { a.addEdge(v, u, w); });
  });
  a.update();
}

/**
 * Transpose a graph.
 * @param x graph to transpose
 * @returns transposed graph
 */
template <class G>
inline auto transpose(const G& x) {
  G a; transposeW(a, x);
  return a;
}


#ifdef OPENMP
/**
 * Transpose a graph in parallel.
 * @param a transposed graph (output)
 * @param x graph to transpose
 */
template <class H, class G>
inline void transposeOmpW(H& a, const G& x) {
  a.reserve(x.span());
  x.forEachVertex([&](auto u, auto d) { a.addVertex(u, d); });
  #pragma omp parallel
  {
    x.forEachVertexKey([&](auto u) {
      x.forEachEdge(u, [&](auto v, auto w) { addEdgeOmpU(a, v, u, w); });
    });
  }
  updateOmpU(a);
}

/**
 * Transpose a graph in parallel.
 * @param x graph to transpose
 * @returns transposed graph
 */
template <class G>
inline auto transposeOmp(const G& x) {
  G a; transposeOmpW(a, x);
  return a;
}
#endif
#pragma endregion




#pragma region TRANSPOSE WITH DEGREE
/**
 * Transpose a graph with degree.
 * @param a transposed graph with degree (output)
 * @param x graph to transpose
 */
template <class H, class G>
inline void transposeWithDegreeW(H& a, const G& x) {
  a.reserve(x.span());
  x.forEachVertexKey([&](auto u) { a.addVertex(u, x.degree(u)); });
  x.forEachVertexKey([&](auto u) {
    x.forEachEdge(u, [&](auto v, auto w) { a.addEdge(v, u, w); });
  });
  a.update();
}

/**
 * Transpose a graph with degree.
 * @param x graph to transpose
 * @returns transposed graph with degree
 */
template <class G>
inline auto transposeWithDegree(const G& x) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  using H = DiGraph<K, K, E>;
  H a; transposeWithDegreeW(a, x);
  return a;
}


#ifdef OPENMP
/**
 * Transpose a graph with degree in parallel.
 * @param a transposed graph with degree (output)
 * @param x graph to transpose
 */
template <class H, class G>
inline void transposeWithDegreeOmpW(H& a, const G& x) {
  a.reserve(x.span());
  x.forEachVertexKey([&](auto u) { a.addVertex(u, x.degree(u)); });
  #pragma omp parallel
  {
    x.forEachVertexKey([&](auto u) {
      x.forEachEdge(u, [&](auto v, auto w) { addEdgeOmpU(a, v, u, w); });
    });
  }
  updateOmpU(a);
}

/**
 * Transpose a graph with degree in parallel.
 * @param x graph to transpose
 * @returns transposed graph with degree
 */
template <class G>
inline auto transposeWithDegreeOmp(const G& x) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  using H = DiGraph<K, K, E>;
  H a; transposeWithDegreeOmpW(a, x);
  return a;
}
#endif
#pragma endregion



#pragma region ARENA-BASED DIGRAPH
template <int PARTITIONS=8, int CHUNK_SIZE=32, class H, class G>
inline void transposeArenaOmpW(H &a, const G& x) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  using O = size_t;
  printf("Transposing graph ...\n");
  size_t S = x.span();
  size_t M = x.size();
  int    T = omp_get_max_threads();
  vector<size_t> buf(T);
  K *degrees[PARTITIONS];
  O *offsets[PARTITIONS];
  K *edgeKeys[PARTITIONS];
  E *edgeValues[PARTITIONS];
  for (int p=0; p<PARTITIONS; ++p) {
    degrees[p]    = new K[S+1];
    offsets[p]    = new O[S+1];
    edgeKeys[p]   = new K[M];
    edgeValues[p] = sizeof(E)? new E[M] : nullptr;
    fillValueOmpU(degrees[p], S+1, K());
  }
  auto t0 = timeNow();
  // Find the per-partition vertex degrees (transposed).
  #pragma omp parallel for schedule(static, CHUNK_SIZE)
  for (K u=0; u<S; ++u) {
    int t = omp_get_thread_num();
    int p = t % PARTITIONS;
    if (!x.hasVertex(u)) continue;
    x.forEachEdgeKey(u, [&](auto v) {
      #pragma omp atomic update
      ++degrees[p][v];
    });
  }
  auto fdeg = [&](auto v) {
    if (PARTITIONS==1) return degrees[0][v];
    else if (PARTITIONS==2) return degrees[0][v] + degrees[1][v];
    else if (PARTITIONS==4) return degrees[0][v] + degrees[1][v] + degrees[2][v] + degrees[3][v];
    else if (PARTITIONS==8) return degrees[0][v] + degrees[1][v] + degrees[2][v] + degrees[3][v] + degrees[4][v] + degrees[5][v] + degrees[6][v] + degrees[7][v];
    else {
      K d = K();
      for (int p=0; p<PARTITIONS; ++p)
        d += degrees[p][v];
      return d;
    }
  };
  // auto t1 = timeNow();
  // Compute per-partition shifted offsets.
  for (int p=0; p<PARTITIONS; ++p) {
    offsets[p][0] = O();
    exclusiveScanOmpW(offsets[p]+1, buf.data(), degrees[p], S);
  }
  // auto t2 = timeNow();
  // Populate per-partition CSR.
  #pragma omp parallel for schedule(static, CHUNK_SIZE)
  for (K u=0; u<S; ++u) {
    int t = omp_get_thread_num();
    int p = t % PARTITIONS;
    if (!x.hasVertex(u)) continue;
    x.forEachEdge(u, [&](auto v, auto w) {
      size_t j = 0;
      #pragma omp atomic capture
      j = offsets[p][v+1]++;
      edgeKeys[p][j] = u;
      if (sizeof(E)) edgeValues[p][j] = w;
      ++j;
    });
  }
  // auto t3 = timeNow();
  // Prepare the output graph.
  a.clearOmp();
  a.reserveOmp(S);
  // auto t4 = timeNow();
  // Add vertices.
  #pragma omp parallel for schedule(static, 2048)
  for (K v=0; v<S; ++v) {
    if (!x.hasVertex(v)) continue;
    a.addVertex(v);
  }
  // auto t5 = timeNow();
  // Reserve space for edges.
  #pragma omp parallel for schedule(dynamic, 1024)
  for (size_t v=0; v<S; ++v) {
    K d = fdeg(v);
    if (d) a.allocateEdges(v, d);
  }
  // auto t6 = timeNow();
  // Populate edges.
  #pragma omp parallel for schedule(dynamic, 1024)
  for (K v=0; v<S; ++v) {
    if (!x.hasVertex(v)) continue;
    for (int p=0; p<PARTITIONS; ++p) {
      size_t i = offsets[p][v];
      size_t I = offsets[p][v+1];
      for (; i<I; ++i)
        a.addEdgeUnsafe(v, edgeKeys[p][i], edgeValues[p][i]);
    }
  }
  // auto t7 = timeNow();
  // Update the graph.
  a.updateOmp(true, false);
  auto t8 = timeNow();
  // Free space.
  for (int i=0; i<PARTITIONS; ++i) {
    delete degrees[i];
    delete offsets[i];
    delete edgeKeys[i];
    delete edgeValues[i];
  }
  printf("{%09.1fms} %s\n", duration(t0, t8), "transposeArenaOmp");
  // printf("{%09.1fms} Measure degrees\n", duration(t0, t1));
  // printf("{%09.1fms} Compute offsets\n", duration(t1, t2));
  // printf("{%09.1fms} Populate CSR\n", duration(t2, t3));
  // printf("{%09.1fms} Prepare graph\n", duration(t3, t4));
  // printf("{%09.1fms} Add vertices\n", duration(t4, t5));
  // printf("{%09.1fms} Allocate edges\n", duration(t5, t6));
  // printf("{%09.1fms} Populate edges\n", duration(t6, t7));
  // printf("{%09.1fms} Update\n", duration(t7, t8));
}
#pragma endregion
#pragma endregion
