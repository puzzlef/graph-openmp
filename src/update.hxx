#pragma once
#include <cstdint>
#include <utility>
#include <vector>
#include <omp.h>

using std::pair;
using std::vector;




// UPDATE
// ------

// Update changes made to a graph.
template <class G>
inline void updateU(G& a) {
  a.update();
}


template <class G>
inline void updateOmpU(G& a) {
  using  K = typename G::key_type;
  using  E = typename G::edge_value_type;
  size_t N = a.order();
  // Create per-thread buffers for update operation.
  int THREADS = omp_get_max_threads();
  vector<vector<pair<K, E>>*> bufs(THREADS);
  for (int i=0; i<THREADS; ++i)
    bufs[i] = new vector<pair<K, E>>();
  // Update edges of each vertex individually.
  #pragma omp parallel for schedule(auto)
  for (K u=0; u<N; ++u) {
    int t = omp_get_thread_num();
    a.update(u, bufs[t]);
  }
  // Update the entire graph, find total edges.
  a.update();
  // Clean up.
  for (int i=0; i<THREADS; ++i)
    delete bufs[i];
}
