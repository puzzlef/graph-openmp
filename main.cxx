#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>
#include <string>
#include <algorithm>
#include <omp.h>
#include "inc/main.hxx"

using namespace std;




#pragma region CONFIGURATION
#ifndef TYPE
/** Type of edge weights. */
#define TYPE float
#endif
#ifndef MAX_THREADS
/** Maximum number of threads to use. */
#define MAX_THREADS 64
#endif
#ifndef REPEAT_BATCH
/** Number of times to repeat each batch. */
#define REPEAT_BATCH 5
#endif
#ifndef REPEAT_METHOD
/** Number of times to repeat each method. */
#define REPEAT_METHOD 5
#endif
#ifndef VISIT_STEPS
/** Number of steps to visit vertices. */
#define VISIT_STEPS 42
#endif
#pragma endregion




#pragma region METHODS
#pragma region PERFORM EXPERIMENT
template <class G, class GD>
inline void testSubractGraphInplace(const G& x, const GD& ydel) {
  using K = typename G::key_type;
  using V = typename G::vertex_value_type;
  using E = typename G::edge_value_type;
  ArenaDiGraph<K, V, E> y;
  y.setAllocator(x.allocator());
  printf("Applying edge deletions in-place ...\n");
  float t0 = measureDuration([&]() {
    duplicateArenaOmpW(y, x);
  });
  float t1 = measureDuration([&]() {
    subtractGraphOmpU(y, ydel);
  });
  println(y);
  printf("{%09.1fms; %09.1fms duplicate} %s\n", t1, t0, "subtractGraphInplace");
  testVisitCountBfs(y, VISIT_STEPS);
  y.clearOmp();
}


template <class G, class GI>
inline void testAddGraphInplace(const G& x, const GI& yins) {
  using K = typename G::key_type;
  using V = typename G::vertex_value_type;
  using E = typename G::edge_value_type;
  ArenaDiGraph<K, V, E> y;
  y.setAllocator(x.allocator());
  printf("Applying edge insertions in-place ...\n");
  float t0 = measureDuration([&]() {
    duplicateArenaOmpW(y, x);
  });
  float t1 = measureDuration([&]() {
    addGraphOmpU(y, yins);
  });
  println(y);
  printf("{%09.1fms; %09.1fms duplicate} %s\n", t1, t0, "addGraphInplace");
  testVisitCountBfs(y, VISIT_STEPS);
  y.clearOmp();
}


template <class G, class GD>
inline void testSubtractGraphNew(const G& x, const GD& ydel) {
  using K = typename G::key_type;
  using V = typename G::vertex_value_type;
  using E = typename G::edge_value_type;
  ArenaDiGraph<K, V, E> y;
  y.setAllocator(x.allocator());
  printf("Applying edge deletions into a new graph ...\n");
  float t = measureDuration([&]() {
    subtractGraphOmpW(y, x, ydel);
  });
  println(y);
  printf("{%09.1fms; %09.1fms duplicate} %s\n", t, 0.0, "subtractGraphNew");
  testVisitCountBfs(y, VISIT_STEPS);
  y.clearOmp();
}


template <class G, class GI>
inline void testAddGraphNew(const G& x, const GI& yins) {
  using K = typename G::key_type;
  using V = typename G::vertex_value_type;
  using E = typename G::edge_value_type;
  ArenaDiGraph<K, V, E> y;
  y.setAllocator(x.allocator());
  printf("Applying edge insertions into a new graph ...\n");
  float t = measureDuration([&]() {
    addGraphOmpW(y, x, yins);
  });
  println(y);
  printf("{%09.1fms; %09.1fms duplicate} %s\n", t, 0.0, "addGraphNew");
  testVisitCountBfs(y, VISIT_STEPS);
  y.clearOmp();
}


template <class G>
inline void testVisitCountBfs(const G& x, int steps) {
  size_t S = x.span();
  printf("Testing visit count with BFS [%d steps] ...\n", steps);
  vector<size_t> visits0(S, 1);
  vector<size_t> visits1(S, 0);
  auto t0 = timeNow();
  for (int i=0; i<steps; ++i) {
    #pragma omp parallel for schedule(dynamic, 512)
    for (size_t u=0; u<S; ++u) {
      if (!x.hasVertex(u)) continue;
      visits1[u] = 0;
      x.forEachEdgeKey(u, [&](auto v) {
        visits1[u] += visits0[v];
      });
    }
    swap(visits0, visits1);
  }
  auto t1 = timeNow();
  size_t total = 0;
  # pragma omp parallel for reduction(+:total) schedule(dynamic, 2048)
  for (size_t u=0; u<S; ++u)
    total += visits0[u];
  printf("Total visits: %zu\n", total);
  printf("{%09.1fms} %s\n", duration(t0, t1), "visitCountBfs");
}




/**
 * Perform the experiment.
 * @param x input graph
 */
template <class G>
inline void runExperiment(const G& x, int run) {
  using K = typename G::key_type;
  using V = typename G::vertex_value_type;
  using E = typename G::edge_value_type;
  printf("Running experiment %d ...\n", run);
  // Create random number generator.
  random_device dev;
  default_random_engine rnd(dev());
  // Test transpose graph.
  {
    ArenaDiGraph<K, V, E> ytr;
    ytr.setAllocator(x.allocator());
    transposeArenaOmpW(ytr, x);
    ytr.clearOmp();
  }
  // Experiment of various batch fractions.
  for (double frac=1e-7; frac<=1e-1; frac*=10) {
    printf("Batch fraction: %.1e\n", frac);
    // Generate random edge deletions and insertions.
    printf("Generating random edge deletions and insertions ...\n");
    auto deletions  = generateEdgeDeletions (rnd, x, size_t(frac * x.size()), 0, x.span(), false);
    auto insertions = generateEdgeInsertions(rnd, x, size_t(frac * x.size()), 0, x.span(), false, E(1));
    printf("Edge deletions:  %zu\n", deletions.size());
    printf("Edge insertions: %zu\n", insertions.size());
    // Create edge deletions and insertions graph.
    printf("Creating edge deletions and insertions graph ...\n");
    ArenaDiGraph<K, V, E> ydel;
    ArenaDiGraph<K, V, E> yins;
    ydel.setAllocator(x.allocator());
    yins.setAllocator(x.allocator());
    ydel.reserveOmp(x.span());
    yins.reserveOmp(x.span());
    for (auto [u, v, w] : deletions)
      ydel.addEdgeUnchecked(u, v, w);
    for (auto [u, v, w] : insertions)
      yins.addEdgeUnchecked(u, v, w);
    printf("Updating edge counts ...\n");
    ydel.updateOmp();
    yins.updateOmp();
    printf("Edge deletions:  "); println(ydel);
    printf("Edge insertions: "); println(yins);
    // Appy edge deletions/insertions in-place.
    testSubractGraphInplace(x, ydel);
    testAddGraphInplace(x, yins);
    // Apply edge deletions/insertions to a new graph.
    testSubtractGraphNew(x, ydel);
    testAddGraphNew(x, yins);
    // Clear memory.
    ydel.clearOmp();
    yins.clearOmp();
    printf("\n");
  }
}


/**
 * Main function.
 * @param argc argument count
 * @param argv argument values
 * @returns exit code
 */
int main(int argc, char **argv) {
  using K = uint32_t;
  using V = None;
  using E = TYPE;
  char *file = argv[1];
  bool symmetric = argc>2? atoi(argv[2]) : false;
  bool weighted  = argc>3? atoi(argv[3]) : false;
  omp_set_num_threads(MAX_THREADS);
  LOG("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  LOG("Loading graph %s ...\n", file);
  // Read graph as CSR.
  {
    DiGraphCsr<K, V, E> xc;
    MappedFile mf(file);
    size_t size = mf.size();
    string_view data((const char*) mf.data(), mf.size());
    double tr = measureDuration([&]() {
      if (weighted) readMtxFormatToCsrOmpW<true> (xc, data);
      else          readMtxFormatToCsrOmpW<false>(xc, data);
    });
    LOG(""); println(xc);
    printf("{%09.1fms} %s\n", tr, "readMtxFormatToCsrOmpW");

    // DiGraph<K, V, E> x;
    // double ts = measureDuration([&]() {
    //   duplicateOmpW(x, xc);
    // });
    // LOG(""); println(x);
    // printf("{%09.1fms} %s\n", ts, "duplicateOmpW");

    // ArenaDiGraph<K, V, E> xa;
    // double ta = measureDuration([&]() {
    //   duplicateArenaOmpW(xa, xc);
    // });
    // LOG(""); println(xa);
    // printf("{%09.1fms} %s\n", ta, "duplicateArenaOmpW");

    // runExperiment(xa, 1);
    // runExperiment(xa, 2);  // Try again.
  }
  printf("\n");
  return 0;
}
#pragma endregion
#pragma endregion
