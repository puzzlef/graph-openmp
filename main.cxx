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
#pragma endregion




#pragma region METHODS
#pragma region PERFORM EXPERIMENT
template <int CHUNK=2048, class G, class GD>
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
    subtractGraphOmpU<CHUNK>(y, ydel);
  });
  println(y);
  printf("{%09.1fms; %09.1fms duplicate} %s%d\n", t1, t0, "subtractGraphInplace", CHUNK);
}


template <int CHUNK=2048, class G, class GI>
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
    addGraphOmpU<CHUNK>(y, yins);
  });
  println(y);
  printf("{%09.1fms; %09.1fms duplicate} %s%d\n", t1, t0, "addGraphInplace", CHUNK);
}


template <int CHUNK=2048, class G, class GD>
inline void testSubtractGraphNew(const G& x, const GD& ydel) {
  using K = typename G::key_type;
  using V = typename G::vertex_value_type;
  using E = typename G::edge_value_type;
  ArenaDiGraph<K, V, E> y;
  y.setAllocator(x.allocator());
  printf("Applying edge deletions into a new graph ...\n");
  float t = measureDuration([&]() {
    subtractGraphOmpW<CHUNK>(y, x, ydel);
  });
  println(y);
  printf("{%09.1fms; %09.1fms duplicate} %s%d\n", t, 0.0, "subtractGraphNew", CHUNK);
}


template <int CHUNK=2048, class G, class GI>
inline void testAddGraphNew(const G& x, const GI& yins) {
  using K = typename G::key_type;
  using V = typename G::vertex_value_type;
  using E = typename G::edge_value_type;
  ArenaDiGraph<K, V, E> y;
  y.setAllocator(x.allocator());
  printf("Applying edge insertions into a new graph ...\n");
  float t = measureDuration([&]() {
    addGraphOmpW<CHUNK>(y, x, yins);
  });
  println(y);
  printf("{%09.1fms; %09.1fms duplicate} %s%d\n", t, 0.0, "addGraphNew", CHUNK);
}


/**
 * Perform the experiment.
 * @param x input graph
 */
template <class G>
inline void runExperiment(const G& x) {
  using K = typename G::key_type;
  using V = typename G::vertex_value_type;
  using E = typename G::edge_value_type;
  printf("Running experiment ...\n");
  // Create random number generator.
  random_device dev;
  default_random_engine rnd(dev());
  // Experiment of various batch fractions.
  for (double frac=1e-1; frac<=1e-1; frac*=10) {
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
    // Appy edge deletions in-place.
    testSubractGraphInplace<256>(x, ydel);
    testSubractGraphInplace<512>(x, ydel);
    testSubractGraphInplace<1024>(x, ydel);
    testSubractGraphInplace<2048>(x, ydel);
    // Apply edge insertions in-place.
    testAddGraphInplace<256>(x, yins);
    testAddGraphInplace<512>(x, yins);
    testAddGraphInplace<1024>(x, yins);
    testAddGraphInplace<2048>(x, yins);
    // Apply edge deletions to a new graph.
    testSubtractGraphNew<256>(x, ydel);
    testSubtractGraphNew<512>(x, ydel);
    testSubtractGraphNew<1024>(x, ydel);
    testSubtractGraphNew<2048>(x, ydel);
    // Apply edge insertions to a new graph.
    testAddGraphNew<256>(x, yins);
    testAddGraphNew<512>(x, yins);
    testAddGraphNew<1024>(x, yins);
    testAddGraphNew<2048>(x, yins);
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

    ArenaDiGraph<K, V, E> xa;
    double ta = measureDuration([&]() {
      duplicateArenaOmpW(xa, xc);
    });
    LOG(""); println(xa);
    printf("{%09.1fms} %s\n", ta, "duplicateArenaOmpW");

    runExperiment(xa);
  }
  printf("\n");
  return 0;
}
#pragma endregion
#pragma endregion
