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
  }
  // Read graph-ps in parallel.
  {
    // DiGraph<K, V, E> x;
    // MappedFile mf(file);
    // size_t size = mf.size();
    // string_view data((const char*) mf.data(), mf.size());
    // double tr = measureDuration([&]() {
    //   if (weighted) readMtxFormatToGraphPsOmpW<true> (x, data);
    //   else          readMtxFormatToGraphPsOmpW<false>(x, data);
    // });
    // LOG(""); println(x);
    // printf("{%09.1fms} %s\n", tr, "readMtxFormatToGraphFromCsrOmpW");
  }
  // Read graph-ps in parallel.
  {
    // DiGraph<K, V, E> x;
    // MappedFile mf(file);
    // size_t size = mf.size();
    // string_view data((const char*) mf.data(), mf.size());
    // double tr = measureDuration([&]() {
    //   if (weighted) readMtxFormatToGraphSpOmpW<true> (x, data);
    //   else          readMtxFormatToGraphSpOmpW<false>(x, data);
    // });
    // LOG(""); println(x);
    // printf("{%09.1fms} %s\n", tr, "readMtxFormatToGraphFromEdgelistsOmpW");
  }
  // Read graph in parallel.
  {
    // DiGraph<K, V, E> x;
    // double tr = measureDuration([&]() {
    //   readMtxFormatFileToGraphOmpW(x, file, weighted);
    // });
    // LOG(""); println(x);
    // printf("{%09.1fms} %s\n", tr, "readMtxFormatFileToGraphOmpW");
  }
  printf("\n");
  return 0;
}
#pragma endregion
#pragma endregion
