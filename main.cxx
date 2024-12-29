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
int mainTestGraph(int argc, char **argv) {
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
  // Read graph in parallel.
  {
    DiGraph<K, V, E> x;
    double tr = measureDuration([&]() {
      readMtxFormatFileToGraphOmpW(x, file, weighted);
    });
    LOG(""); println(x);
    printf("{%09.1fms} %s\n", tr, "readMtxFormatFileToGraphOmpW");
  }
  printf("\n");
  return 0;
}


int main(int argc, char **argv) {
  omp_set_num_threads(MAX_THREADS);
  LOG("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  constexpr size_t EACH  = 64;
  constexpr size_t COUNT = 1024 * 1024 * 256;
  // Test libc memory allocation/deallocation performance.
  vector<void*> ptrs(COUNT);
  float t0 = measureDuration([&]() {
    #pragma omp parallel for schedule(dynamic, 2048)
    for (int i=0; i<COUNT; ++i)
      ptrs[i] = malloc(EACH);
  });
  printf("Allocated %zu blocks of %zu bytes in %.3fms, using libc\n", COUNT, EACH, t0);
  float t1 = measureDuration([&]() {
    #pragma omp parallel for schedule(dynamic, 2048)
    for (int i=0; i<COUNT; ++i)
      free(ptrs[i]);
  });
  printf("Freed %zu blocks of %zu bytes in %.3fms, using libc\n", COUNT, EACH, t1);
  // Test Arena memory allocation/deallocation performance.
  ptrs.clear();
  auto t2 = timeNow();
  void *pool = mmapAlloc(COUNT * EACH);
  ConcurrentArenaAllocator<EACH> arena(pool, EACH*COUNT);
  #pragma omp parallel for schedule(dynamic, 2048)
  for (int i=0; i<COUNT; ++i)
    ptrs.push_back(arena.allocate());
  auto t3 = timeNow();
  printf("Allocated %zu blocks of %zu bytes in %.3fms, using arena allocator\n", COUNT, EACH, duration(t2, t3));
  #pragma omp parallel for schedule(dynamic, 2048)
  for (int i=0; i<COUNT; ++i)
    arena.free(ptrs[i]);
  auto t4 = timeNow();
  printf("Freed %zu blocks of %zu bytes in %.3fms, using arena allocator\n", COUNT, EACH, duration(t3, t4));
  return 0;
}
#pragma endregion
#pragma endregion
