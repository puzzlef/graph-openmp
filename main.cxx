#include <cstdint>
#include <cstdio>
#include <utility>
#include <random>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>
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
#define REPEAT_METHOD 1
#endif
#pragma endregion




#pragma region METHODS
#pragma region PERFORM EXPERIMENT
/**
 * Main function.
 * @param argc argument count
 * @param argv argument values
 * @returns zero on success, non-zero on failure
 */
int main(int argc, char **argv) {
  using K = uint32_t;
  using V = None;
  using Edge = tuple<uint32_t, uint32_t, float>;
  char *file    = argv[1];
  bool weighted = false;
  omp_set_num_threads(MAX_THREADS);
  LOG("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  DiGraph<K, None, V> x, xt;
  auto  ft = [](auto u) { return true; };
  // Disabled old code, as i am testing mmap() now.
  #if 0
  LOG("Loading graph %s ...\n", file);
  float tr = measureDuration([&]() { readMtxFormatOmpW(x, file, weighted); });
  LOG("{%09.1fms} ", tr); print(x);  printf(" readMtx\n");
  float ts = measureDuration([&]() { addSelfLoopsOmpU(x, V(), ft); });
  LOG("{%09.1fms} ", ts); print(x);  printf(" addSelfLoops\n");
  float tt = measureDuration([&]() { transposeOmpW(xt, x); });
  LOG("{%09.1fms} ", tt); print(xt); printf(" transpose\n");
  #endif
  // mmap() the file.
  struct stat sb;
  int   fd = open(file, O_RDONLY);  // O_DIRECT?
  fstat(fd, &sb);
  void *addr = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED | MAP_NORESERVE, fd, 0);  // MAP_PRIVATE?
  madvise(addr, sb.st_size, MADV_WILLNEED);  // MADV_SEQUENTIAL?
  // Read the header.
  LOG("Loading graph %s ...\n", file);
  string_view data((const char*) addr, sb.st_size);
  // bool symmetric; size_t rows, cols, size;
  // readMtxFormatHeaderU(symmetric, rows, cols, size, data);
  // printf("order=%zu, size=%zu\n", rows, size);
  #if 0
  // Allocate memory for the graph.
  int T = omp_get_max_threads();
  vector<vector<Edge>*> edges(T);
  for (int t=0; t<T; ++t) {
    edges[t] = new vector<Edge>();
    edges[t]->reserve(size);  // Over-allocate to avoid reallocation
  }
  // Read the graph.
  auto fb = [&](auto u, auto v, auto w) {
    int t = omp_get_thread_num();
    edges[t]->push_back({uint32_t(u), uint32_t(v), float(w)});
  };
  symmetric = false;
  float tm = measureDuration([&]() {
    readEdgelistFormatSeparateDoOmpU(data, symmetric, weighted, fb);
  });
  size_t readSize = 0;
  for (int t=0; t<T; ++t)
    readSize += edges[t]->size();
  printf("Read %zu edges\n", readSize);
  LOG("{%09.1fms} readMtxMmap\n", tm);
  #endif
  float tm = measureDuration([&]() {
    auto fv = [](auto u, auto d)         { return true; };
    auto fe = [](auto u, auto v, auto w) { return true; };
    auto fh = [&](auto symmetric, auto rows, auto cols, auto size) {
      addVerticesIfU(x, K(1), K(max(rows, cols)+1), V(), fv);
    };
    auto fb = [&](auto u, auto v, auto w) { if (fe(K(u), K(v), K(w))) addEdgeOmpU(x, K(u), K(v), V(w)); };
    readMtxFormatDoOmp(data, weighted, fh, fb);
  });
  LOG("{%09.1fms} readMtxMmap\n", tm);
  updateOmpU(x);
  println(x);
  // Cleanup.
  // for (int t=0; t<T; ++t)
  //   delete edges[t];
  munmap(addr, sb.st_size);
  close(fd);
  printf("\n");
  return 0;
}
#pragma endregion
#pragma endregion
