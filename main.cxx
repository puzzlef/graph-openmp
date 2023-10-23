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
  using V = TYPE;
  char *file    = argv[1];
  bool weighted = false;
  omp_set_num_threads(MAX_THREADS);
  LOG("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  DiGraph<K, None, V> x, xt;
  auto  ft = [](auto u) { return true; };
  // LOG("Loading graph %s ...\n", file);
  // float tr = measureDuration([&]() { readMtxFormatOmpW(x, file, weighted); });
  // LOG("{%09.1fms} ", tr); print(x);  printf(" readMtx\n");
  // float ts = measureDuration([&]() { addSelfLoopsOmpU(x, V(), ft); });
  // LOG("{%09.1fms} ", ts); print(x);  printf(" addSelfLoops\n");
  // float tt = measureDuration([&]() { transposeOmpW(xt, x); });
  // LOG("{%09.1fms} ", tt); print(xt); printf(" transpose\n");
  // // Test with mmap().
  int fd = open(file, O_RDONLY);  // O_DIRECT?
  if (fd == -1) { perror("open"); exit(1); }
  struct stat sb;
  if (fstat(fd, &sb) == -1) { perror("fstat"); exit(1); }
  void *addr = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED | MAP_NORESERVE, fd, 0);
  if (addr == MAP_FAILED) { perror("mmap"); exit(1); }
  if (madvise(addr, sb.st_size, MADV_WILLNEED) == -1) { perror("madvise"); exit(1); }
  LOG("Loading graph %s ...\n", file);
  string_view data((const char*) addr, sb.st_size);
  string line = string(data.substr(0, data.find('\n')));
  LOG("mmap() success: %s\n", line.c_str());
  float tm = measureDuration([&]() {
    bool symmetric; size_t rows, cols, size;
    readMtxFormatHeaderU(symmetric, rows, cols, size, data);
    vector<tuple<uint32_t, uint32_t, float>> edges;
    edges.reserve(size);
    auto fb = [&](auto u, auto v, auto w) { edges.push_back({uint32_t(u), uint32_t(v), float(w)}); };
    readEdgelistFormatDoU(data, symmetric, weighted, fb);
  });
  LOG("{%09.1fms} ", tm); print(x);  printf(" readMtxMmap\n");
  if (munmap(addr, sb.st_size) == -1) { perror("munmap"); exit(1); }
  if (close(fd) == -1) { perror("close"); exit(1); }
  printf("\n");
  return 0;
}
#pragma endregion
#pragma endregion
