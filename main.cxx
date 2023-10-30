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
  struct stat sb;
  omp_set_num_threads(MAX_THREADS);
  LOG("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  DiGraph<K, None, None> x, xt;
  auto  ft = [](auto u) { return true; };
  // mmap() the file.
  int   fd = open(file, O_RDONLY);  // O_DIRECT?
  fstat(fd, &sb);
  void *addr = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED | MAP_NORESERVE, fd, 0);  // MAP_PRIVATE?
  madvise(addr, sb.st_size, MADV_WILLNEED);  // MADV_SEQUENTIAL?
  LOG("Loading graph %s ...\n", file);
  string_view data((const char*) addr, sb.st_size);
  float tr = measureDuration([&]() { readMtxFormatOmpW(x, file, weighted); });
  LOG("{%09.1fms} ", tr); print(x);  printf(" readMtx\n");
  float ts = measureDuration([&]() { addSelfLoopsOmpU(x, V(), ft); });
  LOG("{%09.1fms} ", ts); print(x);  printf(" addSelfLoops\n");
  float tt = measureDuration([&]() { transposeOmpW(xt, x); });
  LOG("{%09.1fms} ", tt); print(xt); printf(" transpose\n");
  // Cleanup.
  munmap(addr, sb.st_size);
  close(fd);
  printf("\n");
  return 0;
}
#pragma endregion
#pragma endregion
