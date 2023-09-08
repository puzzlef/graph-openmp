#include <cstdint>
#include <cstdio>
#include <utility>
#include <random>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include "src/main.hxx"

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
  LOG("Loading graph %s ...\n", file);
  float tr = measureDuration([&]() { readMtxOmpW(x, file, weighted); });
  LOG("{%09.1fms} ", tr); print(x);  printf(" readMtx\n");
  float ts = measureDuration([&]() { addSelfLoopsOmpU(x, V(), ft); });
  LOG("{%09.1fms} ", ts); print(x);  printf(" addSelfLoops\n");
  float tt = measureDuration([&]() { transposeOmpW(xt, x); });
  LOG("{%09.1fms} ", tt); print(xt); printf(" transpose\n");
  printf("\n");
  return 0;
}
#pragma endregion
#pragma endregion
