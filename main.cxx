#include <cstdint>
#include <cstdio>
#include <utility>
#include <random>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
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
  char *file    = argv[1];
  bool weighted = false;
  omp_set_num_threads(MAX_THREADS);
  LOG("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  DiGraph<K, None, V> x;
  MappedFile fmap(file);
  LOG("Loading graph %s ...\n", file);
  string_view data((const char*) fmap.addr, fmap.size);
  float tr = measureDuration([&]() { readMtxFormatOmpW(x, data, weighted); });
  LOG("{%09.1fms} ", tr); print(x);  printf(" readMtx\n");
  printf("\n");
  return 0;
}
#pragma endregion
#pragma endregion
