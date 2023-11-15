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
// Read floats from string.
inline size_t readFloats(uint32_t *edges, const uint8_t *data, size_t N) {
  size_t i = 0;
  const uint8_t *ib = data;
  const uint8_t *ie = data + N;
  while (ib<ie) {
    ib = findNextDigit(ib, ie);
    if (ib<ie) ib = parseNumberW(edges[i++], ib, ie);
  }
  return i;
}

// Read floats from string, using OpenMP.
template <size_t BLOCK=256*1024>
inline size_t readFloatsOmp(uint32_t **edges, const uint8_t *data, size_t N) {
  size_t i = 0;
  const uint8_t *ib = data;
  const uint8_t *ie = data + N;
  auto fu = [](char c) { return false; };
  auto fw = [](char c) { return false; };
  #pragma omp parallel for schedule(dynamic, 1) reduction(+:i)
  for (size_t b=0; b<N; b+=BLOCK) {
    int t = omp_get_thread_num();
    const uint8_t *bb = ib + b;
    const uint8_t *be = min(bb + BLOCK, ie);
    if (bb!=ib && isDigit(*(bb-1))) bb = findNextNonDigit(bb, ie);
    if (be!=ie && isDigit(*(be-1))) be = findNextNonDigit(be, ie);
    while (bb<be) {
      bb = findNextDigit(bb, be);
      // if (bb<be) { bb = findNextNonDigit(bb, be); i++; }
      bb = parseNumberSimdW(edges[t][i++], bb, be);
    }
  }
  return i;
}


/**
 * Main function.
 * @param argc argument count
 * @param argv argument values
 * @returns zero on success, non-zero on failure
 */
int mainOld(int argc, char **argv) {
  char  *file  = argv[1];
  bool   PAR   = argc>2 ? atoi(argv[2]) : 1;     // 0=serial, 1=parallel
  omp_set_num_threads(MAX_THREADS);
  printf("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  // Map file to memory.
  auto [fd, addr, size] = mmapOpenFile(file);
  // Allocate memory for storing converted numbers (overallocate 64K pages).
  int  T = PAR? MAX_THREADS : 1;
  vector<uint32_t*> numbers(T);
  for (size_t t=0; t<T; ++t)
    numbers[t] = (uint32_t*) mmapAlloc(sizeof(uint32_t) * size / 4);
  printf("Counting and converting numbers in file %s ...\n", file);
  size_t n = 0;
  double tr = measureDuration([&]() {
    if (PAR) n = readFloatsOmp(numbers.data(), (uint8_t*) addr, size);
    else     n = readFloats   (numbers[0],     (uint8_t*) addr, size);
  });
  if (n<100) {
    for (size_t t=0; t<T; ++t) {
      for (size_t i=0; i<n; ++i)
        printf("%e\n", numbers[t][i]);
    }
  }
  printf("{%09.1fms, n=%zu} %s\n", tr, n, PAR? "readNumbersOmp" : "readNumbers");
  // Free memory.
  for (size_t t=0; t<T; ++t)
    mmapFree(numbers[t], sizeof(uint32_t) * size / 4);
  mmapCloseFile(fd, addr, size);
  printf("\n");
  return 0;
}


int main(int argc, char **argv) {
  char *file = argv[1];
  omp_set_num_threads(MAX_THREADS);
  printf("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  MappedFile mf(file);
  size_t size = mf.size();
  int T = MAX_THREADS;
  vector<uint32_t*> sources(T);
  vector<uint32_t*> targets(T);
  vector<float*>    weights(T);
  vector<uint32_t*> degrees(T);
  vector<size_t>    counts;
  for (size_t t=0; t<T; ++t) {
    sources[t] = (uint32_t*) mmapAlloc(sizeof(uint32_t) * size / 4);
    targets[t] = (uint32_t*) mmapAlloc(sizeof(uint32_t) * size / 4);
    weights[t] = (float*)    mmapAlloc(sizeof(float)    * size / 4);
    degrees[t] = (uint32_t*) mmapAlloc(sizeof(uint32_t) * size / 4);
  }
  printf("Reading edgelist in file %s ...\n", file);
  string_view data((const char*) mf.data(), mf.size());
  bool symmetric, weighted = false;
  size_t rows, cols, edges, m = 0;
  size_t head = readMtxFormatHeaderW(symmetric, rows, cols, edges, data);
  data.remove_prefix(head);
  vector<size_t*>   offsets(T);
  vector<uint32_t*> edgeKeys(T);
  vector<float*>    edgeValues(T);
  for (int t=0; t<T; ++t) {
    offsets[t]    = (size_t*)   mmapAlloc(sizeof(size_t)   * (rows+1));
    edgeKeys[t]   = (uint32_t*) mmapAlloc(sizeof(uint32_t) * edges);
    edgeValues[t] = (float*)    mmapAlloc(sizeof(float)    * edges);
  }
  printf("rows=%zu, cols=%zu, edges=%zu\n", rows, cols, edges);
  double tr = measureDuration([&]() {
    counts = readEdgelistFormatOmpU<false, 4>(sources.data(), targets.data(), weights.data(), degrees.data(), data, symmetric, weighted);
    convertToCsrOmpU(offsets.data(), degrees.data(), edgeKeys.data(), edgeValues.data(), sources.data(), targets.data(), (float**) nullptr, counts, rows);
  });
  asm("vzeroupper");  // Avoid AVX-SSE transition penalty
  size_t mm = 0;
  for (size_t t=0; t<T; ++t)
    mm += counts[t];
  assert(mm==edges);
  printf("{%09.1fms, size=%zu} %s\n", tr, edges, "readEdgesOmp");
  // Free memory.
  for (size_t t=0; t<T; ++t) {
    mmapFree(sources[t], sizeof(uint32_t) * size / 4);
    mmapFree(targets[t], sizeof(uint32_t) * size / 4);
    mmapFree(weights[t], sizeof(float)    * size / 4);
    mmapFree(degrees[t], sizeof(uint32_t) * size / 4);
  }
  printf("\n");
  return 0;
}
#pragma endregion
#pragma endregion
