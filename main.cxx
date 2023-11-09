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
  bool  PAR  = argc>2 ? atoi(argv[2]) : 1;  // 0=serial, 1=parallel
  omp_set_num_threads(MAX_THREADS);
  printf("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  MappedFile mf(file);
  size_t size = mf.size();
  int T = PAR? MAX_THREADS : 1;
  vector<uint32_t*> sources(T);
  vector<uint32_t*> targets(T);
  vector<float*>    weights(T);
  vector<uint32_t*> degrees(T);
  // vector<unique_ptr<size_t>> scounts;
  vector<size_t> scounts;
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
  MappedPtr<size_t>   offsets(rows+1);
  MappedPtr<uint32_t> edgeKeys(edges);
  MappedPtr<float>    edgeValues(edges);
  // vector<size_t>   offsets(rows+1);
  // vector<uint32_t> edgeKeys(edges);
  // vector<float>    edgeValues(edges);
  printf("rows=%zu, cols=%zu, edges=%zu\n", rows, cols, edges);
  double tr = measureDuration([&]() {
    if (PAR) scounts = readEdgelistFormatOmpU(sources.data(), targets.data(), weights.data(), degrees.data(), data, symmetric, weighted);
    else     m = readEdgelistFormatU(sources[0], targets[0], weights[0], degrees[0], data, symmetric, weighted);
    // csrCreateFromEdgelistOmpU((size_t*) offsets.data(), (uint32_t*) edgeKeys.data(), (float*) edgeValues.data(), (uint32_t**) degrees.data(), (const uint32_t**) sources.data(), (const uint32_t**) targets.data(), (const float**) weights.data(), (const size_t*) counts.data(), rows);
  });
  vector<size_t> counts(T);
  for (size_t t=0; t<T; ++t)
    counts[t] = scounts[t];
  size_t mm = 0;
  for (size_t t=0; t<T; ++t)
    mm += counts[t];
  printf("edges=%zu, mm=%zu\n", edges, mm);
  assert(mm==edges);
  printf("{%09.1fms, m=%zu, size=%zu} %s\n", tr, m, edges, PAR? "readEdgesOmp" : "readEdges");
  // csrCreateFromEdgelistW((size_t*) offsets.data(), (uint32_t*) edgeKeys.data(), (float*) edgeValues.data(), degrees[0], sources[0], targets[0], weights[0], rows, m);
  csrCreateFromEdgelistOmpU((size_t*) offsets.data(), (uint32_t*) edgeKeys.data(), (float*) edgeValues.data(), (uint32_t**) degrees.data(), (const uint32_t**) sources.data(), (const uint32_t**) targets.data(), (const float**) weights.data(), (const size_t*) counts.data(), rows);
  printf("Created CSR graph with %zu nodes and %zu edges.\n", rows, m);
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
