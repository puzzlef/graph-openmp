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
      bb = scanNumberThrowerW<false>(edges[t][i++], bb, be, fu, fw);
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
int main(int argc, char **argv) {
  char  *file  = argv[1];
  bool   PAR   = argc>2 ? atoi(argv[2]) : 1;     // 0=serial, 1=parallel
  omp_set_num_threads(MAX_THREADS);
  printf("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  // Map file to memory.
  auto [fd, addr, size] = mapFileToMemory(file);
  // Allocate memory for storing converted numbers (overallocate 64K pages).
  int  T = PAR? MAX_THREADS : 1;
  vector<uint32_t*> numbers(T);
  for (size_t t=0; t<T; ++t)
    numbers[t] = (uint32_t*) allocateMemoryMmap(sizeof(uint32_t) * size / 4);
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
    freeMemoryMmap(numbers[t], sizeof(uint32_t) * size / 4);
  unmapFileFromMemory(fd, addr, size);
  printf("\n");
  return 0;
}


int mainNew(int argc, char **argv) {
  char *file = argv[1];
  bool  PAR  = argc>2 ? atoi(argv[2]) : 1;  // 0=serial, 1=parallel
  omp_set_num_threads(MAX_THREADS);
  printf("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  auto [fd, addr, size] = mapFileToMemory(file);
  int T = PAR? MAX_THREADS : 1;
  vector<uint32_t*> sources(T);
  vector<uint32_t*> targets(T);
  vector<size_t*>   indices(T);
  for (size_t t=0; t<T; ++t) {
    sources[t] = (uint32_t*) allocateMemoryMmap(sizeof(uint32_t) * size / 4);
    targets[t] = (uint32_t*) allocateMemoryMmap(sizeof(uint32_t) * size / 4);
    indices[t] = new size_t[1];
  }
  printf("Reading edgelist in file %s ...\n", file);
  string_view data((const char*) addr, size);
  bool symmetric, weighted = false;
  size_t rows, cols;
  auto fbs = [&](auto u, auto v, auto w) {
    size_t i = indices[0][0]++;
    sources[0][i] = u;
    targets[0][i] = v;
  };
  auto fbp = [&](auto u, auto v, auto w) {
    int t = omp_get_thread_num();
    size_t i = indices[t][0]++;
    sources[t][i] = u;
    targets[t][i] = v;
  };
  readMtxFormatHeaderU(symmetric, rows, cols, size, data);
  double tr = measureDuration([&]() {
    if (PAR) readEdgelistFormatDoOmpU(data, symmetric, weighted, fbs);
    else     readEdgelistFormatDoU   (data, symmetric, weighted, fbp);
  });
  printf("{%09.1fms} %s\n", tr, PAR? "readNumbersOmp" : "readNumbers");
  // Free memory.
  for (size_t t=0; t<T; ++t) {
    freeMemoryMmap(sources[t], sizeof(uint32_t) * size / 4);
    freeMemoryMmap(targets[t], sizeof(uint32_t) * size / 4);
    delete indices[t];
  }
  unmapFileFromMemory(fd, addr, size);
  printf("\n");
  return 0;
}
#pragma endregion
#pragma endregion
