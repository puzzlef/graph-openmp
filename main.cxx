#include <utility>
#include <random>
#include <vector>
#include <string>
#include <cstdio>
#include <iostream>
#include "src/main.hxx"

using namespace std;




int main(int argc, char **argv) {
  using K = uint32_t;
  using V = TYPE;
  char *file = argv[1];
  bool sym   = argc>2? stoi(argv[2]) : false;
  int repeat = argc>3? stoi(argv[3]) : 5;
  OutDiGraph<K, None, V> x, y;  // V w = 1;
  LOG("Loading graph %s ...\n", file);
  omp_set_num_threads(MAX_THREADS);
  LOG("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  // float tx = measureDuration([&]() { readMtxW(x, file); });
  // println(x); printf("[%09.3f ms] readMtxW\n", tx); x.clear();
  float ty = measureDuration([&]() { readMtxOmpW(y, file, true); });
  LOG(""); println(y); LOG("[%09.3f ms] readMtxOmpW\n", ty); y.clear();
  printf("\n");
  return 0;
}
