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
  omp_set_num_threads(MAX_THREADS);
  LOG("OMP_NUM_THREADS=%d\n", MAX_THREADS);
  OutDiGraph<K> x;
  auto fl = [](auto u) { return true; };
  LOG("Loading graph %s ...\n", file);
  readMtxW(x, file);              LOG(""); println(x);
  x = selfLoopOmp(x, None(), fl); LOG(""); println(x);
  x.forEachVertexKey([&](auto u) {
    K vd = K();
    if (x.degree(u)!=2) return;
    x.forEachEdgeKey(u, [&](auto v) { if (v!=u) vd = v; });
    x.removeEdge(u, vd);
  });
  updateOmpU(x);
  x.forEachVertexKey([&](auto u) {
    if (x.degree(u)==0) LOG("x.degree(%d)==0\n", u);
  });
  printf("\n");
  return 0;
}
