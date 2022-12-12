#pragma once
#include <cstdlib>
#include <utility>
#include <algorithm>
#include <string>
#include <istream>
#include <sstream>
#include <fstream>
#include <omp.h>
#include "_main.hxx"
#include "Graph.hxx"
#include "update.hxx"

using std::tuple;
using std::string;
using std::istream;
using std::istringstream;
using std::ifstream;
using std::ofstream;
using std::move;
using std::max;
using std::getline;




// READ MTX
// --------

inline size_t readMtxHeader(istream& s, bool& sym, size_t& rows, size_t& cols, size_t& size) {
  string line, h0, h1, h2, h3, h4;
  // Skip past the comments and read the graph type.
  while (true) {
    getline(s, line);
    if (line.find('%')!=0) break;
    if (line.find("%%")!=0) continue;
    istringstream sline(line);
    sline >> h0 >> h1 >> h2 >> h3 >> h4;
  }
  if (h1!="matrix" || h2!="coordinate") return 0;
  sym = h4=="symmetric" || h4=="skew-symmetric";
  // Read rows, cols, size.
  istringstream sline(line);
  sline >> rows >> cols >> size;
  return max(rows, cols);
}


template <class FV, class FE>
inline void readMtxDo(istream& s, FV fv, FE fe) {
  bool sym; size_t rows, cols, size;
  size_t n = readMtxHeader(s, sym, rows, cols, size);
  if (n==0) return;
  // Add all vertices first.
  // Prevent unnecessary respan first.
  fv(n);
  for (size_t u=1; u<=n; ++u)
    fv(u);  // a.addVertex(u);
  // Then we add the edges.
  string line;
  while (getline(s, line)) {
    size_t u, v; double w = 1;
    istringstream sline(line);
    if (!(sline >> u >> v)) break;
    sline >> w;
    fe(u, v, w);           // a.addEdge(u, v);
    if (sym) fe(v, u, w);  // a.addEdge(v, u);
  }
}


template <class G>
inline void readMtxW(G& a, istream& s) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  auto fv = [&](auto u) { a.addVertex(K(u)); };
  auto fe = [&](auto u, auto v, auto w) { a.addEdge(K(u), K(v), E(w)); };
  PERFORMI( auto t0 = timeNow() );
  readMtxDo(s, fv, fe);
  PERFORMI( auto t1 = timeNow() );
  a.update();
  PERFORMI( auto t2 = timeNow() );
  PRINTFI("readMtxOmpW(): read=%.1fms, update=%.1fms\n", duration(t0, t1), duration(t1, t2));
}
template <class G>
inline void readMtxW(G& a, const char *pth) {
  ifstream f(pth);
  readMtxW(a, f);
}




// READ MTX (OPENMP)
// -----------------

template <class FV, class FE>
inline void readMtxDoOmp(istream& s, FV fv, FE fe) {
  bool sym; size_t rows, cols, size;
  size_t n = readMtxHeader(s, sym, rows, cols, size);
  if (n==0) return;
  // Add all vertices first.
  // Prevent unnecessary respan first.
  fv(n);
  for (size_t u=1; u<=n; ++u)
    fv(u);  // a.addVertex(u);
  // Then we add the edges.
  const int THREADS = omp_get_max_threads();
  const int LINES   = 131072;
  const size_t CHUNK_SIZE = 1024;
  vector<string> lines(LINES);
  vector<tuple<size_t, size_t, double>> edges(LINES);
  PERFORMI( float dread  = 0 );
  PERFORMI( float dparse = 0 );
  PERFORMI( float dadd   = 0 );
  while (true) {
    PERFORMI( auto t0 = timeNow() );
    // Read several lines from the stream.
    int READ = 0;
    for (int i=0; i<LINES; ++i, ++READ)
      if (!getline(s, lines[i])) break;
    if (READ==0) break;
    PERFORMI( auto t1 = timeNow() );
    // Parse lines using multiple threads.
    #pragma omp parallel for schedule(dynamic, 1024)
    for (int i=0; i<READ; ++i) {
      char *line = (char*) lines[i].c_str();
      size_t u = strtoull(line, &line, 10);
      size_t v = strtoull(line, &line, 10);
      double w = strtod  (line, &line);
      edges[i] = {u, v, w? w : 1};
    }
    PERFORMI( auto t2 = timeNow() );
    // Add edges to the graph.
    #pragma omp parallel
    {
      int T = omp_get_num_threads();
      int t = omp_get_thread_num();
      for (int i=0; i<READ; ++i) {
        const auto& [u, v, w] = edges[i];
        size_t cu = u / CHUNK_SIZE;
        size_t cv = v / CHUNK_SIZE;
        if (cu % T == t) fe(u, v, w);         // a.addEdge(u, v, w);
        if (sym && cv % T == t) fe(v, u, w);  // a.addEdge(v, u, w);
      }
    }
    PERFORMI( auto t3 = timeNow() );
    PERFORMI( dread  += duration(t0, t1) );
    PERFORMI( dparse += duration(t1, t2) );
    PERFORMI( dadd   += duration(t2, t3) );
  }
  PRINTFI("readMtxDoOmp(): read=%.1fms, parse=%.1fms, add=%.1fms\n", dread, dparse, dadd);
}


template <class G>
inline void readMtxOmpW(G& a, istream& s) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  auto fv = [&](auto u) { a.addVertex(K(u)); };
  auto fe = [&](auto u, auto v, auto w) { a.addEdge(K(u), K(v), E(w)); };
  PERFORMI( auto t0 = timeNow() );
  readMtxDoOmp(s, fv, fe);
  PERFORMI( auto t1 = timeNow() );
  updateOmpU(a);
  PERFORMI( auto t2 = timeNow() );
  PRINTFI("readMtxOmpW(): read=%.1fms, update=%.1fms\n", duration(t0, t1), duration(t1, t2));
}
template <class G>
inline void readMtxOmpW(G& a, const char *pth) {
  ifstream f(pth);
  readMtxOmpW(a, f);
}




// WRITE MTX
// ---------

template <class G>
inline void writeMtx(ostream& a, const G& x) {
  a << "%%MatrixMarket matrix coordinate real asymmetric\n";
  a << x.order() << " " << x.order() << " " << x.size() << "\n";
  x.forEachVertexKey([&](auto u) {
    x.forEachEdge(u, [&](auto v, auto w) {
      a << u << " " << v << " " << w << "\n";
    });
  });
}
template <class G>
inline void writeMtx(string pth, const G& x) {
  ofstream f(pth);
  writeMtx(f, x);
  f.close();
}
