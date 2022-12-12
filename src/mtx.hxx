#pragma once
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
using std::stringstream;
using std::ifstream;
using std::ofstream;
using std::move;
using std::max;
using std::getline;




// READ MTX
// --------

inline size_t readMtxHeader(istream& s, bool& sym, size_t& rows, size_t& cols, size_t& size) {
  string ln, h0, h1, h2, h3, h4;
  // Skip past the comments and read the graph type.
  while (true) {
    getline(s, ln);
    if (ln.find('%')!=0) break;
    if (ln.find("%%")!=0) continue;
    stringstream ls(ln);
    ls >> h0 >> h1 >> h2 >> h3 >> h4;
  }
  if (h1!="matrix" || h2!="coordinate") return 0;
  sym = h4=="symmetric" || h4=="skew-symmetric";
  // Read rows, cols, size.
  stringstream ls(ln);
  ls >> rows >> cols >> size;
  return max(rows, cols);
}


template <class FE>
inline bool readMtxLineDo(const string& ln, bool sym, FE fe) {
  // Line: <from> <to> [weight].
  size_t u, v;
  double w = 1;
  stringstream ls(ln);
  if (!(ls >> u >> v)) return false;
  ls >> w;
  fe(u, v, w);  // a.addEdge(u, v);
  if (sym) fe(v, u, w);  // a.addEdge(v, u);
  return true;
}


template <class FV, class FE>
inline void readMtxDo(istream& s, FV fv, FE fe) {
  bool sym; size_t rows, cols, size;
  size_t n = readMtxHeader(s, sym, rows, cols, size);
  if (n==0) return;
  // Add all vertices first.
  // Prevent unnecessary respan.
  fv(n);
  for (size_t u=1; u<=n; ++u)
    fv(u);
  // Then we add the edges.
  string line;
  while (getline(s, line))
    if (!readMtxLineDo(line, sym, fe)) break;
}


template <class G>
inline void readMtxW(G& a, istream& s) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  auto fv = [&](auto u) { a.addVertex(K(u)); };
  auto fe = [&](auto u, auto v, auto w) { a.addEdge(K(u), K(v), E(w)); };
  readMtxDo(s, fv, fe);
  a.update();
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
  using EDGE = tuple<size_t, size_t, double>;
  bool sym; size_t rows, cols, size;
  size_t n = readMtxHeader(s, sym, rows, cols, size);
  if (n==0) return;
  // Add all vertices first.
  // Prevent unnecessary respan.
  fv(n);
  for (size_t u=1; u<=n; ++u)
    fv(u);
  // Then we add the edges.
  // Using 2d vector (for `edges`) causes false sharing.
  const int THREADS = omp_get_max_threads();
  const int LINES   = 8192;
  vector<string> lines;
  vector<vector<EDGE>*> edges(THREADS);
  for (int i=0; i<THREADS; ++i)
    edges[i] = new vector<EDGE>();
  while (true) {
    // Read several lines from the stream.
    for (int l=0; l<LINES; ++l) {
      string line;
      if (!getline(s, line)) break;
      lines.push_back(move(line));
    }
    if (lines.empty()) break;
    // Parse lines using multiple threads onto personal edge lists.
    int L = lines.size();
    #pragma omp parallel for schedule(auto)
    for (int l=0; l<L; ++l) {
      int t = omp_get_thread_num();
      readMtxLineDo(lines[l], sym, [&](auto u, auto v, auto w) {
        edges[t]->push_back({u, v, w});
      });
    }
    // Scan all edge lists in each thread, and add edge if the
    // source vertex belongs to this thread.
    const size_t CHUNK_SIZE = 2048;
    #pragma omp parallel
    {
      int T = omp_get_num_threads();
      int t = omp_get_thread_num();
      for (int i=0; i<THREADS; ++i) {
        for (const auto& [u, v, w] : *edges[i]) {
          size_t chunk = u / CHUNK_SIZE;
          if (chunk % T == t) fe(u, v, w);
        }
      }
    }
    // Reset buffers.
    lines.clear();
    for (int i=0, I=edges.size(); i<I; ++i)
      edges[i]->clear();
  }
  for (int i=0, I=edges.size(); i<I; ++i)
    delete edges[i];
}


template <class G>
inline void readMtxOmpW(G& a, istream& s) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  auto fv = [&](auto u) { a.addVertex(K(u)); };
  auto fe = [&](auto u, auto v, auto w) { a.addEdge(K(u), K(v), E(w)); };
  readMtxDoOmp(s, fv, fe);
  updateOmpU(a);
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
