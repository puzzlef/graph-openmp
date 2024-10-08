#pragma once
#include <utility>
#include <string>
#include <istream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cstdlib>
#include "_main.hxx"
#include "Graph.hxx"
#include "update.hxx"
#ifdef OPENMP
#include <omp.h>
#endif

using std::tuple;
using std::string;
using std::istream;
using std::istringstream;
using std::ifstream;
using std::ofstream;
using std::move;
using std::max;
using std::getline;




#pragma region READ FROM FILE
#pragma region READ HEADER
/**
 * Read header of an MTX format file (check for errors, handle comments).
 * @param symmetric is graph symmetric (updated)
 * @param rows number of rows (updated)
 * @param cols number of columns (updated)
 * @param size number of lines/edges (updated)
 * @param stream input stream (updated)
 */
inline void readMtxFormatFileHeaderW(bool& symmetric, size_t& rows, size_t& cols, size_t& size, istream& stream) {
  string line, h0, h1, h2, h3, h4;
  // Skip past the comments and read the graph type.
  while (true) {
    getline(stream, line);
    if (line.find('%')!=0) break;
    if (line.find("%%")!=0) continue;
    istringstream sline(line);
    sline >> h0 >> h1 >> h2 >> h3 >> h4;
  }
  // Check the graph type.
  if (h1!="matrix" || h2!="coordinate") throw FormatError("Invalid MTX header (unknown format)");
  symmetric = h4=="symmetric" || h4=="skew-symmetric";
  // Read rows, cols, size.
  istringstream sline(line);
  sline >> rows >> cols >> size;
}


/**
 * Read header of an MTX format file (check for errors, handle comments).
 * @param symmetric is graph symmetric (updated)
 * @param rows number of rows (updated)
 * @param cols number of columns (updated)
 * @param size number of lines/edges (updated)
 * @param pth file path
 */
inline void readMtxFormatFileHeaderW(bool& symmetric, size_t& rows, size_t& cols, size_t& size, const char *pth) {
  ifstream stream(pth);
  return readMtxFormatFileHeaderW(symmetric, rows, cols, size, stream);
}


/**
 * Read order of graph from an MTX format file.
 * @param stream input stream (updated)
 * @returns number of vertices (1 to N)
 */
inline size_t readMtxFormatFileOrder(istream& stream) {
  bool symmetric; size_t rows, cols, size;
  readMtxFormatFileHeaderW(symmetric, rows, cols, size, stream);
  return max(rows, cols);
}


/**
 * Read order of graph from an MTX format file.
 * @param pth file path
 * @returns number of vertices (1 to N)
 */
inline size_t readMtxFormatFileOrder(const char *pth) {
  ifstream stream(pth);
  return readMtxFormatFileOrder(stream);
}


/**
 * Read size of graph from an MTX format file.
 * @param stream input stream (updated)
 * @returns number of edges
 */
inline size_t readMtxFormatFileSize(istream& stream) {
  bool symmetric; size_t rows, cols, size;
  readMtxFormatFileHeaderW(symmetric, rows, cols, size, stream);
  return size;
}


/**
 * Read size of graph from an MTX format file.
 * @param pth file path
 * @returns number of edges
 */
inline size_t readMtxFormatFileSize(const char *pth) {
  ifstream stream(pth);
  return readMtxFormatFileSize(stream);
}


/**
 * Read span of graph from an MTX format file.
 * @param stream input stream (updated)
 * @returns span of graph (max vertex id + 1)
 */
inline size_t readMtxFormatFileSpan(istream& stream) {
  bool symmetric; size_t rows, cols, size;
  readMtxFormatFileHeaderW(symmetric, rows, cols, size, stream);
  return max(rows, cols) + 1;
}


/**
 * Read span of graph from an MTX format file.
 * @param pth file path
 * @returns span of graph (max vertex id + 1)
 */
inline size_t readMtxFormatFileSpan(const char *pth) {
  ifstream stream(pth);
  return readMtxFormatFileSpan(stream);
}
#pragma endregion




#pragma region READ DO
/**
 * Read contents of an MTX format file.
 * @param stream input stream (updated)
 * @param weighted is graph weighted?
 * @param fh on header (symmetric, rows, cols, size)
 * @param fb on body line (u, v, w)
 */
template <class FH, class FB>
inline void readMtxFormatFileDo(istream& stream, bool weighted, FH fh, FB fb) {
  bool symmetric; size_t rows, cols, size;
  readMtxFormatFileHeaderW(symmetric, rows, cols, size, stream);
  fh(symmetric, rows, cols, size);
  size_t N = max(rows, cols);
  if (N==0) return;
  // Process body lines sequentially.
  string line;
  while (getline(stream, line)) {
    size_t u, v; double w = 1;
    istringstream sline(line);
    if (!(sline >> u >> v)) break;
    if (weighted) sline >> w;
    fb(u, v, w);
    if (symmetric) fb(v, u, w);
  }
}


/**
 * Read contents of an MTX format file.
 * @param pth file path
 * @param weighted is graph weighted?
 * @param fh on header (symmetric, rows, cols, size)
 * @param fb on body line (u, v, w)
 */
template <class FH, class FB>
inline void readMtxFormatFileDo(const char *pth, bool weighted, FH fh, FB fb) {
  ifstream stream(pth);
  readMtxFormatFileDo(stream, weighted, fh, fb);
}


#ifdef OPENMP
/**
 * Read contents of an MTX format file.
 * @param stream input stream (updated)
 * @param weighted is graph weighted?
 * @param fh on header (symmetric, rows, cols, size)
 * @param fb on body line (u, v, w)
 */
template <class FH, class FB>
inline void readMtxFormatFileDoOmp(istream& stream, bool weighted, FH fh, FB fb) {
  bool symmetric; size_t rows, cols, size;
  readMtxFormatFileHeaderW(symmetric, rows, cols, size, stream);
  fh(symmetric, rows, cols, size);
  size_t N = max(rows, cols);
  if (N==0) return;
  // Process body lines in parallel.
  constexpr int LINES = 131072;
  vector<string> lines(LINES);
  vector<tuple<size_t, size_t, double>> edges(LINES);
  while (true) {
    // Read several lines from the stream.
    int READ = 0;
    for (int i=0; i<LINES; ++i, ++READ)
      if (!getline(stream, lines[i])) break;
    if (READ==0) break;
    // Parse lines using multiple threads.
    #pragma omp parallel for schedule(dynamic, 1024)
    for (int i=0; i<READ; ++i) {
      char *line = (char*) lines[i].c_str();
      size_t u = strtoull(line, &line, 10);
      size_t v = strtoull(line, &line, 10);
      double w = weighted? strtod(line, &line) : 0;
      edges[i] = {u, v, w? w : 1};
    }
    // Notify parsed lines.
    #pragma omp parallel
    {
      for (int i=0; i<READ; ++i) {
        const auto& [u, v, w] = edges[i];
        fb(u, v, w);
        if (symmetric) fb(v, u, w);
      }
    }
  }
}


/**
 * Read contents of an MTX format file.
 * @param pth file path
 * @param weighted is graph weighted?
 * @param fh on header (symmetric, rows, cols, size)
 * @param fb on body line (u, v, w)
 */
template <class FH, class FB>
inline void readMtxFormatFileDoOmp(const char *pth, bool weighted, FH fh, FB fb) {
  ifstream stream(pth);
  readMtxFormatFileDoOmp(stream, weighted, fh, fb);
}
#endif
#pragma endregion




#pragma region READ CONDITIONAL
/**
 * Read MTX format file as graph, with conditions on vertices and edges.
 * @param a output graph (updated)
 * @param stream input stream (updated)
 * @param weighted is graph weighted?
 * @param fv include vertex? (u, d)
 * @param fe include edge? (u, v, w)
 */
template <class G, class FV, class FE>
inline void readMtxFormatFileAsGraphConditionalW(G& a, istream& stream, bool weighted, FV fv, FE fe) {
  using K = typename G::key_type;
  using V = typename G::vertex_value_type;
  using E = typename G::edge_value_type;
  auto fh = [&](auto symmetric, auto rows, auto cols, auto size) {
    addVerticesIfU(a, K(1), K(max(rows, cols) + 1), V(), fv);
  };
  auto fb = [&](auto u, auto v, auto w) {
    if (fe(K(u), K(v), K(w))) a.addEdge(K(u), K(v), E(w));
  };
  readMtxFormatFileDo(stream, weighted, fh, fb);
  a.update();
}


/**
 * Read MTX format file as graph, with conditions on vertices and edges.
 * @param a output graph (updated)
 * @param pth file path
 * @param weighted is graph weighted?
 * @param fv include vertex? (u, d)
 * @param fe include edge? (u, v, w)
 */
template <class G, class FV, class FE>
inline void readMtxFormatFileAsGraphConditionalW(G& a, const char *pth, bool weighted, FV fv, FE fe) {
  ifstream stream(pth);
  readMtxFormatFileAsGraphConditionalW(a, stream, weighted, fv, fe);
}


#ifdef OPENMP
/**
 * Read MTX format file as graph, with conditions on vertices and edges.
 * @param a output graph (updated)
 * @param stream input stream (updated)
 * @param weighted is graph weighted?
 * @param fv include vertex? (u, d)
 * @param fe include edge? (u, v, w)
 */
template <class G, class FV, class FE>
inline void readMtxFormatFileAsGraphConditionalOmpW(G& a, istream& stream, bool weighted, FV fv, FE fe) {
  using K = typename G::key_type;
  using V = typename G::vertex_value_type;
  using E = typename G::edge_value_type;
  auto fh = [&](auto symmetric, auto rows, auto cols, auto size) {
    addVerticesIfU(a, K(1), K(max(rows, cols)+1), V(), fv);
  };
  auto fb = [&](auto u, auto v, auto w) {
    if (fe(K(u), K(v), K(w))) addEdgeOmpU(a, K(u), K(v), E(w));
  };
  readMtxDoOmp(stream, weighted, fh, fb);
  updateOmpU(a);
}


/**
 * Read MTX format file as graph, with conditions on vertices and edges.
 * @param a output graph (updated)
 * @param pth file path
 * @param weighted is graph weighted?
 * @param fv include vertex? (u, d)
 * @param fe include edge? (u, v, w)
 */
template <class G, class FV, class FE>
inline void readMtxFormatFileAsGraphConditionalOmpW(G& a, const char *pth, bool weighted, FV fv, FE fe) {
  ifstream stream(pth);
  readMtxFormatFileAsGraphConditionalOmpW(a, stream, weighted, fv, fe);
}
#endif
#pragma endregion




#pragma region READ
/**
 * Read MTX format file as graph.
 * @param a output graph (updated)
 * @param stream input stream (updated)
 * @param weighted is graph weighted?
 */
template <class G>
inline void readMtxFormatFileW(G& a, istream& stream, bool weighted=false) {
  auto fv = [](auto u, auto d)         { return true; };
  auto fe = [](auto u, auto v, auto w) { return true; };
  readMtxFormatFileAsGraphConditionalW(a, stream, weighted, fv, fe);
}


/**
 * Read MTX format file as graph.
 * @param a output graph (updated)
 * @param pth file path
 * @param weighted is graph weighted?
 */
template <class G>
inline void readMtxFormatFileW(G& a, const char *pth, bool weighted=false) {
  ifstream stream(pth);
  readMtxFormatFileW(a, stream, weighted);
}


#ifdef OPENMP
/**
 * Read MTX format file as graph.
 * @param a output graph (updated)
 * @param stream input stream (updated)
 * @param weighted is it weighted?
 */
template <class G>
inline void readMtxFormatFileOmpW(G& a, istream& stream, bool weighted=false) {
  auto fv = [](auto u, auto d)         { return true; };
  auto fe = [](auto u, auto v, auto w) { return true; };
  readMtxFormatFileAsGraphConditionalOmpW(a, stream, weighted, fv, fe);
}


/**
 * Read MTX format file as graph.
 * @param a output graph (updated)
 * @param pth file path
 * @param weighted is it weighted?
 */
template <class G>
inline void readMtxFormatFileOmpW(G& a, const char *pth, bool weighted=false) {
  ifstream stream(pth);
  readMtxFormatFileOmpW(a, stream, weighted);
}
#endif
#pragma endregion
#pragma endregion
