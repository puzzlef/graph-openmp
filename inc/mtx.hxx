#pragma once
#include <cstdint>
#include <cstring>
#include <utility>
#include <memory>
#include <tuple>
#include <vector>
#include <string>
#include <string_view>
#include <istream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "_main.hxx"
#include "Graph.hxx"
#include "update.hxx"
#ifdef OPENMP
#include <omp.h>
#endif

using std::unique_ptr;
using std::tuple;
using std::vector;
using std::string;
using std::string_view;
using std::istream;
using std::istringstream;
using std::ifstream;
using std::memcpy;
using std::make_unique;
using std::move;
using std::min;
using std::max;
using std::getline;
using std::strtoull;
using std::strtod;




#pragma region READ FROM FILE
#pragma region READ MTX FORMAT HEADER
/**
 * Read header of an MTX format file (check for errors, handle comments).
 * @param symmetric is graph symmetric (output)
 * @param rows number of rows (output)
 * @param cols number of columns (output)
 * @param size number of lines/edges (output)
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
 * @param symmetric is graph symmetric (output)
 * @param rows number of rows (output)
 * @param cols number of columns (output)
 * @param size number of lines/edges (output)
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




#pragma region READ MTX FORMAT
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




#pragma region READ MTX FORMAT TO GRAPH
/**
 * Read MTX format file as graph, with conditions on vertices and edges.
 * @param a output graph (updated)
 * @param stream input stream (updated)
 * @param weighted is graph weighted?
 * @param fv include vertex? (u, d)
 * @param fe include edge? (u, v, w)
 */
template <class G, class FV, class FE>
inline void readMtxFormatFileToGraphConditionalW(G& a, istream& stream, bool weighted, FV fv, FE fe) {
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
inline void readMtxFormatFileToGraphConditionalW(G& a, const char *pth, bool weighted, FV fv, FE fe) {
  ifstream stream(pth);
  readMtxFormatFileToGraphConditionalW(a, stream, weighted, fv, fe);
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
inline void readMtxFormatFileToGraphConditionalOmpW(G& a, istream& stream, bool weighted, FV fv, FE fe) {
  using K = typename G::key_type;
  using V = typename G::vertex_value_type;
  using E = typename G::edge_value_type;
  auto fh = [&](auto symmetric, auto rows, auto cols, auto size) {
    addVerticesIfU(a, K(1), K(max(rows, cols)+1), V(), fv);
  };
  auto fb = [&](auto u, auto v, auto w) {
    if (fe(K(u), K(v), K(w))) addEdgeOmpU(a, K(u), K(v), E(w));
  };
  readMtxFormatFileDoOmp(stream, weighted, fh, fb);
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
inline void readMtxFormatFileToGraphConditionalOmpW(G& a, const char *pth, bool weighted, FV fv, FE fe) {
  ifstream stream(pth);
  readMtxFormatFileToGraphConditionalOmpW(a, stream, weighted, fv, fe);
}
#endif


/**
 * Read MTX format file as graph.
 * @param a output graph (updated)
 * @param stream input stream (updated)
 * @param weighted is graph weighted?
 */
template <class G>
inline void readMtxFormatFileToGraphW(G& a, istream& stream, bool weighted=false) {
  auto fv = [](auto u, auto d)         { return true; };
  auto fe = [](auto u, auto v, auto w) { return true; };
  readMtxFormatFileToGraphConditionalW(a, stream, weighted, fv, fe);
}


/**
 * Read MTX format file as graph.
 * @param a output graph (updated)
 * @param pth file path
 * @param weighted is graph weighted?
 */
template <class G>
inline void readMtxFormatFileToGraphW(G& a, const char *pth, bool weighted=false) {
  ifstream stream(pth);
  readMtxFormatFileToGraphW(a, stream, weighted);
}


#ifdef OPENMP
/**
 * Read MTX format file as graph.
 * @param a output graph (updated)
 * @param stream input stream (updated)
 * @param weighted is it weighted?
 */
template <class G>
inline void readMtxFormatFileToGraphOmpW(G& a, istream& stream, bool weighted=false) {
  auto fv = [](auto u, auto d)         { return true; };
  auto fe = [](auto u, auto v, auto w) { return true; };
  readMtxFormatFileToGraphConditionalOmpW(a, stream, weighted, fv, fe);
}


/**
 * Read MTX format file as graph.
 * @param a output graph (updated)
 * @param pth file path
 * @param weighted is it weighted?
 */
template <class G>
inline void readMtxFormatFileToGraphOmpW(G& a, const char *pth, bool weighted=false) {
  ifstream stream(pth);
  readMtxFormatFileToGraphOmpW(a, stream, weighted);
}
#endif
#pragma endregion
#pragma endregion




#pragma region READ FROM STRING
#pragma region READ COO FORMAT HEADER
/**
 * Read header of COO format data.
 * @param rows number of rows (output)
 * @param cols number of columns (output)
 * @param size number of lines/edges (output)
 * @param data input data
 * @returns size of header
 */
inline size_t readCooFormatHeaderW(size_t& rows, size_t& cols, size_t& size, string_view data) {
  auto fu = [](char c) { return false; };
  auto fw = [](char c) { return false; };
  auto ib = data.begin(), ie = data.end(), it = ib;
  // Skip past empty lines and comments.
  for (; it!=ie; it = findNextLine(it, ie)) {
    it = findNextNonBlank(it, ie, fu);
    if (*it!='%' || *it!='#' || !isNewline(*it)) break;
  }
  // Read rows, cols, size.
  it = readNumberW<true>(rows, it, ie, fu, fw);  // Number of vertices
  it = readNumberW<true>(cols, it, ie, fu, fw);  // Number of vertices
  it = readNumberW<true>(size, it, ie, fu, fw);  // Number of edges
  // Jump to the next line.
  it = findNextLine(it, ie);
  return it-ib;
}
#pragma endregion




#pragma region READ MTX FORMAT HEADER
/**
 * Read header of MTX format data (check for errors, handle comments).
 * @param symmetric is graph symmetric (output)
 * @param rows number of rows (output)
 * @param cols number of columns (output)
 * @param size number of lines/edges (output)
 * @param data input data
 * @returns size of header
 */
inline size_t readMtxFormatHeaderW(bool& symmetric, size_t& rows, size_t& cols, size_t& size, string_view data) {
  auto fu = [](char c) { return false; };
  auto fw = [](char c) { return false; };
  auto ib = data.begin(), ie = data.end(), it = ib;
  // Skip past the comments and read the graph type.
  string_view h0, h1, h2, h3, h4;
  for (; it!=ie; it = findNextLine(it, ie)) {
    if (*it!='%') break;
    if (data.substr(it-ib, 14)!="%%MatrixMarket") continue;
    it = readTokenW(h0, it, ie, fu, fw);  // %%MatrixMarket
    it = readTokenW(h1, it, ie, fu, fw);  // Graph
    it = readTokenW(h2, it, ie, fu, fw);  // Format
    it = readTokenW(h3, it, ie, fu, fw);  // Field
    it = readTokenW(h4, it, ie, fu, fw);  // Symmetry
  }
  // Check the graph type.
  if (h1!="matrix" || h2!="coordinate") throw FormatError("Invalid MTX header (unknown format)", ib);
  symmetric = h4=="symmetric" || h4=="skew-symmetric";
  // Read rows, cols, size.
  it = readNumberW<true>(rows, it, ie, fu, fw);  // Number of vertices
  it = readNumberW<true>(cols, it, ie, fu, fw);  // Number of vertices
  it = readNumberW<true>(size, it, ie, fu, fw);  // Number of edges
  // Jump to the next line.
  it = findNextLine(it, ie);
  return it-ib;
}
#pragma endregion




#pragma region READ EDGELIST FORMAT
/**
 * Read data in Edgelist format.
 * @tparam WEIGHTED is graph weighted?
 * @tparam BASE base vertex id (0 or 1)
 * @param data input data
 * @param symmetric is graph symmetric?
 * @param fb on body line (u, v, w)
 */
template <bool WEIGHTED=false, int BASE=1, class FB>
inline void readEdgelistFormatDoChecked(string_view data, bool symmetric, FB fb) {
  auto fu = [](char c) { return c==','; };                      // Support CSV
  auto fw = [](char c) { return c==',' || c=='%' || c=='#'; };  // Support CSV, comments
  auto ib = data.begin(), ie = data.end(), it = ib;
  for (; it!=ie; it = findNextLine(it, ie)) {
    // Skip past empty lines and comments.
    it = findNextNonBlank(it, ie, fu);
    if (it==ie || *it=='%' || *it=='#' || isNewline(*it)) continue;
    // Read u, v, w (if weighted).
    int64_t u = 0, v = 0; double w = 1; auto il = it;
    it = readNumberW<true>(u, it, ie, fu, fw);  // Source vertex
    it = readNumberW<true>(v, it, ie, fu, fw);  // Target vertex
    if constexpr (WEIGHTED) {
      it = readNumberW<true>(w, it, ie, fu, fw);  // Edge weight
    }
    if constexpr (BASE) { --u; --v; }  // Convert to zero-based
    if (u<0 || v<0) throw FormatError("Invalid Edgelist body (negative vertex-id)", il);
    fb(u, v, w);
    if (symmetric && u!=v) fb(v, u, w);
  }
}


/**
 * Read EdgeList format data (crazy frog version).
 * @tparam WEIGHTED is graph weighted?
 * @tparam BASE base vertex id (0 or 1)
 * @param data input data
 * @param symmetric is graph symmetric?
 * @param fb on body line (u, v, w)
 */
template <bool WEIGHTED=false, int BASE=1, class FB>
inline void readEdgelistFormatDoUnchecked(string_view data, bool symmetric, FB fb) {
  auto ib = data.begin(), ie = data.end(), it = ib;
  while (true) {
    // Read u, v, w (if weighted).
    uint64_t u = 0, v = 0; double w = 1;
    it = findNextDigit(it, ie);
    if (it==ie) break;  // No more lines
    it = parseWholeNumberW(u, it, ie);  // Source vertex
    it = findNextDigit(it, ie);
    it = parseWholeNumberW(v, it, ie);  // Target vertex
    if constexpr (WEIGHTED) {
      it = findNextDigit(it, ie);
      it = parseFloatW(w, it, ie);  // Edge weight
    }
    if constexpr (BASE) { --u; --v; }  // Convert to zero-based
    fb(u, v, w);
    if (symmetric && u!=v) fb(v, u, w);
  }
}


/**
 * Read EdgeList format data.
 * @tparam WEIGHTED is graph weighted?
 * @tparam BASE base vertex id (0 or 1)
 * @tparam CHECK check for error?
 * @param data input data
 * @param symmetric is graph symmetric?
 * @param fb on body line (u, v, w)
 */
template <bool WEIGHTED=false, int BASE=1, bool CHECK=false, class FB>
inline void readEdgelistFormatDo(string_view data, bool symmetric, FB fb) {
  if constexpr (CHECK) readEdgelistFormatDoChecked<WEIGHTED, BASE>(data, symmetric, fb);
  else readEdgelistFormatDoUnchecked<WEIGHTED, BASE>(data, symmetric, fb);
}


/**
 * Read data in Edgelist format, and record the edges.
 * @tparam WEIGHTED is graph weighted?
 * @tparam BASE base vertex id (0 or 1)
 * @tparam CHECK check for error?
 * @param degrees vertex degrees (updated)
 * @param sources source vertices (output)
 * @param targets target vertices (output)
 * @param weights edge weights (output)
 * @param data input data
 * @param symmetric is graph symmetric?
 */
template <bool WEIGHTED=false, int BASE=1, bool CHECK=false, class IK, class IE>
inline void readEdgelistFormatToListsU(IK degrees, IK sources, IK targets, IE weights, string_view data, bool symmetric) {
  size_t i = 0;
  readEdgelistFormatDo<WEIGHTED, BASE, CHECK>(data, symmetric, [&](auto u, auto v, auto w) {
    sources[i] = u;
    targets[i] = v;
    if constexpr (WEIGHTED) weights[i] = w;
    ++degrees[u];
    ++i;
  });
}


#ifdef OPENMP
/**
 * Get characters to process for an EdgeList format block, skip first partial line [helper function].
 * @param data input data
 * @param b block index
 * @param B block size
 * @returns characters to process for a block
 */
inline string_view readEdgelistFormatBlock(string_view data, size_t b, size_t B) {
  auto db = data.begin(), de = data.end();
  auto bb = db+b, be = min(bb+B, de);
  if (bb!=db && !isNewline(*bb-1)) bb = findNextLine(bb, de);
  if (be!=db && !isNewline(*be-1)) be = findNextLine(be, de);
  return data.substr(bb-db, be-bb);
}


/**
 * Read EdgeList format data, and record the edges and vertex degrees.
 * @tparam WEIGHTED is graph weighted?
 * @tparam BASE base vertex id (0 or 1)
 * @tparam CHECK check for error?
 * @tparam PARTITIONS number of partitions for vertex degrees
 * @param degrees per-partition vertex degrees (updated)
 * @param sources per-thread source vertices (output)
 * @param targets per-thread target vertices (output)
 * @param weights per-thread edge weights (output)
 * @param data input data
 * @param symmetric is graph symmetric
 * @returns per-thread number of edges read
 */
template <bool WEIGHTED=false, int BASE=1, bool CHECK=false, int PARTITIONS=4, class IIK, class IIE>
inline vector<unique_ptr<size_t>> readEdgelistFormatToListsOmpU(IIK degrees, IIK sources, IIK targets, IIE weights, string_view data, bool symmetric) {
  const size_t DATA  = data.size();
  const size_t BLOCK = 256 * 1024;  // Characters per block (256KB)
  const int T = omp_get_max_threads();
  FormatError err;  // Common error
  // Allocate space for per-thread index.
  vector<unique_ptr<size_t>> is(T);
  for (int t=0; t<T; ++t)
    is[t] = make_unique<size_t>();
  // Process a grid in parallel with dynamic scheduling.
  #pragma omp parallel for schedule(dynamic, 1) shared(err)
  for (size_t b=0; b<DATA; b+=BLOCK) {
    int    t = omp_get_thread_num();
    // Get per-thread index.
    size_t i = *is[t];
    // Skip if error occurred.
    if (CHECK && !err.empty()) continue;
    // Read a block of data.
    string_view bdata = readEdgelistFormatBlock(data, b, BLOCK);
    auto fb = [&](auto u, auto v, auto w) {
      // Record the edge.
      sources[t][i] = u;
      targets[t][i] = v;
      if constexpr (WEIGHTED) weights[t][i] = w;
      #pragma omp atomic
      ++degrees[t % PARTITIONS][u];
      ++i;
    };
    if constexpr (CHECK) {
      try { readEdgelistFormatDo<WEIGHTED, BASE, true>(bdata, symmetric, fb); }
      catch (const FormatError& e) { if (err.empty()) err = e; }
    }
    else readEdgelistFormatDo<WEIGHTED, BASE>(bdata, symmetric, fb);
    // Update per-thread index.
    *is[t] = i;
  }
  // Throw error if any.
  if (CHECK && !err.empty()) throw err;
  return is;
}
#endif
#pragma endregion




#pragma region CONVERT EDGELIST TO CSR
/**
 * Convert Edgelist to CSR (lists).
 * @tparam WEIGHTED is graph weighted?
 * @param offsets CSR offsets (output)
 * @param edgeKeys CSR edge keys (output)
 * @param edgeValues CSR edge values (output)
 * @param degrees vertex degrees
 * @param sources source vertices
 * @param targets target vertices
 * @param weights edge weights
 * @param rows number of rows/vertices
 */
template <bool WEIGHTED=false, class IO, class IK, class IE>
inline void convertEdgelistToCsrListsW(IO offsets, IK edgeKeys, IE edgeValues, IK degrees, IK sources, IK targets, IE weights, size_t rows) {
  // Compute offsets.
  exclusiveScanW(&offsets[0], &degrees[0][0], rows+1);
  // Populate CSR format.
  for (size_t i=0; i<rows; ++i) {
    size_t u = sources[i];
    size_t v = targets[i];
    size_t j = offsets[u]++;
    edgeKeys[j] = v;
    if constexpr (WEIGHTED) edgeValues[j] = weights[i];
  }
  // Fix offsets.
  memcpy(&offsets[1], &offsets[0], rows * sizeof(offsets[0]));
  offsets[0] = 0;
}


/**
 * Convert Edgelist to CSR (lists).
 * @tparam WEIGHTED is graph weighted?
 * @tparam PARTITIONS number of partitions for vertex degrees
 * @param offsets CSR offsets (output)
 * @param edgeKeys CSR edge keys (output)
 * @param edgeValues CSR edge values (output)
 * @param poffsets per-partition CSR offsets (output)
 * @param pedgeKeys per-partition CSR edge keys (output)
 * @param pedgeValues per-partition CSR edge values (output)
 * @param degrees per-partition vertex degrees
 * @param sources per-thread source vertices
 * @param targets per-thread target vertices
 * @param weights per-thread edge weights
 * @param counts per-thread number of edges read
 * @param rows number of rows/vertices
 */
template <bool WEIGHTED=false, int PARTITIONS=4, class IO, class IK, class IE, class IIO, class IIK, class IIE>
inline void convertEdgelistToCsrListsOmpW(IO offsets, IK edgeKeys, IE edgeValues, IIO poffsets, IIK pedgeKeys, IIE pedgeValues, IIK degrees, IIK sources, IIK targets, IIE weights, const vector<unique_ptr<size_t>>& counts, size_t rows) {
  int T = omp_get_max_threads();
  vector<size_t> buf(T);
  // Compute offsets.
  if (PARTITIONS==1) exclusiveScanOmpW(&offsets[0], &buf[0], &degrees[0][0], rows+1);
  else {
    for (int t=0; t<PARTITIONS; ++t)
      exclusiveScanOmpW(&poffsets[t][0], &buf[0], &degrees[t][0], rows+1);
  }
  if (PARTITIONS==1) {
    // Populate CSR format.
    #pragma omp parallel
    {
      int t = omp_get_thread_num();
      size_t I = *counts[t];
      for (size_t i=0; i<I; ++i) {
        size_t u = sources[t][i];
        size_t v = targets[t][i];
        size_t j = 0;
        #pragma omp atomic capture
        j = offsets[u]++;
        edgeKeys[j] = v;
        if constexpr (WEIGHTED) edgeValues[j] = weights[t][i];
      }
    }
    // Fix offsets.
    memcpy(&offsets[1], &offsets[0], rows * sizeof(offsets[0]));
    offsets[0] = 0;
  }
  else {
    // Populate per-partition CSR format.
    #pragma omp parallel
    {
      int t = omp_get_thread_num();
      size_t I = *counts[t];
      for (size_t i=0; i<I; ++i) {
        size_t u = sources[t][i];
        size_t v = targets[t][i];
        size_t j = 0;
        #pragma omp atomic capture
        j = poffsets[t % PARTITIONS][u]++;
        pedgeKeys[t % PARTITIONS][j] = v;
        if constexpr (WEIGHTED) pedgeValues[t % PARTITIONS][j] = weights[t][i];
      }
    }
    // Fix per-partition offsets.
    #pragma omp parallel
    {
      int t = omp_get_thread_num();
      if (t<PARTITIONS) memcpy(&poffsets[t][1], &poffsets[t][0], rows * sizeof(poffsets[t][0]));
      if (t<PARTITIONS) poffsets[t][0] = 0;
    }
    // Combine per-partition degrees.
    #pragma omp parallel for schedule(static, 2048)
    for (size_t u=0; u<rows; ++u) {
      if (PARTITIONS==1) {}
      else if (PARTITIONS==2) degrees[0][u] += degrees[1][u];
      else if (PARTITIONS==4) degrees[0][u] += degrees[1][u] + degrees[2][u] + degrees[3][u];
      else if (PARTITIONS==8) degrees[0][u] += degrees[1][u] + degrees[2][u] + degrees[3][u] + degrees[4][u] + degrees[5][u] + degrees[6][u] + degrees[7][u];
      else {
        for (int t=1; t<PARTITIONS; ++t)
          degrees[0][u] += degrees[t][u];
      }
    }
    // Compute global offsets.
    exclusiveScanOmpW(&offsets[0], &buf[0], &degrees[0][0], rows+1);
    // Combine per-partition CSR format.
    #pragma omp parallel for schedule(dynamic, 2048)
    for (size_t u=0; u<rows; ++u) {
      size_t j = offsets[u];
      for (int t=0; t<PARTITIONS; ++t) {
        size_t i = poffsets[t][u];
        size_t I = poffsets[t][u+1];
        for (; i<I; ++i, ++j) {
          edgeKeys[j] = pedgeKeys[t][i];
          if constexpr (WEIGHTED) edgeValues[j] = pedgeValues[t][i];
        }
      }
    }
  }
}
#pragma endregion




#pragma region READ EDGELIST FORMAT TO CSR
/**
 * Read data in MTX format, and convert to CSR.
 * @tparam WEIGHTED is graph weighted?
 * @tparam BASE base vertex id (0 or 1)
 * @tparam CHECK check for error?
 * @param a output csr graph (updated)
 * @param data input data
 */
template <bool WEIGHTED=false, int BASE=1, bool CHECK=false, class G>
inline void readMtxFormatToCsrListsW(G& a, string_view data) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  // Read MTX format header.
  bool symmetric; size_t rows, cols, size;
  size_t head = readMtxFormatHeaderW(symmetric, rows, cols, size, data);
  data.remove_prefix(head);
  // Allocate space for CSR.
  const size_t N = max(rows, cols);
  const size_t M = size;
  a.resize(N, M);
  // Allocate space for sources, targets, and weights.
  const size_t SOURCES = bytesof<K>(M);
  const size_t TARGETS = bytesof<K>(M);
  const size_t WEIGHTS = bytesof<E>(M);
  MappedPtr<char> buf(SOURCES + TARGETS + WEIGHTS);
  K *sources = (K*) buf.data();
  K *targets = (K*) (buf.data() + SOURCES);
  E *weights = (E*) (buf.data() + SOURCES + TARGETS);
  K *degrees = a.degrees.data();
  K *offsets = a.offsets.data();
  K *edgeKeys = a.edgeKeys.data();
  E *edgeValues = a.edgeValues.data();
  // Read Edgelist and convert to CSR.
  readEdgelistFormatToListsU<WEIGHTED, BASE, CHECK>(degrees, sources, targets, weights, data, symmetric);
  convertEdgelistToCsrListsW<WEIGHTED>(offsets, edgeKeys, edgeValues, degrees, sources, targets, weights, N);
}


/**
 * Read data in MTX format, and convert to CSR.
 * @tparam WEIGHTED is graph weighted?
 * @tparam BASE base vertex id (0 or 1)
 * @tparam CHECK check for error?
 * @tparam PARTITIONS number of partitions for vertex degrees
 * @param a output csr graph (updated)
 * @param data input data
 */
template <bool WEIGHTED=false, int BASE=1, bool CHECK=false, int PARTITIONS=4, class G>
inline void readMtxFormatToCsrListsOmpW(G& a, string_view data) {
  using O = typename G::offset_type;
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  // Read MTX format header.
  bool symmetric; size_t rows, cols, size;
  size_t head = readMtxFormatHeaderW(symmetric, rows, cols, size, data);
  data.remove_prefix(head);
  // Allocate space for CSR.
  const size_t N = max(rows, cols);
  const size_t M = size;
  a.resize(N, M);
  // Allocate space for sources, targets, weights, degrees, offsets, edge keys, and edge values.
  const int T = omp_get_max_threads();
  const size_t SOURCES = bytesof<K>(M);
  const size_t TARGETS = bytesof<K>(M);
  const size_t WEIGHTS = bytesof<E>(M);
  const size_t DEGREES     = bytesof<K>(N);
  const size_t OFFSETS     = bytesof<O>(N+1);
  const size_t EDGE_KEYS   = bytesof<K>(M);
  const size_t EDGE_VALUES = bytesof<E>(M);
  MappedPtr<char> buf(T * (SOURCES + TARGETS + WEIGHTS) + (PARTITIONS-1) * (DEGREES + OFFSETS + EDGE_KEYS + EDGE_VALUES));
  vector<K*> sources(T);
  vector<K*> targets(T);
  vector<E*> weights(T);
  for (int i=0; i<T; ++i) {
    sources[i] = (K*) (buf.data() + i * SOURCES);
    targets[i] = (K*) (buf.data() + T * (SOURCES) + i * TARGETS);
    weights[i] = (E*) (buf.data() + T * (SOURCES + TARGETS) + i * WEIGHTS);
  }
  vector<K*> degrees(PARTITIONS);
  vector<O*> offsets(PARTITIONS);
  vector<K*> edgeKeys(PARTITIONS);
  vector<E*> edgeValues(PARTITIONS);
  const char *base = buf.data() + T * (SOURCES + TARGETS + WEIGHTS);
  for (int i=0; i<PARTITIONS; ++i) {
    degrees[i]    = (K*) (base + (i-1) * DEGREES);
    offsets[i]    = (O*) (base + (PARTITIONS-1) * (DEGREES) + (i-1) * OFFSETS);
    edgeKeys[i]   = (K*) (base + (PARTITIONS-1) * (DEGREES + OFFSETS) + (i-1) * EDGE_KEYS);
    edgeValues[i] = (E*) (base + (PARTITIONS-1) * (DEGREES + OFFSETS + EDGE_KEYS) + (i-1) * EDGE_VALUES);
  }
  // Read Edgelist and convert to CSR.
  vector<unique_ptr<size_t>> counts = readEdgelistFormatToListsOmpU<WEIGHTED, BASE, CHECK, PARTITIONS>(degrees, sources, targets, weights, data, symmetric);
  convertEdgelistToCsrListsOmpW<WEIGHTED, PARTITIONS>(a.offsets.data(), a.edgeKeys.data(), a.edgeValues.data(), offsets, edgeKeys, edgeValues, degrees, sources, targets, weights, counts, N);
}
#pragma endregion
#pragma endregion
