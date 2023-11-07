#pragma once
#include <memory>
#include <tuple>
#include <string>
#include <string_view>
#include <vector>
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
using std::string;
using std::string_view;
using std::vector;
using std::istream;
using std::istringstream;
using std::ifstream;
using std::make_unique;
using std::min;
using std::max;




#pragma region METHODS
#pragma region READ COO FORMAT HEADER
/**
 * Read header of a COO format file.
 * @param rows number of rows (output)
 * @param cols number of columns (output)
 * @param size number of lines/edges (output)
 * @param data input file data
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
 * Read header of an MTX format file.
 * @param symmetric is graph symmetric (output)
 * @param rows number of rows (output)
 * @param cols number of columns (output)
 * @param size number of lines/edges (output)
 * @param data input file data
 * @returns size of header
 */
inline auto readMtxFormatHeaderW(bool& symmetric, size_t& rows, size_t& cols, size_t& size, string_view data) {
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
 * Read an EdgeList format file (check for errors, handle comments).
 * @param data input file data (updated)
 * @param symmetric is graph symmetric
 * @param weighted is graph weighted
 * @param fb on body line (u, v, w)
 */
template <class FB>
inline void readEdgelistFormatDoChecked(string_view data, bool symmetric, bool weighted, FB fb) {
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
    if (weighted) {
      it = readNumberW<true>(w, it, ie, fu, fw);  // Edge weight
    }
    if (u<0 || v<0) throw FormatError("Invalid Edgelist body (negative vertex-id)", il);
    fb(u, v, w);
    if (symmetric && u!=v) fb(v, u, w);
  }
}


/**
 * Read an EdgeList format file (crazy frog version).
 * @param data input file data (updated)
 * @param symmetric is graph symmetric
 * @param weighted is graph weighted
 * @param fb on body line (u, v, w)
 */
template <class FB>
inline void readEdgelistFormatDoUnchecked(string_view data, bool symmetric, bool weighted, FB fb) {
  auto ib = data.begin(), ie = data.end(), it = ib;
  while (true) {
    // Read u, v, w (if weighted).
    uint64_t u = 0, v = 0; double w = 1;
    it = findNextDigit(it, ie);
    if (it==ie) break;  // No more lines
    it = parseWholeNumberSimdW(u, it, ie);  // Source vertex
    it = findNextDigit(it, ie);
    it = parseWholeNumberSimdW(v, it, ie);  // Target vertex
    if (weighted) {
      it = findNextDigit(it, ie);
      it = parseFloatSimdW(w, it, ie);  // Edge weight
    }
    fb(u, v, w);
    if (symmetric && u!=v) fb(v, u, w);
  }
}


/**
 * Read an EdgeList format file.
 * @tparam CHECK check for error?
 * @param data input file data (updated)
 * @param symmetric is graph symmetric
 * @param weighted is graph weighted
 * @param fb on body line (u, v, w)
 */
template <bool CHECK=false, class FB>
inline void readEdgelistFormatDo(string_view data, bool symmetric, bool weighted, FB fb) {
  if constexpr (CHECK) readEdgelistFormatDoChecked(data, symmetric, weighted, fb);
  else readEdgelistFormatDoUnchecked(data, symmetric, weighted, fb);
}


/**
 * Read an EdgeList format file, and record the edges and vertex degrees.
 * @tparam CHECK check for error?
 * @param sources source vertices (output)
 * @param targets target vertices (output)
 * @param weights edge weights (output)
 * @param degrees vertex degrees (updated, must be initialized, optional)
 * @param data input file data
 * @param symmetric is graph symmetric
 * @param weighted is graph weighted
 * @returns number of edges read
 */
template <bool CHECK=false, class K, class E>
inline size_t readEdgelistFormatU(K *sources, K *targets, E *weights, K *degrees, string_view data, bool symmetric, bool weighted) {
  size_t i = 0;
  readEdgelistFormatDo<CHECK>(data, symmetric, weighted, [&](auto u, auto v, auto w) {
    // Record the edge.
    sources[i] = u;
    targets[i] = v;
    if (weighted) weights[i] = w;
    // Update degree of source vertex.
    if (degrees) ++degrees[u];
    ++i;
  });
  // Return number of edges.
  return i;
}


#ifdef OPENMP
/**
 * Get characters to process for an EdgeList format block, skip first partial line [helper function].
 * @param data input file data
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
 * Read an EdgeList format file, and record the edges and vertex degrees.
 * @tparam CHECK check for error?
 * @param sources per-thread source vertices (output)
 * @param targets per-thread target vertices (output)
 * @param weights per-thread edge weights (output)
 * @param degrees per-thread vertex degrees (updated, must be initialized, optional)
 * @param data input file data
 * @param symmetric is graph symmetric
 * @param weighted is graph weighted
 * @returns total number of edges read
 */
template <bool CHECK=false, class K, class E>
inline vector<unique_ptr<size_t>> readEdgelistFormatOmpU(K **sources, K **targets, E **weights, K **degrees, string_view data, bool symmetric, bool weighted) {
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
      if (weighted) weights[t][i] = w;
      // Update degree of source vertex.
      if (degrees) ++degrees[t][u];
      ++i;
    };
    if constexpr (CHECK) {
      try { readEdgelistFormatDoChecked(bdata, symmetric, weighted, fb); }
      catch (const FormatError& e) { if (err.empty()) err = e; }
    }
    else readEdgelistFormatDoUnchecked(bdata, symmetric, weighted, fb);
    // Update per-thread index.
    *is[t] = i;
  }
  // Update counts, and get total number of edges.
  // size_t m = 0;
  // for (int t=0; t<T; ++t) {
  //   size_t i = *is[t];
  //   if (counts) counts[t] = i;
  //   m += i;
  // }
  // Throw error if any.
  if (CHECK && !err.empty()) throw err;
  // Return number of edges.
  return is;
}
#endif
#pragma endregion




#pragma region READ COO FORMAT AND PERFORM
/**
 * Read a COO format file.
 * @param data input file data
 * @param symmetric is graph symmetric
 * @param weighted is graph weighted
 * @param fh on header (rows, cols, size)
 * @param fb on body line (u, v, w)
 * @returns true if error occurred
 */
template <class FH, class FB>
inline bool readCooFormatDo(string_view data, bool symmetric, bool weighted, FH fh, FB fb) {
  bool err = false;
  // size_t rows, cols, size;
  // err |= readCooFormatHeaderU(rows, cols, size, data);
  // if (err) return err;
  // fh(rows, cols, size);
  // size_t n = max(rows, cols);
  // if (n==0) return err;  // Empty graph
  // err |= readEdgelistFormatDoU(data, symmetric, weighted, fb);
  return err;
}


#ifdef OPENMP
/**
 * Read a COO format file.
 * @param data input file data
 * @param symmetric is graph symmetric
 * @param weighted is graph weighted
 * @param fh on header (rows, cols, size)
 * @param fb on body line (u, v, w)
 * @returns true if error occurred
 */
template <class FH, class FB>
inline bool readCooFormatDoOmp(string_view data, bool symmetric, bool weighted, FH fh, FB fb) {
  bool err = false;
  // size_t rows, cols, size;
  // err |= readCooFormatHeaderU(rows, cols, size, data);
  // if (err) return err;
  // fh(rows, cols, size);
  // size_t n = max(rows, cols);
  // if (n==0) return err;  // Empty graph
  // err |= readEdgelistFormatDoOmpU(data, symmetric, weighted, fb);
  return err;
}
#endif
#pragma endregion




#pragma region READ MTX FORMAT AND PERFORM
/**
 * Read an MTX format file.
 * @param data input file data
 * @param weighted is graph weighted
 * @param fh on header (symmetric, rows, cols, size)
 * @param fb on body line (u, v, w)
 * @returns true if error occurred
 */
template <class FH, class FB>
inline bool readMtxFormatDo(string_view data, bool weighted, FH fh, FB fb) {
  bool err = false;
  // bool symmetric; size_t rows, cols, size;
  // err |= readMtxFormatHeaderU(symmetric, rows, cols, size, data);
  // if (err) return err;
  // fh(symmetric, rows, cols, size);
  // size_t n = max(rows, cols);
  // if (n==0) return err;  // Empty graph
  // err |= readEdgelistFormatDoU(data, symmetric, weighted, fb);
  return err;
}


#ifdef OPENMP
/**
 * Read an MTX format file.
 * @param data input file data
 * @param weighted is graph weighted
 * @param fh on header (symmetric, rows, cols, size)
 * @param fb on body line (u, v, w)
 * @returns true if error occurred
 */
template <class FH, class FB>
inline bool readMtxFormatDoOmp(string_view data, bool weighted, FH fh, FB fb) {
  bool err = false;
  // bool symmetric; size_t rows, cols, size;
  // err |= readMtxFormatHeaderU(symmetric, rows, cols, size, data);
  // if (err) return err;
  // fh(symmetric, rows, cols, size);
  // size_t n = max(rows, cols);
  // if (n==0) return err;
  // err |= readEdgelistFormatDoOmpU(data, symmetric, weighted, fb);
  return err;
}
#endif
#pragma endregion




#pragma region READ COO FORMAT CONDITIONALLY
/**
 * Read a COO format file as graph if test passes.
 * @param a output graph (output)
 * @param data input file data
 * @param symmetric is graph symmetric
 * @param weighted is graph weighted
 * @param fv include vertex? (u, d)
 * @param fe include edge? (u, v, w)
 */
template <class G, class FV, class FE>
inline void readCooFormatIfW(G& a, string_view data, bool symmetric, bool weighted, FV fv, FE fe) {
  using K = typename G::key_type;
  using V = typename G::vertex_value_type;
  using E = typename G::edge_value_type;
  a.clear();  // Ensure that the graph is empty
  auto fh = [&](auto rows, auto cols, auto size) { addVerticesIfU(a, K(1), K(max(rows, cols)+1), V(), fv); };
  auto fb = [&](auto u, auto v, auto w) { if (fe(K(u), K(v), K(w))) a.addEdge(K(u), K(v), E(w)); };
  readCooFormatDo(data, symmetric, weighted, fh, fb);
  a.update();
}


#ifdef OPENMP
/**
 * Read a COO format file as graph if test passes.
 * @param a output graph (output)
 * @param data input file data
 * @param symmetric is graph symmetric
 * @param weighted is graph weighted
 * @param fv include vertex? (u, d)
 * @param fe include edge? (u, v, w)
 */
template <class G, class FV, class FE>
inline void readCooFormatIfOmpW(G& a, string_view data, bool symmetric, bool weighted, FV fv, FE fe) {
  using K = typename G::key_type;
  using V = typename G::vertex_value_type;
  using E = typename G::edge_value_type;
  a.clear();  // Ensure that the graph is empty
  auto fh = [&](auto rows, auto cols, auto size) { addVerticesIfU(a, K(1), K(max(rows, cols)+1), V(), fv); };
  auto fb = [&](auto u, auto v, auto w) { if (fe(K(u), K(v), K(w))) addEdgeOmpU(a, K(u), K(v), E(w)); };
  readCooFormatDoOmp(data, symmetric, weighted, fh, fb);
  a.update();
}
#endif
#pragma endregion




#pragma region READ MTX FORMAT CONDITIONALLY
/**
 * Read an MTX format file as graph if test passes.
 * @param a output graph (output)
 * @param data input file data
 * @param weighted is graph weighted
 * @param fv include vertex? (u, d)
 * @param fe include edge? (u, v, w)
 */
template <class G, class FV, class FE>
inline void readMtxFormatIfW(G& a, string_view data, bool weighted, FV fv, FE fe) {
  using K = typename G::key_type;
  using V = typename G::vertex_value_type;
  using E = typename G::edge_value_type;
  a.clear();  // Ensure that the graph is empty
  auto fh = [&](auto symmetric, auto rows, auto cols, auto size) { addVerticesIfU(a, K(1), K(max(rows, cols)+1), V(), fv); };
  auto fb = [&](auto u, auto v, auto w) { if (fe(K(u), K(v), K(w))) a.addEdge(K(u), K(v), E(w)); };
  readMtxFormatDo(data, weighted, fh, fb);
  a.update();
}


#ifdef OPENMP
/**
 * Read an MTX format file as graph if test passes.
 * @param a output graph (output)
 * @param data input file data
 * @param weighted is graph weighted
 * @param fv include vertex? (u, d)
 * @param fe include edge? (u, v, w)
 */
template <class G, class FV, class FE>
inline void readMtxFormatIfOmpW(G& a, string_view data, bool weighted, FV fv, FE fe) {
  using K = typename G::key_type;
  using V = typename G::vertex_value_type;
  using E = typename G::edge_value_type;
  a.clear();  // Ensure that the graph is empty
  auto fh = [&](auto symmetric, auto rows, auto cols, auto size) { addVerticesIfU(a, K(1), K(max(rows, cols)+1), V(), fv); };
  auto fb = [&](auto u, auto v, auto w) { if (fe(K(u), K(v), K(w))) addEdgeOmpU(a, K(u), K(v), E(w)); };
  readMtxFormatDoOmp(data, weighted, fh, fb);
  updateOmpU(a);
}
#endif
#pragma endregion




#pragma region READ COO FORMAT
/**
 * Read a COO format file as graph.
 * @param a output graph (output)
 * @param data input file data
 * @param symmetric is graph symmetric
 * @param weighted is graph weighted
 */
template <class G>
inline void readCooFormatW(G& a, string_view data, bool symmetric, bool weighted=false) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  a.clear();  // Ensure that the graph is empty
  auto fh = [&](auto rows, auto cols, auto size) {
    addVerticesU(a, K(1), K(max(rows, cols)+1));
  };
  auto fb = [&](auto u, auto v, auto w) { addEdgeU(a, K(u), K(v), E(w)); };
  readCooFormatDo(data, symmetric, weighted, fh, fb);
  a.update();
}


#ifdef OPENMP
/**
 * Read a COO format file as graph.
 * @param a output graph (output)
 * @param data input file data
 * @param symmetric is graph symmetric
 * @param weighted is graph weighted
 */
template <class G>
inline void readCooFormatOmpW(G& a, string_view data, bool symmetric, bool weighted=false) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  a.clear();  // Ensure that the graph is empty
  auto fh = [&](auto rows, auto cols, auto size) {
    addVerticesU(a, K(1), K(max(rows, cols)+1));
  };
  auto fb = [&](auto u, auto v, auto w) { addEdgeOmpU(a, K(u), K(v), E(w)); };
  readCooFormatDoOmp(data, symmetric, weighted, fh, fb);
  updateOmpU(a);
}
#endif
#pragma endregion




#pragma region READ MTX FORMAT
/**
 * Read an MTX format file as graph.
 * @param a output graph (output)
 * @param data input file data
 * @param weighted is graph weighted
 */
template <class G>
inline void readMtxFormatW(G& a, string_view data, bool weighted=false) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  a.clear();  // Ensure that the graph is empty
  auto fh = [&](auto symmetric, auto rows, auto cols, auto size) {
    addVerticesU(a, K(1), K(max(rows, cols)+1));
  };
  auto fb = [&](auto u, auto v, auto w) { addEdgeU(a, K(u), K(v), E(w)); };
  readMtxFormatDo(data, weighted, fh, fb);
  a.update();
}


#ifdef OPENMP
/**
 * Read an MTX format file as graph.
 * @param a output graph (output)
 * @param data input file data
 * @param weighted is graph weighted
 */
template <class G>
inline void readMtxFormatOmpW(G& a, string_view data, bool weighted=false) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  a.clear();  // Ensure that the graph is empty
  auto fh = [&](auto symmetric, auto rows, auto cols, auto size) {
    addVerticesU(a, K(1), K(max(rows, cols)+1));
  };
  auto fb = [&](auto u, auto v, auto w) { addEdgeOmpU(a, K(u), K(v), E(w)); };
  readMtxFormatDoOmp(data, weighted, fh, fb);
  updateOmpU(a);
}
#endif
#pragma endregion
#pragma endregion
