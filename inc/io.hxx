#pragma once
#include <tuple>
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

using std::tuple;
using std::string;
using std::string_view;
using std::istream;
using std::istringstream;
using std::ifstream;
using std::min;
using std::max;
using std::getline;




#pragma region METHODS
#pragma region READ COO FORMAT HEADER
/**
 * Read header of a COO format file.
 * @param rows number of rows (output)
 * @param cols number of columns (output)
 * @param size number of lines/edges (output)
 * @param data input file data (updated)
 * @returns iterator to error in data
 */
inline auto readCooFormatHeaderU(size_t& rows, size_t& cols, size_t& size, string_view& data) {
  using I = decltype(data.begin());
  auto fu = [](char c) { return false; };
  auto fw = [](char c) { return false; };
  auto ib = data.begin(), ie = data.end(), it = ib;
  // Skip past empty lines and comments.
  for (; it!=ie; it = findNextLine(it, ie)) {
    it = findNextNonBlank(it, ie, fu);
    if (*it!='%' || *it!='#' || !isNewline(*it)) break;
  }
  // Read rows, cols, size.
  I err = I();
  do {
    err = readNumberU<true>(rows, it, ie, fu, fw); if (err) break;  // Number of vertices
    err = readNumberU<true>(cols, it, ie, fu, fw); if (err) break;  // Number of vertices
    err = readNumberU<true>(size, it, ie, fu, fw);                  // Number of edges
  } while (false);
  if (err) return err;  // Error: Invalid COO format header
  // Jump to the next line.
  it = findNextLine(it, ie);
  data.remove_prefix(it-ib);
  return I();
}
#pragma endregion




#pragma region READ MTX FORMAT HEADER
/**
 * Read header of an MTX format file.
 * @param symmetric is graph symmetric (output)
 * @param rows number of rows (output)
 * @param cols number of columns (output)
 * @param size number of lines/edges (output)
 * @param data input file data (updated)
 * @returns iterator to error in data
 */
inline auto readMtxFormatHeaderU(bool& symmetric, size_t& rows, size_t& cols, size_t& size, string_view& data) {
  using I = decltype(data.begin());
  auto fu = [](char c) { return false; };
  auto fw = [](char c) { return false; };
  auto ib = data.begin(), ie = data.end(), it = ib;
  // Skip past the comments and read the graph type.
  string_view h0, h1, h2, h3, h4;
  for (; it!=ie; it = findNextLine(it, ie)) {
    if (*ib!='%') break;
    if (data.substr(it-ib, 14)!="%%MatrixMarket") continue;
    readTokenU(h0, it, ie, fu, fw);  // %%MatrixMarket
    readTokenU(h1, it, ie, fu, fw);  // Graph
    readTokenU(h2, it, ie, fu, fw);  // Format
    readTokenU(h3, it, ie, fu, fw);  // Field
    readTokenU(h4, it, ie, fu, fw);  // Symmetry
  }
  // Check the graph type.
  if (h1!="matrix" || h2!="coordinate") return ib;  // Error: Invalid MTX format header
  symmetric = h4=="symmetric" || h4=="skew-symmetric";
  // Read rows, cols, size.
  I err = I();
  do {
    err = readNumberU<true>(rows, it, ie, fu, fw); if (err) break;  // Number of vertices
    err = readNumberU<true>(cols, it, ie, fu, fw); if (err) break;  // Number of vertices
    err = readNumberU<true>(size, it, ie, fu, fw);                  // Number of edges
  } while (false);
  if (err) return err;  // Error: Invalid MTX format header
  // Jump to the next line.
  it = findNextLine(it, ie);
  data.remove_prefix(it-ib);
  return I();
}
#pragma endregion




#pragma region READ EDGELIST FORMAT
/**
 * Read an EdgeList format file.
 * @tparam CHECK check for error?
 * @param data input file data (updated)
 * @param symmetric is graph symmetric
 * @param weighted is graph weighted
 * @param fb on body line (u, v, w)
 * @returns iterator to error in data
 */
template <bool CHECK=false, class FB>
inline auto readEdgelistFormatDoU(string_view& data, bool symmetric, bool weighted, FB fb) {
  using I = decltype(data.begin());
  auto fu = [](char c) { return c==','; };                      // Support CSV
  auto fw = [](char c) { return c==',' || c=='%' || c=='#'; };  // Support CSV, comments
  auto ib = data.begin(), ie = data.end(), it = ib;
  I err = I();
  for (; it!=ie; it = findNextLine(it, ie)) {
    // Skip past empty lines and comments.
    it = findNextNonBlank(it, ie, fu);
    if (*it=='%' || *it=='#' || isNewline(*it)) continue;
    // Read u, v, w (if weighted).
    size_t u = 0, v = 0; double w = 1;
    err = readNumberU<CHECK>(u, it, ie, fu, fw); if (CHECK && err) break;  // Source vertex
    err = readNumberU<CHECK>(v, it, ie, fu, fw); if (CHECK && err) break;  // Target vertex
    if (weighted) { err = readNumberU<CHECK>(w, it, ie, fu, fw); if (CHECK && err) break; }  // Edge weight
    fb(u, v, w);
    if (symmetric) fb(v, u, w);
  }
  if (err) return err;  // Error: Invalid EdgeList format body
  data.remove_prefix(it-ib);
  return I();
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
 * Read an EdgeList format file.
 * @tparam CHECK check for error?
 * @param data input file data (updated)
 * @param symmetric is graph symmetric
 * @param weighted is graph weighted
 * @param fb on body line (u, v, w)
 * @returns iterator to error in data
 */
template <bool CHECK=false, class FB>
inline void readEdgelistFormatDoOmpU(string_view& data, bool symmetric, bool weighted, FB fb) {
  using I = decltype(data.begin());
  const int T = omp_get_max_threads();
  const size_t DATA  = data.size();
  const size_t BLOCK = 4096;             // Characters per block (1 page)
  const size_t GRID  = 128 * T * BLOCK;  // Characters per grid (128T pages)
  I err = I();
  // Process COO file in grids.
  for (size_t g=0; g<DATA; g+=GRID) {
    size_t B = min(g+GRID, DATA);
    // Process a grid in parallel with dynamic scheduling.
    #pragma omp parallel for schedule(dynamic, 1) shared(err)
    for (size_t b=g; b<B; b+=BLOCK) {
      string_view bdata = readEdgelistFormatBlock(data, b, BLOCK);
      I berr = readEdgelistFormatDoU<CHECK>(bdata, weighted, symmetric, fb);
      if (CHECK && berr && !err) err = berr;  // Error: Invalid EdgeList format body
    }
  }
  data.remove_prefix(DATA);
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
