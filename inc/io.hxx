#pragma once
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <vector>
#include <algorithm>
#include "_main.hxx"
#include "Graph.hxx"
#include "update.hxx"
#ifdef OPENMP
#include <omp.h>
#endif

using std::unique_ptr;
using std::remove_reference_t;
using std::string;
using std::string_view;
using std::vector;
using std::min;




#pragma region READ FROM STRING
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
 * Read header of a MTX format file (check for errors, handle comments).
 * @param symmetric is graph symmetric (output)
 * @param rows number of rows (output)
 * @param cols number of columns (output)
 * @param size number of lines/edges (output)
 * @param data input file data
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
 * Read a file in Edgelist format, using mmap and sscanf.
 * @tparam WEIGHTED is graph weighted?
 * @tparam BASE base vertex id (0 or 1)
 * @param data input file data
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
 * Read an EdgeList format file (crazy frog version).
 * @tparam WEIGHTED is graph weighted?
 * @tparam BASE base vertex id (0 or 1)
 * @param data input file data (updated)
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
 * Read an EdgeList format file.
 * @tparam WEIGHTED is graph weighted?
 * @tparam BASE base vertex id (0 or 1)
 * @tparam CHECK check for error?
 * @param data input file data (updated)
 * @param symmetric is graph symmetric?
 * @param fb on body line (u, v, w)
 */
template <bool WEIGHTED=false, int BASE=1, bool CHECK=false, class FB>
inline void readEdgelistFormatDo(string_view data, bool symmetric, FB fb) {
  if constexpr (CHECK) readEdgelistFormatDoChecked<WEIGHTED, BASE>(data, symmetric, fb);
  else readEdgelistFormatDoUnchecked<WEIGHTED, BASE>(data, symmetric, fb);
}


/**
 * Read a file in Edgelist format, and record the edges.
 * @tparam WEIGHTED is graph weighted?
 * @tparam BASE base vertex id (0 or 1)
 * @tparam CHECK check for error?
 * @param degrees vertex degrees (updated)
 * @param sources source vertices (output)
 * @param targets target vertices (output)
 * @param weights edge weights (output)
 * @param data input file data
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
inline vector<size_t> readEdgelistFormatToListsOmpU(IIK degrees, IIK sources, IIK targets, IIE weights, string_view data, bool symmetric) {
  const size_t DATA  = data.size();
  const size_t BLOCK = 256 * 1024;  // Characters per block (256KB)
  const int T = omp_get_max_threads();
  FormatError err;       // Common error
  vector<size_t> is(T);  // Per-thread index
  // Process a grid in parallel with dynamic scheduling.
  #pragma omp parallel shared(err)
  {
    int    t = omp_get_thread_num();
    size_t i = 0;
    #pragma omp for schedule(dynamic) nowait
    for (size_t b=0; b<DATA; b+=BLOCK) {
      if (CHECK && !err.empty()) continue;
      // Read a block of data, and process it.
      string_view bdata = readEdgelistFormatBlock(data, b, BLOCK);
      auto fb = [&](auto u, auto v, auto w) {
        const int p = t % PARTITIONS;
        sources[t][i] = u;
        targets[t][i] = v;
        if constexpr (WEIGHTED) weights[t][i] = w;
        #pragma omp atomic
        ++degrees[p][u];
        ++i;
      };
      if constexpr (CHECK) {
        try { readEdgelistFormatDo<WEIGHTED, BASE, true>(bdata, symmetric, fb); }
        catch (const FormatError& e) { if (err.empty()) err = e; }
      }
      else readEdgelistFormatDo<WEIGHTED, BASE>(bdata, symmetric, fb);
    }
    // Update per-thread index.
    is[t] = i;
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
 * @returns number of edges in the Edgelist
 */
template <bool WEIGHTED=false, class IO, class IK, class IE>
inline size_t convertEdgelistToCsrListsW(IO offsets, IK edgeKeys, IE edgeValues, IK degrees, IK sources, IK targets, IE weights, size_t rows) {
  using O = remove_reference_t<decltype(offsets[0])>;
  // Compute shifted offsets.
  offsets[0] = O();
  size_t   M = exclusiveScanW(offsets+1, degrees, rows);
  // Populate CSR.
  for (size_t i=0; i<M; ++i) {
    size_t u = sources[i];
    size_t v = targets[i];
    size_t j = offsets[u+1]++;
    edgeKeys[j] = v;
    if constexpr (WEIGHTED) edgeValues[j] = weights[i];
  }
  return M;
}


/**
 * Convert Edgelist to CSR (lists).
 * @tparam WEIGHTED is graph weighted?
 * @tparam PARTITIONS number of partitions for vertex degrees
 * @param offsets per-partition CSR offsets (output)
 * @param edgeKeys per-partition CSR edge keys (output)
 * @param edgeValues per-partition CSR edge values (output)
 * @param degrees per-partition vertex degrees
 * @param sources per-thread source vertices
 * @param targets per-thread target vertices
 * @param weights per-thread edge weights
 * @param counts per-thread number of edges read
 * @param rows number of rows/vertices
 * @returns number of edges in the Edgelist
 * @note offsets[0], edgeKeys[0], and edgeValues[0] are special.
 * @note They are used to store global offsets, edge keys, and edge values.
 */
template <bool WEIGHTED=false, int PARTITIONS=4, class IIO, class IIK, class IIE>
inline size_t convertEdgelistToCsrListsOmpW(IIO offsets, IIK edgeKeys, IIE edgeValues, IIK degrees, IIK sources, IIK targets, IIE weights, const vector<size_t>& counts, size_t rows) {
  using  O = remove_reference_t<decltype(offsets[0][0])>;
  int    T = omp_get_max_threads();
  size_t M = 0;
  vector<size_t> buf(T);
  if (PARTITIONS==1) {
    // Compute shifted global offsets at offsets[0].
    offsets[0][0] = O();
    M = exclusiveScanOmpW(offsets[0]+1, buf.data(), degrees[0], rows);
    // Populate global CSR at edgeKeys[0] and edgeValues[0].
    #pragma omp parallel
    {
      int t = omp_get_thread_num();
      size_t I = counts[t];
      for (size_t i=0; i<I; ++i) {
        size_t u = sources[t][i];
        size_t v = targets[t][i];
        size_t j = 0;
        #pragma omp atomic capture
        j = offsets[0][u+1]++;
        edgeKeys[0][j] = v;
        if constexpr (WEIGHTED) edgeValues[0][j] = weights[t][i];
      }
    }
    return M;
  }
  // Compute global degrees at degrees[0].
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
  // Compute per-partition shifted offsets at offsets[p].
  for (int p=0; p<PARTITIONS; ++p) {
    offsets[p][0] = O();
    size_t MP = exclusiveScanOmpW(offsets[p]+1, buf.data(), degrees[p], rows);
    if (p==0) M = MP;
  }
  // Populate per-partition CSR at edgeKeys[p] and edgeValues[p].
  #pragma omp parallel
  {
    int t = omp_get_thread_num();
    size_t I = counts[t];
    for (size_t i=0; i<I; ++i) {
      const int p = t % PARTITIONS;
      size_t u = sources[t][i];
      size_t v = targets[t][i];
      size_t j = 0;
      #pragma omp atomic capture
      j = offsets[p][u+1]++;
      edgeKeys[p][j] = v;
      if constexpr (WEIGHTED) edgeValues[p][j] = weights[t][i];
    }
  }
  // Populate global CSR at edgeKeys[0] and edgeValues[0], from per-partition CSRs.
  #pragma omp parallel for schedule(dynamic, 2048)
  for (size_t u=0; u<rows; ++u) {
    size_t j = offsets[0][u+1];
    for (int p=1; p<PARTITIONS; ++p) {
      size_t i = offsets[p][u];
      size_t I = offsets[p][u+1];
      for (; i<I; ++i, ++j) {
        edgeKeys[0][j] = edgeKeys[p][i];
        if constexpr (WEIGHTED) edgeValues[0][j] = edgeValues[p][i];
      }
    }
    offsets[0][u+1] = j;
  }
  return M;
}
#pragma endregion




#pragma region READ MTX FORMAT TO CSR
/**
 * Read data in MTX format, and convert to CSR.
 * @tparam WEIGHTED is graph weighted?
 * @tparam BASE base vertex id (0 or 1)
 * @tparam CHECK check for error?
 * @param a output csr graph (updated)
 * @param data input data
 */
template <bool WEIGHTED=false, int BASE=1, bool CHECK=false, class G>
inline void readMtxFormatToCsrW(G& a, string_view data) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  // Read MTX format header.
  bool symmetric; size_t rows, cols, size;
  size_t head = readMtxFormatHeaderW(symmetric, rows, cols, size, data);
  data.remove_prefix(head);
  // Allocate space for CSR.
  const size_t N = max(rows, cols);
  const size_t M = symmetric? 2 * size : size;
  a.resize(N, M);
  // Allocate space for sources, targets, and weights.
  K *sources = new K[M];
  K *targets = new K[M];
  E *weights = WEIGHTED? new E[M] : nullptr;
  K *degrees    = a.degrees;
  K *offsets    = a.offsets;
  K *edgeKeys   = a.edgeKeys;
  E *edgeValues = WEIGHTED? a.edgeValues : nullptr;
  fillValueU(degrees, N, K());
  // Read Edgelist and convert to CSR.
  readEdgelistFormatToListsU<WEIGHTED, BASE, CHECK>(degrees, sources, targets, weights, data, symmetric);
  size_t MA = convertEdgelistToCsrListsW<WEIGHTED>(offsets, edgeKeys, edgeValues, degrees, sources, targets, weights, N);
  if (symmetric) a.resize(N, MA);
  // Free space for sources, targets, and weights.
  delete sources;
  delete targets;
  if (WEIGHTED) delete weights;
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
inline void readMtxFormatToCsrOmpW(G& a, string_view data) {
  using O = typename G::offset_type;
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  // Read MTX format header.
  bool symmetric; size_t rows, cols, size;
  size_t head = readMtxFormatHeaderW(symmetric, rows, cols, size, data);
  data.remove_prefix(head);
  auto t0 = timeNow();
  // Allocate space for CSR.
  const size_t N = max(rows, cols);
  const size_t M = symmetric? 2 * size : size;
  a.resize(N, M);
  // Allocate space for sources, targets, weights.
  const int T = omp_get_max_threads();
  vector<K*> sources(T);
  vector<K*> targets(T);
  vector<E*> weights(T);
  for (int i=0; i<T; ++i) {
    sources[i] = new K[M];
    targets[i] = new K[M];
    weights[i] = WEIGHTED? new E[M] : nullptr;
  }
  // Allocate space for degrees, offsets, edge keys, and edge values
  // Note that degrees[0], offsets[0], edgeKeys[0], and edgeValues[0] are special.
  // They point to the global degrees, offsets, edge keys, and edge values in the CSR.
  vector<K*> degrees(PARTITIONS);
  vector<O*> offsets(PARTITIONS);
  vector<K*> edgeKeys(PARTITIONS);
  vector<E*> edgeValues(PARTITIONS);
  degrees[0]    = a.degrees;  // NOTE: Assuming that a.degrees is zero-initialized.
  offsets[0]    = a.offsets;
  edgeKeys[0]   = a.edgeKeys;
  edgeValues[0] = WEIGHTED? a.edgeValues : nullptr;
  fillValueOmpU(degrees[0], N, K());
  for (int i=1; i<PARTITIONS; ++i) {
    degrees[i]    = new K[N+1];
    offsets[i]    = new O[N+1];
    edgeKeys[i]   = new K[M];
    edgeValues[i] = WEIGHTED? new E[M] : nullptr;
    fillValueOmpU(degrees[i], N+1, K());
  }
  auto t1 = timeNow();
  // Read Edgelist and convert to CSR.
  vector<size_t> counts = readEdgelistFormatToListsOmpU<WEIGHTED, BASE, CHECK, PARTITIONS>(degrees, sources, targets, weights, data, symmetric);
  auto t2 = timeNow();
  size_t MA = convertEdgelistToCsrListsOmpW<WEIGHTED, PARTITIONS>(offsets, edgeKeys, edgeValues, degrees, sources, targets, weights, counts, N);
  auto t3 = timeNow();
  if (symmetric) a.resize(N, MA);
  // Free space for sources, targets, and weights.
  for (int i=0; i<T; ++i) {
    delete sources[i];
    delete targets[i];
    if (WEIGHTED) delete weights[i];
  }
  // Free space for degrees, offsets, edge keys, and edge values.
  for (int i=1; i<PARTITIONS; ++i) {
    delete degrees[i];
    delete offsets[i];
    delete edgeKeys[i];
    if (WEIGHTED) delete edgeValues[i];
  }
  auto t4 = timeNow();
  // Print time taken.
  printf("readMtxFormatToCsrOmpW: {%09.3fms} Allocate memory\n", duration(t0, t1));
  printf("readMtxFormatToCsrOmpW: {%09.3fms} Read Edgelist\n", duration(t1, t2));
  printf("readMtxFormatToCsrOmpW: {%09.3fms} Convert to CSR\n", duration(t2, t3));
  printf("readMtxFormatToCsrOmpW: {%09.3fms} Free memory\n", duration(t3, t4));
}
#pragma endregion




#pragma region CONVERT EDGELIST TO GRAPH
/**
 * Convert Edgelist to Graph (Arena-allocator based).
 * @tparam WEIGHTED is graph weighted?
 * @param a output graph (output)
 * @param degrees vertex degrees
 * @param sources source vertices
 * @param targets target vertices
 * @param weights edge weights
 * @param rows number of rows/vertices
 * @returns number of edges in the Edgelist
 */
template <bool WEIGHTED=false, class G, class IO, class IK, class IE>
inline size_t convertEdgelistToGraphW(G& a, IK degrees, IK sources, IK targets, IE weights, size_t rows) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  a.clear();
  a.reserve(rows);
  // Add vertices.
  for (K u=0; u<rows; ++u)
    a.addVertex(u);
  // Allocate space for edges.
  size_t M = 0;
  for (K u=0; u<rows; ++u) {
    a.allocateEdges(u, degrees[u]);
    M += degrees[u];
  }
  // Populate edges.
  for (size_t i=0; i<M; ++i)
    a.addEdge(sources[i], targets[i], WEIGHTED? weights[i] : E(1));
  return M;
}


/**
 * Convert Edgelist to Graph (Arena-allocator based).
 * @tparam WEIGHTED is graph weighted?
 * @tparam PARTITIONS number of partitions for vertex degrees
 * @param a output graph (output)
 * @param offsets per-partition CSR offsets (scratch)
 * @param edgeKeys per-partition CSR edge keys (scratch)
 * @param edgeValues per-partition CSR edge values (scratch)
 * @param degrees per-partition vertex degrees (updated)
 * @param sources per-thread source vertices
 * @param targets per-thread target vertices
 * @param weights per-thread edge weights
 * @param counts per-thread number of edges read
 * @param rows number of rows/vertices
 * @returns number of edges in the Edgelist
 * @note offsets[0], edgeKeys[0], and edgeValues[0] are special.
 * @note They are used to store global offsets, edge keys, and edge values.
 */
template <bool WEIGHTED=false, int PARTITIONS=4, class G, class IIO, class IIK, class IIE>
inline size_t convertEdgelistToGraphOmpW(G &a, IIO offsets, IIK edgeKeys, IIE edgeValues, IIK degrees, IIK sources, IIK targets, IIE weights, const vector<size_t>& counts, size_t rows) {
  using  O = remove_reference_t<decltype(offsets[0][0])>;
  using  K = typename G::key_type;
  using  E = typename G::edge_value_type;
  int    T = omp_get_max_threads();
  size_t M = 0;
  vector<size_t> buf(T);
  // Compute per-partition shifted offsets at offsets[p].
  for (int p=0; p<PARTITIONS; ++p) {
    offsets[p][0] = O();
    size_t MP = exclusiveScanOmpW(offsets[p]+1, buf.data(), degrees[p], rows);
    M += MP;
  }
  // Populate per-partition CSR at edgeKeys[p] and edgeValues[p].
  #pragma omp parallel
  {
    int t = omp_get_thread_num();
    size_t I = counts[t];
    for (size_t i=0; i<I; ++i) {
      const int p = t % PARTITIONS;
      size_t u = sources[t][i];
      size_t v = targets[t][i];
      size_t j = 0;
      #pragma omp atomic capture
      j = offsets[p][u+1]++;
      edgeKeys[p][j] = v;
      if constexpr (WEIGHTED) edgeValues[p][j] = weights[t][i];
    }
  }
  auto fdeg = [&](auto u) {
    if (PARTITIONS==1) return degrees[0][u];
    else if (PARTITIONS==2) return degrees[0][u] + degrees[1][u];
    else if (PARTITIONS==4) return degrees[0][u] + degrees[1][u] + degrees[2][u] + degrees[3][u];
    else if (PARTITIONS==8) return degrees[0][u] + degrees[1][u] + degrees[2][u] + degrees[3][u] + degrees[4][u] + degrees[5][u] + degrees[6][u] + degrees[7][u];
    else {
      K d = K();
      for (int t=0; t<PARTITIONS; ++t)
        d += degrees[t][u];
      return d;
    }
  };
  a.clearOmp();
  a.reserveOmp(rows);
  // Allocate space for vertices.
  #pragma omp parallel for schedule(static, 2048)
  for (size_t u=0; u<rows; ++u)
    a.addVertex(u);
  // Allocate space for edges.
  #pragma omp parallel for schedule(dynamic, 2048)
  for (size_t u=0; u<rows; ++u)
    a.allocateEdges(u, fdeg(u));
  // Populate edges.
  #pragma omp parallel for schedule(dynamic, 1024)
  for (size_t u=0; u<rows; ++u) {
    for (int p=0; p<PARTITIONS; ++p) {
      size_t i = offsets[p][u];
      size_t I = offsets[p][u+1];
      for (; i<I; ++i)
        a.addEdgeUnsafe(u, edgeKeys[p][i], WEIGHTED? edgeValues[p][i] : E(1));
    }
  }
  a.updateOmp(true, false);
  return M;
}
#pragma endregion




#pragma region READ MTX FORMAT TO GRAPH
/**
 * Read data in MTX format, and convert to Graph (Arena-allocator based).
 * @tparam WEIGHTED is graph weighted?
 * @tparam BASE base vertex id (0 or 1)
 * @tparam CHECK check for error?
 * @param a output csr graph (updated)
 * @param data input data
 */
template <bool WEIGHTED=false, int BASE=1, bool CHECK=false, class G>
inline void readMtxFormatToGraphW(G& a, string_view data) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  // Read MTX format header.
  bool symmetric; size_t rows, cols, size;
  size_t head = readMtxFormatHeaderW(symmetric, rows, cols, size, data);
  data.remove_prefix(head);
  // Allocate space for sources, targets, and weights.
  size_t N = max(rows, cols);
  size_t M = symmetric? 2 * size : size;
  K *sources = new K[M];
  K *targets = new K[M];
  E *weights = WEIGHTED? new E[M] : nullptr;
  K *degrees = new K[N+1];
  fillValueU(degrees, N+1, K());
  // Read Edgelist and convert to CSR.
  readEdgelistFormatToListsU<WEIGHTED, BASE, CHECK>(degrees, sources, targets, weights, data, symmetric);
  size_t MA = convertEdgelistToGraphW<WEIGHTED>(a, degrees, sources, targets, weights, N);
  // Free space for sources, targets, and weights.
  delete sources;
  delete targets;
  if (WEIGHTED) delete weights;
  delete degrees;
}


/**
 * Read data in MTX format, and convert to Graph (Arena-allocator based).
 * @tparam WEIGHTED is graph weighted?
 * @tparam BASE base vertex id (0 or 1)
 * @tparam CHECK check for error?
 * @tparam PARTITIONS number of partitions for vertex degrees
 * @param a output csr graph (updated)
 * @param data input data
 */
template <bool WEIGHTED=false, int BASE=1, bool CHECK=false, int PARTITIONS=4, class G>
inline void readMtxFormatToGraphOmpW(G& a, string_view data) {
  using O = typename G::offset_type;
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  // Read MTX format header.
  bool symmetric; size_t rows, cols, size;
  size_t head = readMtxFormatHeaderW(symmetric, rows, cols, size, data);
  data.remove_prefix(head);
  // Allocate space for sources, targets, weights.
  const int T = omp_get_max_threads();
  size_t N = max(rows, cols);
  size_t M = symmetric? 2 * size : size;
  vector<K*> sources(T);
  vector<K*> targets(T);
  vector<E*> weights(T);
  for (int i=0; i<T; ++i) {
    sources[i] = new K[M];
    targets[i] = new K[M];
    weights[i] = WEIGHTED? new E[M] : nullptr;
  }
  // Allocate space for degrees, offsets, edge keys, and edge values
  // Note that degrees[0], offsets[0], edgeKeys[0], and edgeValues[0] are special.
  // They point to the global degrees, offsets, edge keys, and edge values in the CSR.
  vector<K*> degrees(PARTITIONS);
  vector<O*> offsets(PARTITIONS);
  vector<K*> edgeKeys(PARTITIONS);
  vector<E*> edgeValues(PARTITIONS);
  for (int i=0; i<PARTITIONS; ++i) {
    degrees[i]    = new K[N+1];
    offsets[i]    = new O[N+1];
    edgeKeys[i]   = new K[M];
    edgeValues[i] = WEIGHTED? new E[M] : nullptr;
    fillValueOmpU(degrees[i], N+1, K());
  }
  // Read Edgelist and convert to CSR.
  vector<size_t> counts = readEdgelistFormatToListsOmpU<WEIGHTED, BASE, CHECK, PARTITIONS>(degrees, sources, targets, weights, data, symmetric);
  size_t MA = convertEdgelistToGraphOmpW<WEIGHTED, PARTITIONS>(a, offsets, edgeKeys, edgeValues, degrees, sources, targets, weights, counts, N);
  // Free space for sources, targets, and weights.
  for (int i=0; i<T; ++i) {
    delete sources[i];
    delete targets[i];
    if (WEIGHTED) delete weights[i];
  }
  // Free space for degrees, offsets, edge keys, and edge values.
  for (int i=0; i<PARTITIONS; ++i) {
    delete degrees[i];
    delete offsets[i];
    delete edgeKeys[i];
    if (WEIGHTED) delete edgeValues[i];
  }
}
#pragma endregion
#pragma endregion
