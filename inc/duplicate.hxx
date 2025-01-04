#pragma once
#include "update.hxx"




#pragma region METHODS
#pragma region DUPLICATE IF
/**
 * Duplicate vertices/edges of a graph if test passes.
 * @param a output graph (output)
 * @param x input graph
 * @param fv include vertex? (u, d)
 * @param fe include edge? (u, v, w)
 */
template <class H, class G, class FV, class FE>
inline void duplicateIfW(H& a, const G& x, FV fv, FE fe) {
  size_t S = x.span();
  // Delete existing data.
  a.clear();
  // Add vertices and reserve space for edges.
  a.reserve(S);
  x.forEachVertex([&](auto u, auto d) {
    if (fv(u, d)) a.addVertex(u, d);
    a.reserveEdges(u, x.degree(u));
  });
  // Populate the edges.
  x.forEachVertexKey([&](auto u) {
    x.forEachEdge(u, [&](auto v, auto w) {
      if (fe(u, v, w)) a.addEdge(u, v, w);
    });
  });
  a.update();
}

/**
 * Duplicate vertices/edges of a graph if test passes.
 * @param x a graph
 * @param fv include vertex? (u, d)
 * @param fe include edge? (u, v, w)
 * @returns duplicate of graph
 */
template <class G, class FV, class FE>
inline G duplicateIf(const G& x, FV fv, FE fe) {
  G a; duplicateIfW(a, x, fv, fe);
  return a;
}


#ifdef OPENMP
/**
 * Duplicate vertices/edges of a graph if test passes.
 * @param a output graph (output)
 * @param x input graph
 * @param fv include vertex? (u, d)
 * @param fe include edge? (u, v, w)
 */
template <class H, class G, class FV, class FE>
inline void duplicateIfOmpW(H& a, const G& x, FV fv, FE fe) {
  using  K = typename G::key_type;
  size_t S = x.span();
  // Delete existing data.
  a.clear();
  // Add vertices and reserve space for edges.
  a.reserve(S);
  #pragma omp parallel for schedule(dynamic, 2048)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    a.addVertex(u, x.vertexValue(u));
    a.reserveEdges(u, x.degree(u));
  }
  // Populate the edges.
  #pragma omp parallel for schedule(dynamic, 2048)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    x.forEachEdge(u, [&](auto v, auto w) {
      if (fe(u, v, w)) a.addEdge(u, v, w);
    });
  }
  updateOmpU(a);
}

/**
 * Duplicate vertices/edges of a graph if test passes.
 * @param x a graph
 * @param fv include vertex? (u, d)
 * @param fe include edge? (u, v, w)
 * @returns duplicate of graph
 */
template <class G, class FV, class FE>
inline G duplicateIfOmp(const G& x, FV fv, FE fe) {
  G a; duplicateIfOmpW(a, x, fv, fe);
  return a;
}
#endif
#pragma endregion




#pragma region DUPLICATE
/**
 * Duplicate vertices/edges of a graph.
 * @param a output graph (updated)
 * @param x input graph
 */
template <class H, class G>
inline void duplicateW(H& a, const G& x) {
  auto fv = [](auto u, auto d)         { return true; };
  auto fe = [](auto u, auto v, auto w) { return true; };
  duplicateIfW(a, x, fv, fe);
}

/**
 * Duplicate a graph.
 * @param x a graph
 * @returns duplicate of graph
 */
template <class G>
inline G duplicate(const G& x) {
  G a = x;  // Just use the copy constructor.
  return a;
}


#ifdef OPENMP
/**
 * Duplicate vertices/edges of a graph.
 * @param a output graph (updated)
 * @param x input graph
 */
template <class H, class G>
inline void duplicateOmpW(H& a, const G& x) {
  auto fv = [](auto u, auto d)         { return true; };
  auto fe = [](auto u, auto v, auto w) { return true; };
  duplicateIfOmpW(a, x, fv, fe);
}

/**
 * Duplicate a graph.
 * @param x a graph
 * @returns duplicate of graph
 */
template <class G>
inline G duplicateOmp(const G& x) {
  G a; duplicateOmpW(a, x);
  return a;
}
#endif
#pragma endregion
#pragma endregion
