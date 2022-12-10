#pragma once
#include <cstdint>
#include <utility>
#include <vector>
#include <ostream>
#include <iostream>
#include "_main.hxx"

using std::pair;
using std::vector;
using std::ostream;
using std::cout;




// HELPER MACROS
// -------------
// Helps create graphs.

#ifndef GRAPH_TYPES
#define GRAPH_TYPES(K, V, E) \
  using key_type = K; \
  using vertex_key_type   = K; \
  using vertex_value_type = V; \
  using vertex_pair_type  = pair<K, V>; \
  using edge_key_type     = K; \
  using edge_value_type   = E; \
  using edge_pair_type    = pair<K, E>;

#define GRAPH_SHORT_TYPES_FROM(G) \
  using K = typename G::key_type; \
  using V = typename G::vertex_value_type; \
  using E = typename G::edge_value_type;
#endif


#ifndef GRAPH_SIZE
#define GRAPH_SIZE(K, V, E, N, M, vexists)  \
  inline size_t span()  const noexcept { return vexists.size(); } \
  inline size_t order() const noexcept { return N; } \
  inline size_t size()  const noexcept { return M; } \
  inline bool   empty() const noexcept { return N==0; }

#define GRAPH_SIZE_FROM(K, V, E, x) \
  inline size_t span()  const noexcept { return x.span(); } \
  inline size_t order() const noexcept { return x.order(); } \
  inline size_t size()  const noexcept { return x.size(); } \
  inline bool   empty() const noexcept { return x.empty(); }
#endif


#ifndef GRAPH_DIRECTED
#define GRAPH_DIRECTED(K, V, E, de) \
  inline bool directed() const noexcept { return de; }

#define GRAPH_DIRECTED_FROM(K, V, E, x) \
  inline bool directed() const noexcept { return x.directed(); }
#endif


#ifndef GRAPH_ENTRIES
#define GRAPH_VERTICES(K, V, E, vexists, vvalues) \
  inline auto vertexKeys() const noexcept { \
    auto vkeys = rangeIterable(span()); \
    return conditionalIterable(vkeys, vexists); \
  } \
  inline auto vertexValues() const noexcept { \
    return conditionalIterable(vvalues, vexists); \
  } \
  inline auto vertices() const noexcept { \
    auto vkeys = rangeIterable(span()); \
    auto pairs = pairIterable(vkeys, vvalues); \
    return conditionalIterable(pairs, vexists); \
  }

#define GRAPH_EDGES(K, V, E, eto, enone) \
  inline auto edgeKeys(K u) const noexcept { \
    return u<span()? eto[u].keys()   : enone.keys(); \
  } \
  inline auto edgeValues(K u) const noexcept { \
    return u<span()? eto[u].values() : enone.values(); \
  } \
  inline auto edges(K u) const noexcept { \
    return u<span()? eto[u].pairs()  : enone.pairs(); \
  }

#define GRAPH_INEDGES(K, V, E, efrom, enone) \
  inline auto inEdgeKeys(K v) const noexcept { \
    return v<span()? efrom[v].keys()   : enone.keys(); \
  } \
  inline auto inEdgeValues(K v) const noexcept { \
    return v<span()? efrom[v].values() : enone.values(); \
  } \
  inline auto inEdges(K v) const noexcept { \
    return v<span()? efrom[v].pairs()  : enone.pairs(); \
  }

#define GRAPH_INEDGES_SEARCH(K, V, E, eto) \
  inline auto inEdgeKeys(K v) const noexcept { \
    auto vkeys = rangeIterable(span()); \
    auto fedge = [&](K u) { return eto[u].has(v); }; \
    return filterIterable(vkeys, fedge); \
  } \
  inline auto inEdgeValues(K v) const noexcept { \
    auto fvals = [&](K u) { return eto[u].get(v); }; \
    return transformIterable(inEdgeKeys(v), fvals); \
  } \
  inline auto inEdges(K v) const noexcept { \
    return pairIterable(inEdgeKeys(v), inEdgeValues(v)); \
  }

#define GRAPH_VERTICES_FROM(K, V, E, x) \
  inline auto vertexKeys()   const noexcept { return x.vertexKeys(); } \
  inline auto vertexValues() const noexcept { return x.vertexValues(); } \
  inline auto vertices()     const noexcept { return x.vertices(); }
#define GRAPH_EDGES_FROM(K, V, E, x, name) \
  inline auto edgeKeys(K u)   const noexcept { return x.##name##Keys(u); } \
  inline auto edgeValues(K u) const noexcept { return x.##name##Values(u); } \
  inline auto edges(K u)      const noexcept { return x.##name##s(u); }
#define GRAPH_INEDGES_FROM(K, V, E, x, name) \
  inline auto inEdgeKeys(K u)   const noexcept { return x.##name##Keys(u); } \
  inline auto inEdgeValues(K u) const noexcept { return x.##name##Values(u); } \
  inline auto inEdges(K u)      const noexcept { return x.##name##s(u); }

#define GRAPH_ENTRIES(K, V, E, vexists, vvalues, eto, efrom, enone) \
  GRAPH_VERTICES(K, V, E, vexists, vvalues) \
  GRAPH_EDGES(K, V, E, eto, enone) \
  GRAPH_INEDGES(K, V, E, efrom, enone)

#define GRAPH_ENTRIES_FROM(K, V, E, x, ename, iname) \
  GRAPH_VERTICES_FROM(K, V, E, x) \
  GRAPH_EDGES_FROM(K, V, E, x, ename) \
  GRAPH_INEDGES_FROM(K, V, E, x, iname)
#endif


#ifndef GRAPH_FOREACH
#define GRAPH_FOREACH_VERTEX(K, V, E, vexists, vvalues) \
  template <class F> \
  inline void forEachVertexKey(F fn) const noexcept { \
    for (K u=0; u<span(); ++u) \
      if (vexists[u]) fn(u); \
  } \
  template <class F> \
  inline void forEachVertexValue(F fn) const noexcept { \
    for (K u=0; u<span(); ++u) \
      if (vexists[u]) fn(vvalues[u]); \
  } \
  template <class F> \
  inline void forEachVertex(F fn) const noexcept { \
    for (K u=0; u<span(); ++u) \
      if (vexists[u]) fn(u, vvalues[u]); \
  }

#define GRAPH_FOREACH_XEDGE(K, V, E, name, u, fn, eto) \
  template <class F> \
  inline void forEach##name##Key(K u, F fn) const noexcept { \
    if (u<span()) eto[u].forEachKey(fn); \
  } \
  template <class F> \
  inline void forEach##name##Value(K u, F fn) const noexcept { \
    if (u<span()) eto[u].forEachValue(fn); \
  } \
  template <class F> \
  inline void forEach##name(K u, F fn) const noexcept { \
    if (u<span()) eto[u].forEach(fn); \
  }

#define GRAPH_FOREACH_INEDGE_SEARCH(K, V, E, eto) \
  template <class F> \
  inline void forEachInEdgeKey(K v, F fn) const noexcept { \
    for (K u=0; u<span(); ++u) \
      if (eto[u].has(v)) fn(u); \
  } \
  template <class F> \
  inline void forEachInEdgeValue(K v, F fn) const noexcept { \
    for (K u=0; u<span(); ++u) \
      if (eto[u].has(v)) fn(eto[u].get(v)); \
  } \
  template <class F> \
  inline void forEachInEdge(K v, F fn) const noexcept { \
    for (K u=0; u<span(); ++u) \
      if (eto[u].has(v)) fn(u, eto[u].get(v)); \
  }

#define GRAPH_FOREACH_EDGE(K, V, E, eto) \
  GRAPH_FOREACH_XEDGE(K, V, E, Edge, u, fn, eto)
#define GRAPH_FOREACH_INEDGE(K, V, E, efrom) \
  GRAPH_FOREACH_XEDGE(K, V, E, InEdge, v, fn, efrom)

#define GRAPH_FOREACH_VERTEX_FROM(K, V, E, x) \
  template <class F> \
  inline void forEachVertexKey(F fn)   const noexcept { x.forEachVertexKey(fn); } \
  template <class F> \
  inline void forEachVertexValue(F fn) const noexcept { x.forEachVertexValue(fn); } \
  template <class F> \
  inline void forEachVertex(F fn)      const noexcept { x.forEachVertex(fn); }

#define GRAPH_FOREACH_XEDGE_FROM(K, V, E, name, x, ename) \
  template <class F> \
  inline void forEach##name##Key(K u, F fn)   const noexcept { x.forEach##ename##Key(u, fn); } \
  template <class F> \
  inline void forEach##name##Value(K u, F fn) const noexcept { x.forEach##ename##Value(u, fn); } \
  template <class F> \
  inline void forEach##name(K u, F fn)        const noexcept { x.forEach##ename(u, fn); }

#define GRAPH_FOREACH_EDGE_FROM(K, V, E, x, name) \
  GRAPH_FOREACH_XEDGE_FROM(K, V, E, Edge, x, name)
#define GRAPH_FOREACH_INEDGE_FROM(K, V, E, x, name) \
  GRAPH_FOREACH_XEDGE_FROM(K, V, E, InEdge, x, name)

#define GRAPH_FOREACH(K, V, E, vexists, vvalues, eto, efrom) \
  GRAPH_FOREACH_VERTEX(K, V, E, vexists, vvalues) \
  GRAPH_FOREACH_EDGE(K, V, E, eto) \
  GRAPH_FOREACH_INEDGE(K, V, E, efrom)
#define GRAPH_FOREACH_FROM(K, V, E, x, ename, iname) \
  GRAPH_FOREACH_VERTEX_FROM(K, V, E, x) \
  GRAPH_FOREACH_EDGE_FROM(K, V, E, x, ename) \
  GRAPH_FOREACH_INEDGE_FROM(K, V, E, x, iname)
#endif


#ifndef GRAPH_HAS
#define GRAPH_HAS(K, V, E, vexists, eto) \
  inline bool hasVertex(K u) const noexcept { \
    return u<span() && vexists[u]; \
  } \
  inline bool hasEdge(K u, K v) const noexcept { \
    return u<span() && eto[u].has(v); \
  }

#define GRAPH_HAS_FROM(K, V, E, x, u, v, ve, ee) \
  inline bool hasVertex(K u)    const noexcept { return x.ve; } \
  inline bool hasEdge(K u, K v) const noexcept { return x.ee; }
#endif


#ifndef GRAPH_DEGREES
#define GRAPH_XDEGREE(K, V, E, name, eto) \
  inline K name(K u) const noexcept { \
    return u<span()? K(eto[u].size()) : 0; \
  }
#define GRAPH_INDEGREE_SEARCH(K, V, E, eto) \
  inline K inDegree(K v) const noexcept { \
    auto fedge = [&](K u) { return eto[u].has(v); }; \
    return countIf(rangeIterable(span()), fedge); \
  }

#define GRAPH_DEGREE(K, V, E, eto) \
  GRAPH_XDEGREE(K, V, E, degree, eto)
#define GRAPH_INDEGREE(K, V, E, efrom) \
  GRAPH_XDEGREE(K, V, E, inDegree, efrom)

#define GRAPH_DEGREES(K, V, E, eto, efrom) \
  GRAPH_DEGREE(K, V, E, eto) \
  GRAPH_INDEGREE(K, V, E, efrom)
#define GRAPH_DEGREES_SEARCH(K, V, E, eto) \
  GRAPH_DEGREE(K, V, E, eto) \
  GRAPH_INDEGREE_SEARCH(K, V, E, eto)

#define GRAPH_DEGREES_FROM(K, V, E, x, u, v, de, ie) \
  inline K degree(K u)   const noexcept { return x.de; } \
  inline K inDegree(K v) const noexcept { return x.ie; }
#endif


#ifndef GRAPH_VALUES
#define GRAPH_VERTEX_VALUE(K, V, E, vvalues) \
  inline V vertexValue(K u) const noexcept { \
    return u<span()? vvalues[u] : V(); \
  }
#define GRAPH_EDGE_VALUE(K, V, E, eto) \
  inline E edgeValue(K u, K v) const noexcept { \
    return u<span()? eto[u].get(v) : E(); \
  }
#define GRAPH_VALUES(K, V, E, vvalues, eto) \
  GRAPH_VERTEX_VALUE(K, V, E, vvalues) \
  GRAPH_EDGE_VALUE(K, V, E, eto)

#define GRAPH_SET_VERTEX_VALUE(K, V, E, vvalues) \
  inline void setVertexValue(K u, V d) noexcept { \
    if (!hasVertex(u)) return; \
    vvalues[u] = d; \
  }
#define GRAPH_SET_EDGE_VALUE_X(K, V, E, e0, e1) \
  inline void setEdgeValue(K u, K v, E w) noexcept { \
    if (!hasVertex(u) || !hasVertex(v)) return; \
    e0; \
    e1; \
  }
#define GRAPH_SET_EDGE_VALUE(K, V, E, eto, efrom) \
  GRAPH_SET_EDGE_VALUE_X(K, V, E, eto[u].set(v, w), efrom[v].set(u, w))
#define GRAPH_SET_EDGE_VALUE_SEARCH(K, V, E, eto) \
  GRAPH_SET_EDGE_VALUE_X(K, V, E, eto[u].set(v, w), false)
#define GRAPH_SET_VALUES(K, V, E, vvalues, eto, efrom) \
  GRAPH_SET_VERTEX_VALUE(K, V, E, vvalues) \
  GRAPH_SET_EDGE_VALUE(K, V, E, eto, efrom)
#define GRAPH_SET_VALUES_SEARCH(K, V, E, vvalues, eto) \
  GRAPH_SET_VERTEX_VALUE(K, V, E, vvalues) \
  GRAPH_SET_EDGE_VALUE_SEARCH(K, V, E, eto)
#endif


#ifndef GRAPH_UPDATE
#define GRAPH_UPDATE(K, V, E, M, u, buf, e0, e1) \
  inline void update(K u, vector<pair<K, E>> *buf=nullptr) { \
    e0; \
    e1; \
  } \
  inline void update() { \
    M = 0; \
    vector<pair<K, E>> buf; \
    forEachVertexKey([&](K u) { \
      update(u, &buf); \
      M += degree(u); \
    }); \
  }
#endif


#ifndef GRAPH_RESPAN
#define GRAPH_RESPAN_X(K, V, E, n, vexists, vvalues, eto, extra) \
  inline void respan(size_t n) { \
    vexists.resize(n); \
    vvalues.resize(n); \
    eto.resize(n); \
    extra; \
  }
#define GRAPH_RESPAN(K, V, E, vexists, vvalues, eto, efrom) \
  GRAPH_RESPAN_X(K, V, E, n, vexists, vvalues, eto, efrom.resize(n))
#define GRAPH_RESPAN_SEARCH(K, V, E, vexists, vvalues, eto) \
  GRAPH_RESPAN_X(K, V, E, n, vexists, vvalues, eto,)
#endif


#ifndef GRAPH_CLEAR
#define GRAPH_CLEAR_X(K, V, E, N, M, vexists, vvalues, eto, extra) \
  inline void clear() noexcept { \
    N = 0; M = 0; \
    vexists.clear(); \
    vvalues.clear(); \
    eto.clear(); \
    extra; \
  }
#define GRAPH_CLEAR(K, V, E, N, M, vexists, vvalues, eto, efrom) \
  GRAPH_CLEAR_X(K, V, E, N, M, vexists, vvalues, eto, efrom.clear())
#define GRAPH_CLEAR_SEARCH(K, V, E, N, M, vexists, vvalues, eto) \
  GRAPH_CLEAR_X(K, V, E, N, M, vexists, vvalues, eto,)
#endif


#ifndef GRAPH_ADD_VERTEX
#define GRAPH_ADD_VERTEX(K, V, E, N, vexists, vvalues) \
  inline void addVertex(K u, V d=V()) { \
    if (hasVertex(u)) return; \
    if (u>=span()) respan(u+1); \
    vexists[u] = true; \
    vvalues[u] = d; \
    ++N; \
  }
#endif


#ifndef GRAPH_ADD_EDGE
#define GRAPH_ADD_EDGE_X(K, V, E, u, v, w, M, e0, e1) \
  inline void addEdge(K u, K v, E w=E()) { \
    addVertex(u); addVertex(v); \
    e0; \
    e1; \
  }
#define GRAPH_ADD_EDGE(K, V, E, M, eto, efrom) \
  GRAPH_ADD_EDGE_X(K, V, E, u, v, w, M, eto[u].add(v, w), efrom[v].add(u, w))
#define GRAPH_ADD_EDGE_SEARCH(K, V, E, M, eto) \
  GRAPH_ADD_EDGE_X(K, V, E, u, v, w, M, eto[u].add(v, w), false)
#endif


#ifndef GRAPH_REMOVE_EDGE
#define GRAPH_REMOVE_EDGE_X(K, V, E, u, v, M, e0, e1) \
  inline void removeEdge(K u, K v) { \
    if (!hasVertex(u) || !hasVertex(v)) return; \
    e0; \
    e1; \
  }

#define GRAPH_REMOVE_EDGE(K, V, E, M, eto, efrom) \
  GRAPH_REMOVE_EDGE_X(K, V, E, u, v, M, eto[u].remove(v), efrom[v].remove(u))
#define GRAPH_REMOVE_EDGE_SEARCH(K, V, E, M, eto) \
  GRAPH_REMOVE_EDGE_X(K, V, E, u, v, M, eto[u].remove(v), false)
#endif


#ifndef GRAPH_REMOVE_EDGES
#define GRAPH_REMOVE_EDGES(K, V, E, M, eto, efrom) \
  inline void removeEdges(K u) { \
    if (!hasVertex(u)) return; \
    eto[u].forEachKey([&](K v) { efrom[v].remove(u); }); \
    eto[u].clear(); \
  }

#define GRAPH_REMOVE_EDGES_SEARCH(K, V, E, M, eto) \
  inline void removeEdges(K u) { \
    if (!hasVertex(u)) return; \
    eto[u].clear(); \
  }
#endif


#ifndef GRAPH_REMOVE_INEDGES
#define GRAPH_REMOVE_INEDGES(K, V, E, M, eto, efrom) \
  inline void removeInEdges(K v) { \
    if (!hasVertex(v)) return; \
    efrom[v].forEachKey([&](K u) { eto[u].remove(v); }); \
    efrom[v].clear(); \
  }

#define GRAPH_REMOVE_INEDGES_SEARCH(K, V, E, M, eto) \
  inline void removeInEdges(K v) { \
    if (!hasVertex(v)) return; \
    forEachVertexKey([&](K u) { eto[u].remove(v); }); \
  }
#endif


#ifndef GRAPH_REMOVE_VERTEX
#define GRAPH_REMOVE_VERTEX(K, V, E, N, vexists, vvalues) \
  inline void removeVertex(K u) { \
    if (!hasVertex(u)) return; \
    removeEdges(u); \
    removeInEdges(u); \
    vexists[u] = false; \
    vvalues[u] = V(); \
    --N; \
  }
#endif


#ifndef GRAPH_WRITE
#define GRAPH_WRITE(K, V, E, Bitset, Graph) \
  template <class K, class V, class E, tclass2 Bitset> \
  inline void write(ostream& a, const Graph<K, V, E, Bitset>& x, bool detailed=false) { writeGraph(a, x, detailed); } \
  template <class K, class V, class E, tclass2 Bitset> \
  inline ostream& operator<<(ostream& a, const Graph<K, V, E, Bitset>& x) { write(a, x); return a; }

#define GRAPH_WRITE_VIEW(G, Graph) \
  template <class G> \
  inline void write(ostream& a, const Graph<G>& x, bool detailed=false) { writeGraph(a, x, detailed); } \
  template <class G> \
  inline ostream& operator<<(ostream& a, const Graph<G>& x) { write(a, x); return a; }
#endif




// DI-GRAPH
// --------
// Directed graph that memorizes in- and out-edges for each vertex.

template <class K=uint32_t, class V=NONE, class E=NONE, tclass2 Bitset=LazyBitset>
class DiGraph {
  // Data.
  protected:
  size_t N = 0, M = 0;
  vector<bool> vexists;
  vector<V>    vvalues;
  vector<Bitset<K, E>> eto;
  vector<Bitset<K, E>> efrom;
  Bitset<K, E> enone;

  // Types.
  public:
  GRAPH_TYPES(K, V, E)

  // Property operations.
  public:
  GRAPH_SIZES(K, V, E, N, M, vexists)
  GRAPH_DIRECTED(K, V, E, true)

  // Scan operations.
  public:
  GRAPH_VERTICES(K, V, E, vexists, vvalues)
  GRAPH_EDGES(K, V, E, eto, enone)
  GRAPH_INEDGES(K, V, E, efrom, enone)
  GRAPH_FOREACH_VERTEX(K, V, E, vexists, vvalues)
  GRAPH_FOREACH_EDGE(K, V, E, eto)
  GRAPH_FOREACH_INEDGE(K, V, E, efrom)

  // Access operations.
  public:
  GRAPH_HAS(K, V, E, vexists, eto)
  GRAPH_DEGREES(K, V, E, eto, efrom)
  GRAPH_VALUES(K, V, E, vvalues, eto)
  GRAPH_SET_VALUES(K, V, E, vvalues, eto, efrom)

  // Update operations.
  public:
  GRAPH_UPDATE(K, V, E, M, buf, u, eto[u].update(buf), efrom[u].update(buf))
  GRAPH_CLEAR(K, V, E, N, M, vexists, vvalues, eto, efrom)
  GRAPH_RESPAN(K, V, E, vexists, vvalues, eto, efrom)
  GRAPH_ADD_VERTEX(K, V, E, N, vexists, vvalues)
  GRAPH_ADD_EDGE(K, V, E, M, eto, efrom)
  GRAPH_REMOVE_EDGE(K, V, E, M, eto, efrom)
  GRAPH_REMOVE_EDGES(K, V, E, M, eto, efrom)
  GRAPH_REMOVE_INEDGES(K, V, E, M, eto, efrom)
  GRAPH_REMOVE_VERTEX(K, V, E, N, vexists, vvalues)
};

template <class K=uint32_t, class V=NONE, class E=NONE>
using UnorderedDiGraph = DiGraph<K, V, E, LazyBitset;




// OUT DI-GRAPH
// ------------
// Directed graph that memorizes only out-edges for each vertex.

template <class K=uint32_t, class V=NONE, class E=NONE, tclass2 Bitset=LazyBitset>
class OutDiGraph {
  // Data.
  protected:
  size_t N = 0, M = 0;
  vector<bool> vexists;
  vector<V>    vvalues;
  vector<Bitset<K, E>> eto;
  Bitset<K, E> enone;

  // Types.
  public:
  GRAPH_TYPES(K, V, E)

  // Property operations.
  public:
  GRAPH_SIZE(K, V, E, N, M, vexists)
  GRAPH_DIRECTED(K, V, E, true)

  // Scan operations.
  public:
  GRAPH_VERTICES(K, V, E, vexists, vvalues)
  GRAPH_EDGES(K, V, E, eto, enone)
  GRAPH_INEDGES_SEARCH(K, V, E, eto)
  GRAPH_FOREACH_VERTEX(K, V, E, vexists, vvalues)
  GRAPH_FOREACH_EDGE(K, V, E, eto)
  GRAPH_FOREACH_INEDGE_SEARCH(K, V, E, eto)

  // Access operations.
  public:
  GRAPH_HAS(K, V, E, vexists, eto)
  GRAPH_DEGREES_SEARCH(K, V, E, eto)
  GRAPH_VALUES(K, V, E, vvalues, eto)
  GRAPH_SET_VALUES_SEARCH(K, V, E, vvalues, eto)

  // Update operations.
  public:
  GRAPH_UPDATE(K, V, E, M, buf, u, eto[u].update(buf), false)
  GRAPH_CLEAR_SEARCH(K, V, E, N, M, vexists, vvalues, eto)
  GRAPH_RESPAN_SEARCH(K, V, E, vexists, vvalues, eto)
  GRAPH_ADD_VERTEX(K, V, E, N, vexists, vvalues)
  GRAPH_ADD_EDGE_SEARCH(K, V, E, M, eto)
  GRAPH_REMOVE_EDGE_SEARCH(K, V, E, M, eto)
  GRAPH_REMOVE_EDGES_SEARCH(K, V, E, M, eto)
  GRAPH_REMOVE_INEDGES_SEARCH(K, V, E, M, eto)
  GRAPH_REMOVE_VERTEX(K, V, E, N, vexists, vvalues)
};

template <class K=uint32_t, class V=NONE, class E=NONE>
using LazyOutDiGraph = OutDiGraph<K, V, E, LazyBitset>;




// GRAPH
// -----
// Undirected graph.

template <class K=uint32_t, class V=NONE, class E=NONE, tclass2 Bitset=LazyBitset>
class Graph : public OutDiGraph<K, V, E, Bitset> {
  using G = OutDiGraph<K, V, E, Bitset>;

  // Property operations.
  public:
  inline size_t size() const noexcept { return G::size()/2; }
  GRAPH_DIRECTED(K, V, E, false)

  // Scan operations.
  public:
  GRAPH_INEDGES_FROM(K, V, E, (*this), edge)
  GRAPH_FOREACH_INEDGE_FROM(K, V, E, (*this), Edge)

  // Access operations.
  public:
  inline K    inDegree(K v) const noexcept { return degree(v); }
  inline void setEdgeValue(K u, K v, E w) noexcept {
    G::setEdgeValue(u, v, w);
    G::setEdgeValue(v, u, w);
  }

  // Update operations.
  public:
  inline void addEdge(K u, K v, E w=E()) {
    G::addEdge(u, v, w);
    G::addEdge(v, u, w);
  }
  inline void removeEdge(K u, K v) {
    G::removeEdge(u, v);
    G::removeEdge(v, u);
  }
  inline void removeEdges(K u) {
    forEachEdgeKey(u, [&](K v) { G::removeEdge(v, u); });
    G::removeEdges(u);
  }
  inline void removeInEdges(K v) {
    removeEdges(v);
  }
};

template <class K=uint32_t, class V=NONE, class E=NONE>
using LazyGraph = Graph<K, V, E, LazyBitset>;




// RETYPE
// ------

template <class K, class V, class E, tclass2 B, class KA=K, class VA=V, class EA=E>
constexpr auto retype(const DiGraph<K, V, E, B>& x, KA _k=KA(), VA _v=VA(), EA _e=E()) {
  return DiGraph<KA, VA, EA, B>();
}
template <class K, class V, class E, tclass2 B, class KA=K, class VA=V, class EA=E>
constexpr auto retype(const OutDiGraph<K, V, E, B>& x, KA _k=KA(), VA _v=VA(), EA _e=E()) {
  return OutDiGraph<KA, VA, EA, B>();
}
template <class K, class V, class E, tclass2 B, class KA=K, class VA=V, class EA=E>
constexpr auto retype(const Graph<K, V, E, B>& x, KA _k=KA(), VA _v=VA(), EA _e=E()) {
  return Graph<KA, VA, EA, B>();
}




// WTITE
// -----

template <class G>
void writeGraphSizes(ostream& a, const G& x) {
  a << "order: " << x.order() << " size: " << x.size();
  a << (x.directed()? " [directed]" : " [undirected]") << " {}";
}

template <class G>
void writeGraphDetailed(ostream& a, const G& x) {
  a << "order: " << x.order() << " size: " << x.size();
  a << (x.directed()? " [directed]" : " [undirected]") << " {\n";
  x.forEachVertex([&](auto u, auto d) {
    a << u << ":" << d << " ->";
    x.forEachEdge(u, [&](auto v, auto w) {
      a << " " << v << ":" << w;
    });
    a << "\n";
  });
  a << "}";
}

template <class G>
inline void writeGraph(ostream& a, const G& x, bool detailed=false) {
  if (detailed) writeGraphDetailed(a, x);
  else writeGraphSizes(a, x);
}

GRAPH_WRITE(K, V, E, Bitset, DiGraph)
GRAPH_WRITE(K, V, E, Bitset, OutDiGraph)
GRAPH_WRITE(K, V, E, Bitset, Graph)
