#pragma once
#include <cstring>
#include <utility>
#include <iterator>
#include <vector>
#include <ostream>
#include <algorithm>
#include "_main.hxx"
#ifdef OPENMP
#include <omp.h>
#endif

using std::pair;
using std::vector;
using std::ostream;
using std::memcpy;
using std::distance;
using std::max;
using std::find_if;
using std::lower_bound;
using std::unique;
using std::sort;
using std::set_union;




#pragma region CLASSES
/**
 * Directed graph that memorizes only out-edges for each vertex (and uses arena allocation).
 * @tparam K key type (vertex id)
 * @tparam V vertex value type (vertex data)
 * @tparam E edge value type (edge weight)
 */
template <class K=uint32_t, class V=None, class E=None>
class ArenaDiGraph {
  #pragma region TYPES
  protected:
  /** Memory allocator for the graph. */
  using Allocator = ConcurrentPow2Allocator<>;
  public:
  /** Key type (vertex id). */
  using key_type = K;
  /** Vertex value type (vertex data). */
  using vertex_value_type = V;
  /** Edge value type (edge weight). */
  using edge_value_type = E;
  #pragma endregion


  #pragma region CONSTANTS
  protected:
  /** Threshold for using binary search. */
  static constexpr size_t BSEARCH = 64;
  /** Size of each edge in bytes. */
  static constexpr size_t EDGE = sizeof(pair<K, E>);
  #pragma endregion


  #pragma region DATA
  protected:
  /** Number of vertices. */
  size_t N = 0;
  /** Number of edges. */
  size_t M = 0;
  /** Vertex existence flags. */
  vector<bool> exists;
  /** Vertex values. */
  vector<V> values;
  /** Out-degree of each vertex. */
  vector<K> degrees;
  /** Edge capacity of each vertex. */
  vector<K> capacities;
  /** Outgoing edges for each vertex (including edge weights). */
  vector<pair<K, E>*> edges;
  /** Memory allocator. */
  Allocator *mx;
  #pragma endregion


  #pragma region METHODS
  #pragma region PROPERTIES
  public:
  /**
   * Get the size of buffer required to store data associated with each vertex
   * in the graph, indexed by its vertex-id.
   * @returns size of buffer required
   */
  inline size_t span() const noexcept {
    return exists.size();
  }

  /**
   * Get the number of vertices in the graph.
   * @returns |V|
   */
  inline size_t order() const noexcept {
    return N;
  }

  /**
   * Get the number of edges in the graph.
   * @returns |E|
   */
  inline size_t size() const noexcept {
    return M;
  }

  /**
   * Check if the graph is empty.
   * @returns is the graph empty?
   */
  inline bool empty() const noexcept {
    return N == 0;
  }

  /**
   * Check if the graph is directed.
   * @returns is the graph directed?
   */
  inline bool directed() const noexcept {
    return true;
  }
  #pragma endregion


  #pragma region FOREACH
  public:
  /**
   * Iterate over the vertices in the graph.
   * @param fp process function (vertex id, vertex data)
   */
  template <class FP>
  inline void forEachVertex(FP fp) const noexcept {
    for (K u=0; u<span(); ++u)
      if (exists[u]) fp(u, values[u]);
  }

  /**
   * Iterate over the vertex ids in the graph.
   * @param fp process function (vertex id)
   */
  template <class FP>
  inline void forEachVertexKey(FP fp) const noexcept {
    for (K u=0; u<span(); ++u)
      if (exists[u]) fp(u);
  }

  /**
   * Iterate over the outgoing edges of a source vertex in the graph.
   * @param u source vertex id
   * @param fp process function (target vertex id, edge weight)
   */
  template <class FP>
  inline void forEachEdge(K u, FP fp) const noexcept {
    K d = degrees[u];
    auto ib = edges[u], ie = edges[u] + d;
    for (auto it=ib; it!=ie; ++it)
      fp((*it).first, (*it).second);
  }

  /**
   * Iterate over the target vertex ids of a source vertex in the graph.
   * @param u source vertex id
   * @param fp process function (target vertex id)
   */
  template <class FP>
  inline void forEachEdgeKey(K u, FP fp) const noexcept {
    K d = degrees[u];
    auto ib = edges[u], ie = edges[u] + d;
    for (auto it=ib; it!=ie; ++it)
      fp((*it).first);
  }
  #pragma endregion


  #pragma region ITERATORS
  public:
  /**
   * Get an iterator to the begin of edges of a vertex.
   * @param u vertex id
   * @returns begin iterator of edges
   */
  inline const pair<K, V>* beginEdges(K u) const noexcept {
    return edges[u];
  }


  /**
   * Get an iterator to the end of edges of a vertex.
   * @param u vertex id
   * @returns end iterator of edges
   */
  inline const pair<K, V>* endEdges(K u) const noexcept {
    return edges[u] + degrees[u];
  }


  /**
   * Get an iterator to the begin of edges of a vertex.
   * @param u vertex id
   * @returns begin iterator of edges
   */
  inline pair<K, V>* beginEdges(K u) noexcept {
    return edges[u];
  }


  /**
   * Get an iterator to the end of edges of a vertex.
   * @param u vertex id
   * @returns end iterator of edges
   */
  inline pair<K, V>* endEdges(K u) noexcept {
    return edges[u] + degrees[u];
  }
  #pragma endregion


  #pragma region ACCESS
  protected:
  /**
   * Find an entry in a range of pairs.
   * @param ib begin iterator
   * @param ie end iterator
   * @param k key to find
   * @returns iterator to the entry, or end iterator if not found
   */
  template <class I>
  static inline auto findEntry(I ib, I ie, K k) noexcept {
    auto fe = [&](const auto& p)      { return p.first == k; };
    auto fl = [ ](const auto& p, K k) { return p.first  < k; };
    // Use linear search if the range is small.
    if (distance(ib, ie) < BSEARCH) return find_if(ib, ie, fe);
    // Else, use binary search.
    auto it = lower_bound(ib, ie, k, fl);
    return it==ie || (*it).first!=k? ie : it;
  }

  public:
  /**
   * Check if a vertex exists in the graph.
   * @param u vertex id
   * @returns does the vertex exist?
   */
  inline bool hasVertex(K u) const noexcept {
    return u < span() && exists[u];
  }

  /**
   * Check if an edge exists in the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @returns does the edge exist?
   */
  inline bool hasEdge(K u, K v) const noexcept {
    if (u >= span()) return false;
    auto ib = edges[u], ie = edges[u] + degrees[u];
    return findEntry(ib, ie, v) != ie;
  }

  /**
   * Get the number of outgoing edges of a vertex in the graph.
   * @param u vertex id
   * @returns number of outgoing edges of the vertex
   */
  inline K degree(K u) const noexcept {
    return u < span()? degrees[u] : 0;
  }

  /**
   * Get the capacity of a vertex in the graph.
   * @param u vertex id
   * @returns current edge capacity of the vertex
   */
  inline K capacity(K u) const noexcept {
    return u < span()? capacities[u] : 0;
  }

  /**
   * Get the vertex data of a vertex in the graph.
   * @param u vertex id
   * @returns associated data of the vertex
   */
  inline V vertexValue(K u) const noexcept {
    return u < span()? values[u] : V();
  }

  /**
   * Set the vertex data of a vertex in the graph.
   * @param u vertex id
   * @param d associated data of the vertex
   * @returns success?
   */
  inline bool setVertexValue(K u, V d) noexcept {
    if (!hasVertex(u)) return false;
    values[u] = d;
    return true;
  }

  /**
   * Get the edge weight of an edge in the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @returns associated weight of the edge
   */
  inline E edgeValue(K u, K v) const noexcept {
    if (u >= span()) return E();
    auto ib = edges[u], ie = edges[u] + degrees[u];
    auto it = findEntry(ib, ie, v);
    return it!=ie? (*it).second : E();
  }

  /**
   * Set the edge weight of an edge in the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @param w associated weight of the edge
   * @returns success?
   */
  inline bool setEdgeValue(K u, K v, E w) noexcept {
    if (!hasVertex(u) || !hasVertex(v)) return false;
    auto ib = edges[u], ie = edges[u] + degrees[u];
    auto it = findEntry(ib, ie, v);
    if (it == ie) return false;
    (*it).second = w;
  }
  #pragma endregion


  #pragma region UPDATE
  protected:
  /**
   * Get the allocation capacity for a number of elements.
   * @param n number of elements
   * @returns allocation capacity
   */
  static inline constexpr K allocationCapacity(K n) noexcept {
    return Allocator::allocationCapacity(n*EDGE) / EDGE;
  }


  /**
   * Allocate memory for a number of edges.
   * @param c allocation capacity
   * @returns pointer to the allocated memory
   */
  inline void* allocate(K c) {
    return mx->allocate(c*EDGE);
  }


  /**
   * Deallocate memory for a number of edges.
   * @param ptr pointer to the memory
   * @param c allocation capacity
   */
  inline void deallocate(void *ptr, K c) {
    mx->deallocate(ptr, c*EDGE);
  }


  /**
   * Reset all arena allocators.
   */
  inline void resetAllocators() {
    mx->reset();
  }


  /**
   * Resize arrays to specified size.
   * @param n new size
   */
  inline void resizeArrays(size_t n) {
    exists.resize(n);
    values.resize(n);
    degrees.resize(n);
    capacities.resize(n);
    edges.resize(n);
  }


  /**
   * Sort the outgoing edges of a vertex in the graph.
   * @param u source vertex id
   */
  inline void sortEdges(K u) {
    if (u >= span()) return;
    auto fl = [](const auto& a, const auto& b) { return a.first < b.first; };
    auto ib = edges[u], ie = edges[u] + degrees[u];
    sort(ib, ie, fl);
  }


  /**
   * Ensure that the outgoing edges of a vertex in the graph are unique.
   * @param u source vertex id
   */
  inline void uniqueEdges(K u) {
    if (u >= span()) return;
    auto fe = [](const auto& a, const auto& b) { return a.first == b.first; };
    auto ib = edges[u], ie = edges[u] + degrees[u];
    auto it = unique(ib, ie, fe);
    degrees[u] = it - ib;
  }


  public:
  /**
   * Remove all vertices and edges from the graph.
   */
  inline void clear() noexcept {
    N = 0; M = 0;
    exists.clear();
    values.clear();
    degrees.clear();
    capacities.clear();
    edges.clear();
    resetAllocators();
  }


  /**
   * Clear the outgoing edges of a vertex in the graph.
   * @param u source vertex id
   */
  inline void clearEdges(K u) {
    if (u >= span() || !edges[u]) return;
    deallocate(edges[u], capacities[u]);
    edges[u] = nullptr;
    degrees[u] = 0;
    capacities[u] = 0;
  }


  /**
   * Allocate space for outgoing edges of a vertex in the graph.
   * @param u source vertex id
   * @param deg expected degree of the vertex
   * @note This works only if the vertex has no edges yet.
   */
  inline void allocateEdges(K u, K deg) {
    if (u >= span() || edges[u]) return;
    K cap = allocationCapacity(deg);
    edges[u] = (pair<K, E>*) allocate(cap);
    capacities[u] = cap;
  }


  /**
   * Reserve space for outgoing edges of a vertex in the graph.
   * @param u source vertex id
   * @param deg expected degree of the vertex
   * @note Also supports shrinking the capacity.
   */
  inline void reserveEdges(K u, K deg) {
    if (u >= span()) return;
    // Deallocate if no edges are expected.
    if (deg==0 && degrees[u]==0) {
      deallocate(edges[u], capacities[u]);
      edges[u] = nullptr;
      capacities[u] = 0;
      return;
    }
    // Skip if no change in capacity.
    K cap = allocationCapacity(deg);
    if (cap == capacities[u]) return;
    // Allocate new memory and copy old data.
    void *ptr = allocate(cap);
    memcpy(ptr, edges[u], degrees[u] * EDGE);
    deallocate(edges[u], capacities[u]);
    // Update pointer and capacities.
    edges[u] = (pair<K, E>*) ptr;
    capacities[u] = cap;
  }


  /**
   * Reserve space for a number of vertices and edges in the graph.
   * @param n number of vertices to reserve space for
   * @param deg expected average degree of vertices
   */
  inline void reserve(size_t n, size_t deg=0) {
    size_t S = max(n, span());
    resizeArrays(S);
    if (deg==0) return;
    for (K u=0; u<S; ++u)
      reserveEdges(u, deg);
  }


  /**
   * Reserve space for a number of vertices and edges in the graph [parallel].
   * @param n number of vertices to reserve space for
   * @param deg expected average degree of vertices
   */
  inline void reserveOmp(size_t n, size_t deg=0) {
    size_t S = max(n, span());
    resizeArrays(S);
    if (deg==0) return;
    #pragma omp parallel for schedule(dynamic, 2048)
    for (K u=0; u<S; ++u)
      reserveEdges(u, deg);
  }


  /**
   * Adjust the span of the graph.
   * @param n new span
   */
  inline void respan(size_t n) {
    size_t  S = span();
    size_t dN = 0, dM = 0;
    for (K u=n; u<S; ++u) {
      if (!exists[u]) continue;
      ++dN; dM += degrees[u];
      clearEdges(u);
    }
    N -= dN; M -= dM;
    resizeArrays(n);
  }


  /**
   * Adjust the span of the graph [parallel].
   * @param n new span
   */
  inline void respanOmp(size_t n) {
    size_t S = span();
    size_t dN = 0, dM = 0;
    #pragma omp parallel for schedule(dynamic, 2048) reduction(+:dN,dM)
    for (K u=n; u<S; ++u) {
      if (!exists[u]) continue;
      ++dN; dM += degrees[u];
      clearEdges(u);
    }
    N -= dN; M -= dM;
    resizeArrays(n);
  }


  /**
   * Update the count of vertices and edges in the graph.
   */
  inline void updateCounts() {
    N = 0; M = 0;
    for (K u=0; u < span(); ++u) {
      if (!exists[u]) continue;
      ++N; M += degrees[u];
    }
  }


  /**
   * Update the count of vertices and edges in the graph [parallel].
   */
  inline void updateCountsOmp() {
    N = 0; M = 0;
    #pragma omp parallel for schedule(auto) reduction(+:N,M)
    for (K u=0; u < span(); ++u) {
      if (!exists[u]) continue;
      ++N; M += degrees[u];
    }
  }


  /**
   * Update the graph after changes.
   * @param isUnique are the edges unique?
   * @param isSorted are the edges sorted?
   */
  inline void update(bool isUnique=false, bool isSorted=false) {
    if (!isSorted) {
      for (K u=0; u < span(); ++u)
        sortEdges(u);
    }
    if (!isUnique) {
      for (K u=0; u < span(); ++u)
        uniqueEdges(u);
    }
    updateCounts();
  }


  /**
   * Update the graph after changes [parallel].
   * @param isUnique are the edges unique?
   * @param isSorted are the edges sorted?
   */
  inline void updateOmp(bool isUnique=false, bool isSorted=false) {
    if (!isSorted) {
      #pragma omp parallel for schedule(dynamic, 2048)
      for (K u=0; u < span(); ++u)
        sortEdges(u);
    }
    if (!isUnique) {
      #pragma omp parallel for schedule(dynamic, 2048)
      for (K u=0; u < span(); ++u)
        uniqueEdges(u);
    }
    updateCountsOmp();
  }


  /**
   * Add a vertex to the graph.
   * @param u vertex id
   * @note `update()` must be called after all vertices are added.
   */
  inline void addVertex(K u) {
    if (u >= span()) respan(u+1);
    exists[u] = true;
  }


  /**
   * Add a vertex to the graph.
   * @param u vertex id
   * @param d associated data of the vertex
   * @note `update()` must be called after all vertices are added.
   */
  inline void addVertex(K u, V d) {
    if (u >= span()) respan(u+1);
    exists[u] = true;
    values[u] = d;
  }


  /**
   * Remove a vertex from the graph.
   * @param u vertex id
   * @note `update()` must be called after all vertices are removed.
   */
  inline void removeVertex(K u) {
    if (!hasVertex(u)) return;
    exists[u] = false;
    values[u] = V();
    clearEdges(u);
  }


  /**
   * Add an outgoing edge to the graph, without "any" checks.
   * @param u source vertex id
   * @param v target vertex id
   * @param w associated weight of the edge
   */
  inline void addEdgeUnchecked(K u, K v, E w=E()) {
    auto *ptr = edges[u];
    ptr[degrees[u]++] = {v, w};
  }


  /**
   * Add an outgoing edge to the graph, without "any" checks [parallel].
   * @param u source vertex id
   * @param v target vertex id
   * @param w associated weight of the edge
   */
  inline void addEdgeUncheckedOmp(K u, K v, E w=E()) {
    auto *ptr = edges[u];
    K i = K();
    #pragma omp atomic capture
    i = degrees[u]++;
    ptr[i] = {v, w};
  }


  /**
   * Remove outgoing edges from a vertex in the graph.
   * @param u source vertex id
   * @param ib begin iterator of edge keys to remove
   * @param ie end iterator of edge keys to remove
   * @note [ib, ie) must be sorted and unique.
   */
  template <class I>
  inline void removeEdges(K u, I ib, I ie) {
    if (!hasVertex(u)) return;
    auto *eb = edges[u], *ee = edges[u] + degrees[u];
    auto  fl = [](const auto& a, const auto& b) { return a.first <  b; };
    auto  fe = [](const auto& a, const auto& b) { return a.first == b; };
    auto  it = set_difference_inplace(eb, ee, ib, ie, fl, fe);
    degrees[u] = it - eb;
  }


  /**
   * Add outgoing edges to a vertex in the graph.
   * @param u source vertex id
   * @param ib begin iterator of edges to add
   * @param ie end iterator of edges to add
   * @param buf scratch buffer for the update (of size at least 3 + distance(ib, ie))
   * @note [ib, ie) must be sorted and unique.
   */
  template <class I>
  inline void addEdges(K u, I ib, I ie) {
    if (!hasVertex(u) || ib==ie) return;
    auto *eb = edges[u], *ee = edges[u] + degrees[u];
    size_t deg = degrees[u] + distance(ib, ie);
    size_t cap = allocationCapacity(deg);
    void  *ptr = allocate(cap);
    auto fl = [](const auto& a, const auto& b) { return a.first <  b.first; };
    auto fe = [](const auto& a, const auto& b) { return a.first == b.first; };
    auto it = set_union(eb, ee, ib, ie, (pair<K, E>*) ptr, fl);
    degrees[u] = it - eb;
  }
  #pragma endregion


  #pragma region CONSTRUCTORS
  public:
  /**
   * Create an empty graph.
   */
  ArenaDiGraph() {
    mx = new Allocator();
  }


  /**
   * Destroy the Arena DiGraph.
   */
  ~ArenaDiGraph() {
    delete mx;
  }
  #pragma endregion
  #pragma endregion
};




/**
 * Directed graph that memorizes only out-edges for each vertex.
 * @tparam K key type (vertex id)
 * @tparam V vertex value type (vertex data)
 * @tparam E edge value type (edge weight)
 */
template <class K=uint32_t, class V=None, class E=None>
class DiGraph {
  #pragma region TYPES
  public:
  /** Key type (vertex id). */
  using key_type = K;
  /** Vertex value type (vertex data). */
  using vertex_value_type = V;
  /** Edge value type (edge weight). */
  using edge_value_type   = E;
  #pragma endregion


  #pragma region DATA
  protected:
  /** Number of vertices. */
  size_t N = 0;
  /** Number of edges. */
  size_t M = 0;
  /** Vertex existence flags. */
  vector<bool> exists;
  /** Vertex values. */
  vector<V> values;
  /** Outgoing edges for each vertex (including edge weights). */
  vector<LazyBitset<K, E>> edges;
  #pragma endregion


  #pragma region METHODS
  #pragma region PROPERTIES
  public:
  /**
   * Get the size of buffer required to store data associated with each vertex
   * in the graph, indexed by its vertex-id.
   * @returns size of buffer required
   */
  inline size_t span() const noexcept {
    return exists.size();
  }

  /**
   * Get the number of vertices in the graph.
   * @returns |V|
   */
  inline size_t order() const noexcept {
    return N;
  }

  /**
   * Get the number of edges in the graph.
   * @returns |E|
   */
  inline size_t size() const noexcept {
    return M;
  }

  /**
   * Check if the graph is empty.
   * @returns is the graph empty?
   */
  inline bool empty() const noexcept {
    return N == 0;
  }

  /**
   * Check if the graph is directed.
   * @returns is the graph directed?
   */
  inline bool directed() const noexcept {
    return true;
  }
  #pragma endregion


  #pragma region FOREACH
  public:
  /**
   * Iterate over the vertices in the graph.
   * @param fp process function (vertex id, vertex data)
   */
  template <class FP>
  inline void forEachVertex(FP fp) const noexcept {
    for (K u=0; u<span(); ++u)
      if (exists[u]) fp(u, values[u]);
  }

  /**
   * Iterate over the vertex ids in the graph.
   * @param fp process function (vertex id)
   */
  template <class FP>
  inline void forEachVertexKey(FP fp) const noexcept {
    for (K u=0; u<span(); ++u)
      if (exists[u]) fp(u);
  }

  /**
   * Iterate over the outgoing edges of a source vertex in the graph.
   * @param u source vertex id
   * @param fp process function (target vertex id, edge weight)
   */
  template <class FP>
  inline void forEachEdge(K u, FP fp) const noexcept {
    edges[u].forEach(fp);
  }

  /**
   * Iterate over the target vertex ids of a source vertex in the graph.
   * @param u source vertex id
   * @param fp process function (target vertex id)
   */
  template <class FP>
  inline void forEachEdgeKey(K u, FP fp) const noexcept {
    edges[u].forEachKey(fp);
  }
  #pragma endregion


  #pragma region ACCESS
  public:
  /**
   * Check if a vertex exists in the graph.
   * @param u vertex id
   * @returns does the vertex exist?
   */
  inline bool hasVertex(K u) const noexcept {
    return u < span() && exists[u];
  }

  /**
   * Check if an edge exists in the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @returns does the edge exist?
   */
  inline bool hasEdge(K u, K v) const noexcept {
    return u < span() && edges[u].has(v);
  }

  /**
   * Get the number of outgoing edges of a vertex in the graph.
   * @param u vertex id
   * @returns number of outgoing edges of the vertex
   */
  inline size_t degree(K u) const noexcept {
    return u < span()? edges[u].size() : 0;
  }

  /**
   * Get the vertex data of a vertex in the graph.
   * @param u vertex id
   * @returns associated data of the vertex
   */
  inline V vertexValue(K u) const noexcept {
    return u < span()? values[u] : V();
  }

  /**
   * Set the vertex data of a vertex in the graph.
   * @param u vertex id
   * @param d associated data of the vertex
   * @returns success?
   */
  inline bool setVertexValue(K u, V d) noexcept {
    if (!hasVertex(u)) return false;
    values[u] = d;
    return true;
  }

  /**
   * Get the edge weight of an edge in the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @returns associated weight of the edge
   */
  inline E edgeValue(K u, K v) const noexcept {
    return u < span()? edges[u].get(v) : E();
  }

  /**
   * Set the edge weight of an edge in the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @param w associated weight of the edge
   * @returns success?
   */
  inline bool setEdgeValue(K u, K v, E w) noexcept {
    if (!hasVertex(u) || !hasVertex(v)) return false;
    return edges[u].set(v, w);
  }
  #pragma endregion


  #pragma region UPDATE
  public:
  /**
   * Remove all vertices and edges from the graph.
   */
  inline void clear() noexcept {
    N = 0; M = 0;
    exists.clear();
    values.clear();
    edges.clear();
  }

  /**
   * Reserve space for outgoing edges of a vertex in the graph.
   * @param u source vertex id
   * @param deg expected degree of the vertex
   */
  inline void reserveEdges(K u, size_t deg) {
    if (u < span()) edges[u].reserve(deg);
  }


  /**
   * Reserve space for a number of vertices and edges in the graph.
   * @param n number of vertices to reserve space for
   * @param deg expected average degree of vertices
   */
  inline void reserve(size_t n, size_t deg=0) {
    size_t S = max(n, span());
    exists.resize(S);
    values.resize(S);
    edges.resize(S);
    if (deg==0) return;
    for (K u=0; u<S; ++u)
      edges[u].reserve(deg);
  }

  /**
   * Adjust the span of the graph.
   * @param n new span
   * @note This operation is lazy.
   */
  inline void respan(size_t n) {
    exists.resize(n);
    values.resize(n);
    edges.resize(n);
  }

  /**
   * Update the outgoing edges of a vertex in the graph to reflect the changes.
   * @param u source vertex id
   * @param buf scratch buffer for the update
   */
  inline void updateEdges(K u, vector<pair<K, E>> *buf=nullptr) {
    if (u < span()) edges[u].update(buf);
  }

  /**
   * Update the graph to reflect the changes.
   * @note This is an expensive operation.
   */
  inline void update() {
    vector<pair<K, E>> buf;
    N = 0; M = 0;
    forEachVertexKey([&](K u) {
      edges[u].update(&buf);
      M += degree(u); ++N;
    });
  }

  /**
   * Add a vertex to the graph.
   * @param u vertex id
   * @note This operation is lazy.
   */
  inline void addVertex(K u) {
    if (hasVertex(u)) return;
    if (u >= span()) respan(u+1);
    exists[u] = true;
  }

  /**
   * Add a vertex to the graph.
   * @param u vertex id
   * @param d associated data of the vertex
   * @note This operation is lazy.
   */
  inline void addVertex(K u, V d) {
    if (hasVertex(u)) { values[u] = d; return; }
    if (u >= span()) respan(u+1);
    exists[u] = true;
    values[u] = d;
  }

  /**
   * Add an outgoing edge to the graph if a condition is met.
   * @param u source vertex id
   * @param v target vertex id
   * @param w associated weight of the edge
   * @param ft test function (source vertex id)
   */
  template <class FT>
  inline void addEdgeIf(K u, K v, E w, FT ft) {
    addVertex(u);
    addVertex(v);
    if (ft(u)) edges[u].add(v, w);
  }

  /**
   * Add an outgoing edge to the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @param w associated weight of the edge
   * @note This operation is lazy.
   */
  inline void addEdge(K u, K v, E w=E()) {
    auto ft = [](K u) { return true; };
    addEdgeIf(u, v, w, ft);
  }

  /**
   * Remove an outgoing edge from the graph if a condition is met.
   * @param u source vertex id
   * @param v target vertex id
   * @param ft test function (source vertex id)
   */
  template <class FT>
  inline void removeEdgeIf(K u, K v, FT ft) {
    if (!hasVertex(u) || !hasVertex(v)) return;
    if (ft(u)) edges[u].remove(v);
  }

  /**
   * Remove an outgoing edge from the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @note This operation is lazy.
   */
  inline void removeEdge(K u, K v) {
    auto ft = [](K u) { return true; };
    removeEdgeIf(u, v, ft);
  }

  /**
   * Remove a vertex from the graph.
   * @param u vertex id
   * @note This operation is lazy.
   */
  inline void removeVertex(K u) {
    if (!hasVertex(u)) return;
    exists[u] = false;
    values[u] = V();
    edges[u].clear();
  }
  #pragma endregion
  #pragma endregion
};



/**
 * A directed graph with CSR representation.
 * @tparam K key type (vertex id)
 * @tparam V vertex value type (vertex data)
 * @tparam E edge value type (edge weight)
 * @tparam O offset type
 */
template <class K=uint32_t, class V=None, class E=None, class O=size_t>
class DiGraphCsr {
  #pragma region TYPES
  public:
  /** Key type (vertex id). */
  using key_type = K;
  /** Vertex value type (vertex data). */
  using vertex_value_type = V;
  /** Edge value type (edge weight). */
  using edge_value_type   = E;
  /** Offset type (edge offset). */
  using offset_type       = O;
  #pragma endregion


  #pragma region DATA
  public:
  /** Offsets of the outgoing edges of vertices. */
  vector<O> offsets;
  /** Degree of each vertex. */
  vector<K> degrees;
  /** Vertex values. */
  vector<V> values;
  /** Vertex ids of the outgoing edges of each vertex (lookup using offsets). */
  vector<K> edgeKeys;
  /** Edge weights of the outgoing edges of each vertex (lookup using offsets). */
  vector<E> edgeValues;
  #pragma endregion


  #pragma region METHODS
  #pragma region PROPERTIES
  public:
  /**
   * Get the size of buffer required to store data associated with each vertex
   * in the graph, indexed by its vertex-id.
   * @returns size of buffer required
   */
  inline size_t span() const noexcept {
    return degrees.size();
  }

  /**
   * Get the number of vertices in the graph.
   * @returns |V|
   */
  inline size_t order() const noexcept {
    return degrees.size();
  }

  /**
   * Obtain the number of edges in the graph.
   * @returns |E|
   */
  inline size_t size() const noexcept {
    size_t M = 0;
    for (auto d : degrees)
      M += d;
    return M;
  }

  /**
   * Check if the graph is empty.
   * @returns is the graph empty?
   */
  inline bool empty() const noexcept {
    return degrees.empty();
  }

  /**
   * Check if the graph is directed.
   * @returns is the graph directed?
   */
  inline bool directed() const noexcept {
    return true;
  }
  #pragma endregion


  #pragma region FOREACH
  public:
  /**
   * Iterate over the vertices in the graph.
   * @param fp process function (vertex id, vertex data)
   */
  template <class FP>
  inline void forEachVertex(FP fp) const noexcept {
    for (K u=0; u<span(); ++u)
      fp(u, values[u]);
  }

  /**
   * Iterate over the vertex ids in the graph.
   * @param fp process function (vertex id)
   */
  template <class FP>
  inline void forEachVertexKey(FP fp) const noexcept {
    for (K u=0; u<span(); ++u)
      fp(u);
  }

  /**
   * Iterate over the outgoing edges of a source vertex in the graph.
   * @param u source vertex id
   * @param fp process function (target vertex id, edge weight)
   */
  template <class FP>
  inline void forEachEdge(K u, FP fp) const noexcept {
    size_t i = offsets[u];
    size_t d = degrees[u];
    for (size_t I=i+d; i<I; ++i)
      fp(edgeKeys[i], edgeValues[i]);
  }

  /**
   * Iterate over the target vertex ids of a source vertex in the graph.
   * @param u source vertex id
   * @param fp process function (target vertex id)
   */
  template <class FP>
  inline void forEachEdgeKey(K u, FP fp) const noexcept {
    size_t i = offsets[u];
    size_t d = degrees[u];
    for (size_t I=i+d; i<I; ++i)
      fp(edgeKeys[i]);
  }
  #pragma endregion


  #pragma region OFFSET
  public:
  /**
   * Get the offset of an edge in the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @returns offset of the edge, or -1 if it does not exist
   */
  inline size_t edgeOffset(K u, K v) const noexcept {
    if (!hasVertex(u) || !hasVertex(v)) return size_t(-1);
    size_t  i = offsets[u];
    size_t  d = degrees[u];
    auto   ib = edgeKeys.begin() + i;
    auto   ie = edgeKeys.begin() + i + d;
    auto   it = find(ib, ie, v);
    return it!=ie? it - edgeKeys.begin() : size_t(-1);
  }
  #pragma endregion


  #pragma region ACCESS
  public:
  /**
   * Check if a vertex exists in the graph.
   * @param u vertex id
   * @returns does the vertex exist?
   */
  inline bool hasVertex(K u) const noexcept {
    return u < span();
  }

  /**
   * Check if an edge exists in the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @returns does the edge exist?
   */
  inline bool hasEdge(K u, K v) const noexcept {
    size_t o = edgeOffset(u, v);
    return o != size_t(-1);
  }

  /**
   * Get the number of outgoing edges of a vertex in the graph.
   * @param u vertex id
   * @returns number of outgoing edges of the vertex
   */
  inline size_t degree(K u) const noexcept {
    return u < span()? degrees[u] : 0;
  }

  /**
   * Get the vertex data of a vertex in the graph.
   * @param u vertex id
   * @returns associated data of the vertex
   */
  inline V vertexValue(K u) const noexcept {
    return u < span()? values[u] : V();
  }

  /**
   * Set the vertex data of a vertex in the graph.
   * @param u vertex id
   * @param d associated data of the vertex
   * @returns success?
   */
  inline bool setVertexValue(K u, V d) noexcept {
    if (!hasVertex(u)) return false;
    values[u] = d;
    return true;
  }

  /**
   * Get the edge weight of an edge in the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @returns associated weight of the edge
   */
  inline E edgeValue(K u, K v) const noexcept {
    size_t o = edgeOffset(u, v);
    return o != size_t(-1)? edgeValues[o] : E();
  }

  /**
   * Set the edge weight of an edge in the graph.
   * @param u source vertex id
   * @param v target vertex id
   * @param w associated weight of the edge
   * @returns success?
   */
  inline bool setEdgeValue(K u, K v, E w) noexcept {
    size_t o = edgeOffset(u, v);
    if (o == size_t(-1)) return false;
    edgeValues[o] = w;
    return true;
  }
  #pragma endregion


  #pragma region UPDATE
  public:
  /**
   * Adjust the order of the graph (or the number of vertices).
   * @param n new order, or number of vertices
   */
  inline void resize(size_t n) {
    offsets.resize(n+1);
    degrees.resize(n);
    values.resize(n);
  }


  /**
   * Adjust the order and size of the graph (or the number of vertices and edges).
   * @param n new order, or number of vertices
   * @param m new size, or number of edges
   */
  inline void resize(size_t n, size_t m) {
    offsets.resize(n+1);
    degrees.resize(n);
    values.resize(n);
    edgeKeys.resize(m);
    edgeValues.resize(m);
  }


  /**
   * Adjust the span of the graph (or the number of vertices).
   * @param n new span
   */
  inline void respan(size_t n) {
    resize(n);
  }
  #pragma endregion
  #pragma endregion


  #pragma region CONSTRUCTORS
  public:
  /**
   * Create an empty CSR representation of a directed graph.
   */
  DiGraphCsr() {
    offsets.push_back(0);
  }


  /**
   * Allocate space for CSR representation of a directed graph.
   * @param n number of vertices
   * @param m number of edges
   */
  DiGraphCsr(size_t n, size_t m) {
    offsets.resize(n+1);
    degrees.resize(n);
    values.resize(n);
    edgeKeys.resize(m);
    edgeValues.resize(m);
  }
  #pragma endregion
};
#pragma endregion




#pragma region METHODS
#pragma region SET OPERATIONS
/**
 * Subtract a graph's edges from another graph.
 * @param a graph to subtract from (updated)
 * @param x graph to subtract
 */
template <class H, class G>
inline void subtractGraphEdgesU(H& a, const G& x) {
  using K = typename G::key_type;
  using E = typename G::edge_value_type;
  a.forEachVertexKey([&](auto u) {
    if (!x.hasVertex(u)) return;
    auto ib = static_transform_iterator(x.beginEdges(u), ConstPairFirst<K, E>());
    auto ie = static_transform_iterator(x.endEdges(u),   ConstPairFirst<K, E>());
    a.removeEdges(u, ib, ie);
  });
  a.update(true, true);
}


/**
 * Subtract a graph's edges from another graph [parallel].
 * @param a graph to subtract from (updated)
 * @param x graph to subtract
 */
template <class H, class G>
inline void subtractGraphOmpU(H& a, const G& x) {
  using  K = typename G::key_type;
  using  E = typename G::edge_value_type;
  size_t S = a.span();
  #pragma omp parallel for schedule(dynamic, 2048)
  for (K u=0; u<S; ++u) {
    if (!a.hasVertex(u) || !x.hasVertex(u)) continue;
    auto ib = static_transform_iterator(x.beginEdges(u), ConstPairFirst<K, E>());
    auto ie = static_transform_iterator(x.endEdges(u),   ConstPairFirst<K, E>());
    a.removeEdges(u, ib, ie);
  }
  a.updateOmp(true, true);
}


/**
 * Add a graph's edges to another graph.
 * @param a graph to add to (updated)
 * @param x graph to add
 */
template <class H, class G>
inline void addGraphU(H& a, const G& x) {
  using  K = typename G::key_type;
  using  E = typename G::edge_value_type;
  size_t A = a.span();
  size_t X = x.span();
  size_t S = max(A, X);
  a.respan(S);
  x.forEachVertex([&](auto u, auto d) {
    a.addVertex(u, d);
    auto ib = x.beginEdges(u);
    auto ie = x.endEdges(u);
    a.addEdges(u, ib, ie);
  });
  a.update(true, true);
}


#ifdef OPENMP
/**
 * Add a graph's edges to another graph [parallel].
 * @param a graph to add to (updated)
 * @param x graph to add
 */
template <class H, class G>
inline void addGraphOmpU(H& a, const G& x) {
  using  K = typename G::key_type;
  using  E = typename G::edge_value_type;
  size_t A = a.span();
  size_t X = x.span();
  size_t S = max(A, X);
  a.respan(S);
  #pragma omp parallel for schedule(dynamic, 2048)
  for (K u=0; u<X; ++u) {
    int t = omp_get_thread_num();
    if (!x.hasVertex(u)) continue;
    a.addVertex(u, x.vertexValue(u));
    auto ib = x.beginEdges(u);
    auto ie = x.endEdges(u);
    a.addEdges(u, ib, ie);
  }
  a.updateOmp(true, true);
}
#endif
#pragma endregion




#pragma region WRITE
/**
 * Write the only the sizes of a graph to an output stream.
 * @param a output stream
 * @param x graph
 */
template <class G>
inline void writeGraphSizes(ostream& a, const G& x) {
  a << "order: " << x.order() << " size: " << x.size();
  a << (x.directed()? " [directed]" : " [undirected]") << " {}";
}

/**
 * @brief Write the full details of a graph to an output stream.
 * @param a output stream
 * @param x graph
 */
template <class G>
inline void writeGraphDetailed(ostream& a, const G& x) {
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

/**
 * Write a graph to an output stream.
 * @param a output stream
 * @param x graph
 * @param detailed write detailed information?
 */
template <class G>
inline void writeGraph(ostream& a, const G& x, bool detailed=false) {
  if (detailed) writeGraphDetailed(a, x);
  else writeGraphSizes(a, x);
}


/**
 * Write a graph to an output stream.
 * @param a output stream
 * @param x graph
 * @param detailed write detailed information?
 */
template <class K, class V, class E>
inline void write(ostream& a, const ArenaDiGraph<K, V, E>& x, bool detailed=false) {
  writeGraph(a, x, detailed);
}


/**
 * Write a graph to an output stream.
 * @param a output stream
 * @param x graph
 * @param detailed write detailed information?
 */
template <class K, class V, class E>
inline void write(ostream& a, const DiGraph<K, V, E>& x, bool detailed=false) {
  writeGraph(a, x, detailed);
}


/**
 * Write a graph to an output stream.
 * @param a output stream
 * @param x csr graph
 * @param detailed write detailed information?
 */
template <class K, class V, class E, class O>
inline void write(ostream& a, const DiGraphCsr<K, V, E, O>& x, bool detailed=false) {
  writeGraph(a, x, detailed);
}


/**
 * Write only the sizes of a graph to an output stream.
 * @param a output stream
 * @param x graph
 */
template <class K, class V, class E>
inline ostream& operator<<(ostream& a, const ArenaDiGraph<K, V, E>& x) {
  write(a, x);
  return a;
}


/**
 * Write only the sizes of a graph to an output stream.
 * @param a output stream
 * @param x graph
 */
template <class K, class V, class E>
inline ostream& operator<<(ostream& a, const DiGraph<K, V, E>& x) {
  write(a, x);
  return a;
}


/**
 * Write only the sizes of a graph to an output stream.
 * @param a output stream
 * @param x csr graph
 */
template <class K, class V, class E, class O>
inline ostream& operator<<(ostream& a, const DiGraphCsr<K, V, E, O>& x) {
  write(a, x);
  return a;
}
#pragma endregion
#pragma endregion
