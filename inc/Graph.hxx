#pragma once
#include <cstring>
#include <utility>
#include <iterator>
#include <memory>
#include <vector>
#include <ostream>
#include <algorithm>
#include "_main.hxx"
#ifdef OPENMP
#include <omp.h>
#endif

using std::pair;
using std::shared_ptr;
using std::vector;
using std::ostream;
using std::make_shared;
using std::memcpy;
using std::distance;
using std::max;
using std::find_if;
using std::lower_bound;
using std::copy;
using std::sort;
using std::unique;
using std::set_union;
using std::set_difference;




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
  /** Type for storing boolean flags. */
  using bool_type = uint64_t;
  public:
  /** Key type (vertex id). */
  using key_type = K;
  /** Vertex value type (vertex data). */
  using vertex_value_type = V;
  /** Edge value type (edge weight). */
  using edge_value_type = E;
  /** Offset type (edge offset). */
  using offset_type = size_t;
  /** Memory allocator type. */
  using allocator_type = ConcurrentPow2Allocator<>;
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
  /** Vertex existence flags. */
  bool_type *exists = nullptr;
  /** Outgoing edges for each vertex (including edge weights). */
  pair<K, E> **edges = nullptr;
  /** Out-degree of each vertex. */
  K *degrees = nullptr;
  /** Edge capacity of each vertex. */
  K *capacities = nullptr;
  /** Memory allocator. */
  shared_ptr<allocator_type> mx;
  /** Span of the graph. */
  size_t SPAN = 0;
  /** Vertex values. */
  V *values = nullptr;
  /** Number of vertices. */
  size_t N = 0;
  /** Number of edges. */
  size_t M = 0;
  /** Memory reserved for vertices. */
  size_t RESV = 0;
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
    return SPAN;
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

  /**
   * Set the number of vertices in the graph.
   * @param n new number of vertices
   * @note This is used for internal purposes only.
   */
  inline void setOrder(size_t n) noexcept {
    N = n;
  }

  /**
   * Set the number of edges in the graph.
   * @param m new number of edges
   * @note This is used for internal purposes only.
   */
  inline void setSize(size_t m) noexcept {
    M = m;
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
      if (getBit(exists, u)) fp(u, values[u]);
  }

  /**
   * Iterate over the vertex ids in the graph.
   * @param fp process function (vertex id)
   */
  template <class FP>
  inline void forEachVertexKey(FP fp) const noexcept {
    for (K u=0; u<span(); ++u)
      if (getBit(exists, u)) fp(u);
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
  inline const pair<K, E>* beginEdges(K u) const noexcept {
    return u < span()? edges[u] : nullptr;
  }


  /**
   * Get an iterator to the end of edges of a vertex.
   * @param u vertex id
   * @returns end iterator of edges
   */
  inline const pair<K, E>* endEdges(K u) const noexcept {
    return u < span()? edges[u] + degrees[u] : nullptr;
  }


  /**
   * Get an iterator to the begin of edges of a vertex.
   * @param u vertex id
   * @returns begin iterator of edges
   */
  inline pair<K, E>* beginEdges(K u) noexcept {
    return u < span()? edges[u] : nullptr;
  }


  /**
   * Get an iterator to the end of edges of a vertex.
   * @param u vertex id
   * @returns end iterator of edges
   */
  inline pair<K, E>* endEdges(K u) noexcept {
    return u < span()? edges[u] + degrees[u] : nullptr;
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
    return u < span() && getBit(exists, u);
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
    if (u >= span()) return false;
    auto ib = edges[u], ie = edges[u] + degrees[u];
    auto it = findEntry(ib, ie, v);
    if (it == ie) return false;
    (*it).second = w;
    return true;
  }


  /**
   * Get the memory allocator in use.
   * @returns memory allocator
   */
  inline shared_ptr<allocator_type> allocator() const noexcept {
    return mx;
  }
  #pragma endregion


  #pragma region UPDATE
  protected:
  /**
   * Get the allocation capacity for a number of elements.
   * @param n number of elements
   * @returns allocation capacity
   */
  static inline constexpr K allocationCapacity(size_t n) noexcept {
    return allocator_type::allocationCapacity(n*EDGE) / EDGE;
  }


  /**
   * Setup the memory allocator, if not already set.
   */
  inline void setupAllocator() {
    if (!mx) mx = make_shared<allocator_type>();
  }


  /**
   * Allocate memory for a number of edges.
   * @param c allocation capacity
   * @returns pointer to the allocated memory
   */
  inline pair<K, E>* allocate(K c) {
    return (pair<K, E>*) mx->allocate(c*EDGE);
  }


  /**
   * Deallocate memory for a number of edges.
   * @param ptr pointer to the memory
   * @param c allocation capacity
   */
  inline void deallocate(pair<K, E> *ptr, K c) {
    mx->deallocate(ptr, c*EDGE);
  }


  /**
   * Resize arrays to specified size.
   * @param n new size
   * @note The edges are not cleared!
   */
  inline void resizeArrays(size_t n) {
    constexpr size_t B = 8 * sizeof(bool_type);
    // Compute new reserved size (round up to page size).
    size_t resv = ceilDiv(n, size_t(PAGE_SIZE)) * PAGE_SIZE;
    if (n <= SPAN && resv == RESV) return;
    // Allocate new memory.
    edges      = reallocateValues(edges, SPAN, RESV, n, resv);
    degrees    = reallocateValues(degrees, SPAN, RESV, n, resv);
    capacities = reallocateValues(capacities, SPAN, RESV, n, resv);
    values     = reallocateValues(values, SPAN, RESV, n, resv);
    exists     = reallocateValues(exists, ceilDiv(SPAN, B), ceilDiv(RESV, B), ceilDiv(n, B), ceilDiv(resv, B));
    // Update span and reserved size.
    SPAN = n;
    RESV = resv;
  }


  /**
   * Resize arrays to specified size [parallel].
   * @param n new size
   * @note The edges are not cleared!
   */
  inline void resizeArraysOmp(size_t n) {
    constexpr size_t B = 8 * sizeof(bool_type);
    // Compute new reserved size (round up to page size).
    size_t resv = ceilDiv(n, size_t(PAGE_SIZE)) * PAGE_SIZE;
    if (n <= SPAN && resv == RESV) return;
    // Allocate new memory.
    edges      = reallocateValuesOmp(edges, SPAN, RESV, n, resv);
    degrees    = reallocateValuesOmp(degrees, SPAN, RESV, n, resv);
    capacities = reallocateValuesOmp(capacities, SPAN, RESV, n, resv);
    values     = reallocateValuesOmp(values, SPAN, RESV, n, resv);
    exists     = reallocateValuesOmp(exists, ceilDiv(SPAN, B), ceilDiv(RESV, B), ceilDiv(n, B), ceilDiv(resv, B));
    // Update span and reserved size.
    SPAN = n;
    RESV = resv;
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
   * Remove all vertices and edges from the graph.
   */
  inline void clear() {
    if (mx == nullptr) return;
    // Clear all edges first, if the allocator is shared.
    bool isShared = mx.use_count() > 1;
    if (!isShared)  mx->reset();
    else {
      for (K u=0; u<SPAN; ++u)
        clearEdges(u);
    }
    // Update counts.
    SPAN = 0; N = 0; M = 0;
  }


  /**
   * Remove all vertices and edges from the graph [parallel].
   */
  inline void clearOmp() {
    if (mx == nullptr) return;
    // Clear all edges first, if the allocator is shared.
    bool isShared = mx.use_count() > 1;
    if (!isShared)  mx->reset();
    else {
      #pragma omp parallel for schedule(dynamic, 2048)
      for (K u=0; u<SPAN; ++u)
        clearEdges(u);
    }
    // Update counts.
    SPAN = 0; N = 0; M = 0;
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
    edges[u] = allocate(cap);
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
    pair<K, E> *ptr = allocate(cap);
    memcpy(ptr, edges[u], degrees[u] * EDGE);
    deallocate(edges[u], capacities[u]);
    // Update pointer and capacities.
    edges[u] = ptr;
    capacities[u] = cap;
  }


  /**
   * Reserve space for a number of vertices and edges in the graph.
   * @param n number of vertices to reserve space for
   * @param deg expected average degree of vertices
   */
  inline void reserve(size_t n, size_t deg=0) {
    setupAllocator();
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
    setupAllocator();
    size_t S = max(n, span());
    resizeArraysOmp(S);
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
    setupAllocator();
    size_t  S = span();
    size_t dN = 0, dM = 0;
    for (K u=n; u<S; ++u) {
      if (!getBit(exists, u)) continue;
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
    setupAllocator();
    size_t S = span();
    size_t dN = 0, dM = 0;
    #pragma omp parallel for schedule(dynamic, 2048) reduction(+:dN,dM)
    for (K u=n; u<S; ++u) {
      if (!getBit(exists, u)) continue;
      ++dN; dM += degrees[u];
      clearEdges(u);
    }
    N -= dN; M -= dM;
    resizeArraysOmp(n);
  }


  /**
   * Update the count of vertices and edges in the graph.
   */
  inline void updateCounts() {
    N = 0; M = 0;
    for (K u=0; u < span(); ++u) {
      if (!getBit(exists, u)) continue;
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
      if (!getBit(exists, u)) continue;
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
      for (K u=0; u < span(); ++u) {
        #pragma omp task if(degree(u) > 2048)
        sortEdges(u);
      }
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
    setBit(exists, u);
  }


  /**
   * Add a vertex to the graph.
   * @param u vertex id
   * @param d associated data of the vertex
   * @note `update()` must be called after all vertices are added.
   */
  inline void addVertex(K u, V d) {
    if (u >= span()) respan(u+1);
    setBit(exists, u);
    values[u] = d;
  }


  /**
   * Remove a vertex from the graph.
   * @param u vertex id
   * @note `update()` must be called after all vertices are removed.
   */
  inline void removeVertex(K u) {
    if (!hasVertex(u)) return;
    clearBit(exists, u);
    values[u] = V();
    clearEdges(u);
  }


  /**
   * Set the degree of a vertex in the graph, without "any" checks.
   * @param u vertex id
   * @param d new degree of the vertex
   * @note Ensure that the vertex exists, and it has at least `d` edges.
   */
  inline void setDegreeUnsafe(K u, K d) {
    degrees[u] = d;
  }


  /**
   * Add an outgoing edge to the graph, without "any" checks.
   * @param u source vertex id
   * @param v target vertex id
   * @param w associated weight of the edge
   * @note Ensure that the vertex exists, and enough capacity is reserved.
   */
  inline void addEdgeUnsafe(K u, K v, E w=E()) {
    auto *ptr = edges[u];
    K i = degrees[u]++;
    ptr[i] = {v, w};
  }


  /**
   * Add an outgoing edge to the graph, without "any" checks [parallel].
   * @param u source vertex id
   * @param v target vertex id
   * @param w associated weight of the edge
   * @note Ensure that the vertex exists, and enough capacity is reserved.
   */
  inline void addEdgeUnsafeOmp(K u, K v, E w=E()) {
    auto *ptr = edges[u];
    K i = K();
    #pragma omp atomic capture
    i = degrees[u]++;
    ptr[i] = {v, w};
  }


  /**
   * Add an outgoing edge to the graph, without uniqueness check.
   * @param u source vertex id
   * @param v target vertex id
   * @param w associated weight of the edge
   * @note Ensure that the span of the graph is sufficient.
   */
  inline void addEdgeUnchecked(K u, K v, E w=E()) {
    if (!getBit(exists, u)) setBit(exists, u);
    if (!getBit(exists, v)) setBit(exists, v);
    K i = degrees[u]++;
    if (i >= capacities[u]) {
      K cap = allocationCapacity(degrees[u]);
      pair<K, E> *tmp = allocate(cap);
      if (i>0) memcpy(tmp, edges[u], i*EDGE);
      if (i>0) deallocate(edges[u], capacities[u]);
      edges[u] = tmp;
      capacities[u] = cap;
    }
    pair<K, E> *ptr = edges[u];
    ptr[i] = {v, w};
  }


  /**
   * Remove outgoing edges from a vertex in the graph.
   * @param u source vertex id
   * @param ib begin iterator of edge keys to remove
   * @param ie end iterator of edge keys to remove
   * @param fl comparison function for edge keys
   * @note [ib, ie) must be sorted and unique.
   * @returns number of edges removed
   */
  template <class I, class FL>
  inline size_t removeEdges(K u, I ib, I ie, FL fl) {
    if (!hasVertex(u)) return 0;
    auto *eb = edges[u], *ee = edges[u] + degrees[u];
    auto  it = set_difference(eb, ee, ib, ie, eb, fl);
    degrees[u] = it - eb;
    return ee - it;
  }


  /**
   * Remove outgoing edges from a vertex in the graph.
   * @param u source vertex id
   * @param ib begin iterator of edge keys to remove
   * @param ie end iterator of edge keys to remove
   * @note [ib, ie) must be sorted and unique.
   * @returns number of edges removed
   */
  template <class I>
  inline size_t removeEdges(K u, I ib, I ie) {
    auto fl = [](const auto& a, const auto& b) { return a.first <  b; };
    return removeEdges(u, ib, ie, fl);
  }


  /**
   * Add outgoing edges to a vertex in the graph.
   * @param u source vertex id
   * @param ib begin iterator of edges to add
   * @param ie end iterator of edges to add
   * @note [ib, ie) must be sorted and unique.
   * @returns number of edges added
   */
  template <class I>
  inline size_t addEdges(K u, I ib, I ie) {
    if (!hasVertex(u) || ib==ie) return 0;
    auto *eb = edges[u], *ee = edges[u] + degrees[u];
    size_t deg = degrees[u] + distance(ib, ie);
    size_t cap = allocationCapacity(deg);
    pair<K, E> *ptr = allocate(cap);
    auto fl = [](const auto& a, const auto& b) { return a.first <  b.first; };
    auto it = set_union(eb, ee, ib, ie, ptr, fl);
    deallocate(eb, capacities[u]);
    edges[u]   = ptr;
    degrees[u] = it - ptr;
    capacities[u] = cap;
    return (it - ptr) - (ee - eb);
  }


  /**
   * Set the memory allocator in use.
   * @param alloc memory allocator
   */
  inline void setAllocator(shared_ptr<allocator_type> alloc) {
    mx = alloc;
  }
  #pragma endregion


  #pragma region CONSTRUCTORS
  public:
  /**
   * Create an empty graph.
   */
  ArenaDiGraph() {}


  /**
   * Destroy the Arena DiGraph.
   */
  ~ArenaDiGraph() {
    clear();
    delete[] exists;
    delete[] edges;
    delete[] degrees;
    delete[] capacities;
    delete[] values;
    mx = nullptr;
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
  using edge_value_type = E;
  /** Offset type (edge offset). */
  using offset_type = O;
  #pragma endregion


  #pragma region DATA
  public:
  /** Offsets of the outgoing edges of vertices. */
  O *offsets = nullptr;
  /** Degree of each vertex. */
  K *degrees = nullptr;
  /** Vertex ids of the outgoing edges of each vertex (lookup using offsets). */
  K *edgeKeys = nullptr;
  /** Edge weights of the outgoing edges of each vertex (lookup using offsets). */
  E *edgeValues = nullptr;
  /** Span of the graph. */
  size_t SPAN = 0;
  /** Vertex values. */
  V *values = nullptr;
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
    return SPAN;
  }

  /**
   * Get the number of vertices in the graph.
   * @returns |V|
   */
  inline size_t order() const noexcept {
    return SPAN;
  }

  /**
   * Obtain the number of edges in the graph.
   * @returns |E|
   */
  inline size_t size() const noexcept {
    if (SPAN == 0) return 0;
    size_t M = 0;
    for (size_t u=0; u<SPAN; ++u)
      M += degrees[u];
    return M;
  }

  /**
   * Check if the graph is empty.
   * @returns is the graph empty?
   */
  inline bool empty() const noexcept {
    return SPAN == 0;
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
  protected:
  /**
   * Free the memory allocated for the CSR representation of the graph.
   */
  inline void freeArrays() {
    delete offsets;
    delete degrees;
    delete edgeKeys;
    delete edgeValues;
    delete values;
  }


  public:
  /**
   * Adjust the order of the graph (or the number of vertices).
   * @param n new order, or number of vertices
   */
  inline void resize(size_t n) {
    if (n <= SPAN) return;
    freeArrays();
    offsets  = new O[n+1];
    degrees  = new K[n];
    edgeKeys = new K[n];
  }


  /**
   * Adjust the order and size of the graph (or the number of vertices and edges).
   * @param n new order, or number of vertices
   * @param m new size, or number of edges
   */
  inline void resize(size_t n, size_t m) {
    if (n <= SPAN && m <= size()) return;
    freeArrays();
    offsets  = new O[n+1];
    degrees  = new K[n];
    values   = new V[n];
    edgeKeys = new K[m];
    edgeValues = new E[m];
    SPAN = n;
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
  DiGraphCsr() {}


  /**
   * Allocate space for CSR representation of a directed graph.
   * @param n number of vertices
   * @param m number of edges
   */
  DiGraphCsr(size_t n, size_t m) {
    resize(n, m);
  }
  #pragma endregion
};
#pragma endregion




#pragma region METHODS
#pragma region SET OPERATIONS
/**
 * Subtract a graph's edges from another graph.
 * @param a graph to subtract from (updated)
 * @param y graph to subtract
 * @returns number of edges removed
 */
template <class H, class G>
inline size_t subtractGraphEdgesU(H& a, const G& y) {
  size_t dM = 0;
  y.forEachVertexKey([&](auto u) {
    if (!a.hasVertex(u)) return;
    auto yb = y.beginEdges(u), ye = y.endEdges(u);
    auto fl = [](const auto& a, const auto& b) { return a.first < b.first; };
    dM += a.removeEdges(u, yb, ye, fl);
  });
  // Update the number of edges.
  a.setSize(a.size() - dM);
  return dM;
}


/**
 * Subtract a graph's edges from another graph [parallel].
 * @param a graph to subtract from (updated)
 * @param y graph to subtract
 * @returns number of edges removed
 */
template <int CHUNK=1024, class H, class G>
inline size_t subtractGraphOmpU(H& a, const G& y) {
  using  K = typename H::key_type;
  size_t S = y.span(), dM = 0;
  #pragma omp parallel for schedule(dynamic, CHUNK) reduction(+:dM)
  for (K u=0; u<S; ++u) {
    if (!y.hasVertex(u) || !a.hasVertex(u)) continue;
    auto yb = y.beginEdges(u), ye = y.endEdges(u);
    auto fl = [](const auto& a, const auto& b) { return a.first < b.first; };
    dM += a.removeEdges(u, yb, ye, fl);
  }
  // Update the number of edges.
  a.setSize(a.size() - dM);
  return dM;
}


/**
 * Subtract a graph's edges from another graph.
 * @param a output graph (output)
 * @param x graph to subtract from
 * @param y graph to subtract
 * @returns number of edges removed
 */
template <class H, class GX, class GY>
inline size_t subtractGraphW(H& a, const GX& x, const GY& y) {
  size_t S = x.span(), dM = 0;
  a.clear();
  a.reserve(S);
  // Add the vertices.
  x.forEachVertex([&](auto u, auto d) {
    a.addVertex(u, d);
  });
  // Reserve space for the edges.
  x.forEachVertexKey([&](auto u) {
    a.allocateEdges(u, x.degree(u));
  });
  // Now add edges of vertices that are untouched.
  x.forEachVertexKey([&](auto u) {
    if (y.hasVertex(u)) return;
    auto xb = x.beginEdges(u), xe = x.endEdges(u);
    auto ab = a.beginEdges(u);
    auto it = copy(xb, xe, ab);
    a.setDegreeUnsafe(u, it - ab);
  });
  // Now add edges of vertices that are touched.
  y.forEachVertexKey([&](auto u) {
    if (!x.hasVertex(u)) return;
    auto xb = x.beginEdges(u), xe = x.endEdges(u);
    auto yb = y.beginEdges(u), ye = y.endEdges(u);
    auto ab = a.beginEdges(u);
    auto fl = [](const auto& a, const auto& b) { return a.first < b.first; };
    auto it = set_difference(xb, xe, yb, ye, ab, fl);
    a.setDegreeUnsafe(u, it - ab);
    dM += (xe - xb) - (it - ab);
  });
  // Update the number of edges.
  a.update(true, true);
  return dM;
}


#ifdef OPENMP
/**
 * Subtract a graph's edges from another graph [parallel].
 * @param a output graph (output)
 * @param x graph to subtract from
 * @param y graph to subtract
 * @returns number of edges removed
 */
template <int CHUNK=1024, class H, class GX, class GY>
inline size_t subtractGraphOmpW(H& a, const GX& x, const GY& y) {
  using  K = typename H::key_type;
  size_t S = x.span(), dM = 0;
  a.clearOmp();
  a.reserveOmp(S);
  // Add the vertices.
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    a.addVertex(u, x.vertexValue(u));
  }
  // Reserve space for the edges.
  #pragma omp parallel for schedule(dynamic, 2048)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    a.allocateEdges(u, x.degree(u));
  }
  // Now add edges of vertices that are untouched.
  #pragma omp parallel for schedule(dynamic, CHUNK)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u) || y.hasVertex(u)) continue;
    auto xb = x.beginEdges(u), xe = x.endEdges(u);
    auto ab = a.beginEdges(u);
    auto it = copy(xb, xe, ab);
    a.setDegreeUnsafe(u, it - ab);
  }
  // Now add edges of vertices that are touched.
  #pragma omp parallel for schedule(dynamic, CHUNK) reduction(+:dM)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u) || !y.hasVertex(u)) continue;
    auto xb = x.beginEdges(u), xe = x.endEdges(u);
    auto yb = y.beginEdges(u), ye = y.endEdges(u);
    auto ab = a.beginEdges(u);
    auto fl = [](const auto& a, const auto& b) { return a.first < b.first; };
    auto it = set_difference(xb, xe, yb, ye, ab, fl);
    a.setDegreeUnsafe(u, it - ab);
    dM += (xe - xb) - (it - ab);
  }
  // Update the number of edges.
  a.updateOmp(true, true);
  return dM;
}
#endif


/**
 * Add a graph's edges to another graph.
 * @param a graph to add to (updated)
 * @param y graph to add
 * @returns number of edges added
 */
template <class H, class G>
inline size_t addGraphU(H& a, const G& y) {
  size_t A = a.span(),   Y = y.span();
  size_t S = max(A, Y), dM = 0;
  if (S!=A) a.respan(S);
  // Add new vertices.
  y.forEachVertex([&](auto u, auto d) {
    a.addVertex(u, d);
  });
  // Add new edges.
  y.forEachVertex([&](auto u, auto d) {
    auto yb = y.beginEdges(u), ye = y.endEdges(u);
    dM += a.addEdges(u, yb, ye);
  });
  // Update the number of edges.
  a.setSize(a.size() + dM);
  return dM;
}


#ifdef OPENMP
/**
 * Add a graph's edges to another graph [parallel].
 * @param a graph to add to (updated)
 * @param y graph to add
 * @returns number of edges added
 */
template <int CHUNK=512, class H, class G>
inline size_t addGraphOmpU(H& a, const G& y) {
  using  K = typename H::key_type;
  size_t A = a.span(),   Y = y.span();
  size_t S = max(A, Y), dM = 0;
  if (S!=A) a.respan(S);
  // Add new vertices.
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    if (!y.hasVertex(u)) continue;
    a.addVertex(u, y.vertexValue(u));
  }
  // Add new edges.
  #pragma omp parallel for schedule(dynamic, CHUNK) reduction(+:dM)
  for (K u=0; u<Y; ++u) {
    if (!y.hasVertex(u)) continue;
    auto yb = y.beginEdges(u), ye = y.endEdges(u);
    dM += a.addEdges(u, yb, ye);
  }
  // Update the number of edges.
  a.setSize(a.size() + dM);
  return dM;
}
#endif


/**
 * Add a graph's edges to another graph.
 * @param a output graph (output)
 * @param x graph to add from
 * @param y graph to add
 * @returns number of edges added
 */
template <class H, class GX, class GY>
inline size_t addGraphW(H& a, const GX& x, const GY& y) {
  using  K = typename H::key_type;
  size_t X = x.span(),   Y = y.span();
  size_t S = max(X, Y), dM = 0;
  a.clear();
  a.reserve(S);
  // Add the vertices.
  for (K u=0; u<S; ++u) {
    if (y.hasVertex(u))      a.addVertex(u, y.vertexValue(u));
    else if (x.hasVertex(u)) a.addVertex(u, x.vertexValue(u));
  }
  // Reserve space for the edges.
  for (K u=0; u<S; ++u)
    a.allocateEdges(u, x.degree(u) + y.degree(u));
  // Now add edges of vertices.
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u) && !y.hasVertex(u)) continue;
    auto xb = x.beginEdges(u), xe = x.endEdges(u);
    auto yb = y.beginEdges(u), ye = y.endEdges(u);
    auto ab = a.beginEdges(u);
    auto fl = [](const auto& a, const auto& b) { return a.first < b.first; };
    auto it = set_union(xb, xe, yb, ye, ab, fl);
    a.setDegreeUnsafe(u, it - ab);
    dM += (it - ab) - (xe - xb);
  }
  // Update the number of edges.
  a.update(true, true);
  return dM;
}


#ifdef OPENMP
/**
 * Add a graph's edges to another graph [parallel].
 * @param a output graph (output)
 * @param x graph to add from
 * @param y graph to add
 */
template <int CHUNK=512, class H, class GX, class GY>
inline size_t addGraphOmpW(H& a, const GX& x, const GY& y) {
  using  K = typename H::key_type;
  size_t X = x.span(),   Y = y.span();
  size_t S = max(X, Y), dM = 0;
  a.clearOmp();
  a.reserveOmp(S);
  // Add the vertices.
  #pragma omp parallel for schedule(static, 2048)
  for (K u=0; u<S; ++u) {
    if (y.hasVertex(u))      a.addVertex(u, y.vertexValue(u));
    else if (x.hasVertex(u)) a.addVertex(u, x.vertexValue(u));
  }
  // Reserve space for the edges.
  #pragma omp parallel for schedule(dynamic, 2048)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u)) continue;
    a.allocateEdges(u, x.degree(u) + y.degree(u));
  }
  // Now add edges of vertices that are touched.
  #pragma omp parallel for schedule(dynamic, CHUNK) reduction(+:dM)
  for (K u=0; u<S; ++u) {
    if (!x.hasVertex(u) && !y.hasVertex(u)) continue;
    auto xb = x.beginEdges(u), xe = x.endEdges(u);
    auto yb = y.beginEdges(u), ye = y.endEdges(u);
    auto ab = a.beginEdges(u);
    auto fl = [](const auto& a, const auto& b) { return a.first < b.first; };
    auto it = set_union(xb, xe, yb, ye, ab, fl);
    a.setDegreeUnsafe(u, it - ab);
    dM += (it - ab) - (xe - xb);
  }
  // Update the number of edges.
  a.updateOmp(true, true);
  return dM;
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
