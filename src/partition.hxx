#pragma once
#include <vector>
#include "_main.hxx"
#include "bfs.hxx"

using std::vector;




#pragma region METHODS
#pragma region PARTITION BY ID
/**
 * Partition vertices contiguously by id.
 * @param a vertex keys grouped by partition (output)
 * @param x original graph
 * @param p start partition
 * @param P stop partition
 */
template <class G>
inline void partitionByIdW(vector<int>& a, const G& x, int p, int P) {
  using  K = typename G::key_type;
  size_t N = x.order(), i = 0;
  int    Q = P - p;
  x.forEachVertexKey([&](auto u) {
    a[u] = p + int(i*Q/N); ++i;
  });
}


/**
 * Partition vertices contiguously by id.
 * @param x original graph
 * @param p start partition
 * @param P stop partition
 * @returns vertex keys grouped by partition
 */
template <class G>
inline vector<int> partitionById(const G& x, int p, int P) {
  vector<int> a(x.span());
  partitionByIdW(a, x, p, P);
  return a;
}
#pragma endregion




#pragma region PARTITION BY BFS
/**
 * Partition vertices by a custom cost function using BFS.
 * @tparam FLAG visited flag type
 * @param a vertex keys grouped by partition (output)
 * @param x original graph
 * @param p start partition
 * @param P stop partition
 * @param fc cost function
 * @param ft cost limit test function
 */
template <class FLAG=bool, class G, class FC, class FT>
inline void partitionByBfsW(vector<int>& a, const G& x, int p, int P, FC fc, FT ft) {
  using  K = typename G::key_type;
  using  C = decltype(fc(K()));
  using  B = FLAG;
  size_t N = x.order();
  int    Q = P - p, i = 0;
  C      c = C();
  vector<B> vis(x.span());
  vector<K> us, vs;
  auto fp = [&](auto v, auto d) { a[v] = p + i; c += fc(v); };
  auto fu = [&](auto v, auto d) { return ft(c); };
  x.forEachVertexKey([&](auto u) {
    if (vis[u]) return;
    us.clear(); us.push_back(u);
    bfsVisitedForEachU(vis, us, vs, x, fu, fp);
    if (!ft(c)) { i = (i + 1) % Q; c = C(); }
  });
}


/**
 * Partition vertices by a custom cost function using BFS.
 * @tparam FLAG visited flag type
 * @param x original graph
 * @param p start partition
 * @param P stop partition
 * @param fc cost function
 * @param ft cost limit test function
 * @returns vertex keys grouped by partition
 */
template <class FLAG=bool, class G, class FC, class FT>
inline vector<int> partitionByBfs(const G& x, int p, int P, FC fc, FT ft) {
  vector<int> a(x.span());
  partitionByBfsW<FLAG>(a, x, p, P, fc, ft);
  return a;
}
#pragma endregion
#pragma endregion
