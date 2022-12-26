#pragma once
#include <vector>

using std::vector;




// DFS
// ---
// Depth First Search (DFS) graph traversal.

/**
 * Find vertices visited with DFS.
 * @param vis vertex visited? (updated)
 * @param x original graph
 * @param u start vertex
 * @param fn action to perform on every visited vertex
 */
template <class G, class K, class F>
inline void dfsVisitedForEachW(vector<bool>& vis, const G& x, K u, F fn) {
  if (vis[u]) return;
  vis[u] = true; fn(u);
  x.forEachEdgeKey(u, [&](K v) {
    if (!vis[v]) dfsVisitedForEachW(vis, x, v, fn);
  });
}
template <class G, class K, class F>
inline vector<bool> dfsVisitedForEach(const G& x, K u, F fn) {
  vector<bool> vis(x.span());
  dfsVisitedForEachW(vis, x, u, fn);
  return vis;
}
