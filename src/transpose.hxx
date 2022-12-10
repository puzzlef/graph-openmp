#pragma once
#include "Graph.hxx"




// TRANSPOSE
// ---------

template <class H, class G>
void transposeW(H& a, const G& x) {
  x.forEachVertex([&](auto u, auto d) { a.addVertex(u, d); });
  x.forEachVertexKey([&](auto u) {
    x.forEachEdge(u, [&](auto v, auto w) { a.addEdge(v, u, w); });
  });
  a.update();
}

template <class G>
auto transpose(const G& x) {
  G a; transposeW(a, x);
  return a;
}




// TRANSPOSE WITH DEGREE
// ---------------------

template <class H, class G>
void transposeWithDegreeW(H& a, const G& x) {
  x.forEachVertexKey([&](auto u) { a.addVertex(u, x.degree(u)); });
  x.forEachVertexKey([&](auto u) {
    x.forEachEdge(u, [&](auto v, auto w) { a.addEdge(v, u, w); });
  });
  a.update();
}

template <class G>
auto transposeWithDegree(const G& x) {
  using K = typename G::key_type;
  using H = decltype(retype(x, K(), K()));
  H a; transposeWithDegreeW(a, x);
  return a;
}
