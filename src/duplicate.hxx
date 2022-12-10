#pragma once




// DUPLICATE
// ---------

template <class H, class G, class FV, class FE>
void duplicateW(H& a, const G& x, FV fv, FE fe) {
  x.forEachVertex([&](auto u, auto d) { if (fv(u)) a.addVertex(u, d); });
  x.forEachVertex([&](auto u, auto d) {
    if (fv(u)) x.forEachEdge(u, [&](auto v, auto w) {
      if (fv(v) && fe(u, v)) a.addEdge(u, v, w);
    });
  });
  a.update();
}
template <class H, class G, class FV>
inline void duplicateW(H& a, const G& x, FV fv) {
  duplicateW(a, x, fv, [](auto u, auto v) { return true; });
}
template <class H, class G>
inline void duplicateW(H& a, const G& x) {
  duplicateW(a, x, [](auto u) { return true; });
}

template <class G, class FV, class FE>
inline auto duplicate(const G& x, FV fv, FE fe) {
  G a; duplicateW(a, x, fv, fe);
  return a;
}
template <class G, class FV>
inline auto duplicate(const G& x, FV fv) {
  G a; duplicateW(a, x, fv);
  return a;
}
template <class G>
inline auto duplicate(const G& x) {
  G a; duplicateW(a, x);
  return a;
}
