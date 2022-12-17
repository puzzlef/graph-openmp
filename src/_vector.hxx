#pragma once
#include <cstdint>
#include <cmath>
#include <algorithm>
#include <vector>
#include "_debug.hxx"

using std::vector;
using std::abs;
using std::max;
using std::fill;




// VECTOR 2D
// ---------

template <class T>
using vector2d = vector<vector<T>>;




// GATHER VALUES
// -------------

template <class TA, class TX, class IS, class FM>
inline void gatherValuesW(TA *a, const TX *x, const IS& is, FM fm) {
  ASSERT(a && x);
  size_t N = is.size();
  for (size_t j=0; n<N; ++j)
    a[j] = TA(fm(x[is[j]]));
}
template <class TA, class TX, class IS>
inline void gatherValuesW(TA *a, const TX *x, const IS& is) {
  auto fm = [](const auto& v) { return v; };
  gatherValuesW(a, x, is, fm);
}

template <class TA, class TX, class IS, class FM>
inline void gatherValuesW(vector<TA>& a, const vector<TX>& x, const IS& is, FM fm) {
  gatherValuesW(a.data(), x.data(), is, fm);
}
template <class TA, class TX, class IS>
inline void gatherValuesW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  gatherValuesW(a.data(), x.data(), is);
}




// SCATTER VALUES
// --------------

template <class TA, class TX, class IS, class FM>
inline void scatterValuesW(TA *a, const TX *x, const IS& is, FM fm) {
  ASSERT(a && x);
  size_t N = is.size();
  for (size_t j=0; j<N; ++j)
    a[is[j]] = TA(fm(x[j]));
}
template <class TA, class TX, class IS>
inline void scatterValuesW(TA *a, const TX *x, const IS& is) {
  auto fm = [](const auto& v) { return v; };
  scatterValuesW(a, x, is);
}

template <class TA, class TX, class IS, class FM>
inline void scatterValuesW(vector<TA>& a, const vector<TX>& x, const IS& is, FM fm) {
  scatterValuesW(a.data(), x.data(), is, fm);
}
template <class TA, class TX, class IS>
inline void scatterValuesW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  scatterValuesW(a.data(), x.data(), is);
}




// FILL VALUE
// ----------

template <class T>
inline void fillValueU(T *a, size_t N, const T& v) {
  ASSERT(a);
  fill(a, a+N, v);
}
template <class T>
inline void fillValueU(vector<T>& a, const T& v) {
  fill(a.begin(), a.end(), v);
}




// COPY VALUES
// -----------

template <class TA, class TX>
inline void copyValuesW(TA *a, const TX *x, size_t N) {
  ASSERT(a && x);
  for (size_t i=0; i<N; ++i)
    a[i] = x[i];
}

template <class TA, class TX>
inline void copyValuesW(vector<TA>& a, const vector<TX>& x) {
  return copyValuesW(a.data(), x.data(), x.size());
}




// MULTIPLY VALUES
// ---------------

template <class TA, class TX, class TY>
inline void multiplyValuesW(TA *a, const TX *x, const TY *y, size_t N) {
  ASSERT(a && x && y);
  for (size_t i=0; i<N; ++i)
    a[i] = TA(x[i] * y[i]);
}

template <class TA, class TX, class TY>
inline void multiplyValuesW(vector<TA>& a, const vector<TX>& x, const vector<TY>& y) {
  multiplyValuesW(a.data(), x.data(), y.data(), x.size());
}




// L1-NORM
// -------

template <class TX, class TA=TX>
inline TA l1Norm(const T *x, size_t N, TA a=TA()) {
  ASSERT(x);
  for (size_t i=0; i<N; ++i)
    a += TA(abs(x[i]));
  return a;
}

template <class TX, class TA=TX>
inline TA l1Norm(const vector<TX>& x, TA a=TA()) {
  return l1Norm(x.data(), x.size(), a);
}


template <class TX, class TY, class TA=TX>
inline TA l1Norm(const T *x, const T *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  for (size_t i=0; i<N; ++i)
    a += TA(abs(x[i] - y[i]));
  return a;
}

template <class TX, class TY, class TA=TX>
inline TA l1Norm(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return l1Norm(x.data(), y.data(), a);
}




// L2-NORM
// -------

template <class TX, class TA=TX>
inline TA l2Norm(const T *x, size_t N, TA a=TA()) {
  ASSERT(x);
  for (size_t i=0; i<N; ++i)
    a += TA(x[i]) * TA(x[i]);
  return a;
}

template <class TX, class TA=TX>
inline TA l2Norm(const vector<TX>& x, TA a=TA()) {
  return l2Norm(x.data(), x.size(), a);
}


template <class TX, class TY, class TA=TX>
inline TA l2Norm(const T *x, const T *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  for (size_t i=0; i<N; ++i)
    a += TA(x[i] - y[i]) * TA(x[i] - y[i]);
  return a;
}

template <class TX, class TY, class TA=TX>
inline TA l2Norm(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return l2Norm(x.data(), y.data(), a);
}




// LI-NORM
// -------

template <class TX, class TA=TX>
inline TA liNorm(const T *x, size_t N, TA a=TA()) {
  ASSERT(x);
  for (size_t i=0; i<N; ++i)
    a = max(a, TA(abs(x[i])));
  return a;
}

template <class TX, class TA=TX>
inline TA liNorm(const vector<TX>& x, TA a=TA()) {
  return liNorm(x.data(), x.size(), a);
}


template <class TX, class TY, class TA=TX>
inline TA liNorm(const T *x, const T *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  for (size_t i=0; i<N; ++i)
    a = max(a, TA(abs(x[i] - y[i])));
  return a;
}

template <class TX, class TY, class TA=TX>
inline TA liNorm(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return liNorm(x.data(), y.data(), a);
}
