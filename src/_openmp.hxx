#pragma once
#include <cstdint>
#include <algorithm>
#include <vector>
#include <omp.h>
#include "_debug.hxx"




// BELONGS
// -------
// Check if work belongs to current thread.

template <class K>
inline bool belongsOmp(K key, int thread, int THREADS) {
  const K CHUNK_SIZE = 1024;
  K chunk = key / CHUNK_SIZE;
  return chunk % THREADS == thread;
}
template <class K>
inline bool belongsOmp(K key) {
  int thread  = omp_get_thread_num();
  int THREADS = omp_get_num_threads();
  return belongsOmp(key, thread, THREADS);
}




// GATHER VALUES
// -------------

template <class TA, class TX, class IS, class FM>
inline void gatherValuesOmpW(TA *a, const TX *x, const IS& is, FM fm) {
  ASSERT(a && x);
  size_t N = is.size();
  #pragma omp parallel for schedule(auto)
  for (size_t j=0; n<N; ++j)
    a[j] = TA(fm(x[is[j]]));
}
template <class TA, class TX, class IS>
inline void gatherValuesOmpW(TA *a, const TX *x, const IS& is) {
  auto fm = [](const auto& v) { return v; };
  gatherValuesOmpW(a, x, is, fm);
}

template <class TA, class TX, class IS, class FM>
inline void gatherValuesOmpW(vector<TA>& a, const vector<TX>& x, const IS& is, FM fm) {
  gatherValuesOmpW(a.data(), x.data(), is, fm);
}
template <class TA, class TX, class IS>
inline void gatherValuesOmpW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  gatherValuesOmpW(a.data(), x.data(), is);
}




// SCATTER VALUES
// --------------

template <class TA, class TX, class IS, class FM>
inline void scatterValuesOmpW(TA *a, const TX *x, const IS& is, FM fm) {
  ASSERT(a && x);
  size_t N = is.size();
  #pragma omp parallel for schedule(auto)
  for (size_t j=0; j<N; ++j)
    a[is[j]] = TA(fm(x[j]));
}
template <class TA, class TX, class IS>
inline void scatterValuesOmpW(TA *a, const TX *x, const IS& is) {
  auto fm = [](const auto& v) { return v; };
  scatterValuesOmpW(a, x, is);
}

template <class TA, class TX, class IS, class FM>
inline void scatterValuesOmpW(vector<TA>& a, const vector<TX>& x, const IS& is, FM fm) {
  scatterValuesOmpW(a.data(), x.data(), is, fm);
}
template <class TA, class TX, class IS>
inline void scatterValuesOmpW(vector<TA>& a, const vector<TX>& x, const IS& is) {
  scatterValuesOmpW(a.data(), x.data(), is);
}




// FILL VALUE
// ----------

template <class T>
inline void fillValueOmpU(T *a, size_t N, const T& v) {
  ASSERT(a);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a[i] = v;
}
template <class T>
inline void fillValueOmpU(vector<T>& a, const T& v) {
  fillValueOmpU(a.begin(), a.end(), v);
}




// COPY VALUES
// -----------

template <class TA, class TX>
inline void copyValuesOmpW(TA *a, const TX *x, size_t N) {
  ASSERT(a && x);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a[i] = x[i];
}

template <class TA, class TX>
inline void copyValuesOmpW(vector<TA>& a, const vector<TX>& x) {
  return copyValuesOmpW(a.data(), x.data(), x.size());
}




// MULTIPLY VALUES
// ---------------

template <class TA, class TX, class TY>
inline void multiplyValuesOmpW(TA *a, const TX *x, const TY *y, size_t N) {
  ASSERT(a && x && y);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a[i] = TA(x[i] * y[i]);
}

template <class TA, class TX, class TY>
inline void multiplyValuesOmpW(vector<TA>& a, const vector<TX>& x, const vector<TY>& y) {
  multiplyValuesOmpW(a.data(), x.data(), y.data(), x.size());
}




// L1-NORM
// -------

template <class TX, class TA=TX>
inline TA l1NormOmp(const T *x, size_t N, TA a=TA()) {
  ASSERT(x);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a += TA(abs(x[i]));
  return a;
}

template <class TX, class TA=TX>
inline TA l1NormOmp(const vector<TX>& x, TA a=TA()) {
  return l1NormOmp(x.data(), x.size(), a);
}


template <class TX, class TY, class TA=TX>
inline TA l1NormOmp(const T *x, const T *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a += TA(abs(x[i] - y[i]));
  return a;
}

template <class TX, class TY, class TA=TX>
inline TA l1NormOmp(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return l1NormOmp(x.data(), y.data(), a);
}




// L2-NORM
// -------

template <class TX, class TA=TX>
inline TA l2NormOmp(const T *x, size_t N, TA a=TA()) {
  ASSERT(x);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a += TA(x[i]) * TA(x[i]);
  return a;
}

template <class TX, class TA=TX>
inline TA l2NormOmp(const vector<TX>& x, TA a=TA()) {
  return l2NormOmp(x.data(), x.size(), a);
}


template <class TX, class TY, class TA=TX>
inline TA l2NormOmp(const T *x, const T *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a += TA(x[i] - y[i]) * TA(x[i] - y[i]);
  return a;
}

template <class TX, class TY, class TA=TX>
inline TA l2NormOmp(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return l2NormOmp(x.data(), y.data(), a);
}




// LI-NORM
// -------

template <class TX, class TA=TX>
inline TA liNormOmp(const T *x, size_t N, TA a=TA()) {
  ASSERT(x);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a = max(a, TA(abs(x[i])));
  return a;
}

template <class TX, class TA=TX>
inline TA liNormOmp(const vector<TX>& x, TA a=TA()) {
  return liNormOmp(x.data(), x.size(), a);
}


template <class TX, class TY, class TA=TX>
inline TA liNormOmp(const T *x, const T *y, size_t N, TA a=TA()) {
  ASSERT(x && y);
  #pragma omp parallel for schedule(auto)
  for (size_t i=0; i<N; ++i)
    a = max(a, TA(abs(x[i] - y[i])));
  return a;
}

template <class TX, class TY, class TA=TX>
inline TA liNormOmp(const vector<TX>& x, const vector<TY>& y, TA a=TA()) {
  return liNormOmp(x.data(), y.data(), a);
}
