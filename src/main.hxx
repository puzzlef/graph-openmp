#pragma once
#define BUILD  0  // 5 == BUILD_TRACE
#define OPENMP 1

#ifndef TYPE
#define TYPE double
#endif
#ifndef MAX_THREADS
#define MAX_THREADS 32
#endif
#ifndef REPEAT_BATCH
#define REPEAT_BATCH 5
#endif
#ifndef REPEAT_METHOD
#define REPEAT_METHOD 1
#endif

#include "_main.hxx"
#include "Graph.hxx"
#include "mtx.hxx"
#include "snap.hxx"
#include "vertices.hxx"
#include "edges.hxx"
#include "csr.hxx"
#include "duplicate.hxx"
#include "transpose.hxx"
#include "symmetricize.hxx"
#include "selfLoop.hxx"
#include "deadEnds.hxx"
#include "properties.hxx"
#include "modularity.hxx"
#include "dfs.hxx"
#include "bfs.hxx"
#include "partition.hxx"
#include "random.hxx"
