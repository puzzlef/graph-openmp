#pragma once
#include <cstdint>
#include <atomic>
#include <vector>
#include <thread>
#include <algorithm>
#include "_mman.hxx"

using std::memory_order_relaxed;
using std::memory_order_acquire;
using std::memory_order_release;
using std::atomic_flag;
using std::atomic;
using std::vector;
using std::this_thread::yield;
using std::min;




#pragma region CLASSES
/**
 * A thread-safe Arena Allocator, that supports fast allocation and deallocation.
 * @tparam SIZE size of each allocation
 * @tparam MIN_FREED minimum number of freed allocations to keep [128]
 * @note Use for allocations larger than cache line size to avoid false sharing.
 */
template <size_t SIZE, size_t MIN_FREED=128>
class ConcurrentArenaAllocator {
  #pragma region CONSTANTS
  public:
  /** Size of each allocation. */
  static constexpr size_t allocation_size = SIZE;
  #pragma endregion


  #pragma region DATA
  /** These allocations have been freed, and can be reused. */
  vector<void*> freed;
  /** Mutex to protect the freed list. */
  atomic_flag mutex {};
  /** The memory pool. */
  void  *pool;
  /** Number of bytes used in the memory pool. */
  atomic<size_t> used = 0;
  /** Size of the memory pool. */
  size_t capacity;
  #pragma endregion


  #pragma region METHODS
  public:
  /**
   * Allocate a memory block.
   * @returns allocated memory, or nullptr if out of memory
   */
  inline void* allocate() noexcept {
    // Allocate from freed list, if available, and no other thread is using it.
    if (!freed.empty() && !mutex.test_and_set(memory_order_acquire)) {
      void *ptr = freed.back();
      freed.pop_back();
      mutex.clear(memory_order_release);
      return ptr;
    }
    // Else, try to allocate from pool.
    if (used < capacity) return (char*) pool + used.fetch_add(SIZE, memory_order_relaxed);
    // Else, must try from freed list.
    while (!freed.empty()) {
      if (mutex.test_and_set(memory_order_acquire)) { yield(); continue; }
      void *ptr = freed.back();
      freed.pop_back();
      mutex.clear(memory_order_release);
      return ptr;
    }
    return nullptr;
  }


  /**
   * Free a memory block.
   * @param ptr memory to free
   */
  inline void free(void *ptr) {
    // Add to freed list.
    while (mutex.test_and_set(memory_order_acquire))
      yield();
    freed.push_back(ptr);
    mutex.clear(memory_order_release);
  }


  /**
   * Reset the allocator, freeing all memory.
   */
  inline void reset() noexcept {
    // Reset freed list.
    while (mutex.test_and_set(memory_order_acquire))
      yield();
    freed.clear();
    mutex.clear(memory_order_release);
    // Reset used bytes.
    used.store(0);
  }
  #pragma endregion


  #pragma region CONSTRUCTORS
  public:
  /**
   * Create a new Concurrent Arena Allocator.
   * @param pool memory pool
   * @param capacity size of memory pool
   */
  ConcurrentArenaAllocator(void *pool, size_t capacity)
  : pool(pool), capacity(capacity) {
    // Reserve space for freed list.
    freed.reserve(min(capacity/SIZE, MIN_FREED));
  }
  #pragma endregion
};




/**
 * A Recursive Arena Allocator, that supports fast allocation and deallocation.
 * @tparam SIZE size of each allocation
 * @tparam MIN_FREED minimum number of freed allocations to keep [32]
 * @note This allocator is not thread-safe. Use for allocations smaller than cache line size, on each thread.
 */
template <size_t SIZE, size_t MIN_FREED=32>
class RecursiveArenaAllocator {
  #pragma region CONSTANTS
  public:
  /** Size of each allocation. */
  static constexpr size_t allocation_size = SIZE;
  #pragma endregion


  #pragma region DATA
  /** These allocations have been freed, and can be reused. */
  vector<void*> freed;
  /** The memory pools. */
  vector<void*> pools;
  /** Number of bytes used in the last memory pool. */
  size_t used = 0;
  /** Size of each memory pool. */
  size_t capacity;
  #pragma endregion


  #pragma region METHODS
  public:
  /**
   * Allocate a memory block.
   * @param fallocate function to allocate a memory pool, e.g., fallocate(capacity)
   * @returns allocated memory, or nullptr if out of memory
   */
  template <class FA>
  inline void* allocate(FA fallocate) {
    // Allocate from freed list, if available.
    if (!freed.empty()) {
      void *ptr = freed.back();
      freed.pop_back();
      return ptr;
    }
    // Else, try to allocate from last pool.
    if (!pools.empty() && used < capacity) {
      void *ptr = (char*) pools.back() + used;
      used += SIZE;
      return ptr;
    }
    // Else, must allocate a new pool.
    void *pool = fallocate(capacity);
    pools.push_back(pool);
    used = SIZE;
    return pool;
  }


  /**
   * Free a memory block.
   * @param ptr memory to free
   */
  inline void free(void *ptr) {
    freed.push_back(ptr);
  }


  /**
   * Reset the allocator, freeing all memory.
   * @param ffree function to free a memory pool, e.g., ffree(pool)
   */
  template <class FF>
  inline void reset(FF ffree) noexcept {
    // Reset freed list, used bytes.
    freed.clear();
    used = 0;
    // Free all memory pools.
    for (void *pool : pools)
      ffree(pool);
    pools.clear();
  }
  #pragma endregion


  #pragma region CONSTRUCTORS
  public:
  /**
   * Create a new Recursive Arena Allocator.
   * @param capacity size of each memory pool
   */
  RecursiveArenaAllocator(size_t capacity)
  : capacity(capacity) {
    // Reserve space for freed list.
    freed.reserve(min(capacity/SIZE, MIN_FREED));
  }
  #pragma endregion
};
#pragma endregion




#pragma region METHODS
/**
 * Get the number of bytes required to store N elements of type T.
 * @tparam ALIGN memory alignment size [page-aligned]
 * @param N number of elements to store
 * @returns number of bytes required
 */
template <class T, size_t ALIGN=4096>
inline constexpr size_t bytesof(size_t N) {
  return ((N * sizeof(T) + ALIGN-1) / ALIGN) * ALIGN;
}
#pragma endregion
