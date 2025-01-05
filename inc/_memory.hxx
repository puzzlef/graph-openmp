#pragma once
#include <cstdint>
#include <cstdlib>
#include <vector>
#include <algorithm>

using std::vector;
using std::min;




#pragma region CONSTANTS
/** Expected page size. */
#define PAGE_SIZE 4096
/** Expected cache line size. */
#define CACHE_LINE_SIZE 128
/** Minimum number of freed allocations expected. */
#define MIN_FREED 128
/** Minimum number of memory pools expected. */
#define MIN_POOLS 16
#pragma endregion




#pragma region CLASSES
/**
 * A fixed-capacity Arena Allocator.
 * @tparam SIZE size of each allocation
 * @tparam CAPACITY capacity of the memory pool
 */
template <size_t SIZE, size_t CAPACITY>
class FixedArenaAllocator {
  #pragma region CONSTANTS
  public:
  /** Size of each allocation. */
  static constexpr size_t allocation_size = SIZE;
  /** Size of the memory pool. */
  static constexpr size_t pool_size = CAPACITY * SIZE;
  #pragma endregion


  #pragma region DATA
  /** These allocations have been freed, and can be reused. */
  vector<void*> freed {};
  /** Number of bytes used in the memory pool. */
  size_t used = 0;
  /** The memory pool. */
  void *pool;
  #pragma endregion


  #pragma region METHODS
  public:
  /**
   * Allocate memory.
   * @returns allocated memory, or nullptr if out of memory
   */
  inline void* allocate() noexcept {
    // Allocate from free list, if available.
    if (!freed.empty()) {
      void *ptr = freed.back();
      freed.pop_back();
      return ptr;
    }
    // Allocate from pool.
    if (used < pool_size) {
      void *ptr = (char*) pool + used;
      used += SIZE;
      return ptr;
    }
    return nullptr;
  }


  /**
   * Free memory.
   * @param ptr memory to free
   */
  inline void deallocate(void *ptr) {
    freed.push_back(ptr);
  }


  /**
   * Free all memory.
   */
  inline void reset() noexcept {
    used = 0;
    freed.clear();
  }
  #pragma endregion


  #pragma region CONSTRUCTORS
  public:
  /**
   * Create a fixed-capacity Arena Allocator.
   * @param pool memory pool
   */
  FixedArenaAllocator(void *pool)
  : pool(pool) {
    freed.reserve(min(CAPACITY, size_t(MIN_FREED)));
  }
  #pragma endregion
};




/**
 * A variable-capacity Arena Allocator.
 * @tparam SIZE size of each allocation
 * @tparam CAPACITY capacity of each memory pool
 */
template <size_t SIZE, size_t CAPACITY>
class ArenaAllocator {
  #pragma region CONSTANTS
  public:
  /** Size of each allocation. */
  static constexpr size_t allocation_size = SIZE;
  /** Size of each memory pool. */
  static constexpr size_t pool_size = CAPACITY * SIZE;
  #pragma endregion


  #pragma region DATA
  /** These allocations have been freed, and can be reused. */
  vector<void*> freed {};
  /** Number of bytes used in the last memory pool. */
  size_t used = pool_size;
  /** Memory pools. */
  vector<void*> pools {};
  #pragma endregion


  #pragma region METHODS
  public:
  /**
   * Allocate memory.
   * @returns allocated memory, or nullptr if out of memory
   */
  inline void* allocate() {
    // Allocate from free list, if available.
    if (!freed.empty()) {
      void *ptr = freed.back();
      freed.pop_back();
      return ptr;
    }
    // Allocate from pool.
    if (used < pool_size) {
      void *ptr = (char*) pools.back() + used;
      used += SIZE;
      return ptr;
    }
    // Allocate a new pool.
    void *pool = new char[pool_size];
    if (pool) {
      pools.push_back(pool);
      used = SIZE;
      return pool;
    }
    return nullptr;
  }


  /**
   * Free memory.
   * @param ptr memory to free
   */
  inline void deallocate(void *ptr) {
    freed.push_back(ptr);
  }


  /**
   * Free all memory.
   */
  inline void reset() noexcept {
    for (void *pool : pools)
      delete[] (char*) pool;
    pools.clear();
    used = pool_size;
    freed.clear();
  }
  #pragma endregion


  #pragma region CONSTRUCTORS
  public:
  /**
   * Create a variable-capacity Arena Allocator.
   */
  ArenaAllocator() {
    pools.reserve(MIN_POOLS);
    freed.reserve(min(CAPACITY, size_t(MIN_FREED)));
  }


  /**
   * Destroy the Arena Allocator.
   */
  ~ArenaAllocator() {
    reset();
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
template <class T, size_t ALIGN=PAGE_SIZE>
inline constexpr size_t bytesof(size_t N) {
  return ((N * sizeof(T) + ALIGN-1) / ALIGN) * ALIGN;
}
#pragma endregion
