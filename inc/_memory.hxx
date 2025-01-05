#pragma once
#include <cstdint>
#include <vector>
#include <algorithm>
#include "_cmath.hxx"
#ifdef OPENMP
#include <omp.h>
#endif

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
/** Capacity of the power-of-two memory pool, in bytes. */
#define POW2_POOL_SIZE 65536
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
  protected:
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
  protected:
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




/**
 * A fast power-of-two size allocator, based on variable-capacity Arena Allocator.
 * @tparam CAPACITY capacity of each memory pool, in bytes
 */
template <size_t CAPACITY=POW2_POOL_SIZE>
class Pow2Allocator {
  #pragma region CONSTANTS
  public:
  /** Size of each allocation (variable). */
  static constexpr size_t allocation_size = 1;
  /** Size of each memory pool. */
  static constexpr size_t pool_size = CAPACITY;
  #pragma endregion


  #pragma region DATA
  protected:
  /** Arena allocator for 16-byte allocations. */
  ArenaAllocator<16, CAPACITY/16> a16;
  /** Arena allocator for 32-byte allocations. */
  ArenaAllocator<32, CAPACITY/32> a32;
  /** Arena allocator for 64-byte allocations. */
  ArenaAllocator<64, CAPACITY/64> a64;
  /** Arena allocator for 128-byte allocations. */
  ArenaAllocator<128, CAPACITY/128> a128;
  /** Arena allocator for 256-byte allocations. */
  ArenaAllocator<256, CAPACITY/256> a256;
  /** Arena allocator for 512-byte allocations. */
  ArenaAllocator<512, CAPACITY/512> a512;
  /** Arena allocator for 1024-byte allocations. */
  ArenaAllocator<1024, CAPACITY/1024> a1024;
  /** Arena allocator for 2048-byte allocations. */
  ArenaAllocator<2048, CAPACITY/2048> a2048;
  #pragma endregion


  #pragma region METHODS
  public:
  /**
   * Get recommended allocation bytes for a given size.
   * @param n number of bytes needed
   * @returns recommended allocation bytes (power of two or multiple of pages, and >= 16)
   */
  static inline size_t allocationCapacity(size_t n) noexcept {
    if (n <= 16) return 16;
    if (n <= PAGE_SIZE) return nextPow2(n);
    return ceilDiv(n, size_t(PAGE_SIZE)) * PAGE_SIZE;
  }


  /**
   * Allocate memory.
   * @param n number of bytes to allocate (must be a power of two or multiple of pages, and >= 16)
   * @returns allocated memory, or nullptr if out of memory
   */
  inline void* allocate(size_t n) {
    switch (n) {
      case 16:   return a16.allocate();
      case 32:   return a32.allocate();
      case 64:   return a64.allocate();
      case 128:  return a128.allocate();
      case 256:  return a256.allocate();
      case 512:  return a512.allocate();
      case 1024: return a1024.allocate();
      case 2048: return a2048.allocate();
      default:   return new char[n];
    }
  }


  /**
   * Free memory.
   * @param ptr memory to free
   * @param n number of bytes allocated
   */
  inline void deallocate(void *ptr, size_t n) {
    switch (n) {
      case 16:   a16.deallocate(ptr); break;
      case 32:   a32.deallocate(ptr); break;
      case 64:   a64.deallocate(ptr); break;
      case 128:  a128.deallocate(ptr); break;
      case 256:  a256.deallocate(ptr); break;
      case 512:  a512.deallocate(ptr); break;
      case 1024: a1024.deallocate(ptr); break;
      case 2048: a2048.deallocate(ptr); break;
      default:   delete[] (char*) ptr;
    }
  }


  /**
   * Free all memory.
   */
  inline void reset() noexcept {
    a16.reset();
    a32.reset();
    a64.reset();
    a128.reset();
    a256.reset();
    a512.reset();
    a1024.reset();
    a2048.reset();
  }
  #pragma endregion


  #pragma region CONSTRUCTORS
  public:
  /**
   * Create a fast power-of-two size allocator.
   */
  Pow2Allocator() = default;


  /**
   * Destroy the Pow2 Allocator.
   */
  ~Pow2Allocator() {
    reset();
  }
  #pragma endregion
};




/**
 * A fast thread-safe power-of-two size allocator, which supports multiple threads.
 * @tparam CAPACITY capacity of each memory pool, in bytes
 */
template <size_t CAPACITY=POW2_POOL_SIZE>
class ConcurrentPow2Allocator {
  #pragma region CONSTANTS
  public:
  /** Size of each allocation (variable). */
  static constexpr size_t allocation_size = 1;
  /** Size of each memory pool. */
  static constexpr size_t pool_size = CAPACITY;
  #pragma endregion


  #pragma region DATA
  protected:
  vector<Pow2Allocator<CAPACITY>*> a;
  #pragma endregion


  #pragma region METHODS
  public:
  /**
   * Get recommended allocation bytes for a given size.
   * @param n number of bytes needed
   * @returns recommended allocation bytes (power of two or multiple of pages, and >= 16)
   */
  static inline size_t allocationCapacity(size_t n) noexcept {
    return Pow2Allocator<CAPACITY>::allocationCapacity(n);
  }


  /**
   * Allocate memory.
   * @param n number of bytes to allocate (must be a power of two or multiple of pages, and >= 16)
   * @returns allocated memory, or nullptr if out of memory
   */
  inline void* allocate(size_t n) {
    int t = omp_get_thread_num();
    return a[t]->allocate(n);
  }


  /**
   * Free memory.
   * @param ptr memory to free
   * @param n number of bytes allocated
   */
  inline void deallocate(void *ptr, size_t n) {
    int t = omp_get_thread_num();
    a[t]->deallocate(ptr, n);
  }


  /**
   * Free all memory.
   */
  inline void reset() noexcept {
    int T = a.size();
    for (int t=0; t<T; ++t)
      a[t]->reset();
  }
  #pragma endregion


  #pragma region CONSTRUCTORS
  public:
  /**
   * Create a fast thread-safe power-of-two size allocator.
   */
  ConcurrentPow2Allocator() {
    int T = omp_get_max_threads();
    a.resize(T);
    for (int t=0; t<T; ++t)
      a[t] = new Pow2Allocator<CAPACITY>();
  }

  /**
   * Destroy the Concurrent Pow2 Allocator.
   */
  ~ConcurrentPow2Allocator() {
    reset();
    int T = a.size();
    for (int t=0; t<T; ++t)
      delete a[t];
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
