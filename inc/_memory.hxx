#pragma once
#include <cstdint>
#include <vector>
#include <algorithm>
#include "_mman.hxx"

using std::vector;
using std::min;




#pragma region CLASSES
/**
 * A fixed-capacity Arena Allocator.
 * @tparam SIZE size of each allocation
 */
template <size_t SIZE>
class FixedArenaAllocator {
  #pragma region CONSTANTS
  public:
  /** Size of each allocation. */
  static constexpr size_t allocation_size = SIZE;
  #pragma endregion


  #pragma region DATA
  /** The memory pool. */
  void  *pool;
  /** Size of the memory pool. */
  size_t capacity;
  /** Number of bytes used in the memory pool. */
  size_t used;
  /** These allocations have been freed, and can be reused. */
  vector<void*> freed;
  #pragma endregion


  #pragma region METHODS
  public:
  /**
   * Allocate memory.
   * @returns allocated memory, or nullptr if out of memory
   */
  inline void* alloc() {
    // Allocate from free list, if available.
    if (!freed.empty()) {
      void *ptr = freed.back();
      freed.pop_back();
      return ptr;
    }
    // Allocate from pool.
    if (used < capacity) {
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
  inline void free(void *ptr) {
    freed.push_back(ptr);
  }
  #pragma endregion


  #pragma region CONSTRUCTORS
  public:
  /**
   * Create a fixed-capacity Arena Allocator.
   * @param pool memory pool
   * @param capacity size of memory pool
   */
  FixedArenaAllocator(void *pool, size_t capacity) :
  pool(pool), capacity(capacity), used(0) {
    freed.reserve(min(capacity/SIZE, size_t(128)));
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
