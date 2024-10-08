#pragma once
#include <cstdint>
#include <vector>
#include "_mman.hxx"

using std::vector;




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
  vector<void*> free;
  #pragma endregion


  #pragma region METHODS
  public:
  /**
   * Allocate memory.
   * @returns allocated memory, or nullptr if out of memory
   */
  inline void* alloc() {
    // Allocate from free list, if available.
    if (!free.empty()) {
      void *ptr = free.back();
      free.pop_back();
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
    free.push_back(ptr);
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
    free.reserve(min(capacity/SIZE, size_t(128)));
  }
  #pragma endregion
};
#pragma endregion
