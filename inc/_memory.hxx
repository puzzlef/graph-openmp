#pragma once
#include <cstdint>
#include <memory>
#include <atomic>
#include <vector>
#include <thread>
#include <algorithm>
#include <omp.h>
#include "_mman.hxx"
#include "_vector.hxx"

using std::memory_order_relaxed;
using std::memory_order_acquire;
using std::memory_order_release;
using std::memory_order_acq_rel;
using std::unique_ptr;
using std::atomic_flag;
using std::atomic;
using std::vector;
using std::this_thread::yield;
using std::make_unique;
using std::min;




#pragma region CONSTANTS
/** Size of a cache line. */
#define CACHE_LINE_SIZE 128
#pragma endregion




#pragma region CLASSES
template <class T, uint32_t CAPACITY>
struct MultiConsumerRingBuffer {
  #pragma region TYPES, CONSTANTS
  public:
  /** Element type. */
  using value_type = T;
  /** Capacity of the buffer. */
  static constexpr uint32_t capacity = CAPACITY;
  #pragma endregion


  #pragma region DATA
  public:
  /** Head of the buffer (multiple consumers). */
  atomic<uint32_t> head {};
  /** Tail of the buffer (single producer). */
  uint32_t tail = 0;
  /** Data buffer. */
  T data[CAPACITY];
  #pragma endregion


  #pragma region METHODS
  public:
  /**
   * Check if the buffer is empty.
   * @returns true if empty, else false
   */
  inline bool empty() const noexcept {
    return tail == head.load(memory_order_acquire);
  }

  /**
   * Check if the buffer is full.
   * @returns true if full, else false
   */
  inline bool full() const noexcept {
    return (tail + 1) % CAPACITY == head.load(memory_order_acquire);
  }

  /**
   * Get the number of elements in the buffer.
   * @returns number of elements
   */
  inline uint32_t size() const noexcept {
    return (CAPACITY + tail - head.load(memory_order_acquire)) % CAPACITY;
  }

  /**
   * Clear the buffer.
   */
  inline void clear() noexcept {
    head.store(0);
    tail = 0;
  }

  /**
   * Try to push an element to the buffer.
   * @param value element to push
   * @returns 0 if successful, 1 if full
   */
  inline int tryPush(const T& value) noexcept {
    if (full()) return 1;
    data[tail] = value;
    tail = (tail + 1) % CAPACITY;
    return 0;
  }

  /**
   * Try to pop an element from the buffer.
   * @param value element to pop (output)
   * @returns 0 if successful, 1 if empty, -1 if failed
   */
  inline int tryPop(T& value) noexcept {
    uint32_t h = head.load(memory_order_acquire);
    if (h == tail) return -100;
    if (!head.compare_exchange_weak(h, h+1, memory_order_acq_rel)) return -1;
    value = data[h % CAPACITY];
    return 0;
  }

  /**
   * Lock the buffer from consumers.
   * @returns head of the buffer
   */
  inline uint32_t lockHead() noexcept {
    while (1) {
      uint32_t h = head.load(memory_order_acquire);
      if (head.compare_exchange_weak(h, tail, memory_order_acq_rel)) return h;
    }
  }
  #pragma endregion
};




/**
 * A thread-safe Arena Allocator, that supports fast allocation and deallocation.
 * @tparam SIZE size of each allocation
 * @tparam MIN_FREED minimum number of freed allocations to keep [128]
 * @note Use for allocations larger than cache line size to avoid false sharing.
 */
template <size_t SIZE, size_t CAPACITY, uint32_t LOCAL_FREED=128, uint32_t GLOBAL_FREED=4096>
class ConcurrentArenaAllocator {
  #pragma region CONSTANTS
  public:
  /** Size of each allocation. */
  static constexpr size_t allocation_size = SIZE;
  /** Size of the memory pool. */
  static constexpr size_t pool_size = CAPACITY;
  /** Maximum number of thread-local freed allocations. */
  static constexpr size_t max_local_freed = LOCAL_FREED;
  #pragma endregion


  #pragma region DATA
  /** The memory pool. */
  void *pool;
  /** Number of bytes used in the memory pool. */
  atomic<size_t> used = 0;
  /** Mutex to protect the global freed list. */
  atomic_flag gmutex = {};
  /** Global list of freed allocations. */
  vector<void*> gfreed {};
  /** Per-thread list of freed allocations. */
  MultiConsumerRingBuffer<void*, LOCAL_FREED> *lfreed;
  #pragma endregion


  #pragma region METHODS
  private:
  /**
   * Move n freed elements to the global freed list.
   * @param t thread index
   * @param n number of elements to move
   */
  inline void moveFreedGlobal(int t, uint32_t n) noexcept {
    // Lock consumers from using the thread-local freed list.
    uint32_t b = lfreed[t].lockHead();
    uint32_t e = b + n;
    // Push elements to global freed list.
    while (!gmutex.test_and_set(memory_order_acquire));
    for (uint32_t i=b; i<e; ++i)
      gfreed.push_back(lfreed[t].data[i % LOCAL_FREED]);
    gmutex.clear(memory_order_release);
    // Unlock consumers to use the thread-local freed list.
    lfreed[t].head.store(e % LOCAL_FREED, memory_order_release);
  }

  public:
  /**
   * Reset the allocator, freeing all memory.
   */
  inline void reset() noexcept {
    used.store(0);
    gfreed.clear();
    for (auto& f : lfreed)
      f.clear();
  }


  /**
   * Free a memory block.
   * @param ptr memory to free
   */
  inline void free(void *ptr) noexcept {
    int t = omp_get_thread_num();
    // Try adding to thread-local freed list.
    if (lfreed[t].tryPush(ptr) == 0) return;
    // Else, move 50% of elements to global freed list.
    moveFreedGlobal(t, (lfreed[t].size() + 1) / 2);
    // And then, add to thread-local freed list.
    lfreed[t].tryPush(ptr);
  }


  /**
   * Allocate a memory block.
   * @returns allocated memory, or nullptr if out of memory
   */
  inline void* allocate() noexcept {
    int t = omp_get_thread_num();
    // Allocate from thread-local freed list, if available.
    while (1) {
      void *ptr = nullptr;
      // Try to allocate from thread-local freed list.
      int err = lfreed[t].tryPop(ptr);
      if (err==0) return ptr;
      // Else, try to allocate from global freed list.
      if (!gfreed.empty() && !gmutex.test_and_set(memory_order_acquire)) {
        void *ptr = gfreed.back();
        gfreed.pop_back();
        gmutex.clear(memory_order_release);
        return ptr;
      }
      // Else, try to allocate from pool.
      if (used < CAPACITY) return (char*) pool + used.fetch_add(SIZE, memory_order_relaxed);
      // Abort if no memory available.
      if (err==1 && gfreed.empty() && used>=CAPACITY) return nullptr;
    }
  }
  #pragma endregion


  #pragma region CONSTRUCTORS
  public:
  /**
   * Create a new Concurrent Arena Allocator.
   * @param pool memory pool
   * @param capacity size of memory pool
   */
  ConcurrentArenaAllocator(void *pool)
  : pool(pool) {
    int T = omp_get_max_threads();
    // Reserve space for global and thread-local freed lists.
    gfreed.reserve(min(CAPACITY/SIZE, size_t(GLOBAL_FREED)));
    lfreed = new MultiConsumerRingBuffer<void*, LOCAL_FREED>[T];
  }

  /**
   * Destroy the Concurrent Arena Allocator.
   */
  ~ConcurrentArenaAllocator() {
    // delete[] lfreed;
  }
  #pragma endregion
};




// template <size_t ALLOC_SIZE, size_t POOL_SIZE, size_t MIN_FREED=128>
// class MultiArenaAllocator {
//   #pragma region CONSTANTS
//   public:
//   /** Size of each allocation. */
//   static constexpr size_t allocation_size = SIZE;
//   /** Size of the memory pool. */
//   static constexpr size_t capacity = CAPACITY;
//   #pragma endregion


//   #pragma region DATA
//   /** These allocations have been freed, and can be reused. */
//   vector<void*> freed;
//   /** The memory pool. */
//   void  *pool;
//   /** Number of bytes used in the memory pool. */
//   size_t used = 0;
//   #pragma endregion


//   #pragma region METHODS
//   public:
//   /**
//    * Allocate a memory block.
//    * @returns allocated memory, or nullptr if out of memory
//    */
//   inline void* allocate() noexcept {
//     // Allocate from freed list, if available.
//     if (!freed.empty()) {
//       void *ptr = freed.back();
//       freed.pop_back();
//       return ptr;
//     }
//     // Else, try to allocate from pool.
//     if (used < capacity) return (char*) pool + used;
//     return nullptr;
//   }


//   /**
//    * Free a memory block.
//    * @param ptr memory to free
//    */
//   inline void free(void *ptr) {
//     freed.push_back(ptr);
//   }


//   /**
//    * Reset the allocator, freeing all memory.
//    */
//   inline void reset() noexcept {
//     // Reset freed list.
//     freed.clear();
//     // Reset used bytes.
//     used = 0;
//   }
//   #pragma endregion


//   #pragma region CONSTRUCTORS
//   public:
//   /**
//    * Create a new Arena Allocator.
//    */
//   ArenaAllocator() {
//     // Allocate memory pool.
//     pool = mmapAlloc(CAPACITY);
//   }
//   #pragma endregion
// };






// /**
//  * A Recursive Arena Allocator, that supports fast allocation and deallocation.
//  * @tparam SIZE size of each allocation
//  * @tparam MIN_FREED minimum number of freed allocations to keep [32]
//  * @note This allocator is not thread-safe. Use for allocations smaller than cache line size, on each thread.
//  */
// template <size_t SIZE, size_t MIN_FREED=32>
// class RecursiveArenaAllocator {
//   #pragma region CONSTANTS
//   public:
//   /** Size of each allocation. */
//   static constexpr size_t allocation_size = SIZE;
//   #pragma endregion


//   #pragma region DATA
//   /** These allocations have been freed, and can be reused. */
//   vector<void*> freed;
//   /** The memory pools. */
//   vector<void*> pools;
//   /** Number of bytes used in the last memory pool. */
//   size_t used = 0;
//   /** Size of each memory pool. */
//   size_t capacity;
//   #pragma endregion


//   #pragma region METHODS
//   public:
//   /**
//    * Allocate a memory block.
//    * @param fallocate function to allocate a memory pool, e.g., fallocate(capacity)
//    * @returns allocated memory, or nullptr if out of memory
//    */
//   template <class FA>
//   inline void* allocate(FA fallocate) {
//     // Allocate from freed list, if available.
//     if (!freed.empty()) {
//       void *ptr = freed.back();
//       freed.pop_back();
//       return ptr;
//     }
//     // Else, try to allocate from last pool.
//     if (!pools.empty() && used < capacity) {
//       void *ptr = (char*) pools.back() + used;
//       used += SIZE;
//       return ptr;
//     }
//     // Else, must allocate a new pool.
//     void *pool = fallocate(capacity);
//     pools.push_back(pool);
//     used = SIZE;
//     return pool;
//   }


//   /**
//    * Free a memory block.
//    * @param ptr memory to free
//    */
//   inline void free(void *ptr) {
//     freed.push_back(ptr);
//   }


//   /**
//    * Reset the allocator, freeing all memory.
//    * @param ffree function to free a memory pool, e.g., ffree(pool)
//    */
//   template <class FF>
//   inline void reset(FF ffree) noexcept {
//     // Reset freed list, used bytes.
//     freed.clear();
//     used = 0;
//     // Free all memory pools.
//     for (void *pool : pools)
//       ffree(pool);
//     pools.clear();
//   }
//   #pragma endregion


//   #pragma region CONSTRUCTORS
//   public:
//   /**
//    * Create a new Recursive Arena Allocator.
//    * @param capacity size of each memory pool
//    */
//   RecursiveArenaAllocator(size_t capacity)
//   : capacity(capacity) {
//     // Reserve space for freed list.
//     freed.reserve(min(capacity/SIZE, MIN_FREED));
//   }
//   #pragma endregion
// };
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
