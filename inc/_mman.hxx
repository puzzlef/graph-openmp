#pragma once
#include <tuple>
#include <cstdint>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>

using std::tuple;
using std::tie;




#pragma region MEMORY MAPPED FILE
/**
 * Map a file to memory.
 * @param pth file path
 * @returns file descriptor, mapped data, and file size
 */
inline tuple<int, void*, size_t> mmapOpenFile(const char *pth) {
  // Open file as read-only.
  int fd = open(pth, O_RDONLY);
  if (fd==-1) return {-1, nullptr, 0};
  // Get file size.
  struct stat sb;
  if (fstat(fd, &sb)==-1) return {-1, nullptr, 0};
  // Map file to memory.
  void *data = mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE | MAP_NORESERVE, fd, 0);  // MAP_SHARED?
  if (data==MAP_FAILED) return {-1, nullptr, 0};
  madvise(data, sb.st_size, MADV_WILLNEED);  // MADV_SEQUENTIAL?
  // Return file descriptor, mapped data, and file size.
  return {fd, data, sb.st_size};
}


/**
 * Unmap a file from memory.
 * @param fd file descriptor
 * @param data mapped data
 * @param size file size
 */
inline void mmapCloseFile(int fd, void *data, size_t size) {
  munmap(data, size);
  close(fd);
}


/**
 * Represents a memory mapped file.
 */
struct MappedFile {
  #pragma region DATA
  private:
  /** Mapped data. */
  void  *_data;
  /** File size. */
  size_t _size;
  /** File descriptor. */
  int    _fd;
  #pragma endregion


  #pragma region METHODS
  public:
  /**
   * Get mapped data (implicit conversion).
   */
  inline operator void*() const { return _data; }

  /**
   * Get mapped data.
   * @returns mapped data
   */
  inline void* data() const { return _data; }

  /**
   * Get file descriptor.
   * @returns file descriptor
   */
  inline int fd() const { return _fd; }

  /**
   * Get file size.
   * @returns file size
   */
  inline size_t size() const { return _size; }

  /**
   * Unmap file from memory.
   */
  inline void close() {
    if (_fd<=0) return;
    mmapCloseFile(_fd, _data, _size);
    _fd   = -1;
    _data = nullptr;
    _size = 0;
  }
  #pragma endregion


  #pragma region CONSTRUCTORS / DESTRUCTORS
  public:
  /**
   * Map file to memory.
   * @param pth file path
   */
  MappedFile(const char *pth) {
    tie(_fd, _data, _size) = mmapOpenFile(pth);
  }

  /**
   * Unmap file from memory.
   */
  ~MappedFile() { close(); }
  #pragma endregion
};
#pragma endregion




#pragma region ALLOCATE MEMORY
/**
 * Allocate memory using mmap.
 * @param size memory size
 * @returns allocated memory
 */
inline void* mmapAlloc(size_t size) {
  return mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
}


/**
 * Free memory allocated using mmap.
 * @param addr memory address
 * @param size memory size
 */
inline void mmapFree(void *addr, size_t size) {
  munmap(addr, size);
}


/**
 * Represents a memory allocated with mmap().
 */
template <class T>
struct MappedPtr {
  #pragma region DATA
  private:
  /** Base address of allocated memory. */
  T     *_data;
  /** Size of allocated memory. */
  size_t _size;
  #pragma endregion


  #pragma region METHODS
  public:
  /**
   * Get base address of allocated memory (implicit conversion).
   * @returns base address
   */
  inline operator T*() const { return _data; }

  /**
   * Get base address of allocated memory.
   * @returns base address
   */
  inline T* data() const { return _data; }

  /**
   * Get size of allocated memory.
   * @returns allocation size
   */
  inline size_t size() const { return _size; }

  /**
   * Unmap file from memory.
   */
  inline void release() {
    if (_data == nullptr) return;
    mmapFree(_data, _size);
    _data = nullptr;
    _size = 0;
  }
  #pragma endregion


  #pragma region CONSTRUCTORS / DESTRUCTORS
  public:
  /**
   * Create an empty allocation.
   */
  MappedPtr() : _data(nullptr), _size(0) {}

  /**
   * Allocate memory using mmap().
   * @param size size of memory to allocate
   */
  MappedPtr(size_t size) {
    _data = (T*) mmapAlloc(size);
    _size = size;
  }

  /**
   * Free memory allocated using mmap().
   */
  ~MappedPtr() { release(); }
  #pragma endregion
};
#pragma endregion
