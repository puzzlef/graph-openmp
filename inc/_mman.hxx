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
 * Map file to memory.
 * @param pth file path
 * @returns file descriptor, mapped data, and file size
 */
inline tuple<int, void*, size_t> mapFileToMemory(const char *pth) {
  // Open file as read-only.
  int fd = open(pth, O_RDONLY);
  if (fd==-1) return {-1, nullptr, 0};
  // Get file size.
  struct stat sb;
  if (fstat(fd, &sb)==-1) return {-1, nullptr, 0};
  // Map file to memory.
  void *addr = mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE | MAP_NORESERVE, fd, 0);  // MAP_SHARED?
  if (addr==MAP_FAILED) return {-1, nullptr, 0};
  madvise(addr, sb.st_size, MADV_WILLNEED);  // MADV_SEQUENTIAL?
  // Return file descriptor, mapped data, and file size.
  return {fd, addr, sb.st_size};
}


/**
 * Unmap file from memory.
 * @param fd file descriptor
 * @param addr mapped data
 * @param size file size
 */
inline void unmapFileFromMemory(int fd, void *addr, size_t size) {
  munmap(addr, size);
  close(fd);
}


/**
 * Represents a memory mapped file.
 */
struct MappedFile {
  #pragma region DATA
  /** File descriptor. */
  int    fd;
  /** Mapped data. */
  void  *addr;
  /** File size. */
  size_t size;
  #pragma endregion


  #pragma region CONSTRUCTORS / DESTRUCTORS
  /**
   * Map file to memory.
   * @param pth file path
   */
  MappedFile(const char *pth) {
    tie(fd, addr, size) = mapFileToMemory(pth);
  }

  /**
   * Unmap file from memory.
   */
  ~MappedFile() {
    if (fd<=0) return;
    unmapFileFromMemory(fd, addr, size);
    fd = -1;
  }
  #pragma endregion
};
#pragma endregion




#pragma region ALLOCATE MEMORY
/**
 * Allocate memory using mmap.
 * @param size memory size
 * @returns allocated memory
 */
inline void* allocateMemoryMmap(size_t size) {
  return mmap(NULL, size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
}


/**
 * Free memory allocated using mmap.
 * @param addr memory address
 * @param size memory size
 */
inline void freeMemoryMmap(void *addr, size_t size) {
  munmap(addr, size);
}
#pragma endregion
