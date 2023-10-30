#pragma once
#include <utility>
#include <type_traits>
#include <iterator>
#include <tuple>
#include <array>
#include <vector>
#include <ostream>
#include <iostream>
#include <chrono>
#include <ctime>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>

using std::pair;
using std::tuple;
using std::array;
using std::vector;
using std::ostream;
using std::is_fundamental;
using std::iterator_traits;
using std::time_t;
using std::tm;
using std::localtime;
using std::chrono::time_point;
using std::chrono::system_clock;
using std::tie;
using std::cout;




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




#pragma region WRITE
/**
 * Write values to a stream.
 * @param a the stream
 * @param ib begin iterator of values
 * @param ie end iterator of values
 */
template <class I>
inline void write_values(ostream& a, I ib, I ie) {
  using T = typename iterator_traits<I>::value_type;
  if (is_fundamental<T>::value) {
    a << "{";
    for (; ib < ie; ++ib)
      a << " " << *ib;
    a << " }";
  }
  else {
    a << "{\n";
    for (; ib < ie; ++ib)
      a << "  " << *ib << "\n";
    a << "}";
  }
}


/**
 * Write values to a stream.
 * @param a the stream
 * @param x the container of values
 */
template <class J>
inline void writeValues(ostream& a, const J& x) {
  write_values(a, x.begin(), x.end());
}


/**
 * Write a pair to a stream.
 * @param a the stream
 * @param x the pair
 */
template <class K, class V>
inline void write(ostream& a, const pair<K, V>& x) {
  a << x.first << ": " << x.second;
}


/**
 * Write an array to a stream.
 * @param a the stream
 * @param x the array
 */
template <class T, size_t N>
inline void write(ostream& a, const array<T, N>& x) {
  writeValues(a, x);
}


/**
 * Write a vector to a stream.
 * @param a the stream
 * @param x the vector
 */
template <class T>
inline void write(ostream& a, const vector<T>& x) {
  writeValues(a, x);
}


/**
 * Write a pair to a stream.
 * @param a the stream
 * @param x the pair
 */
template <class K, class V>
inline ostream& operator<<(ostream& a, const pair<K, V>& x) {
  write(a, x); return a;
}


/**
 * Write an array to a stream.
 * @param a the stream
 * @param x the array
 */
template <class T, size_t N>
inline ostream& operator<<(ostream& a, const array<T, N>& x) {
  write(a, x); return a;
}


/**
 * Write a vector to a stream.
 * @param a the stream
 * @param x the vector
 * @returns the stream
 */
template <class T>
inline ostream& operator<<(ostream& a, const vector<T>& x) {
  write(a, x); return a;
}
#pragma endregion




#pragma region WRITE TIME
/**
 * Write a time to a stream.
 * @param a the stream
 * @param x the time
 */
inline void writeTime(ostream& a, const time_t& x) {
  const int BUF = 64;
  char  buf[BUF];
  tm* t = localtime(&x);
  sprintf(buf, "%04d-%02d-%02d %02d:%02d:%02d",
    t->tm_year + 1900,
    t->tm_mon  + 1,
    t->tm_mday,
    t->tm_hour,
    t->tm_min,
    t->tm_sec
  );
  a << buf;
}


/**
 * Write a time point to a stream.
 * @param a the stream
 * @param x the time point
 */
inline void writeTimePoint(ostream& a, const time_point<system_clock>& x) {
  writeTime(a, system_clock::to_time_t(x));
}


/**
 * Write a time to a stream.
 * @param a the stream
 * @param x the time
 * @returns the stream
 */
inline ostream& operator<<(ostream& a, const time_t& x) {
  writeTime(a, x); return a;
}


/**
 * Write a time point to a stream.
 * @param a the stream
 * @param x the time point
 * @returns the stream
 */
inline ostream& operator<<(ostream& a, const time_point<system_clock>& x) {
  writeTimePoint(a, x); return a;
}
#pragma endregion




#pragma region PRINT
/**
 * Print an object.
 * @param x the object to print
 */
template <class T>
inline void print(const T& x) {
  cout << x;
}


/**
 * Print an object followed by a newline.
 * @param x the object to print
 */
template <class T>
inline void println(const T& x) {
  cout << x << "\n";
}


/**
 * Print a newline.
 */
inline void println() {
  cout << "\n";
}
#pragma endregion
