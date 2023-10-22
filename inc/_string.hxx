#pragma once
#include <string>
#include <string_view>
#include <charconv>
#include <cstdint>
#include "_debug.hxx"
#include "_cctype.hxx"

using std::string;
using std::string_view;
using std::from_chars;




#pragma region METHODS
#pragma region FIND NEXT LINE
/**
 * Find the next line in a string.
 * @tparam EOL end of line character
 * @param ib begin iterator
 * @param ie end iterator
 * @returns iterator to the next line
 */
template <char EOL='\n'>
inline const char* findNextLine(const char *ib, const char *ie) {
  ASSERT(ib && ie);
  for (; ib<ie; ++ib)
    if (*ib==EOL) return ++ib;
  return ie;
}


/**
 * Find the next line in a string.
 * @tparam EOL end of line character
 * @param x string
 * @param i begin index
 * @returns index to the next line
 */
template <char EOL='\n'>
inline size_t findNextLine(const string_view& x, size_t i) {
  for (; i < x.size(); ++i)
    if (x[i]==EOL) return ++i;
  return x.size();
}
#pragma endregion




#pragma region FIND NEXT WHITESPACE
/**
 * Find the next whitespace in a string.
 * @tparam COMMA comma is considered whitespace?
 * @param ib begin iterator
 * @param ie end iterator
 * @returns iterator to the next whitespace
 */
template <bool COMMA=false>
inline const char* findNextWhitespace(const char *ib, const char *ie) {
  ASSERT(ib && ie);
  for (; ib<ie; ++ib)
    if (isWhitespace(*ib) || (COMMA && *ib==',')) return ib;
  return ie;
}


/**
 * Find the next whitespace in a string.
 * @tparam COMMA comma is considered whitespace?
 * @param x string
 * @param i begin index
 * @returns index to the next whitespace
 */
template <bool COMMA=false>
inline size_t findNextWhitespace(const string_view& x, size_t i) {
  for (; i < x.size(); ++i)
    if (isWhitespace(x[i]) || (COMMA && x[i]==',')) return i;
  return x.size();
}
#pragma endregion




#pragma region FIND NEXT NON-WHITESPACE
/**
 * Find the next non-whitespace in a string.
 * @tparam COMMA comma is considered whitespace?
 * @param ib begin iterator
 * @param ie end iterator
 * @returns iterator to the next non-whitespace
 */
template <bool COMMA=false>
inline const char* findNextNonWhitespace(const char *ib, const char *ie) {
  ASSERT(ib && ie);
  for (; ib<ie; ++ib)
    if (!isWhitespace(*ib) && (!COMMA || *ib!=',')) return ib;
  return ie;
}


/**
 * Find the next non-whitespace in a string.
 * @tparam COMMA comma is considered whitespace?
 * @param x string
 * @param i begin index
 * @returns index to the next non-whitespace
 */
template <bool COMMA=false>
inline size_t findNextNonWhitespace(const string_view& x, size_t i) {
  for (; i < x.size(); ++i)
    if (!isWhitespace(x[i]) && (!COMMA || x[i]!=',')) return i;
  return x.size();
}
#pragma endregion




#pragma region FETCH NEXT
/**
 * Fetch the next word in a string.
 * @tparam COMMA comma is considered whitespace?
 * @param x string
 * @param i begin index (updated to end of word)
 * @returns next word
 */
template <bool COMMA=false>
inline string_view fetchNextWordU(const string_view& x, size_t& i) {
  size_t ib = findNextNonWhitespace<COMMA>(x, i);
  size_t ie = findNextWhitespace<COMMA>(x, ib+1);
  i = ie;  // Update index to end of word
  return x.substr(ib, ie-ib);
}


/**
 * Fetch the next integer/float in a string.
 * @tparam COMMA comma is considered whitespace?
 * @param x string
 * @param i begin index (updated to end of word)
 * @param def default value
 * @returns next integer/float
 */
template <bool COMMA=false, class T>
inline T fetchNextU(const string_view& x, size_t& i, const T& def) {
  string_view w = fetchNextWordU<COMMA>(x, i);
  T v = def;  // Initialize to default value
  from_chars(w.begin(), w.end(), v);
  return v;
}
#pragma endregion




#pragma region COUNT LINES
/**
 * Count the number of lines in a string.
 * @param x string
 * @returns number of lines
 */
inline size_t countLines(const char* x) {
  ASSERT(x);
  size_t a = 1;
  for (; *x; x++) {
    if (*x == '\r' || *x == '\n') ++a;
    else if (*x == '\r' && *(x+1) == '\n') ++x;
  }
  return a;
}


/**
 * Count the number of lines in a string.
 * @param x string
 * @returns number of lines
 */
inline size_t countLines(const string& x) {
  return countLines(x.c_str());
}
#pragma endregion
#pragma endregion
