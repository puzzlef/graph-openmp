#pragma once
#include <cstdint>
#include <cmath>
#include <type_traits>
#include <utility>
#include <string>
#include <string_view>
#include "_debug.hxx"
#include "_cctype.hxx"
#include "_exception.hxx"

using std::is_integral;
using std::is_floating_point;
using std::pair;
using std::string_view;
using std::string;
using std::pow;




#pragma region METHODS
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




#pragma region FIND NEXT
/**
 * Find the next line in a string.
 * @tparam EOL end of line character
 * @param ib begin iterator
 * @param ie end iterator
 * @returns iterator to next line
 */
template <char EOL='\n', class I>
inline I findNextLine(I ib, I ie) {
  for (; ib<ie; ++ib)
    if (*ib==EOL) return ++ib;
  return ie;
}


/**
 * Find the next whitespace in a string.
 * @param ib begin iterator
 * @param ie end iterator
 * @param fw is special whitespace, e.g. comma? (c)
 * @returns iterator to next whitespace
 */
template <class I, class FW>
inline I findNextWhitespace(I ib, I ie, FW fw) {
  for (; ib<ie; ++ib)
    if (isBlank(*ib) || isNewline(*ib) || fw(*ib)) return ib;
  return ie;
}


/**
 * Find the next non-whitespace in a string.
 * @param ib begin iterator
 * @param ie end iterator
 * @param fw is special whitespace, e.g. comma? (c)
 * @returns iterator to next non-whitespace
 */
template <class I, class FW>
inline I findNextNonWhitespace(I ib, I ie, FW fw) {
  for (; ib<ie; ++ib)
    if (!isBlank(*ib) && !isNewline(*ib) && !fw(*ib)) return ib;
  return ie;
}


/**
 * Find the next blank in a string.
 * @param ib begin iterator
 * @param ie end iterator
 * @param fu is special blank, e.g. comma? (c)
 * @returns iterator to next blank
 */
template <class I, class FU>
inline I findNextBlank(I ib, I ie, FU fu) {
  for (; ib<ie; ++ib)
    if (isBlank(*ib) || fu(*ib)) return ib;
  return ie;
}


/**
 * Find the next non-blank in a string.
 * @param ib begin iterator
 * @param ie end iterator
 * @param fu is special blank, e.g. comma? (c)
 * @returns iterator to next non-blank
 */
template <class I, class FU>
inline I findNextNonBlank(I ib, I ie, FU fu) {
  for (; ib<ie; ++ib)
    if (!isBlank(*ib) && !fu(*ib)) return ib;
  return ie;
}


/**
 * Find the next digit in a string.
 * @param ib begin iterator
 * @param ie end iterator
 * @returns iterator to next digit
 */
template <class I>
inline I findNextDigit(I ib, I ie) {
  for (; ib!=ie && !isDigit(*ib); ++ib);
  return ib;
}


/**
 * Find the next non-digit in a string.
 * @param ib begin iterator
 * @param ie end iterator
 * @returns iterator to next non-digit
 */
template <class I>
inline I findNextNonDigit(I ib, I ie) {
  for (; ib!=ie && isDigit(*ib); ++ib);
  return ib;
}


/**
 * Find the next token in a string.
 * @param ib begin iterator
 * @param ie end iterator
 * @param fu is special blank, e.g. comma? (c)
 * @param fw is special whitespace, e.g. comma? (c)
 * @returns [begin, end) iterators to next token
 */
template <class I, class FU, class FW>
inline pair<I, I> findNextToken(I ib, I ie, FU fu, FW fw) {
  auto tb = findNextNonBlank(ib, ie, fu);
  auto te = findNextWhitespace(tb+1, ie, fw);
  return {tb, te};
}
#pragma endregion




#pragma region PARSE NUMBER
/**
 * Parse a whole number from a string.
 * @tparam FULL 0=partial, 1=full range, 2=full range with check
 * @param a obtained number (output)
 * @param ib begin iterator
 * @param ie end iterator
 * @returns iterator to end of number, or error if FULL==2
 */
template <int FULL=0, class T, class I>
inline I parseWholeNumberW(T &a, I ib, I ie) {
  a = T();
  for (; ib!=ie && (FULL==1 || isDigit(*ib)); ++ib)
    a = a*10 + (*ib - '0');
  return ib;
}


/**
 * Parse an integer from a string.
 * @tparam FULL 0=partial, 1=full range, 2=full range with check
 * @param a obtained number (output)
 * @param ib begin iterator
 * @param ie end iterator
 * @returns iterator to end of number, or error if FULL==2
 */
template <int FULL=0, class T, class I>
inline I parseIntegerW(T &a, I ib, I ie) {
  // Skip if empty.
  if (ib==ie) return ib;
  // Handle sign.
  bool neg = *ib=='-';
  if (*ib=='-' || *ib=='+') ++ib;
  // Scan whole number.
  ib = parseWholeNumberW<FULL>(a, ib, ie);
  // Apply sign.
  if (neg) a = -a;
  return ib;
}


/**
 * Parse a floating-point number from a string.
 * @tparam FULL 0=partial, 1=full range, 2=full range with check
 * @param a obtained number (output)
 * @param ib begin iterator
 * @param ie end iterator
 * @returns iterator to end of number, or error if FULL==2
 */
template <int FULL=0, class T, class I>
inline I parseFloatW(T &a, I ib, I ie) {
  // Skip if empty.
  if (ib==ie) return ib;
  // Initialize variables.
  uint64_t u = 0, v = 0;  // Whole and fractional parts.
  int      d = 0, e = 0;  // Fractional and exponent digits.
  // Handle sign.
  bool neg = *ib=='-';
  if (*ib=='-' || *ib=='+') ++ib;
  // Parse whole part, fractional part, and exponent.
  ib = parseWholeNumberW(u, ib, ie);
  if (ib!=ie && *ib=='.') { I id = ++ib; ib = parseWholeNumberW(v, ib, ie); d = int(ib-id); }
  if (ib!=ie && (*ib=='e' || *ib=='E'))  ib = parseIntegerW<FULL>(e, ib+1, ie);
  // Compute number, and apply sign.
  a = (T(u) + T(v) * pow(10, -d)) * T(pow(10, e));
  if (neg) a = -a;
  return ib;
}


/**
 * Parse a number from a string.
 * @tparam FULL 0=partial, 1=full range, 2=full range with check
 * @param a obtained number (output)
 * @param ib begin iterator
 * @param ie end iterator
 * @returns iterator to end of number, or error if FULL==2
 */
template <int FULL=0, class T, class I>
inline I parseNumberW(T& a, I ib, I ie) {
  if constexpr      (is_integral<T>::value)       return parseIntegerW<FULL>(a, ib, ie);
  else if constexpr (is_floating_point<T>::value) return parseFloatW<FULL>  (a, ib, ie);
  return ib;
}
#pragma endregion




#pragma region READ TOKEN
/**
 * Obtain the next token from a string.
 * @tparam CHECK check for error?
 * @param a obtained token (output)
 * @param ib begin iterator (updated)
 * @param ie end iterator
 * @param fu is special blank, e.g. comma? (c)
 * @param fw is special whitespace, e.g. comma? (c)
 * @returns iterator to end of token
 */
template <bool CHECK=false, class I, class FU, class FW>
inline I readTokenW(string_view& a, I ib, I ie, FU fu, FW fw) {
  auto [tb, te] = findNextToken(ib, ie, fu, fw);
  if constexpr (CHECK) { if (tb==te) throw FormatError("Failed to read token (empty)", tb); }
  a  = string_view(&*tb, te-tb);
  return te;
}


/**
 * Obtain the next number from a string.
 * @tparam CHECK check for error?
 * @param a obtained number (output)
 * @param ib begin iterator (updated)
 * @param ie end iterator
 * @param fu is special blank, e.g. comma? (c)
 * @param fw is special whitespace, e.g. comma? (c)
 * @returns iterator to end of number
 */
template <bool CHECK=false, class T, class I, class FU, class FW>
inline I readNumberW(T& a, I ib, I ie, FU fu, FW fw) {
  auto [tb, te] = findNextToken(ib, ie, fu, fw);
  if constexpr (CHECK) { if (tb==te) throw FormatError("Failed to read number (empty)", tb); }
  tb = parseNumberW<CHECK? 2:1>(a, tb, te);
  if constexpr (CHECK) { if (tb!=te) throw FormatError("Failed to read number (bad format)", tb); }
  return te;
}
#pragma endregion
#pragma endregion
