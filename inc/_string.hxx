#pragma once
#include <type_traits>
#include <string>
#include <string_view>
#include <charconv>
#include <cstdint>
#include "_debug.hxx"
#include "_cctype.hxx"

using std::is_integral;
using std::is_floating_point;
using std::string;
using std::string_view;
using std::from_chars;




#pragma region METHODS
#pragma region LINE OPERATIONS
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
 * Remove a line from string view.
 * @tparam EOL end of line character
 * @param x string view (updated)
 * @returns index to next line
 */
template <char EOL='\n'>
inline void dropLineU(string_view& x) {
  auto xb = x.begin(), xe = x.end();
  auto le = findNextLine<EOL>(xb, xe);
  x.remove_prefix(le-xb);
}
#pragma endregion




#pragma region BLANK/WHITESPACE OPERATIONS
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
    if (isWhitespace(*ib) || fw(*ib)) return ib;
  return ie;
}


/**
 * Remove non-whitespace characters from string view.
 * @param x string view (updated)
 * @param fw is special whitespace, e.g. comma? (c)
 */
template <class FW>
inline void dropNonWhitespacesU(string_view& x, FW fw) {
  auto xb = x.begin(), xe = x.end();
  auto wb = findNextWhitespace(xb, xe, fw);
  x.remove_prefix(wb-xb);
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
 * Remove blank characters from string view.
 * @param x string view (updated)
 * @param fu is special blank, e.g. comma? (c)
 */
template <class FU>
inline void dropBlanksU(string_view& x, FU fu) {
  auto xb = x.begin(), xe = x.end();
  auto we = findNextNonBlank(xb, xe, fu);
  x.remove_prefix(we-xb);
}
#pragma endregion




#pragma region TOKEN OPERATIONS
/**
 * Obtain a token from string view.
 * @param a obtained token (output)
 * @param x string view (updated)
 * @param fu is special blank, e.g. comma? (c)
 * @param fw is special whitespace, e.g. comma? (c)
 */
template <class FU, class FW>
inline void readTokenU(string_view& a, string_view& x, FU fu, FW fw) {
  auto xb = x.begin(), xe = x.end();
  auto tb = findNextNonBlank(xb, xe, fu);
  auto te = findNextWhitespace(tb+1, xe, fw);
  a = x.substr(tb-xb, te-tb);
  x.remove_prefix(te-xb);
}


/**
 * Obtain an integer from string view.
 * @param a obtained integer (output)
 * @param x string view (updated)
 * @param fu is special blank, e.g. comma? (c)
 * @param fw is special whitespace, e.g. comma? (c)
 */
template <bool CARELESS=false, class T, class FU, class FW>
inline bool readIntegerU(T& a, string_view& x, FU fu, FW fw) {
  string_view w;
  readTokenU(w, x, fu, fw);
  if (!CARELESS && w.empty()) return true;
  T v = T();
  bool neg = w[0]=='-';
  if (w[0]=='-' || w[0]=='+') w.remove_prefix(1);
  for (auto c : w) {
    if (CARELESS || isDigit(c)) v = v*10 + (c-'0');
    else return true;
  }
  a = neg? -v : v;
  return false;
}


/**
  * Obtain a floating-point value from string view.
  * @param a obtained floating-point value (output)
  * @param x string view (updated)
  * @param fu is special blank, e.g. comma? (c)
  * @param fw is special whitespace, e.g. comma? (c)
  * @returns true if error occurred
  */
template <bool CARELESS=false, class T, class FU, class FW>
inline bool readFloatU(T& a, string_view& x, FU fu, FW fw) {
  string_view w;
  readTokenU(w, x, fu, fw);
  if (!CARELESS && w.empty()) return true;
  int64_t v = 0;
  bool neg = w[0]=='-';
  bool dot = false;
  if (w[0]=='-' || w[0]=='+') w.remove_prefix(1);
  int exp = 0, dec = 0, i = 0;
  for (auto c : w) {
    if (dot) --dec;
    if (isDigit(c)) v = v*10 + (c-'0');
    else if (c=='.') dot = true;
    else if (c=='e' || c=='E') {
      w.remove_prefix(i+1);
      if (readIntegerU<CARELESS>(exp, w, fu, fw)) return true;
      break;
    }
    else return true;
    ++i;
  }
  a = pow(10, exp+dec) * (neg? -v : v);
  return false;
}


/**
 * Obtain a value from string view.
 * @param a obtained value (output)
 * @param x string view (updated)
 * @param fu is special blank, e.g. comma? (c)
 * @param fw is special whitespace, e.g. comma? (c)
 * @returns true if error occurred
 */
template <bool CARELESS=false, class T, class FU, class FW>
inline bool readValueU(T& a, string_view& x, FU fu, FW fw) {
  if constexpr (is_integral<T>::value) return readIntegerU<CARELESS>(a, x, fu, fw);
  else if constexpr (is_floating_point<T>::value) return readFloatU<CARELESS>(a, x, fu, fw);
  else return true;
}
#pragma endregion
#pragma endregion
