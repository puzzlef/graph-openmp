#pragma once
#include <type_traits>
#include <utility>
#include <string_view>
#include <cstdint>
#include <x86intrin.h>
#include "_debug.hxx"
#include "_cctype.hxx"

using std::is_integral;
using std::is_floating_point;
using std::pair;
using std::string_view;




#pragma region METHODS
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
 * Parse a whole number from a string, using SIMD instructions.
 * @tparam FULL 0=partial, 1=full range, 2=full range with check
 * @param a obtained number (output)
 * @param ib begin iterator
 * @param ie end iterator
 * @returns iterator to end of number, or error if FULL==2
 */
template <int FULL=0, class T, class I>
inline I parseWholeNumberSimdW(T &a, I ib, I ie) {
  // Initialize constants.
  const __m256i C0 = _mm256_set1_epi8('0');
  const __m256i D9 = _mm256_set1_epi8(9);
  const __m256i P1 = _mm256_set_epi8(
    1, 10, 1, 10, 1, 10, 1, 10, 1, 10, 1, 10, 1, 10, 1, 10,
    1, 10, 1, 10, 1, 10, 1, 10, 1, 10, 1, 10, 1, 10, 1, 10
  );
  const __m128i P2 = _mm_set_epi8(1, 100, 1, 100, 1, 100, 1, 100, 1, 100, 1, 100, 1, 100, 1, 100);
  const __m128i P4 = _mm_set_epi16(1, 10000, 1, 10000, 1, 10000, 1, 10000);
  // Find the length of integer in text.
  if constexpr (FULL==2) { I it = findNextNonDigit(ib, ie); if (it!=ie) return it; }
  if constexpr (FULL==0) ie = findNextNonDigit(ib, ie);
  size_t n = size_t(ie - ib);
  // Load 32 bytes with proper mask (zero-fill the rest).
  uint32_t mask = uint32_t(0xFFFFFFFF) << (ib - ie + 32);
  auto xc = _mm256_maskz_loadu_epi8(mask, ie - 32);
  // Subtract character '0' per byte as per mask.
  auto xd = _mm256_maskz_sub_epi8(mask, xc, C0);
  // Now do a 1, 10, 1, 10, ... per byte multiply and add to 16bit.
  auto x2_16 = _mm256_maddubs_epi16(xd, P1);
  // Convert the 16bit values to 8bit values (as 99 < 256).
  auto x2_08 = _mm256_cvtepi16_epi8(x2_16);
  // Now do a 1, 100, 1, 100, ... per byte multiple and add to 16 bit.
  auto x4_16 = _mm_maddubs_epi16(x2_08, P2);
  // Now do a 1, 10000, 1, 10000, ... per word (16bit) and add to 32-bit.
  auto x8_32 = _mm_madd_epi16(x4_16, P4);
  // Finally extract 3 32-bit values and calculate the total value.
  uint64_t u = _mm_extract_epi32(x8_32, 3);
  uint64_t v = _mm_extract_epi32(x8_32, 2);
  uint64_t w = _mm_extract_epi32(x8_32, 1);
  if constexpr (sizeof(T)<=4) a = u + v*100000000;
  else a = u + v*100000000 + w*10000000000000000;
  return ie;
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
 * Parse an integer from a string, using SIMD instructions.
 * @tparam FULL 0=partial, 1=full range, 2=full range with check
 * @param a obtained number (output)
 * @param ib begin iterator
 * @param ie end iterator
 * @returns iterator to end of number, or error if FULL==2
 */
template <int FULL=0, class T, class I>
inline I parseIntegerSimdW(T &a, I ib, I ie) {
  // Skip if empty.
  if (ib==ie) return ib;
  // Parse whole number, and apply sign.
  if (*ib=='-') { ib = parseWholeNumberSimdW<FULL>(a, ib+1, ie); a = -a; }
  else          { ib = parseWholeNumberSimdW<FULL>(a, *ib=='+'? ib+1 : ib, ie); }
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
 * Parse a floating-point number from a string, using SIMD instructions.
 * @tparam FULL 0=partial, 1=full range, 2=full range with check
 * @param a obtained number (output)
 * @param ib begin iterator
 * @param ie end iterator
 * @returns iterator to end of number, or error if FULL==2
 */
template <int FULL=0, class T, class I>
inline I parseFloatSimdW(T &a, I ib, I ie) {
  // Initialize constants.
  static constexpr double DV[] = {  // Exponents of 10 for fractional part
    1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10,
    1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17, 1e-18, 1e-19, 1e-20
  };
  // Skip if empty.
  if (ib==ie) return ib;
  // Initialize variables.
  uint64_t u = 0, v = 0;  // Whole and fractional parts.
  int      d = 0, e = 0;  // Fractional and exponent digits.
  // Handle sign.
  bool neg = *ib=='-';
  if (*ib=='-' || *ib=='+') ++ib;
  // Parse whole part, fractional part, and exponent.
  ib = parseWholeNumberSimdW(u, ib, ie);
  if (ib!=ie && *ib=='.') { I id = ++ib; ib = parseWholeNumberSimdW(v, ib, ie); d = int(ib-id); }
  if (ib!=ie && (*ib=='e' || *ib=='E'))  ib = parseIntegerSimdW<FULL>(e, ib+1, ie);
  // Compute number, and apply sign.
  a = (T(u) + (v? T(v) * T(DV[d]) : T())) * (e? T(pow(10, e)) : T(1));
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


/**
 * Parse a number from a string, using SIMD instructions.
 * @tparam FULL 0=partial, 1=full range, 2=full range with check
 * @param a obtained number (output)
 * @param ib begin iterator
 * @param ie end iterator
 * @returns iterator to end of number, or error if FULL==2
 */
template <int FULL=0, class T, class I>
inline I parseNumberSimdW(T& a, I ib, I ie) {
  if constexpr      (is_integral<T>::value)       return parseIntegerSimdW<FULL>(a, ib, ie);
  else if constexpr (is_floating_point<T>::value) return parseFloatSimdW<FULL>  (a, ib, ie);
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
 * @returns iterator to error
 */
template <bool CHECK=false, class I, class FU, class FW>
inline I readTokenU(string_view& a, I& ib, I ie, FU fu, FW fw) {
  auto [tb, te] = findNextToken(ib, ie, fu, fw);
  if constexpr (CHECK) { if (tb==te) return tb; }
  a  = string_view(&*tb, te-tb);
  ib = te;  // Update begin iterator
  return I();
}


/**
 * Obtain the next number from a string.
 * @tparam CHECK check for error?
 * @param a obtained number (output)
 * @param ib begin iterator (updated)
 * @param ie end iterator
 * @param fu is special blank, e.g. comma? (c)
 * @param fw is special whitespace, e.g. comma? (c)
 * @returns iterator to error
 */
template <bool CHECK=false, class T, class I, class FU, class FW>
inline I readNumberU(T& a, I& ib, I ie, FU fu, FW fw) {
  auto [tb, te] = findNextToken(ib, ie, fu, fw);
  if constexpr (CHECK) { if (tb==te) return tb; }
  tb = parseNumberSimdW<CHECK? 2:1>(a, tb, te);
  if constexpr (CHECK) { if (tb!=te) return tb; }
  ib = te;  // Update begin iterator
  return I();
}
#pragma endregion
#pragma endregion
