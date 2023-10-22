#pragma once




#pragma region METHODS
#pragma region CHARACTER CLASSIFICATION
/**
 * Check if a character is a blank.
 * @param c character to check
 * @returns true if c is a blank [ \t]
 */
inline bool isBlank(char c) {
  return c==' ' || c=='\t';
}


/**
 * Check if a character is a newline.
 * @param c character to check
 * @returns true if c is a newline [\r\n]
 */
inline bool isNewline(char c) {
  return c=='\r' || c=='\n';
}


/**
 * Check if a character is a whitespace.
 * @param c character to check
 * @returns true if c is a whitespace [ \t\r\n\v\f]
 */
inline bool isWhitespace(char c) {
  return isBlank(c) || isNewline(c) || c=='\v' || c=='\f';
}


/**
 * Check if a character is a digit.
 * @param c character to check
 * @returns true if c is a digit [0-9]
 */
inline bool isDigit(char c) {
  return c>='0' && c<='9';
}


/**
 * Check if a character is a hex digit.
 * @param c character to check
 * @returns true if c is a hex digit [0-9A-Fa-f]
 */
inline bool isHexDigit(char c) {
  return (c>='0' && c<='9') || (c>='A' && c<='F') || (c>='a' && c<='f');
}


/**
 * Check if a character is an uppercase alphabet.
 * @param c character to check
 * @returns true if c is an uppercase alphabet [A-Z]
 */
inline bool isUppercaseAlphabet(char c) {
  return c>='A' && c<='Z';
}


/**
 * Check if a character is a lowercase alphabet.
 * @param c character to check
 * @returns true if c is a lowercase alphabet [a-z]
 */
inline bool isLowercaseAlphabet(char c) {
  return c>='a' && c<='z';
}


/**
 * Check if a character is an alphabet.
 * @param c character to check
 * @returns true if c is an alphabet [A-Za-z]
 */
inline bool isAlphabet(char c) {
  return isUppercaseAlphabet(c) || isLowercaseAlphabet(c);
}


/**
 * Check if a character is an alphabet or digit.
 * @param c character to check
 * @returns true if c is an alphabet or digit [A-Za-z0-9]
 */
inline bool isAlphabetOrDigit(char c) {
  return isAlphabet(c) || isDigit(c);
}


/**
 * Check if a character is a control character.
 * @param c character to check
 * @returns true if c is a control character [\0-\x1F\x7F]
 */
inline bool isControlCharacter(char c) {
  return c<0x20 || c==0x7F;
}


/**
 * Check if a character is a printable character.
 * @param c character to check
 * @returns true if c is a printable character [\x20-\x7E]
 */
inline bool isPrintableCharacter(char c) {
  return c>=0x20 && c<0x7F;
}


/**
 * Check if a character is a graphical character.
 * @param c character to check
 * @returns true if c is a graphical character [\x21-\x7E]
 */
inline bool isGraphicalCharacter(char c) {
  return c>0x20 && c<0x7F;
}


/**
 * Check if a character is a punctuation character.
 * @param c character to check
 * @returns true if c is a punctuation character
 */
inline bool isPunctuationCharacter(char c) {
  return isGraphicalCharacter(c) && !isAlphabetOrDigit(c);
}
#pragma endregion




#pragma region CHARACTER CONVERSION
/**
 * Convert a character to uppercase.
 * @param c character to convert
 * @returns uppercase character
 */
inline char toUppercaseAlphabet(char c) {
  return isLowercaseAlphabet(c)? c-'a'+'A' : c;
}


/**
 * Convert a character to lowercase.
 * @param c character to convert
 * @returns lowercase character
 */
inline char toLowercaseAlphabet(char c) {
  return isUppercaseAlphabet(c)? c-'A'+'a' : c;
}
#pragma endregion
#pragma endregion
