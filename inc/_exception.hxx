#pragma once
#include <string>
#include <string_view>
#include <exception>

using std::string;
using std::string_view;
using std::exception;




#pragma region TYPES
/**
 * Exception thrown when a string is not in the expected format.
 */
class FormatException : public exception {
  #pragma region DATA
  /** Pointer to the character where format check fails. */
  const char *it;
  #pragma endregion

  public:
  #pragma region CONSTRUCTORS
  /**
   * Create a format exception.
   * @param msg error message
   * @param it iterator to the character where format check fails
   */
  template <class I>
  FormatException(const char *msg, I it) :
  exception(msg), it(&*it) {}
  #pragma endregion

  #pragma region METHODS
  /**
   * Get pointer to the character where format check fails.
   * @returns pointer to the character
   */
  inline const char* where() const noexcept {
    return it;
  }
  #pragma endregion
};
#pragma endregion
