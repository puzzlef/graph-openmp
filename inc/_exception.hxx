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
  /** Error message. */
  const char *msg;
  /** Pointer to the character where format check fails. */
  const void *it;
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
  msg(msg), it(&*it) {}
  #pragma endregion

  #pragma region METHODS
  /**
   * Get error message.
   * @returns error message
   */
  inline const char* what() const noexcept override {
    return msg;
  }

  /**
   * Get pointer to the character where format check fails.
   * @returns pointer to the character
   */
  inline const void* where() const noexcept {
    return it;
  }
  #pragma endregion
};
#pragma endregion
