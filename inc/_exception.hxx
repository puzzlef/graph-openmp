#pragma once
#include <exception>

using std::exception;




#pragma region TYPES
/**
 * Error thrown when a string is not in the expected format.
 */
class FormatError : public exception {
  #pragma region DATA
  /** Error message. */
  const char *msg;
  /** Pointer to the character where format check fails. */
  const void *it;
  #pragma endregion


  #pragma region CONSTRUCTORS
  public:
  /**
   * Create a format error.
   * @param msg error message
   * @param it iterator to the character where format check fails
   */
  template <class I>
  FormatError(const char *msg, I it) :
  msg(msg), it(&*it) {}

  /**
   * Create a format error, without iterator.
   * @param msg error message
   */
  FormatError(const char *msg) :
  msg(msg), it(nullptr) {}

  /**
   * Create an empty format error.
   */
  FormatError() :
  msg(nullptr), it(nullptr) {}
  #pragma endregion


  #pragma region METHODS
  public:
  inline bool empty() const noexcept {
    return msg == nullptr;
  }

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
