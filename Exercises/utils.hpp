#pragma once

#include <cmath>
#include <csignal>
#include <iostream>
#include <ostream>
#include <random>
#include <string>
#include <vector>

namespace utils {

/**
 * @brief Slice a vector
 *
 * @tparam T
 * @param v
 * @param slice_range
 * @return std::vector<T>
 */
template <typename T>
std::vector<T> slice(std::vector<T> const &v, int start, int end) {
  auto first = v.cbegin() + start;
  auto last = v.cbegin() + end + 1;
  return std::vector<T>(first, last);
}

enum LogLevel { INFO, WARN, ERROR };

#ifndef NDEBUG

// TODO: implement stack trace
inline void UtilsLog(const std::string &message, int type, int line,
                     const char *file) {
  // Set prefix
  std::string prefix = "";

  switch (type) {
    case INFO:
      // green
      prefix += "\033[1;32m";
      prefix += "INFO";
      break;
    case WARN:
      // yellow
      prefix += "\033[1;33m";
      prefix += "WARN";
      break;
    case ERROR:
      // red
      prefix += "\033[1;31m";
      prefix += "ERROR";
      break;
  }

  // reset
  prefix += "\033[0m";

  std::cout << "[" << prefix << "] " << message;

    std::cout << " at " << file << ":" << line << std::endl;
}

#define LOG_ASSERT(expr, message, type)          \
  if (!(expr)) {                                 \
    UtilsLog(message, type, __LINE__, __FILE__); \
    if (type == utils::ERROR) {                  \
      raise(SIGABRT);                            \
    }                                            \
  }

#define LOG_INFO(message) UtilsLog(message, utils::INFO, __LINE__, __FILE__)
#define LOG_WARN(message) UtilsLog(message, utils::WARN, __LINE__, __FILE__)
#define LOG_ERROR(message) UtilsLog(message, utils::ERROR, __LINE__, __FILE__)

#else

inline void UtilsLog(const std::string &message, int type) {}

#define LOG_ASSERT(expr, message, type)
#define LOG_INFO(message)
#define LOG_WARN(message)
#define LOG_ERROR(message)

#endif

template <typename T>
T RandomValue(T a, T b) {
  static std::random_device rd;
  static std::mt19937 gen(rd());

  if constexpr (std::is_integral<T>::value) {
    std::uniform_int_distribution<T> dis(a, b);
    return dis(gen);
  } else if constexpr (std::is_floating_point<T>::value) {
    std::uniform_real_distribution<T> dis(a, b);
    return dis(gen);
  } else {
    static_assert(std::is_arithmetic<T>::value,
                  "RandomValue requires an arithmetic type");
  }
}

}  // namespace utils
