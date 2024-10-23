#pragma once

#include <cmath>
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

}  // namespace utils
