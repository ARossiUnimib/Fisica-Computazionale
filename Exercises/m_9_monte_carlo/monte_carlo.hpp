#pragma once

#include <functional>
#include <type_traits>
#include <utility>

#include "../m_3_interpolation/range.hpp"
#include "../utils.hpp"

namespace montecarlo {

template <typename T>
T UniformRandomValue(T a, T b) {
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

template <typename T,
          typename F = std::function<T(std::vector<T> const&)> const&,
          typename V = std::vector<std::pair<T, T>> const&>

std::pair<T, T> UniformSampling(F f, V ranges, int N) {
  int dim = ranges.size();
  T sum = 0;
  T sum2 = 0;
  for (int i = 0; i < N; i++) {
    std::vector<T> x;
    for (int j = 0; j < dim; j++) {
      x.push_back(UniformRandomValue(ranges[j].first, ranges[j].second));
    }
    T y = f(x);
    sum += y;
    sum2 += y * y;
  }

  T mean = (sum) / N;
  T variance = (sum2 / N) - mean * mean;
  T error = std::sqrt(variance / N);
  return {mean, error};
}

template <typename T, typename F = std::function<T(T const)> const&>
std::pair<T, T> ImportanteSampling(F f, F g, F generator, std::pair<T, T> range,
                                   int N) {
  T sum = 0;
  T sum2 = 0;
  for (int i = 0; i < N; i++) {
    T x = generator(UniformRandomValue(0, 1));
    T y = f(x) / g(x);
    sum += y;
    sum2 += y * y;
  }

  T mean = (sum) / N;
  T variance = (sum2 / N) - mean * mean;
  T error = std::sqrt(variance / N);
  return {mean, error};
}

template <typename T, typename R = func::Range<T> const&,
          typename F = std::function<T(T const)> const&>
std::vector<T> AccumulateValues(R range_gen, R range_values, F generator,
                                int N) {
  LOG_ASSERT(std::is_floating_point<T>(),
             "Frequencies must be floating numbers", utils::ERROR);

  // Vector stores the number of times the values was generated in the interval
  // (range.End() - range.Start()) / range.nodes().size() + i*range.Step()
  auto vec = std::vector<T>(range_values.Nodes().size());

  for (int i = 0; i < N; i++) {
    double value =
        generator(UniformRandomValue(range_gen.Start(), range_gen.End()));

    // if its in the nth interval increment
    if (value >= range_values.Start() && value <= range_values.End()) {
      // We get the index from removing the decimal part of
      // (value - start) / step
      int index =
          std::floor((value - range_values.Start()) / range_values.Step());
      vec[index]++;
    }
  }

  // Normalize frequencies
  for (int i = 0; i < vec.size(); i++) {
    vec[i] /= N;
  }

  return vec;
}

}  // namespace montecarlo
