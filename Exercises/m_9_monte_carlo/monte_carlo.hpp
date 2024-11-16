#pragma once

#include <functional>
#include <type_traits>
#include <utility>

#include "../m_3_interpolation/range.hpp"
#include "../utils.hpp"

namespace montecarlo {

template <typename T>
std::pair<T, T> UniformSampling(std::function<T(std::vector<T>)> f,
                                std::vector<std::pair<T, T>> ranges, int N) {
  int dim = ranges.size();
  T sum = 0;
  T sum2 = 0;
  for (int i = 0; i < N; i++) {
    std::vector<T> x;
    for (int j = 0; j < dim; j++) {
      x.push_back(utils::RandomValue(ranges[j].first, ranges[j].second));
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

template <typename T>
std::pair<T, T> ImportanteSampling(std::function<T(T)> f, std::function<T(T)> g,
                                   std::function<T(T)> generator,
                                   std::pair<T, T> range, int N) {
  T sum = 0;
  T sum2 = 0;
  for (int i = 0; i < N; i++) {
    T x = generator(utils::RandomValue(0, 1));
    T y = f(x) / g(x);
    sum += y;
    sum2 += y * y;
  }

  T mean = (sum) / N;
  T variance = (sum2 / N) - mean * mean;
  T error = std::sqrt(variance / N);
  return {mean, error};
}

template <typename T>
std::vector<T> AccumulateValues(int n_iter, func::Range<T> range_gen,
                                func::Range<T> range_values,
                                T (*generator)(T x)) {
  LOG_ASSERT(std::is_floating_point<T>(),
             "Frequencies must be floating numbers", utils::ERROR);

  // Vector stores the number of times the values was generated in the interval
  // (range.End() - range.Start()) / range.nodes().size() + i*range.Step()
  auto vec = std::vector<T>(range_values.Nodes().size());

  for (int i = 0; i < n_iter; i++) {
    double value =
        generator(utils::RandomValue(range_gen.Start(), range_gen.End()));

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
    vec[i] /= n_iter;
  }

  return vec;
}

}  // namespace montecarlo
