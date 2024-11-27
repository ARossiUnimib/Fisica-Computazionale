#pragma once

#include <cmath>
#include <cstdlib>
#include <utility>
#include <vector>

#include "../m_2_matrices/tensor.hpp"
#include "../m_3_interpolation/function_data.hpp"
#include "../m_3_interpolation/range.hpp"
#include "../m_5_roots/roots.hpp"
#include "../utils.hpp"

namespace integration {

template <typename T>
T MaxSecondDerivative(T a, T b, T (*f)(T), T dx = 0.0001) {
  T max = 0;
  for (T i = a; i < b; i += dx) {
    T temp = std::abs(f(i + 2 * dx) - 2 * f(i + dx) + f(i));
    if (temp > max) {
      max = temp;
    }
  }
  return max;
}

template <typename T>
T MaxFourthDerivative(T a, T b, T (*f)(T), T dx = 0.0001) {
  T max = 0;
  for (T i = a; i < b; i += dx) {
    T temp = std::abs(f(i + 4 * dx) - 4 * f(i + 3 * dx) + 6 * f(i + 2 * dx) -
                      4 * f(i + dx) + f(i));
    if (temp > max) {
      max = temp;
    }
  }
  return max;
}

template <typename T>
std::pair<T, T> Trapezoid(T a, T b, T (*f)(T), int n) {
  T dx = (b - a) / (double)n;

  T sum = 0.5 * (f(a) + f(b));

  for (int i = 1; i < n - 2; i++) {
    sum += f(a + i * dx);
  }

  T M = MaxSecondDerivative(a, b, f, dx);
  T error = 1.0 / 12.0 * M * (b - a) * dx * dx;

  return {dx * sum, error};
}

template <typename T>
std::pair<T, T> Simpson(T a, T b, T (*f)(T), int n) {
  T dx = (b - a) / n;
  T sum = f(a) + f(b);

  for (int i = 1; i < n; i++) {
    if (i % 2 == 0) {
      sum += 2 * f(a + i * dx);
    } else {
      sum += 4 * f(a + i * dx);
    }
  }

  T M = MaxFourthDerivative(a, b, f, dx);
  T error = 1.0 / 180.0 * M * (b - a) * dx * dx * dx * dx;

  return {dx / 3 * sum, error};
}

template <typename T>
std::pair<T, T> GaussLegendreQuadrature(T a, T b, std::function<T(T)> f,
                                        int n) {
  T sum = 0;
  T dx = (b - a) / n;
  for (int i = 0; i < n; i++) {
    T x = a + i * dx;
    // c-
    T x1 = x + dx / 2 - dx / 2 / std::sqrt(3);
    // c+
    T x2 = x + dx / 2 + dx / 2 / std::sqrt(3);
    sum += f(x1) + f(x2);
  }
  // TODO: calculate error
  return {dx / 2 * sum, 0};
}

// Calculate nth derivative of a function
template <typename T>
T Derivative(double (*f)(T), T x, int n, double dx = 0.0001) {
  if (n == 0) {
    return f(x);
  } else {
    return (Derivative(f, x + dx, n - 1) - Derivative(f, x, n - 1)) / dx;
  }
}

// Calculate Hermite polynomial of orden n
// Using explicit sum representation
// https://en.wikipedia.org/wiki/Hermite_polynomials
// NOTE: Hermite polynomials scales as n! ... very inefficient to compute and also it is easy to exceed maximum floating type size 
template <typename T>
tensor::Tensor<T> HermitePolynomial(int n) {
  std::vector<T> coefficients(n + 1);

  for (int m = 0; m <= floor(n / 2.0); m++) {
    using namespace utils;
    int sign = m % 2 == 0 ? 1 : -1;
    coefficients[n - 2 * m] =
        static_cast<T>(Factorial(n) * sign * pow(2, n - 2 * m) /
                       (Factorial(m) * Factorial(n - 2 * m)));
  }

  return tensor::Tensor<T>::FromData(coefficients);
}

template <typename T>
T LagrangePolynomial(T x, int i, std::vector<T> j_data) {
  T product = 1;
  for (int j = 0; j < j_data.size(); j++) {
    if (i != j) {
      product *= (x - j_data[j]) / (j_data[i] - j_data[j]);
    }
  }
  return product;
}

template <typename T>
std::vector<std::pair<T, T>> GenerateGaussHermiteCoeffs(T a, T b, int n_poly,
                                                        int n_int,
                                                        int n_eigen) {
  tensor::Tensor<double> coeffs =
      integration::HermitePolynomial<double>(n_poly);
  std::vector<T> zeros = func::PolynomialRoots(coeffs, n_eigen, 0.1);

  std::vector<std::pair<T, T>> weights(zeros.size());

  for (int i = 0; i < zeros.size(); i++) {
    std::function<T(T)> f = [&zeros, &i](T x) {
      return exp(-x * x) * LagrangePolynomial(x, i, zeros);
    };

    weights[i].first = zeros[i];
    weights[i].second = GaussLegendreQuadrature(a, b, f, n_int).first;
  }

  return weights;
}

//  FIXME: It clearly needs a class to store the precomputed Coefficients
template <typename T>
T GaussHermiteQuadrature(T a, T b, std::function<T(T)> f, int n_poly = 10,
                         int n_int = 1000, int n_eigen = 1000) {
  auto weights = GenerateGaussHermiteCoeffs(a, b, n_poly, n_int, n_eigen);
  T sum = 0;
  for (int i = 0; i < weights.size(); i++) {
    sum += weights[i].second * f(weights[i].first);
  }
  return sum;
}

}  // namespace integration
