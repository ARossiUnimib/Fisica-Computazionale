#pragma once

#include <cmath>
#include <cstdlib>
#include <utility>
#include <vector>

#include "../m_3_interpolation/function_data.hpp"
#include "../m_3_interpolation/range.hpp"

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
std::pair<T, T> GaussLegendreQuadrature(T a, T b, T (*f)(T), int n) {
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
// H_n(x) = (-1)^n e^(x^2) (d^n/dx^n e^(-x^2))
template <typename T>
func::FunctionData<T> HermitePolynomial(int n, func::Range<T> x) {
  func::FunctionData<T> p;
  for (auto x_ : x) {
    p.Add(x_, std::pow(-1, n) * std::exp(x_ * x_) *
                  Derivative(std::exp, -x_ * x_, n));
  }
  return p;
}

}  // namespace integration
