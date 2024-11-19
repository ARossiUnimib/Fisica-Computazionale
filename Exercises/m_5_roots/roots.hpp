#pragma once

#include <cassert>
#include <cstdlib>
#include <functional>

#include "../m_2_matrices/tensor.hpp"
#include "../m_4_eigenvalues/eigenvalues.hpp"

namespace func {

template <typename T>
T Bisection(std::function<T(T)> f, T a, T b, T tol = 1e-6, int n_tol = 1000) {
  assert(f(a) * f(b) < 0);

  T c = 0.0;
  int n = 0;
  while (true) {
    if (n > n_tol) {
      // NOTE: Probably we shoud log error
      LOG_WARN("Bisection did not converge");
      break;
    }

    if (std::abs(b - a) < tol) {
      break;
    }

    c = (a + b) / 2.0;

    if (f(c) == 0.0) {
      break;
    } else if (f(c) * f(a) < 0) {
      b = c;
    } else {
      a = c;
    }

    n++;
  }
  return c;
}

template <typename T, typename F = std::function<T(T)> const&>
T NewtonRaphson(F f, F f_prime, T x0, T tol = 1e-16, int n_tol = 1000) {
  double x1;
  int n = 0;

  while (n < n_tol) {
    x1 = x0 - f(x0) / f_prime(x0);

    if (std::abs(x1 - x0) < tol) {
      LOG_WARN("Newton-Raphson did not converge");
      return x1;
    }

    x0 = x1;
    n++;
  }

  LOG_WARN("Newton-Raphson did not converge");
  return x1;
}

template <typename T>
T Secant(std::function<T(T)> const& f, T x0, T x1, T tol = 1e-6,
         int n_tol = 1000) {
  double x2;
  int n = 0;

  while (n < n_tol) {
    x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));

    if (std::abs(x2 - x1) < tol) {
      LOG_WARN("Secant did not converge");
      return x2;
    }

    x0 = x1;
    x1 = x2;
    n++;
  }

  LOG_WARN("Secant did not converge");
  return x2;
}

template <typename T>
std::vector<T> PolynomialRoots(tensor::Tensor<T> coefficients, int n,
                               int n_eigenvalues) {
  LOG_ASSERT(coefficients.Rows() > 1,
             "Polynomial must have at least 2 coefficients", utils::ERROR);

  // Number of normalized coeffcients
  int n_coeffs = coefficients.Rows() - 1;

  auto coeffs_matrix = tensor::Tensor<T>::SMatrix(n_coeffs);

  // Fill the n-1 x n-1 matrix
  for (int i = 0; i < n_coeffs; i++) {
    if (i < n_coeffs - 1) {
      coeffs_matrix(i, i + 1) = 1.0;
    }

    coeffs_matrix(n_coeffs - 1, i) = -coefficients(i) / coefficients(n_coeffs);
  }

  auto solution = eigen::PowerMethodDeflation(coeffs_matrix, n);

  std::vector<T> roots;
  for (auto&& [eigenvalue, eigenvector] : solution) {
    roots.push_back(eigenvalue);
  }

  return roots;
}

}  // namespace func
