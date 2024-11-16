#pragma once

#include <cassert>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <ostream>

#include "../m_2_matrices/tensor.hpp"
#include "../m_4_eigenvalues/eigenvalues.hpp"

namespace func {

template <typename T>
T Bisection(std::function<T(T)> f, T a, T b, T tol = 1e-6) {
  assert(f(a) * f(b) < 0);

  T c = 0.0;
  while (b - a > tol) {
    c = (a + b) / 2.0;

    if (f(c) == 0.0) {
      break;
    } else if (f(c) * f(a) < 0) {
      b = c;
    } else {
      a = c;
    }
  }
  return c;
}

template <typename T>
T NewtonRaphson(std::function<T(T)> f, std::function<T(T)> f_prime, T x0,
                T tol = 1e-16, int n_tol = 1000) {
  T x = x0;
  T x_new = x0;
  int n = 0;
  while (true) {
    x = x_new;
    x_new = x - f(x) / f_prime(x);

    if (n > n_tol) {
      std::cerr << "Newton-Raphson did not converge for value: " << x0
                << std::endl;
      break;
    }
    if (std::abs(x_new - x) < tol) {
      break;
    }
    n++;
  }
  return x_new;
}

template <typename T>
T Secant(std::function<T(T)> f, T x0, T x1, T tol = 1e-6) {
  T x = x0;
  T x_new = x1;
  T x_old = x0;

  T f_old = f(x0);
  T f_new = f(x1);

  while (true) {
    x_old = x;
    x = x_new;
    f_new = f(x);
    x_new = x - f_new * (x - x_old) / (f_new - f_old);
    f_old = f_new;

    if (std::abs(x_new - x) < tol) {
      break;
    }
  }
  return x_new;
}

template <typename T>
std::vector<T> PolynomialRoots(tensor::Tensor<T> coefficients, int n,
                               int n_eigenvalues) {
  assert(coefficients.Rows() > 1);

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

  // Inverse is too unstable
  // auto solution = eigen::PowerMethodDeflation(coeffs_matrix, n,
  // n_eigenvalues);
  auto solution = eigen::PowerMethodDeflation(coeffs_matrix, n);

  // convert in vector
  std::vector<T> roots;
  for (auto &&[eigenvalue, eigenvector] : solution) {
    roots.push_back(eigenvalue);
  }

  return roots;
}

}  // namespace func
