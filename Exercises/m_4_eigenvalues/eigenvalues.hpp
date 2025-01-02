#pragma once

#include <stdexcept>
#include <utility>

#include "../m_2_matrices/tensor.hpp"
#include "../utils.hpp"

namespace eigen {

template <typename T>
T RayleighQuotient(tensor::Tensor<T> const& A, tensor::Tensor<T>& x) {
  return x.Dagger().Dot(A).Dot(x)(0);
}

template <typename T>
std::pair<T, tensor::Tensor<T>> PowerMethod(tensor::Tensor<T> const& A,
                                            int n_max, T alpha = 0,
                                            T tol = 1e-6) {
  tensor::Tensor<T> x = tensor::Tensor<T>::RandomVector(A.Rows());
  x = x / x.Norm();

  tensor::Tensor<T> I = tensor::Tensor<T>::Identity(A.Rows());
  tensor::Tensor<T> shifted_A = A + alpha * I;

  for (int i = 0; i < n_max; i++) {
    tensor::Tensor<T> y = shifted_A.Dot(x);

    // Check the convergence tolerance
    auto delta = std::abs(RayleighQuotient(shifted_A, x) -
                          RayleighQuotient(shifted_A, y));

    if (std::abs(delta) < tol) {
      break;
    }

    x = y / y.Norm();
  }

  // Eigenvalue remains the same: (A-aI)v=Av-av=λv-av=(λ-a)v
  return {RayleighQuotient(shifted_A, x) - alpha, x};
}

template <typename T>
std::pair<T, tensor::Tensor<T>> InversePowerMethod(tensor::Tensor<T> const& A,
                                                   int n, T tol = 1e-6) {
  auto x = tensor::Tensor<T>::RandomVector(A.Rows());

  // NOTE: shift operation creates more problems than it solves: disabled at
  // the moment
  // auto I = tensor::Tensor<T>::Identity(A.Rows());
  // T mu = RayleighQuotient(A, x);

  for (int i = 0; i < n; i++) {
    tensor::Tensor<T> y = tensor::Tensor<T>::Vector(A.Rows());

    // Catching singularity
    try {
      y = tensor::GaussianElimination(A, x);
      // y = tensor::GaussianElimination(A - I * mu, x);
    } catch (const std::runtime_error& e) {
      LOG_WARN("Matrix become singular");
      break;
    }

    auto delta = RayleighQuotient(A, x) - RayleighQuotient(A, y);
    if (std::abs(delta) < tol) {
      LOG_INFO("Tolerance was exceeded");
      break;
    }

    x = y / y.Norm();

    // mu = RayleighQuotient(A, y);
  }

  return {RayleighQuotient(A, x), x};
}

template <typename T>
std::vector<std::pair<T, tensor::Tensor<T>>> PowerMethodDeflation(
    tensor::Tensor<T> const& A, int n, T alpha = 0) {
  tensor::Tensor<T> A_deflated = A;
  std::vector<std::pair<T, tensor::Tensor<T>>> eigenpairs;

  for (int i = 0; i < A.Rows(); i++) {
    T eigenvalue;
    tensor::Tensor<T> eigenvector;
    std::tie(eigenvalue, eigenvector) = PowerMethod(A_deflated, n, alpha);

    eigenpairs.push_back({eigenvalue, eigenvector});

    A_deflated =
        A_deflated - eigenvalue * eigenvector.Dot(eigenvector.Dagger());
  }

  return eigenpairs;
}

// NOTE: this is only for testing purpose, obviously the inverse power method
// deflation creates singular matrices
template <typename T>
std::vector<std::pair<T, tensor::Tensor<T>>> InversePowerMethodDeflation(
    tensor::Tensor<T> const& A, int n) {
  LOG_WARN("Inverse power method deflation creates singular matrices");

  tensor::Tensor<T> A_deflated = A;
  std::vector<std::pair<T, tensor::Tensor<T>>> eigenpairs;

  for (int i = 0; i < A.Rows(); i++) {
    T eigenvalue;
    tensor::Tensor<T> eigenvector;
    std::tie(eigenvalue, eigenvector) = InversePowerMethod(A_deflated, n);

    eigenpairs.push_back({eigenvalue, eigenvector});

    A_deflated =
        A_deflated - eigenvalue * eigenvector.Dot(eigenvector.Dagger());
  }

  return eigenpairs;
}

}  // namespace eigen
