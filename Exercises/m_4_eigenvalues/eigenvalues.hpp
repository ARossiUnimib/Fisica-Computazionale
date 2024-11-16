#pragma once

#include <complex>
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
                                            int n_max) {
  tensor::Tensor<T> x = tensor::Tensor<T>::RandomVector(A.Rows());

  for (int i = 0; i < n_max; i++) {
    tensor::Tensor<T> y = A.Dot(x);

    auto delta = std::abs(RayleighQuotient(A, x) - RayleighQuotient(A, y));
    if (std::abs(delta) < 1e-6) {
      break;
    }

    x = y / y.Norm();
  }

  return {RayleighQuotient(A, x), x};
}

template <typename T>
std::pair<T, tensor::Tensor<T>> InversePowerMethod(tensor::Tensor<T> const& A,
                                                   int n) {
  // NOTE: Initial guess can change the estimated eigenvalue
  // FIXME: imaginary part is 0!
  auto x = tensor::Tensor<T>::RandomVector(A.Rows());

  // NOTE: shift operation creates more problems than it solves: disabled at
  // the moment
  // auto I = tensor::Tensor<T>::Identity(A.Rows());
  // T mu = RayleighQuotient(A, x);

  for (int i = 0; i < n; i++) {
    tensor::Tensor<T> y = tensor::Tensor<T>::Vector(A.Rows());

    try {
      y = tensor::GaussianElimination(A, x);
      // y = tensor::GaussianElimination(A - I * mu, x);
    } catch (const std::runtime_error& e) {
      LOG_INFO("Matrix become singular");
      break;
    }

    auto delta = RayleighQuotient(A, x) - RayleighQuotient(A, y);
    if (std::abs(delta) < 1e-6) {
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
    tensor::Tensor<T> const& A, int n) {
  tensor::Tensor<T> A_deflated = A;
  std::vector<std::pair<T, tensor::Tensor<T>>> eigenpairs;

  for (int i = 0; i < A.Rows(); i++) {
    T eigenvalue;
    tensor::Tensor<T> eigenvector;
    std::tie(eigenvalue, eigenvector) = PowerMethod(A_deflated, n);

    eigenpairs.push_back({eigenvalue, eigenvector});

    A_deflated =
        A_deflated - eigenvalue * eigenvector.Dot(eigenvector.Dagger());
  }

  return eigenpairs;
}

template <typename T>
// NOTE: this is only for testing purpose, obviously the inverse power method
// deflation creates singular matrices
std::vector<std::pair<T, tensor::Tensor<T>>> InversePowerMethodDeflation(
    tensor::Tensor<T> const& A, int n) {
  
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
