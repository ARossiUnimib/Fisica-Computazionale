#pragma once

#include <utility>

#include "../m_2_matrices/tensor.hpp"

namespace eigen {

template <typename T>
T RayleighQuotient(tensor::Tensor<T> const& A, tensor::Tensor<T>& x) {
  auto y = x.Dagger().Dot(A).Dot(x);
  auto z = x.Dagger().Dot(x);
  return y(0) / z(0);
}

// get tolerance based on the precision of the type
template <typename T>
T Tolerance() {
  return std::numeric_limits<T>::epsilon();
}

template <typename T>
std::pair<T, tensor::Tensor<T>> PowerMethod(tensor::Tensor<T> const& A,
                                            int n_max) {
  tensor::Tensor<T> x = tensor::Tensor<T>::RandomVector(A.Rows());

  int i = 0;
  while (i < n_max) {
    tensor::Tensor<T> y = A.Dot(x);

    auto delta = std::abs(RayleighQuotient(A, x) - RayleighQuotient(A, y));
    if (delta < Tolerance<T>()) {
      break;
    }

    x = y / y.Norm();
    i++;
  }

  return {RayleighQuotient(A, x), x};
}

template <typename T>
std::pair<T, tensor::Tensor<T>> InversePowerMethod(tensor::Tensor<T> const& A,
                                                   int n) {
  // auto I = tensor::Tensor<T>::Identity(A.Rows());

  // NOTE: Initial guess can change the estimated eigenvalue
  auto x = tensor::Tensor<T>::RandomVector(A.Rows());

  // NOTE: shift operation creates more problems than it solves: disabled at
  // the moment

  // T mu = RayleighQuotient(A, x);

  for (int i = 0; i < n; i++) {
    // Rule of thumb
    LOG_ASSERT(n < A.Rows() * 3,
               "Inverse power method converges fast, avoid large numbers "
               "otherwise a "
               "singular matrix will be produced in gaussian elimination",
               // NOTE: maybe implement a method to avoid the singularity
               utils::WARN);

    auto y = tensor::GaussianElimination(A, x);
    // auto y = tensor::GaussianElimination(A - I * mu, x);

    x = y / y.Norm();

    // mu = RayleighQuotient(A, y);
  }

  return {RayleighQuotient(A, x), x};
}

template <typename T>
std::vector<std::pair<T, tensor::Tensor<T>>> PowerMethodDeflation(
    tensor::Tensor<T> const& A, int n, int num_eigenvalues) {
  tensor::Tensor<T> A_deflated = A;
  std::vector<std::pair<T, tensor::Tensor<T>>> eigenpairs;

  for (int i = 0; i < num_eigenvalues; i++) {
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
// NOTE: the method is not stable it seems, for close eigenvalues at least
std::vector<std::pair<T, tensor::Tensor<T>>> InversePowerMethodDeflation(
    tensor::Tensor<T> const& A, int n, int num_eigenvalues) {
  tensor::Tensor<T> A_deflated = A;
  std::vector<std::pair<T, tensor::Tensor<T>>> eigenpairs;

  for (int i = 0; i < num_eigenvalues; i++) {
    T eigenvalue;
    tensor::Tensor<T> eigenvector;
    std::tie(eigenvalue, eigenvector) = InversePowerMethod(A_deflated, n);

    eigenpairs.push_back({eigenvalue, eigenvector});

    A_deflated =
        A_deflated - eigenvalue * eigenvector.Dot(eigenvector.Dagger());
  }

  return eigenpairs;  // Return all computed eigenvalue/eigenvector pairs
}

}  // namespace eigen
