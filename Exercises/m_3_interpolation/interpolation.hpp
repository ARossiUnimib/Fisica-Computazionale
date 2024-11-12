#pragma once

#include "../m_2_matrices/tensor.hpp"
#include "../m_2_matrices/tensor_utils.hpp"
#include "function_data.hpp"
#include "function_utils.hpp"
#include "range.hpp"

namespace interp {

/**
 * @brief Compute the Vandermonde matrix
 *
 * @tparam T
 * @param values
 * @return Tensor<T>
 */
template <typename T>
tensor::Tensor<T> VandermondeMatrix(std::vector<T> const &values);

/**
 * @brief Retrieve polynomial coefficients from direct interpolation using
 * Vandermonde matrix
 *
 * @tparam T
 * @param f
 * @return std::vector<T> coefficients of size f.Size()
 */
template <typename T>
std::vector<T> DirectCoefficients(func::FunctionData<T> const &f);

/**
 * @brief Direct interpolation using Vandermonde matrix
 *
 * @tparam T
 * @param f
 * @param range Values to be extrapolated from polynomial
 * @return FunctionData<T> polynomial data
 */
template <typename T>
func::FunctionData<T> DirectPolynomial(func::FunctionData<T> const &f,
                                       func::Range<T> const &range);

/**
 * @brief Compute the Newton coefficients of the function
 * @note use macro _SLOW_ALGORITHM to use the algorithm given during the course
 *
 * @tparam T
 * @param f
 * @return std::vector<T>
 */
template <typename T>
std::vector<T> NewtonCoefficients(func::FunctionData<T> const &f);

/**
 * @brief Compute the Newton polynomial of the function
 *
 * @tparam T
 * @param f
 * @param range
 * @return FunctionData<T>
 */
template <typename T>
func::FunctionData<T> NewtonPolynomial(func::FunctionData<T> const &f,
                                       func::Range<T> const &range);

/* ---------------------------------------------------- */

template <typename T>
tensor::Tensor<T> VandermondeMatrix(std::vector<T> const &values) {
  auto mat = tensor::Tensor<T>::SMatrix(values.size());

  for (int i = 0; i < values.size(); i++) {
    for (int j = 0; j < values.size(); j++) {
      mat(i, j) = pow(values[i], j);
    }
  }

  return mat;
}
template <typename T>
std::vector<T> DirectCoefficients(func::FunctionData<T> const &f) {
  auto f_tensor = tensor::Tensor<T>::FromData(f.F());

  tensor::Tensor<T> vande_matrix = VandermondeMatrix(f.X());
  tensor::Tensor<T> values =
      tensor::GaussianElimination(vande_matrix, f_tensor);

  return values.RawData();
}

template <typename T>
func::FunctionData<T> DirectPolynomial(func::FunctionData<T> const &f,
                                       func::Range<T> const &range) {
  std::vector<T> a = DirectCoefficients(f);
  return func::Polynomial(a, range);
}

template <typename T>
std::vector<T> NewtonCoefficients(func::FunctionData<T> const &f) {
  int N = f.Size();

  std::vector<T> a(N);

  auto A = tensor::Tensor<double>::SMatrix(N);

  for (int i = 0; i < N; i++) {
    A(i, 0) = f.F(i);
  }

  for (int j = 1; j < N; j++) {
    for (int i = 0; i < N - j; i++) {
      A(i, j) = (A(i + 1, j - 1) - A(i, j - 1)) / (f.X(i + j) - f.X(i));
    }
  }

  for (int i = 0; i < N; i++) {
    a[i] = A(0, i);
  }

  return a;
}

template <typename T>
func::FunctionData<T> NewtonPolynomial(func::FunctionData<T> const &f,
                                       func::Range<T> const &range) {
  std::vector<T> a = NewtonCoefficients(f);
  int N = f.Size();

  func::FunctionData<T> p;
  T sum;

  for (auto i : range) {
    sum = a[N - 1];

    for (int j = N - 2; j >= 0; j--) {
      sum = a[j] + (i - f.X(j)) * sum;
    }

    p.Add(i, sum);
  }

  return p;
}

}  // namespace interp
