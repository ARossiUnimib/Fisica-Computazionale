#pragma once

#include <cassert>
#include <iostream>
#include <limits>

#include "../utils.hpp"
#include "tensor.hpp"

namespace tensor {

template <typename T>
using TensorPair = std::pair<Tensor<T>, Tensor<T>>;

/**
 * @brief Calculate the determinant of a matrix using LU decomposition
 *
 * @tparam T
 * @param A
 * @return T
 */
template <typename T>
T DeterminantFromLU(Tensor<T> const &A);

/**
 * @brief Calculate the inverse of a matrix using LU decomposition
 *
 * @tparam T
 * @param A
 * @return Tensor<T>
 */
template <typename T>
TensorPair<T> LUDecomposition(Tensor<T> const &A);

/**
 * @brief Calculate the inverse of a matrix using LU decomposition
 *
 * @tparam T
 * @param A
 * @return Tensor<T>
 */
template <typename T>
Tensor<T> GaussianElimination(Tensor<T> A, Tensor<T> b);

/**
 * @brief Perform Back Substitution Algorithm on an upper triangular matrix
 * Ax = b
 *
 * @tparam T precision
 * @param A upper triangular matrix
 * @param b given tensor (also)
 * @return Tensor<T>
 */
template <typename T>
Tensor<T> BackwardSubstitution(Tensor<T> const &A, Tensor<T> const &b);

template <typename T>
void Print(Tensor<T> const &mat);

/**
 * @brief Calculate the inverse of a matrix using LU decomposition
 *
 * @tparam T precision
 * @param A matrix
 * @return Tensor<T>
 */
template <typename T>
Tensor<T> InverseMatrix(Tensor<T> const &A);

template <typename T>
T DeterminantFromLU(Tensor<T> const &A) {
  auto U = std::get<1>(LUDecomposition(A));
  T det = 1;

  for (int i = 0; i < A.Rows(); i++) {
    det *= U(i, i);
  }

  return det;
}

template <typename T>
TensorPair<T> LUDecomposition(Tensor<T> const &A) {
  LOG_ASSERT(A.Cols() == A.Rows(), "Matrix must be square", utils::ERROR);

  T scalar;
  auto L = Tensor<T>::SMatrix(A.Rows());
  auto U = Tensor<T>(A);

  for (int i = 0; i < U.Rows(); i++) {
    L(i, i) = 1;
  }

  for (int j = 0; j < U.Rows() - 1; j++) {
    for (int i = j + 1; i < U.Cols(); i++) {
      scalar = -U(i, j) / U(j, j);
      L(i, j) = -scalar;
      U.LinearCombRows(i, j, scalar, i);
    }
  }

  return std::make_pair(L, U);
}

template <typename T>
Tensor<T> InverseMatrix(Tensor<T> const &A) {
  LOG_ASSERT(A.Cols() == A.Rows(), "Matrix must be square", utils::ERROR);
  LOG_ASSERT(std::abs(DeterminantFromLU(A)) != 0, "Matrix is singular",
             utils::ERROR);

  auto I = Tensor<T>::Identity(A.Cols());

  return GaussianElimination(A, I);
}

template <typename T>
// TODO: probably throw an exception if the matrix is singular or becomes it
Tensor<T> GaussianElimination(Tensor<T> A, Tensor<T> b) {
  LOG_ASSERT(A.Cols() == A.Rows(), "Matrix must be square", utils::ERROR);
  LOG_ASSERT(std::abs(DeterminantFromLU(A)) != 0, "Matrix is singular",
             utils::ERROR);

  int N = A.Rows();

  LOG_ASSERT(N == b.Rows(), "A and b are dimensionally incompatible",
             utils::ERROR);

  T scalar;

  for (int j = 0; j < N - 1; j++) {
    // Partial pivoting: find the maximum element in the current column below or
    // at the diagonal
    int pivotRow = j;
    for (int i = j + 1; i < N; i++) {
      if (std::abs(A(i, j)) > std::abs(A(pivotRow, j))) {
        pivotRow = i;
      }
    }

    // If pivotRow is not the current row, swap rows in A and b
    if (pivotRow != j) {
      A.SwapRows(j, pivotRow);
      b.SwapRows(j, pivotRow);
    }

    // Normal elimination
    for (int i = j + 1; i < N; i++) {
      scalar = -A(i, j) / A(j, j);
      A.LinearCombRows(i, j, scalar, i);
      b.LinearCombRows(i, j, scalar, i);
    }
  }

  // HACK: cut off floating point values that are very close to zero
  // NOTE: this could lead to a singular matrix
  double tol = std::numeric_limits<T>::epsilon();
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (std::abs(A(i, j)) < tol) {
        // Singular matrix
        LOG_ASSERT(i != j,
                   "rank(A) < dim(A) !, backward subst division by zero!",
                   utils::ERROR);

        A(i, j) = 0;
      }
    }
  }

  return BackwardSubstitution(A, b);
}

template <typename T>
bool IsUpperTriangular(Tensor<T> const &A) {
  for (int i = 1; i < A.Rows(); i++) {
    for (int j = 0; j < i; j++) {
      if (std::abs(A(i, j))) {
        return false;
      }
    }
  }
  return true;
}

template <typename T>
Tensor<T> BackwardSubstitution(Tensor<T> const &A, Tensor<T> const &b) {
  LOG_ASSERT(A.Cols() == A.Rows(), "A must be a square matrix", utils::ERROR);
  LOG_ASSERT(A.Cols() == b.Rows(), "A and b are dimensionally incompatible",
             utils::ERROR);

  LOG_ASSERT(IsUpperTriangular(A), "A is not upper triangular", utils::ERROR);

  int N = A.Rows();
  auto solution = Tensor<T>::Matrix(N, b.Cols());

  T sum;

  // Iterate through b columns
  for (int k = 0; k < b.Cols(); k++) {
    for (int i = N - 1; i >= 0; i--) {
      sum = 0;
      for (int j = i + 1; j <= N - 1; j++) {
        sum += A(i, j) * solution(j, k);
      }

      solution(i, k) = (b(i, k) - sum) / A(i, i);
    }
  }

  return solution;
}

template <typename T>
Tensor<T> ForwardSubstitution(Tensor<T> const &A, Tensor<T> const &b) {
  LOG_ASSERT(A.Cols() == A.Rows(), "A must be a square matrix", utils::ERROR);
  LOG_ASSERT(A.Cols() == b.Rows(), "A and b are dimensionally incompatible",
             utils::ERROR);

  // "lower"
  //  assert(IsUpperTriangular(A));

  int N = A.Rows();
  auto solution = Tensor<T>::Matrix(N, b.Cols());

  T sum;

  // Iterate through b columns
  for (int k = 0; k < b.Cols(); k++) {
    for (int i = 0; i < N; i++) {
      sum = 0;
      for (int j = 0; j < i; j++) {
        sum += A(i, j) * solution(j, k);
      }
      solution(i, k) = (b(i, k) - sum) / A(i, i);
    }
  }
  return solution;
}

template <typename T>
void Print(Tensor<T> const &mat) {
  for (int i = 0; i < mat.Rows(); i++) {
    for (int j = 0; j < mat.Cols(); j++) std::cout << " " << mat(i, j);

    std::cout << "\n";
  }

  std::cout << "\n";
}

}  // namespace tensor
