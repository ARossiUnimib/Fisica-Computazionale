#pragma once

#include <cassert>
#include <iostream>
#include <vector>

#include "tensor.hpp"

template <typename T> using TensorPair = std::pair<Tensor<T>, Tensor<T>>;

/**
 * @brief Builder to construct tensors
 *
 * The builder pattern avoids commons mistakes that could happen when declaring
 * tensors preventing access to the inner logic of the Tensor class
 *
 * @tparam T precision
 */
template <typename T> class TensorBuilder;

namespace TensorUtils
{

/**
 * @brief Calculate the determinant of a matrix using LU decomposition
 *
 * @tparam T
 * @param A
 * @return T
 */
template <typename T> T DeterminantFromLU(Tensor<T> const &A);

/**
 * @brief Calculate the inverse of a matrix using LU decomposition
 *
 * @tparam T
 * @param A
 * @return Tensor<T>
 */
template <typename T> TensorPair<T> LUDecomposition(Tensor<T> const &A);

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

template <typename T> void Print(Tensor<T> const &mat);

/**
 * @brief Calculate the inverse of a matrix using LU decomposition
 *
 * @tparam T precision
 * @param A matrix
 * @return Tensor<T>
 */
template <typename T> Tensor<T> InverseMatrix(Tensor<T> const &A);

} // namespace TensorUtils

template <typename T> class TensorBuilder
{
  public:
    static TensorBuilder<T> FromData(std::vector<T>& data, int rows, int cols)
    {
        return TensorBuilder<T>(std::move(data), rows, cols);
    }

    static TensorBuilder<T> FromData(std::vector<T>& data)
    {
        return FromData(std::move(data), data.size(), 1);
    }

    /*
     * Fill an identity tensor in the correct format for each order of the
     * tensor
     */
    TensorBuilder<T> &Identity();

    TensorBuilder<T> &Ones();

    TensorBuilder<T> &Random();

    Tensor<T> Build()
    {
        return m_Temp;
    }

  private:
    TensorBuilder(int rows, int cols) : m_Temp(Tensor<T>({rows, cols}))
    {
    }

    TensorBuilder(std::vector<T> &&data, int rows, int cols)
        : m_Temp(Tensor<T>(data, rows, cols))
    {
    }

  private:
    Tensor<T> m_Temp;
};

// Implementation
#include "tensor_utils.inl"
