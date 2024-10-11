#pragma once

#include "tensor.hpp"

#include <cassert>
#include <iostream>
#include <tuple>
#include <vector>

#define D_TENSOR_COUPLE std::tuple<Tensor<T>, Tensor<T>>

namespace TensorUtils
{

template <typename T> T DeterminantFromLU(Tensor<T> &A);
template <typename T> D_TENSOR_COUPLE LUDecomposition(Tensor<T> &A);

template <typename T> Tensor<T> GaussianElimination(Tensor<T> &A, Tensor<T> &b);

/*
 * Perform Back Substitution Algorithm on an upper triangular matrix
 *
 * Ax = b
 * A: the matrix
 * b: the given vector
 */
template <typename T> Tensor<T> BackwardSubstitution(Tensor<T> A, Tensor<T> b);

template <typename T> void Print(Tensor<T> mat);

} // namespace TensorUtils

/*
 * Builder to construct tensors
 *
 * The builder pattern avoids commons mistakes that could happen when declaring
 * tensors preventing access to the inner logic of the Tensor class
 */
template <typename T> class TensorBuilder
{
  public:
    static TensorBuilder<T> Vector(int dim)
    {
        return TensorBuilder<T>(dim, 1);
    }

    static TensorBuilder<T> Matrix(int rows, int cols)
    {
        return TensorBuilder<T>(rows, cols);
    }

    static TensorBuilder<T> SMatrix(int size)
    {
        return TensorBuilder<T>(size, size);
    }

    static TensorBuilder<T> FromData(std::vector<T> data, int rows, int cols)
    {
        return TensorBuilder<T>(std::move(data), rows, cols);
    }

    static TensorBuilder<T> FromData(std::vector<T> data)
    {
        return FromData(std::move(data), data.size(), 1);
    }

    /*
     * Fill an identity tensor in the correct format for each order of the tensor
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

    TensorBuilder(std::vector<T> &&data, int rows, int cols) : m_Temp(Tensor<T>(data, rows, cols))
    {
    }

  private:
    Tensor<T> m_Temp;
};

namespace TensorUtils
{

template <typename T> T DeterminantFromLU(Tensor<T> &A)
{
    assert(A.Cols() == A.Rows());
    auto U = std::get<1>(LUDecomposition(A));
    T det = 1;

    for (int i = 0; i < A.Rows(); i++)
    {
        det *= U(i, i);
    }

    return det;
}

template <typename T> D_TENSOR_COUPLE LUDecomposition(Tensor<T> &A)
{
    T scalar;
    auto L = TensorBuilder<T>::SMatrix(A.Rows()).Build();
    auto U = Tensor<T>(A);

    for (int j = 0; j < U.Rows() - 1; j++)
    {
        for (int i = j + 1; i < U.Cols(); i++)
        {
            scalar = -U(i, j) / U(j, j);
            L(i, j) = -scalar;
            U.LinearCombRows(i, j, scalar, i);
        }
    }

    return std::make_tuple(L, U);
}


template <typename T> Tensor<T> InverseMatrix(Tensor<T> A)
{
    assert(A.Cols() == A.Rows());
    assert(DeterminantFromLU(A) != 0);

    auto I = TensorBuilder<T>::SMatrix(A.Cols()).Identity().Build();

    return GaussianElimination(A, I);
}

template <typename T> Tensor<T> GaussianElimination(Tensor<T> &A, Tensor<T> &b)
{
    assert(A.Cols() == A.Rows());
    assert(DeterminantFromLU(A) != 0);

    int N = A.Rows();
    assert(N == b.Rows());

    T scalar;

    for (int j = 0; j < A.Rows() - 1; j++)
    {
        for (int i = j + 1; i < A.Cols(); i++)
        {
            scalar = -A(i, j) / A(j, j);
            A.LinearCombRows(i, j, scalar, i);
            b.LinearCombRows(i, j, scalar, i);
        }
    }

    return BackwardSubstitution(A, b);
}

template <typename T> Tensor<T> BackwardSubstitution(Tensor<T> A, Tensor<T> b)
{
    assert(A.Cols() == A.Rows());
    int N = A.Rows();
    assert(N == b.Rows());


    // TODO: Check if A is upper right triangular
    /*
    bool flag = true;

    for (int i = 0; i < N; i++)
    {
        if (!A(i, i))
        {
            flag = false;
            break;
        }

        for (int j = 0; j < N - i; j++)
        {
        }
    }
    */

    auto solution = TensorBuilder<T>::Matrix(b.Rows(), b.Cols()).Build();

    T sum;

    // Iterate through b columns 
    for (int k = 0; k < b.Cols(); k++)
    {
        for (int i = N - 1; i >= 0; i--)
        {
            sum = 0;
            for (int j = i + 1; j <= N - 1; j++)
            {
                sum += A(i, j) * solution(j, k);
            }

            solution(i, k) = (b(i, k) - sum) / A(i, i);
        }
    }

    return solution;
}

template <typename T> void Print(Tensor<T> mat)
{
    for (int i = 0; i < mat.Rows(); i++)
    {
        for (int j = 0; j < mat.Cols(); j++)
            std::cout << " " << mat(i, j);

        std::cout << "\n";
    }

    std::cout << "\n";
}

} // namespace TensorUtils

template <typename T> TensorBuilder<T> &TensorBuilder<T>::Identity()
{
    // Identity for vector is a tensor filled with ones
    if (m_Temp.m_Cols == 1 || m_Temp.m_Rows == 1)
    {
        return Ones();
    }

    // Identity should be a square matrix
    assert(m_Temp.m_Cols == m_Temp.m_Rows);

    for (int i = 0; i < m_Temp.m_Rows; i++)
    {
        m_Temp(i, i) = 1;
    }

    return *this;
}

template <typename T> TensorBuilder<T> &TensorBuilder<T>::Ones()
{
    for (int i = 0; i < m_Temp.m_Rows * m_Temp.m_Cols; i++)
        m_Temp.m_Data[i] = 1.0;

    return *this;
}

template <typename T> TensorBuilder<T> &TensorBuilder<T>::Random()
{
    for (int i = 0; i < m_Temp.m_Rows * m_Temp.m_Cols; i++)
        m_Temp.m_Data[i] = (T)(std::rand() / (double)RAND_MAX);

    return *this;
}
