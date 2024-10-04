#pragma once

#include "tensor.hpp"

#include <cassert>
#include <iostream>

/*
 * Perform Back Substitution Algorithm on an upper triangular matrix
 *
 * Ax = b
 * A: the matrix
 * b: the given vector
 */

namespace TensorHelper
{

template <typename T> Tensor<T> BackwardSubstitution(Tensor<T> A, Tensor<T> b);

template <typename T> void Print(Tensor<T> mat);

} // namespace TensorHelper

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

    /*
     * Fill an identity tensor in the correct format for each order of the tensor
     */
    TensorBuilder<T> &Identity()
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

    TensorBuilder<T> &Ones()
    {
        for (int i = 0; i < m_Temp.m_Rows * m_Temp.m_Cols; i++)
            m_Temp.m_Data[i] = 1.0;

        return *this;
    }

    TensorBuilder<T> &Random()
    {
        for (int i = 0; i < m_Temp.m_Rows * m_Temp.m_Cols; i++)
            m_Temp.m_Data[i] = (T)(std::rand() / (double)RAND_MAX);

        return *this;
    }

    Tensor<T> Build()
    {
        return m_Temp;
    }

  private:
    TensorBuilder(int rows, int cols) : m_Temp(Tensor<T>({rows, cols}))
    {
    }

  private:
    Tensor<T> m_Temp;
};

namespace TensorHelper
{

template <typename T> Tensor<T> BackwardSubstitution(Tensor<T> A, Tensor<T> b)
{
    // Check if A is a square matrix
    assert(A.Rows() == A.Cols());

    int N = A.Rows();

    // TODO: Check if A is upper right triangular

    auto solution = TensorBuilder<double>::Vector(3).Build();

    double sum;

    for (int i = N - 1; i >= 0; i--)
    {
        sum = 0;
        for (int j = i + 1; j <= N - 1; j++)
        {
            sum += A(i, j) * solution(j);
        }

        solution(i) = (b(i) - sum) / A(i, i);
    }

    return solution;
}

template <typename T> void Print(Tensor<T> mat)
{
    for (int i = 0; i < mat.Rows(); i++)
    {
        for (int j = 0; j < mat.Cols(); j++)
            std::cout << mat(i, j);

        std::cout << "\n";
    }
}

} // namespace TensorHelper
