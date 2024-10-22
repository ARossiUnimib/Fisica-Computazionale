#include "tensor_utils.hpp"

namespace TensorUtils
{

template <typename T> T DeterminantFromLU(Tensor<T> const &A)
{
    auto U = std::get<1>(LUDecomposition(A));
    T det = 1;

    for (int i = 0; i < A.Rows(); i++)
    {
        det *= U(i, i);
    }

    return det;
}

template <typename T> TensorPair<T> LUDecomposition(Tensor<T> const &A)
{
    assert(A.Cols() == A.Rows());

    T scalar;
    auto L = Tensor<T>::SMatrix(A.Rows());
    auto U = Tensor<T>(A);

    for (int i = 0; i < U.Rows(); i++)
    {
        L(i, i) = 1;
    }

    for (int j = 0; j < U.Rows() - 1; j++)
    {
        for (int i = j + 1; i < U.Cols(); i++)
        {
            scalar = -U(i, j) / U(j, j);
            L(i, j) = -scalar;
            U.LinearCombRows(i, j, scalar, i);
        }
    }

    return std::make_pair(L, U);
}

template <typename T> Tensor<T> InverseMatrix(Tensor<T> &A)
{
    assert(A.Cols() == A.Rows());
    assert(DeterminantFromLU(A) != 0);

    auto I = Tensor<T>::Identity(A.Cols());

    return GaussianElimination(A, I);
}

template <typename T>
Tensor<T> GaussianElimination(Tensor<T> A, Tensor<T> b)
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

template <typename T> bool IsUpperTriangular(Tensor<T> const &A)
{
    for (int i = 1; i < A.Rows(); i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (A(i, j))
            {
                return false;
            }
        }
    }
    return true;
}

template <typename T>
Tensor<T> BackwardSubstitution(Tensor<T> const &A, Tensor<T> const &b)
{
    assert(A.Cols() == A.Rows());
    assert(A.Rows() == b.Rows());

    Print(A);

    assert(IsUpperTriangular(A));

    int N = A.Rows();
    auto solution = Tensor<T>::Matrix(N, b.Cols());

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

template <typename T> void Print(Tensor<T> const &mat)
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
