#pragma once

#include "../Module 3: Matrices/tensor_utils.hpp"
#include "function_data.hpp"

namespace Interpolation
{

/**
 * @brief Compute the Vandermonde matrix
 *
 * @tparam T
 * @param values
 * @return Tensor<T>
 */
template <typename T> Tensor<T> VandermondeMatrix(const std::vector<T> &values);

/**
 * @brief Retrieve polynomial coefficients from direct interpolation using
 * Vandermonde matrix
 *
 * @tparam T
 * @param f
 * @return std::vector<T> coefficients of size f.Size()
 */
template <typename T> std::vector<T> DirectCoefficients(const FunctionData<T> &f);

/**
 * @brief Direct interpolation using Vandermonde matrix
 *
 * @tparam T
 * @param f
 * @param range Values to be extrapolated from polynomial
 * @return FunctionData<T> polynomial data
 */
template <typename T> FunctionData<T> DirectPolynomial(const FunctionData<T> &f, const Utils::Range<T> &range);

/**
 * @brief Compute the Newton coefficients of the function
 * @note use macro _SLOW_ALGORITHM to use the algorithm given during the course
 *
 * @tparam T
 * @param f
 * @return std::vector<T>
 */
template <typename T> std::vector<T> NewtonCoefficients(const FunctionData<T> &f);

/**
 * @brief Compute the Newton polynomial of the function
 *
 * @tparam T
 * @param f
 * @param range
 * @return FunctionData<T>
 */
template <typename T> FunctionData<T> NewtonPolynomial(const FunctionData<T> &f, const Utils::Range<T> &range);

} // namespace Interpolation

namespace Interpolation
{

template <typename T> Tensor<T> VandermondeMatrix(const std::vector<T> &values)
{
    auto mat = TensorBuilder<T>::SMatrix(values.size()).Build();

    for (int i = 0; i < values.size(); i++)
    {
        for (int j = 0; j < values.size(); j++)
        {
            mat(i, j) = pow(values[i], j);
        }
    }

    return mat;
}

template <typename T> std::vector<T> DirectCoefficients(const FunctionData<T> &f)
{
    auto f_tensor = TensorBuilder<T>::FromData(f.F()).Build();
    auto vande_matrix = VandermondeMatrix(f.X());
    auto values = TensorUtils::GaussianElimination(vande_matrix, f_tensor);
    return values.RawData();
}

template <typename T> FunctionData<T> DirectPolynomial(const FunctionData<T> &f, const Utils::Range<T> &range)
{
    std::vector<T> a = DirectCoefficients(f);
    return FunctionUtils::Polynomial(a, range);
}

template <typename T> std::vector<T> NewtonCoefficients(const FunctionData<T> &f)
{
    int N = f.Size();

    std::vector<T> a(N);

    auto A = TensorBuilder<double>::SMatrix(N).Build();

    for (int i = 0; i < N; i++)
    {
        A(i, 0) = f.F(i);
    }

    for (int j = 1; j < N; j++)
    {
        for (int i = 0; i < N - j; i++)
        {
            A(i, j) = (A(i + 1, j - 1) - A(i, j - 1)) / (f.X(i + j) - f.X(i));
        }
    }

    for (int i = 0; i < N; i++)
    {
        a[i] = A(0, i);
    }

    return a;
}

template <typename T> FunctionData<T> NewtonPolynomial(const FunctionData<T> &f, const Utils::Range<T> &range)
{
    std::vector<T> a = NewtonCoefficients(f);
    int N = f.Size();

    FunctionData<T> p;
    T sum;

    for (auto i : range)
    {
        sum = a[N - 1];

        for (int j = N - 2; j >= 0; j--)
        {
            sum = a[j] + (i - f.X(j)) * sum;
        }

        p.Add(i, sum);
    }

    return p;
}

} // namespace Interpolation
