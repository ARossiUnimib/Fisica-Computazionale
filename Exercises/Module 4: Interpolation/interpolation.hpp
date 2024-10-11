#pragma once

#include "../Module 3: Matrices/tensor_utils.hpp"
#include "function_data.hpp"

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

template <typename T> FunctionData<T> NewtonPolynomial(const FunctionData<T> &f, const Utils::Range<T> &range)
{
    int N = f.Size();

    auto A = TensorBuilder<double>::SMatrix(N).Build();

    for (int i = 0; i < N; i++)
    {
        A(i, 0) = f.F(i);
    }

    for (int j = 1; j < N; j++)
    {
        for (int k = 0; k < N - j; k++)
        {
            A(k, j) = (A(k + 1, j - 1) - A(k, j - 1)) / (f.X(k + j) - f.X(k));
        }
    }

    auto a = std::vector<T>(N);

    for (int i = 0; i < N; i++)
    {
        a[i] = A(0, i);
    }

    FunctionData<T> p;

    for (const auto &&i : range)
    {
        T sum = a[N - 1];

        for (int j = N - 2; j >= 0; j--)
        {
            sum = a[j] + (i - f.X(j)) * sum;
        }

        p.Add(i, sum);
    }

    return p;
}

} // namespace Interpolation
