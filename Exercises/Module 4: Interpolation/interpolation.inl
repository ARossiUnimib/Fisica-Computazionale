#include "interpolation.hpp"

namespace Interpolation
{

template <typename T> Tensor<T> VandermondeMatrix(std::vector<T> const &values)
{
    auto mat = Tensor<T>::SMatrix(values.size());

    for (int i = 0; i < values.size(); i++)
    {
        for (int j = 0; j < values.size(); j++)
        {
            mat(i, j) = pow(values[i], j);
        }
    }

    return mat;
}

template <typename T>
std::vector<T> DirectCoefficients(FunctionData<T> const &f)
{
    auto f_tensor = Tensor<T>::FromData(f.F());
    auto vande_matrix = VandermondeMatrix(f.X());
    auto values = TensorUtils::GaussianElimination(vande_matrix, f_tensor);
    return values.RawData();
}

template <typename T>
FunctionData<T> DirectPolynomial(FunctionData<T> const &f,
                                 Utils::Range<T> const &range)
{
    std::vector<T> a = DirectCoefficients(f);
    return FunctionUtils::Polynomial(a, range);
}

template <typename T>
std::vector<T> NewtonCoefficients(FunctionData<T> const &f)
{
    int N = f.Size();

    std::vector<T> a(N);

    auto A = Tensor<double>::SMatrix(N);

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

template <typename T>
FunctionData<T> NewtonPolynomial(FunctionData<T> const &f,
                                 Utils::Range<T> const &range)
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
