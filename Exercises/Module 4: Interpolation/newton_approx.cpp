#include "../Module 3: Matrices/tensor_utils.hpp"
#include <cassert>
#include <cmath>
#include <utility>
#include <vector>

using NewtonPoly = std::pair<std::vector<double>, std::vector<double>>;

// TODO: template
NewtonPoly CalculateNewtonPoly(std::vector<double> &x, std::vector<double> &f, double start, double end, double steps)
{
    assert(x.size() == f.size());

    int n = x.size();

    auto A = TensorBuilder<double>::SMatrix(n).Build();

    for (int i = 0; i < n; i++)
    {
        A(i, 0) = f[i];
    }

    for (int j = 1; j < n; j++)
    {
        for (int k = 0; k < n - j; k++)
        {
            A(k, j) = (A(k + 1, j - 1) - A(k, j - 1)) / (x[k + j] - x[k]);
        }
    }

    double a[n];

    for (int i = 0; i < n; i++)
    {
        a[i] = A(0, i);
    }

    auto p = NewtonPoly();
    for (double i = start; i <= end; i += steps)
    {

#ifdef NON_OPTIMIZED
        // Non-optimized version
        double sum = a[0];
        double temp = 1;
        for (int j = 1; j < n; j++)
        {
            temp *= (i - x[j - 1]);
            sum += a[j] * temp;
        }
#else
        double sum = a[n - 1];

        for (int j = n - 2; j >= 0; j--)
        {
            sum = a[j] + (i - x[j]) * sum;
        }
#endif

        p.first.push_back(i);
        p.second.push_back(sum);
    }

    return p;
}

template <typename T> std::vector<T> slice(std::vector<T> const &v, int m, int n)
{
    auto first = v.cbegin() + m;
    auto last = v.cbegin() + n + 1;

    return std::vector<T>(first, last);
}

void print_poly(NewtonPoly &poly)
{
    for (int i = 0; i < poly.first.size(); i++)
    {
        std::cout << poly.first[i] << "\t" << poly.second[i] << std::endl;
    }
}

void CalculateVandermondeMatrix(Tensor<double> &mat, std::vector<double> values)
{
    assert(mat.Rows() == mat.Cols());
    int N = mat.Rows();

    // TODO: create iterator for Tensor class
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            mat(i, j) = pow(values[j], i);
        }
    }
}

#define START 0
#define END 30
#define STEPS 0.25

void select_range(std::vector<double> &x, std::vector<double> &f, int start, int end)
{
    // create a slice avoiding copying elements
    auto x_p = slice(x, start, end);
    auto f_p = slice(f, start, end);

    // convert f_p to Tensor
    auto f_tensor = TensorBuilder<double>::Vector(f_p.size()).Build();

    for (int i = 0; i < f_p.size(); i++)
    {
        f_tensor(i) = f_p[i];
    }

    NewtonPoly new_poly = CalculateNewtonPoly(x_p, f_p, START, END, STEPS);
    print_poly(new_poly);

    auto vande_matrix = TensorBuilder<double>::SMatrix(f_tensor.Rows()).Build();
    auto values = TensorHelper::GaussianElimination(vande_matrix, f_tensor);
}

int main()
{
    auto x = std::vector<double>{0, 10, 15, 20, 22.5, 30};
    auto f = std::vector<double>{0, 227.04, 362.78, 517.35, 602.97, 901.67};

#if !(defined(PART_1) || defined(PART_2) || defined(PART_3))
    std::cout << "Please define parameter PART_1, PART_2 or PART_3" << std::endl;
#endif

    // TODO: we need to split apart xs with f for the inversion of the matrix
    // cause we defined as this the gaussian elimination...
    // TensorHelper::GaussianElimination(Tensor<T> &A, Tensor<T> &b)

#ifdef PART_1
    select_range(x, f, 2, 3);
#endif

#ifdef PART_2
    select_range(x, f, 1, 3);
#endif

#ifdef PART_3
    select_range(x, f, 0, 4);
#endif

    return 0;
}
