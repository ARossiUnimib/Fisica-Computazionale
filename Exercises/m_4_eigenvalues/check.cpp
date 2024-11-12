#include "../m_2_matrices/tensor.hpp"
#include "eigenvalues.hpp"
#include <complex>

int main()
{
    auto A = tensor::Tensor<std::complex<double>>::SMatrix(3);
    A(0,1) = A(1,0) = A(0,2) = A(1,2) = A(2,0) = A(2,1) = 1;
    A(0,0) = {4.0, -1.0};
    A(1,1) = 3.0;
    A(2,2) = 2.0;

    auto vector = eigen::PowerMethodDeflation(A, 100, 3);
    for (int i = 0; i < vector.size(); i++)
    {
        std::cout << vector[i].first << std::endl;
        tensor::Print(vector[i].second);
    }
   return 0;
}
