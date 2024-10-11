#include "../tensor.hpp"
#include "../tensor_utils.hpp"
#include <iostream>

int main()
{
    auto mat = TensorBuilder<double>::Matrix(3, 3).Build();
    mat(0, 1) = mat(0, 2) = mat(1, 0) = mat(1, 1) = mat(2, 0) = mat(2, 2) = 1;
    mat(0, 0) = mat(2, 1) = 2;
    mat(1, 2) = -2;

    auto i_mat = Tensor<double>(mat);

    auto vec = TensorBuilder<double>::Vector(3).Build();
    vec(0) = 8;
    vec(1) = -2;
    vec(2) = 2;

    auto i_vec = Tensor<double>(vec);

    {
        using namespace TensorUtils;
        Print(mat);
        Print(vec);

        auto solution = GaussianElimination(mat, vec);
        Print(solution);
        std::cout << "Check? " << (i_vec == i_mat.Dot(solution)) << std::endl;

        auto decomp = LUDecomposition(i_mat);
        Print(std::get<0>(decomp));
        Print(std::get<1>(decomp));
    }
}
