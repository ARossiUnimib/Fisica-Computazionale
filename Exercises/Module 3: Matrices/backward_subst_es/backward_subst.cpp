#include "tensor.hpp"
#include "tensor_utils.hpp"

int main()
{
    auto U_mat = TensorBuilder<double>::Matrix(3, 3).Build();
    auto b_vec = TensorBuilder<double>::Vector(3).Build();

    // 2 1 1
    // 0 1 -2
    // 0 0 1
    U_mat(0, 1) = U_mat(0, 2) = U_mat(1, 1) = U_mat(2, 2) = 1;
    U_mat(0, 0) = 2;
    U_mat(1, 2) = -2;

    b_vec(0) = 1;
    b_vec(1) = -1;
    b_vec(2) = 4;

    auto solution = TensorHelper::BackwardSubstitution(U_mat, b_vec);
    TensorHelper::Print(solution);

    auto check = U_mat.Dot(solution);
    TensorHelper::Print(check);

    std::cout << (check == b_vec) << std::endl;

    return 0;
}
