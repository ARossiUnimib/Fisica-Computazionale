#include "../tensor.hpp"
#include "../tensor_utils.hpp"

int main()
{
#ifndef OPTIONAL_EXERCISE
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

    {
        using namespace TensorHelper;

        auto solution = BackwardSubstitution(U_mat, b_vec);
        Print(solution);

        auto check = U_mat.Dot(solution);
        Print(check);

        std::cout << (check == b_vec) << std::endl;
    }

#else

    {
        for (int n = 2; n <= 500; n += 10) 
        {

        }
    }

    // TODO: exercise 1.1

#endif

    return 0;
}
