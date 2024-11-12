#include "../tensor.hpp"
#include "../tensor_utils.hpp"

int main() {
  using namespace tensor;

#ifndef OPTIONAL_EXERCISE
  auto u_mat = tensor::Tensor<double>::SMatrix(3);
  auto b_vec = tensor::Tensor<double>::Vector(3);

  // 2 1 1
  // 0 1 -2
  // 0 0 1
  u_mat(0, 1) = u_mat(0, 2) = u_mat(1, 1) = u_mat(2, 2) = 1;
  u_mat(0, 0) = 2;
  u_mat(1, 2) = -2;

  // 1 -1 4
  b_vec(0) = 1;
  b_vec(1) = -1;
  b_vec(2) = 4;

  auto solution = BackwardSubstitution(u_mat, b_vec);
  Print(solution);

  auto check = u_mat.Dot(solution);
  Print(check);

  std::cout << (check == b_vec) << std::endl;

#else

  {
    for (int n = 2; n <= 500; n += 10) {
    }
  }

  // TODO: exercise 1.1

#endif

  return 0;
}
