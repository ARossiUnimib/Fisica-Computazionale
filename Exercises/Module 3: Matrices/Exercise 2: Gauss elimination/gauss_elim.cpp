#include <iostream>

#include "../tensor.hpp"
#include "../tensor_utils.hpp"

int main() {
  using namespace tensor;

  auto mat = Tensor<double>::SMatrix(3);
  mat(0, 1) = mat(0, 2) = mat(1, 0) = mat(1, 1) = mat(2, 0) = mat(2, 2) = 1;
  mat(0, 0) = mat(2, 1) = 2;
  mat(1, 2) = -2;

  auto i_mat = Tensor<double>(mat);

  auto vec = Tensor<double>::Vector(3);
  vec(0) = 8;
  vec(1) = -2;
  vec(2) = 2;

  auto vec_cp = Tensor<double>(vec);

  Print(mat);
  Print(vec);

  auto solution = GaussianElimination(mat, vec);
  Print(solution);
  std::cout << "Check? " << (vec_cp == i_mat.Dot(solution)) << std::endl;

  auto decomp = LUDecomposition(i_mat);
  Print(std::get<0>(decomp));
  Print(std::get<1>(decomp));
}
