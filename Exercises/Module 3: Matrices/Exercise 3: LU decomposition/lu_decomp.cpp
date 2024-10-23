#include "../tensor_utils.hpp"

int main() {
  using namespace tensor;

  auto mat = Tensor<double>::SMatrix(3);
  mat(0, 1) = mat(0, 2) = mat(1, 0) = mat(1, 1) = mat(2, 0) = mat(2, 2) = 1;
  mat(0, 0) = mat(2, 1) = 2;
  mat(1, 2) = -2;

  TensorPair<double> decomp = LUDecomposition(mat);

  Print(std::get<0>(decomp));
  Print(std::get<1>(decomp));

  std::cout << (DeterminantFromLU(mat)) << std::endl;

  return 0;
}
