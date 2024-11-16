#include <complex>

#include "../m_2_matrices/tensor_ops.hpp"
#include "../m_2_matrices/tensor_utils.hpp"
#include "../m_3_interpolation/range.hpp"
#include "eigenvalues.hpp"

int main(int argc, char const **argv) {
  int es_point = argv[1] ? std::stoi(argv[1]) : 0;
  auto A = tensor::Tensor<std::complex<double>>::SMatrix(3);

  using namespace std::complex_literals;

  A(0, 0) = 4.0;
  A(0, 1) = -1i;
  A(0, 2) = A(1, 1) = A(2, 0) = 2.0;
  A(1, 0) = 1i;
  A(1, 2) = 2.0 + 7i;
  A(2, 1) = 2.0 - 7i;
  A(2, 2) = -2.0;

  using ComplexSolutionPair =
      std::pair<std::complex<double>, tensor::Tensor<std::complex<double>>>;

  ComplexSolutionPair solution_pair = eigen::PowerMethod(A, 1000);
  ComplexSolutionPair solution_pair_inv = eigen::InversePowerMethod(A, 1000);

  if (!es_point) {
    tensor::Print(solution_pair.second);
    std::cout << "Power Eigenvalue: " << solution_pair.first << std::endl;

    // tensor::Print(solution_pair_inv.second);
    std::cout << "Inverse Power Eigenvalue: " << solution_pair_inv.first
              << std::endl;

    // Check
    auto Ax = (1.0 / solution_pair.first) * A.Dot(solution_pair.second);
    tensor::Print(Ax);

    auto Ax2 =
        (1.0 / solution_pair_inv.first) * A.Dot(solution_pair_inv.second);
    tensor::Print(Ax2);

  } else {
    auto range = func::Range<int>::Fixed(0, 40, 1);
    for (auto n : range) {
      ComplexSolutionPair new_solution_pair = eigen::PowerMethod(A, n);
      ComplexSolutionPair new_solution_pair_inv =
          eigen::InversePowerMethod(A, n);

      std::cout << n << " "
                << std::abs(new_solution_pair.first / solution_pair.first - 1.0)
                << " "
                << std::abs(new_solution_pair_inv.first /
                                solution_pair_inv.first -
                            1.0)
                << std::endl;
    }
  }

  return 0;
}
