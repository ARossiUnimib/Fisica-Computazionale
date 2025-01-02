#include <complex>

#include "../../m_2_matrices/tensor.hpp"
#include "../../m_3_interpolation/range.hpp"
#include "../eigenvalues.hpp"

using ComplexSolutionPair =
    std::pair<std::complex<double>, tensor::Tensor<std::complex<double>>>;

int main(int argc, char const **argv) {
  int es_point = argv[1] ? std::stoi(argv[1]) : 0;
  auto A = tensor::Tensor<std::complex<double>>::SMatrix(3);

  using namespace std::complex_literals;

  A(0, 0) = 4.0; // 8.0 for similar eigenvalues
  A(0, 1) = -1i;
  A(0, 2) = A(1, 1) = A(2, 0) = 2.0;
  A(1, 0) = 1i;
  A(1, 2) = 2.0 + 7i;
  A(2, 1) = 2.0 - 7i;
  A(2, 2) = -2.0;

  ComplexSolutionPair solution_pair = eigen::PowerMethod(A, 1000);
  ComplexSolutionPair solution_pair_inv = eigen::InversePowerMethod(A, 1000);
  switch (es_point) {
    case 1: {
      std::cout << "Power Eigenvalue: " << solution_pair.first << std::endl;
      std::cout << "Power Eigenvector: " << std::endl;
      tensor::Print(solution_pair.second);

      std::cout << "Ax= " << std::endl;
      auto Ax = (1.0 / solution_pair.first) * A.Dot(solution_pair.second);
      tensor::Print(Ax);
      break;
    }
    case 2: {
      std::cout << "Inverse Power Eigenvalue: " << solution_pair_inv.first
                << std::endl;
      std::cout << "Inverse Power Eigenvector: " << std::endl;
      tensor::Print(solution_pair_inv.second);

      std::cout << "Ax= " << std::endl;
      auto Ax =
          (1.0 / solution_pair_inv.first) * A.Dot(solution_pair_inv.second);
      tensor::Print(Ax);
      break;
    }
    case 3: {
      auto range = func::Range<int>::Fixed(0, 40, 1);

      for (auto n : range) {
        ComplexSolutionPair new_solution_pair = eigen::PowerMethod(A, n);
        ComplexSolutionPair new_solution_pair_inv =
            eigen::InversePowerMethod(A, n);

        std::cout << n << " "
                  << std::abs(new_solution_pair.first / solution_pair.first -
                              1.0)
                  << " "
                  << std::abs(new_solution_pair_inv.first /
                                  solution_pair_inv.first -
                              1.0)
                  << std::endl;
      }
      break;
    }
    case 4: {
      std::vector<
          std::pair<std::complex<double>, tensor::Tensor<std::complex<double>>>>
          sol = eigen::PowerMethodDeflation(A, 1000);

      for (auto &s : sol) {
        std::cout << s.first << "\n";
      }
      std::cout << std::endl;

      break;
    }
    default:
      return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
