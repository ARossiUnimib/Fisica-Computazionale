#include <algorithm>
#include <cstdlib>
#include <functional>
#include <ostream>

#include "../../m_2_matrices/tensor.hpp"
#include "../../m_4_eigenvalues/eigenvalues.hpp"

using namespace std::complex_literals;

const double mL = 8;
const double L = 1.0;

tensor::Tensor<double> GenerateKineticMatrix(int N, double A_bound,
                                             double B_bound) {
  auto kinetic_matrix = tensor::Tensor<double>::SMatrix(N);

  for (int i = 0; i < N - 1; i++) {
    kinetic_matrix(i, i) = kinetic_matrix(i + 1, i + 1) = -2;
    kinetic_matrix(i, i + 1) = kinetic_matrix(i + 1, i) = 1;
  }

  return -N * N / (2 * mL) * kinetic_matrix;
}

double V(double x, double V_0) {
  if (x < 1.0 / 4 || x > 3.0 / 4) {
    return 0;
  }

  return -V_0;  // V_0
}

tensor::Tensor<double> GeneratePotentialMatrix(
    int N, double V_0, std::function<double(double, double)> fV) {
  auto pot_matrix = tensor::Tensor<double>::SMatrix(N);

  for (int i = 0; i < N; i++) {
    pot_matrix(i, i) = V(i * L / N, V_0);
  }

  return pot_matrix;
}

void PrintShiftedEigenStates(tensor::Tensor<double> hamiltonian, int n_iters,
                             int N, double V_0) {
  auto eigensol = eigen::PowerMethodDeflation(hamiltonian, n_iters);

  using EigenPair = std::pair<double, tensor::Tensor<double>>;
  std::sort(
      eigensol.begin(), eigensol.end(),
      [](EigenPair const &a, EigenPair const &b) { return a.first < b.first; });

  // NOTE: eigenvalues are ordered from bigger to smaller
  // reverse the order to study bound states
  for (int i = 0; i < eigensol.size(); i++) {
    std::cout << i * 1.0 / N << " " << V(i * 1.0 / N, V_0) << " ";
    for (int j = 0; j < 4; j++) {
      // From last to j
      auto curr = eigensol[j];

      std::cout << curr.first + mL * curr.second(i) << " ";
    }
    std::cout << "\n";
  }
  std::cout << std::endl;
}

int Usage(char const **argv) {
  std::cout << "Usage: " << argv[0] << " <es_point (1,2,3,4)>" << std::endl;
  return EXIT_FAILURE;
}

int main(int argc, char const **argv) {
  int es_point = argv[1] ? std::stoi(argv[1]) : -1;

  if (es_point == -1) {
    return Usage(argv);
  }

  switch (es_point) {
    case 1: {
      int N = 64;
      auto hamiltonian =
          GenerateKineticMatrix(N, 0, 0) + GeneratePotentialMatrix(N, 10, V);

      PrintShiftedEigenStates(hamiltonian, 100, N, 10);
      break;
    }
    case 2: {
      int N = 64;
      auto hamiltonian =
          GenerateKineticMatrix(N, 0, 0) - GeneratePotentialMatrix(N, -10, V);

      PrintShiftedEigenStates(hamiltonian, 100, N, -10);
      break;
    }
    case 3: {
      int N = 32;

      auto hamiltonian =
          GenerateKineticMatrix(N, 0, 0) + GeneratePotentialMatrix(N, 10, V);

      auto [e_10, v_10] = eigen::InversePowerMethod(hamiltonian, 20);

      hamiltonian =
          GenerateKineticMatrix(N, 0, 0) + GeneratePotentialMatrix(N, 40, V);

      auto [e_40, v_40] = eigen::InversePowerMethod(hamiltonian, 20);

      hamiltonian =
          GenerateKineticMatrix(N, 0, 0) + GeneratePotentialMatrix(N, 80, V);

      auto [e_80, v_80] = eigen::InversePowerMethod(hamiltonian, 20);

      for (int i = 0; i < N; i++) {
        std::cout << i * 1.0 / N << " " << v_10(i) << " " << v_40(i) << " "
                  << v_80(i) << "\n";
      }
      std::cout << std::endl;
      break;
    }
    case 4: {
      double V_0 = 10;

      // TODO: plot 1st, 4th, 8th eigenfunctions with N = 16,32,64

      break;
    }
    default: {
      return Usage(argv);
    }
  }
  return 0;
}
