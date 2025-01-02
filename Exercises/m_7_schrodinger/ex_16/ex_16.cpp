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

double V_func(double x, double V_0) {
  if (x < 1.0 / 4 || x > 3.0 / 4) {
    return 0;
  }

  return -V_0;  // V_0
}

tensor::Tensor<double> GeneratePotentialMatrix(
    int N, double V_0, std::function<double(double, double)> fV) {
  auto pot_matrix = tensor::Tensor<double>::SMatrix(N);

  for (int i = 0; i < N; i++) {
    pot_matrix(i, i) = V_func(i * L / N, V_0);
  }

  return pot_matrix;
}

auto SolveSchroedinger(tensor::Tensor<double> hamiltonian, int n_iters, int N,
                       double V_0, bool print = true) {
  auto eigensol = eigen::PowerMethodDeflation(hamiltonian, n_iters);

  using EigenPair = std::pair<double, tensor::Tensor<double>>;
  // NOTE: eigenvalues are ordered from bigger to smaller in absolute value: we
  // need to sort the sign
  std::sort(
      eigensol.begin(), eigensol.end(),
      [](EigenPair const &a, EigenPair const &b) { return a.first < b.first; });

  if (!print) {
    return eigensol;
  }

  // NOTE: eigenvalues are ordered from bigger to smaller
  // reverse the order to study bound states
  for (int i = 0; i < eigensol.size(); i++) {
    std::cout << i * 1.0 / N << " " << V_func(i * 1.0 / N, V_0) << " ";
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
  int kEsPoint = argv[1] ? std::stoi(argv[1]) : -1;

  if (kEsPoint == -1) {
    return Usage(argv);
  }

  switch (kEsPoint) {
    case 1: {
      int N = 64;
      auto hamiltonian = GenerateKineticMatrix(N, 0, 0) +
                         GeneratePotentialMatrix(N, 10, V_func);

      SolveSchroedinger(hamiltonian, 100, N, 10);
      break;
    }
    case 2: {
      int N = 64;
      auto hamiltonian = GenerateKineticMatrix(N, 0, 0) +
                         GeneratePotentialMatrix(N, -10, V_func);

      SolveSchroedinger(hamiltonian, 100, N, -10);
      break;
    }
    case 3: {
      int N = 64;
      std::vector<int> V_0 = {10, 40, 80};

      std::vector<std::pair<double, tensor::Tensor<double>>> ground_states{};

      for (int V : V_0) {
        auto hamiltonian = GenerateKineticMatrix(N, 0, 0) +
                           GeneratePotentialMatrix(N, V, V_func);
        auto eigensol = SolveSchroedinger(hamiltonian, 100, N, V, false);
        ground_states.push_back(eigensol[0]);
      }

      for (int i = 0; i < N; i++) {
        std::cout << i * 1.0 / N << " ";
        for (int j = 0; j < 3; j++) {
          std::cout << ground_states[j].second(i) << " ";
        }
        std::cout << "\n";
      }

      break;
    }
    case 4: {
      double V_0 = 10;
      // Eigenfunction number
      std::vector<int> pos_vec = {0, 3, 7};
      std::vector<int> N_vec = {16, 32, 64};

      for (int N : N_vec) {
        std::cout << N << " ";

        auto hamiltonian = GenerateKineticMatrix(N, 0, 0) +
                           GeneratePotentialMatrix(N, 10, V_func);

        auto eigensol = eigen::PowerMethodDeflation(hamiltonian, 100);

        using EigenPair = std::pair<double, tensor::Tensor<double>>;
        std::sort(eigensol.begin(), eigensol.end(),
                  [](EigenPair const &a, EigenPair const &b) {
                    return a.first < b.first;
                  });

        for (int pos : pos_vec) {
          std::cout << eigensol[pos].first << " ";
        }
        std::cout << "\n";

        /*
        NOTE: prints all the eigenfunction at different steps

        for (int i = 0; i < eigensol.size(); i++) {
            std::cout << i * 1.0 / N << " ";

            for (int pos : pos_vec) {
                std::cout << eigensol[pos].second(i) << " ";
            }

            std::cout << "\n";
        }

        std::cout << std::endl;
        */
      }
      std::cout << std::endl;

      break;
    }
    default: {
      return Usage(argv);
    }
  }
  return 0;
}
