#include <cmath>

#include "../m_3_interpolation/range.hpp"
#include "monte_carlo.hpp"

// e^-x
double GenerateExponential(double x) {
  // log base e
  return -std::log(1 - x);
}

// x*e^(-x^2)
double GenerateOddGaussian(double x) { return std::sqrt(-std::log(1 - 2 * x)); }

int main(int argc, char *argv[]) {
  int es_point = argc > 1 ? std::stoi(argv[1]) : 0;
  int n_iter = argc > 2 ? std::stoi(argv[2]) : 1e6;

  std::vector<double> frequencies{};

  auto range_val = func::Range<double>::FixedNum(0, 10, 20);

  switch (es_point) {
    case 0: {
      std::cout << "Usage: " << argv[0]
                << " <exercise point> <N gen numbers>\n";
      return EXIT_FAILURE;
    }
    case 1: {
      // Domain is [0,1) (using inclusive boundaries since it's very unlikely to
      // happen to be 1)
      auto range_gen = func::Range<double>::Fixed(0, 1);
      frequencies = montecarlo::AccumulateValues(n_iter, range_gen, range_val,
                                                 GenerateExponential);
      break;
    }
    case 2: {
      // Domain is [0,1/2)
      auto range_gen = func::Range<double>::Fixed(0, 1.0 / 2);
      frequencies = montecarlo::AccumulateValues(n_iter, range_gen, range_val,
                                                 GenerateOddGaussian);
      break;
    }
  }

  for (int i = 0; i < frequencies.size() - 1; i++) {
    // Print interval range and frequencies for histogram plot
    std::cout << range_val.Start() + (i + 1) * range_val.Step() << "\t"
              << frequencies[i] << "\n";
  }
  std::cout << std::endl;

  return 0;
}
