#include "../m_3_interpolation/range.hpp"
#include "monte_carlo.hpp"

int main(int argc, char const **argv) {
  int es_point = argc > 1 ? std::stoi(argv[1]) : 0;

  std::function<double(std::vector<double>)> f = [](std::vector<double> x) {
    return x[0] * x[0] + x[1] * x[1] <= 1 ? 1 : 0;
  };

  // TODO: dint understand what it means by "analytic prediction"
  switch (es_point) {
    case 0:
      std::cout << "Usage: " << argv[0] << " <exercise point>\n";
      return EXIT_FAILURE;
    case 1: {
      auto N_range = {200, 500, 1000, 2000, 5000};

      std::cout << "1/N\tDistance\tError" << "\n";

      for (auto N : N_range) {
        auto [mean, error] =
            montecarlo::UniformSampling(f, {{0, 1}, {0, 1}}, N);

        std::cout << 1.0 / N << "\t" << std::abs(mean - M_PI / 4) << "\t"
                  << error << "\t" << 1.0 / N << "\n";
      }
      std::cout << std::endl;
      break;
    }
    case 2: {
      // study the error as a function of 1/N and compare it with the analytic
      // prediction
      auto N_range = func::Range<double>::Fixed(1, 10000, 100);

      std::cout << "1/N\tDistance\tError" << "\n";
      break;
    }
  }
  return 0;
}
