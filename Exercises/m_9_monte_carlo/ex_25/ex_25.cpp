#include <cmath>

#include "../monte_carlo.hpp"

double f(double const x) { return x * std::cos(x); }

double g(double const x) { return std::cos(x); }

double generator(double const x) { return std::asin(x); }

// FIXME: I don't know if ImportanteSampling is working correctly.
int main(int argc, char** argv) {
  for (int i = 1; i <= 1e3; i++) {
    auto [mean, error] = montecarlo::UniformSampling<double>(
        [](std::vector<double> const& x) { return x[0] * std::cos(x[0]); },
        {{0, M_PI / 2}}, i);

    auto [mean2, error2] = montecarlo::ImportanteSampling<double>(
        f, g, generator, {0, M_PI / 2}, i);
    std::cout << i << " " << error << " " << error2 << std::endl;
  }

  return 0;
}
