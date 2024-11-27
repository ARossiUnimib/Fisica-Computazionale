#include "../../m_5_roots/roots.hpp"
#include "../integration.hpp"

double f(double x) { return 1.0 / (1 + x * x); }

double f_complete(double x) { return exp(-x * x) / (1 + x * x); }

int main() {
  auto n = {3, 5, 6, 8, 10, 12};

  // NOTE / TODO: We should analyze all the different parameters, n scales with
  // n!, we have to calculate eigenvalues with a shifted method not quite
  // precise.
  for (auto i : n) {
    std::cout
        << integration::GaussHermiteQuadrature(-10000.0, 10000.0,
                                               std::function<double(double)>(f),
                                               i, 10000, 10000)
        << " "
        << integration::Trapezoid(-10000.0, 10000.0, f_complete, i * 6000).first
        << " "
        << integration::Simpson(-10000.0, 10000.0, f_complete, i * 6000).first
        << "\n";
  }
  std::cout << std::endl;
  return 0;
}
