#include <cmath>
#include <functional>
#include <iostream>
#include <ostream>

#include "../integration.hpp"

double f(double x) { return std::sin(x / 2) * std::sin(x / 2); }

double analytical_sol(double a, double b) {
  return 0.5 * (b - a + std::sin(a) - std::sin(b));
}

int main(int argc, char const *argv[]) {
  double a = 0;

  for (double b = 1; b > 0.001; b /= 1.05) {
    auto [trap_sol, err1] = integration::Trapezoid(a, b, f, 1);
    auto [simp_sol, err2] = integration::Simpson(a, b, f, 1);

    std::cout << b << " " << std::abs(trap_sol - analytical_sol(a, b)) << " "
              << std::abs(simp_sol - analytical_sol(a, b)) << std::endl;
  }

  return 0;
}
