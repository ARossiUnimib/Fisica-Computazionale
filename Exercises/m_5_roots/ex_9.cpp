#include "roots.hpp"

int main() {
  std::function<double(double)> f = [](double x) {
    return 1.0 / 2 + x - x * x;
  };

  std::function<double(double)> f_prime = [](double x) { return 1 - 2 * x; };

  for (double x = 0.000001; x <= 1; x *= 1.5) {
    double roots = func::Bisection(f, 0.8, 1.6, x);
    double new_roots = func::NewtonRaphson(f, f_prime, 0.8, x);

    std::cout << x << " " << std::abs(f(roots) - 1.0 / 2 - sqrt(3) / 2) << " "
              << std::abs(f(new_roots) - 1.0 / 2 - sqrt(3) / 2) << std::endl;
  }

  return 0;
}
