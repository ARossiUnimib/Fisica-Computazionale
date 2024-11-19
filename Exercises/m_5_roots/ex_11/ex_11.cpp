#include <cmath>
#include <iostream>
#include <ostream>

#include "../roots.hpp"

int main(int argc, char const *argv[]) {
  const int kN = argc > 1 ? std::stoi(argv[1]) : 100;

  std::function<double(double)> f = [](double x) {
    return sqrt(x + 1) * cos(x / 2) * cos(x / 2) * cos(x / 2);
  };

  std::function<double(double)> f_prime = [](double x) {
    return (cos(x / 2) * cos(x / 2) * (cos(x / 2) - 3 * (x + 1) * sin(x / 2))) /
           (2 * sqrt(x + 1));
  };

  double prev_bisect_error = 1.0;
  double prev_newton_error = 1.0;
  double prev_secant_error = 1.0;

  for (int i = 0; i < kN; i++) {
#define DISTANCE(f_, f, a, b, c) std::abs(f_(f, a, b, 1e-5, i) - M_PI)
    double bisect_error = DISTANCE(func::Bisection, f, 0.8, 2 * M_PI, i);
    double newton_error = DISTANCE(func::NewtonRaphson, f, f_prime, 0.8, i);
    // NOTE: secant method converges with a power of the golden ratio istead of 2
    double secant_error = DISTANCE(func::Secant, f, 0.8, 2 * M_PI, i);
#undef DISTANCE
    // print convergence rate only
    std::cout << i << " "<< bisect_error / prev_bisect_error
              << " "<< newton_error / prev_newton_error
              << " " << secant_error / prev_secant_error << std::endl;

    prev_bisect_error = bisect_error;
    prev_newton_error = newton_error;
    prev_secant_error = secant_error;
  }
  std::cout << std::endl;

  return 0;
}
