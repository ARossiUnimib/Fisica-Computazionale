#include <cmath>
#include <iostream>
#include <ostream>

#include "../roots.hpp"

int main() {
  std::function<double(double)> f = [](double x) {
    return sqrt(x + 1) * cos(x / 2) * cos(x / 2) * cos(x / 2);
  };

  std::function<double(double)> f_prime = [](double x) {
    return pow(cos(x / 2), 3) / (2 * sqrt(x + 1)) -
           3.0 / 2 * sqrt(x + 1) * sin(x / 2) * cos(x / 2) * cos(x / 2);
  };

  double bisect_root = func::Bisection(f, 0.0, 2 * M_PI);
  std::cout << "Bisection root: " << bisect_root << std::endl;

  double newton_root = func::NewtonRaphson(f, f_prime, 3.0);
  std::cout << "Newton root: " << newton_root << std::endl;

  double secant_root = func::Secant(f, 0.0, 2 * M_PI);
  std::cout << "Secant root: " << secant_root << std::endl;
  return 0;
}
