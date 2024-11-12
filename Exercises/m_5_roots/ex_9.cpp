#include "roots.hpp"

int main() {
  std::function<double(double)> f = [](double x) {
    return 1.0 / 2 + x - x * x;
  };

  double roots = func::Bisection(f, 0.0, 0.5);

  return 0;
}
