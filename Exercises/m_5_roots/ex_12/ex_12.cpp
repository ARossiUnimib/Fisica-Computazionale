#include <cmath>

#include "../roots.hpp"

int main() {
  // legendre poly at order 10

  auto coefficients = tensor::Tensor<double>::FromData(std::vector<double>{
      -63.0, 0, 3465.0, 0, -30030.0, 0, 90090.0, 0, -109395.0, 0, 46189.0});

  auto f_prime = func::PolynomialRoots(coefficients, 1000, 0.01);

  for (int i = 0; i < f_prime.size(); i++) {
    std::cout << f_prime[i] << std::endl;
  }

    std::cout << "Values:" << std::endl;
  // print values of the polynomial at the roots
  for (int i = 0; i < f_prime.size(); i++) {
    // cycle through the coefficients and the power
    double value = 0.0;
    for (int j = 0; j < 10; j++) {
      value += coefficients(j) * std::pow(f_prime[i], j);
    }
    std::cout << value << std::endl;
  }

  return 0;
}
