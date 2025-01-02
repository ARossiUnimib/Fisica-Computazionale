#include <cmath>

#include "../roots.hpp"

int main() {
  auto legendre_coeffs = tensor::Tensor<double>::FromData(std::vector<double>{
      -63.0, 0, 3465.0, 0, -30030.0, 0, 90090.0, 0, -109395.0, 0, 46189.0});

  auto hermite_coeffs = tensor::Tensor<double>::FromData(
      std::vector<double>{-120, 0, 720, 0, -480, 0, 64});

  auto legendre_roots = func::PolynomialRoots(legendre_coeffs, 1000, 0.01);
  std::cout << "Legendre roots:" << std::endl;
  for (int i = 0; i < legendre_roots.size(); i++) {
    std::cout << legendre_roots[i] << std::endl;
  }


  std::cout << "Legendre Check:" << std::endl;
  // print values of the polynomial at the roots
  for (int i = 0; i < legendre_roots.size(); i++) {
    // cycle through the coefficients and the power
    double value = 0.0;
    for (int j = 0; j < legendre_coeffs.Rows(); j++) {
      value += legendre_coeffs(j) * std::pow(legendre_roots[i], j);
    }
    std::cout << value << std::endl;
  }

  auto hermite_roots = func::PolynomialRoots(hermite_coeffs, 1000, 0.3);
  std::cout << "Hermite roots:" << std::endl;
  for (int i = 0; i < legendre_roots.size(); i++) {
    std::cout << hermite_roots[i] << std::endl;
  }


  std::cout << "Hermite Check:" << std::endl;
  // print values of the polynomial at the roots
  for (int i = 0; i < hermite_roots.size(); i++) {
    // cycle through the coefficients and the power
    double value = 0.0;
    for (int j = 0; j < hermite_coeffs.Rows(); j++) {
      value += hermite_coeffs(j) * std::pow(hermite_roots[i], j);
    }
    std::cout << value << std::endl;
  }


  return 0;
}
