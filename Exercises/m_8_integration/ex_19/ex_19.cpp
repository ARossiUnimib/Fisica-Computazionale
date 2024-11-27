#include "../integration.hpp"

double f(double x) { return 4.0 / (1 + x * x); }

int main() {
  int N = 1000;

  for (int n = 0; n < N; n++) {
    double dx = 1.0 / n;
    auto [trapez, _] = integration::Trapezoid(0.0, 1.0, f, n);
    auto [simpson, __] = integration::Simpson(0.0, 1.0, f, n);
    auto [quad, ___] = integration::GaussLegendreQuadrature(0.0, 1.0, f, n);
    double err_trapez = std::abs(M_PI - trapez);
    double err_simpson = std::abs(M_PI - simpson);
    double err_quad = std::abs(M_PI - quad);

    std::cout << dx << " " << err_trapez << " " << err_simpson << " "
              << err_quad << "\n";
  }

  std::cout << std::endl;

  return 0;
}
