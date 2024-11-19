#include <functional>

#include "../roots.hpp"

std::function<double(double, int)> GetMethodFunction(
    std::function<double(double)> f, int method) {
  auto bisection_method = [&f](double tolerance, int n) {
    return func::Bisection(f, 0.8, 1.6, tolerance, n);
  };

  auto newton_method = [&f](double tolerance, int n) {
    std::function<double(double)> f_prime = [](double x) { return 1 - 2 * x; };
    return func::NewtonRaphson(f, f_prime, 1.0, tolerance, n);
  };

  std::function<double(double, int)> method_func;
  switch (method) {
    case 1:
      method_func = bisection_method;
      break;
    case 2:
      method_func = newton_method;
      break;
  }

  return method_func;
}

int main(int argc, char const *argv[]) {
  const int kEsPoint = argc > 1 ? std::stoi(argv[1]) : 0;
  const int kMethod = argc > 2 ? std::stoi(argv[2]) : 0;
  const double kTolerance = argc > 3 ? std::stod(argv[3]) : 1e-10;
  const int kN = argc > 4 ? std::stoi(argv[4]) : 10000;

  if (kEsPoint > 3 || kEsPoint < 1 || kMethod > 2 || kMethod < 1) {
    std::cout << "Usage: " << argv[0]
              << " <es_point (1,2,3)> <method (1,2)> <tolerance: optional> <n: "
                 "optional>"
              << std::endl;
    return EXIT_FAILURE;
  }

  static std::function<double(double)> f = [](double x) {
    return 1.0 / 2 + x - x * x;
  };

  double analytical_root = 1.0 / 2 + sqrt(3) / 2;

  auto method_func = GetMethodFunction(f, kMethod);

  switch (kEsPoint) {
    case 1: {
      std::cout << "Analytical root: " << analytical_root << "\n";
      std::cout << "Bisection method: " << method_func(kTolerance, kN)
                << std::endl;
      break;
    }

    case 2: {
      std::cout << "x f(root) approx_root-root\n";

      for (double x = 0.000001; x <= 1; x *= 1.5) {
        double roots = method_func(x, kN);
        double error_bisection = std::abs(roots - analytical_root);
        std::cout << x << " " << f(roots) << " " << error_bisection << "\n";
      }

      std::cout << std::endl;
      break;
    }

    case 3: {
      std::cout << "tolerance error_i/error_{i-1}\n";

      // NOTE: newton-raphson converges to fast to get a good estimate of the
      // convergence rate
      // meanwhile bisection convergence with some oscillation on the slope
      // it would be wise to se if the mean remains constant
      double prev_error = 0;
      for (int i = 1; i <= kN; i++) {
        double roots = method_func(kTolerance, i);
        double error = std::abs(roots - analytical_root);

        if (i > 1 && prev_error != 0) {
          std::cout << i << " " << error << " " << error / prev_error << "\n";
        } else {
          std::cout << i << " " << error << " "
                    << "0\n";
        }

        prev_error = error;
      }

    }

    break;
  }
  return EXIT_SUCCESS;
}
