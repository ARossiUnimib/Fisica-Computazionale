#include "../roots.hpp"

int main(int argc, char const *argv[]) {
  const int kEsPoint = argc > 1 ? std::stoi(argv[1]) : 0;

  if (!kEsPoint) {
    std::cout << "Usage: " << argv[0] << " <es_point (1,2,3)>" << std::endl;
    return EXIT_FAILURE;
  }

  switch (kEsPoint) {
    case 1: {
      // NOTE: the error displayed is costant 1/2 since x_{n+1} = x_n/2
      // the error is in fact spoiled by the analytic expression of f(x)
      std::function<double(double)> f = [](double x) { return x * x; };
      std::function<double(double)> f_prime = [](double x) { return 2 * x; };

      double prev_error = 1.0;
      for (int i = 0; i < 100; i++) {
        double root = func::NewtonRaphson(f, f_prime, 0.8, 1e-20, i);
        std::cout << i << " " << root << " " << std::abs(root) / prev_error
                  << "\n";

        prev_error = std::abs(root);
      }
      std::cout << std::endl;

      break;
    }
    case 2: {
      std::function<double(double)> f = [](double x) {
        return x * x * x - 2 * x + 2;
      };

      std::function<double(double)> f_prime = [](double x) {
        return 3 * x * x - 2;
      };

      // NOTE: near 0.5 and 1 newton raphson it has oscillating solutions
      // for more info see Newton's Fractal
      // https://en.wikipedia.org/wiki/Newton_fractal#/media/File:Newton_z3-2z+2.png
      // the red part does not converge
      //
      // To solve the issue one can:
      // change the initial guess
      // introduce a parameter x_n+1 = x_n - f(x_n) / f'(x_n) * alpha

      for (int i = 0; i < 100; i++) {
        double root_0 = func::NewtonRaphson(f, f_prime, 0.0, 1e-20, i);
        double root_1 = func::NewtonRaphson(f, f_prime, 1.0, 1e-20, i);
        std::cout << i << " " << root_0 << " " << root_1 << "\n";
      }
      std::cout << std::endl;

      break;
    }
    case 3: {
      // NOTE: look at the images generated with the python script
      std::function<double(double)> f = [](double x) {
        return x * x * x - 2 * x * x - 11 * x + 12;
      };

      std::function<double(double)> f_prime = [](double x) {
        return 3 * x * x - 4 * x - 11;
      };

      for (int i = 0; i < 100; i++) {
        double root_0 = func::NewtonRaphson(f, f_prime, 2.352837350, 1e-20, i);
        double root_1 = func::NewtonRaphson(f, f_prime, 2.352836327, 1e-20, i);
        double root_2 = func::NewtonRaphson(f, f_prime, 2.352836323, 1e-20, i);
        std::cout << i << " " << root_0 << " " << root_1 << " " << root_2
                  << "\n";
      }
      std::cout << std::endl;
      break;
    }
  }
}
