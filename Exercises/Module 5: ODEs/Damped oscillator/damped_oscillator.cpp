#include <cmath>
#include <string>

#include "../odes.hpp"

double Gamma = -1;
double A = -1;

tensor::Tensor<double> RealOscillator(double t,
                                      tensor::Tensor<double> const &y) {
  auto dydt = tensor::Tensor<double>::Vector(2);

  // real armonic oscillator
  dydt(0) = y(1);
  dydt(1) = -sin(y(0));

  return dydt;
}

tensor::Tensor<double> DampedOscillator(double t,
                                        tensor::Tensor<double> const &y) {
  auto dydt = tensor::Tensor<double>::Vector(2);

  // real armonic oscillator
  dydt(0) = y(1);
  dydt(1) = -sin(y(0)) - Gamma * y(1);

  return dydt;
}

tensor::Tensor<double> ForcedOscillator(double t,
                                        tensor::Tensor<double> const &y) {
  auto dydt = tensor::Tensor<double>::Vector(2);

  // real armonic oscillator
  dydt(0) = y(1);
  dydt(1) = -sin(y(0)) - Gamma * y(1) + A * sin((2.0 / 3.0) * t);

  return dydt;
}

void PrintUsage(char const *name) {
  std::cerr << "Usage: " << name
            << " <es_point (1, 2, 3)> <gamma (es 2, 3)> <A (es 3)> \n";
}

int main(int argc, const char *argv[]) {
  int point = argc < 2 ? -1 : std::stoi(argv[1]);
  double gamma = argc < 3 ? -1 : std::stod(argv[2]);
  double A = argc < 4 ? -1 : std::stod(argv[3]);

  if (argc < 2 || point > 3 || point <= 0) {
    PrintUsage(argv[0]);
    return 1;
  }

  auto initial_tensor = tensor::Tensor<double>::Vector(2);
  initial_tensor(0) = 0.0;
  initial_tensor(1) = 1.0;

  auto time_range = func::Range<double>::Fixed(0.0, 100.0, 0.001);

  ode::Function<double> ode_func;

  switch (point) {
    case 1:
      ode_func = RealOscillator;
      break;

    case 2:

      if (gamma >= 2 || gamma <= 0) {
        PrintUsage(argv[0]);
        return 1;
      }

      Gamma = gamma;

      ode_func = DampedOscillator;

      break;

    case 3:

      if (gamma >= 2 || gamma <= 0 || A >= 2 || A <= 0) {
        PrintUsage(argv[0]);
        return 1;
      }

      Gamma = gamma;
      A = A;

      ode_func = ForcedOscillator;

      break;
  }

  tensor::Tensor<double> result =
      ode::Midpoint(initial_tensor, time_range, ode_func);

  ode::Print(result, time_range.Start(), time_range.Step());
}
