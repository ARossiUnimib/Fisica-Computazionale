#include <cmath>
#include <string>

#include "../ode_resolver.hpp"
#include "../odes.hpp"

double kGamma = -1;
double kA = -1;

tensor::Tensor<double> RealOscillator(double t,
                                      tensor::Tensor<double> const &y) {
  auto dydt = tensor::Tensor<double>::Vector(2);

  dydt(0) = y(1);
  dydt(1) = -sin(y(0));

  return dydt;
}

tensor::Tensor<double> DampedOscillator(double t,
                                        tensor::Tensor<double> const &y) {
  auto dydt = tensor::Tensor<double>::Vector(2);

  dydt(0) = y(1);
  dydt(1) = -sin(y(0)) - kGamma * y(1);

  return dydt;
}

tensor::Tensor<double> ForcedOscillator(double t,
                                        tensor::Tensor<double> const &y) {
  auto dydt = tensor::Tensor<double>::Vector(2);

  dydt(0) = y(1);
  dydt(1) = -sin(y(0)) - kGamma * y(1) + kA * sin((2.0 / 3.0) * t);

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

  ode::Function<double> ode_func;

  switch (point) {
    case 1: {
      ode_func = RealOscillator;
      break;
    }

    case 2: {
      if (gamma >= 2 || gamma <= 0) {
        PrintUsage(argv[0]);
        return 1;
      }

      kGamma = gamma;
      ode_func = DampedOscillator;
      break;
    }

    case 3: {
      if (gamma >= 2 || gamma <= 0 || A >= 2 || A <= 0) {
        PrintUsage(argv[0]);
        return 1;
      }

      kGamma = gamma;
      A = A;
      ode_func = ForcedOscillator;
      break;
    }
  }

  auto time_range = func::Range<double>::Fixed(0.0, 100.0, 0.001);

  auto initial_tensor = tensor::Tensor<double>::Vector(2);
  initial_tensor(0) = 0.0;
  initial_tensor(1) = 1;

  tensor::Tensor<double> result = ode::ODESolver<double>::Builder()
                                      .Method(ode::Method::kRK4)
                                      .SystemFunction(ode_func)
                                      .InitialConditions(initial_tensor)
                                      .CoordinatesRange(time_range)
                                      .Build()
                                      ->Solve();

  ode::Print(result, time_range.Start(), time_range.Step());
}
