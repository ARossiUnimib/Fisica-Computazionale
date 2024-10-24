#include <cmath>
#include <string>

#include "../ode_resolver.hpp"
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
  int _point = argc < 2 ? -1 : std::stoi(argv[1]);
  double _gamma = argc < 3 ? -1 : std::stod(argv[2]);
  double _A = argc < 4 ? -1 : std::stod(argv[3]);

  if (argc < 2 || _point > 3 || _point <= 0) {
    PrintUsage(argv[0]);
    return 1;
  }

  ode::Function<double> ode_func;

  switch (_point) {
    case 1: {
      ode_func = RealOscillator;
      break;
    }

    case 2: {
      if (_gamma >= 2 || _gamma <= 0) {
        PrintUsage(argv[0]);
        return 1;
      }

      Gamma = _gamma;
      ode_func = DampedOscillator;
      break;
    }

    case 3: {
      if (_gamma >= 2 || _gamma <= 0 || _A >= 2 || _A <= 0) {
        PrintUsage(argv[0]);
        return 1;
      }

      Gamma = _gamma;
      _A = _A;
      ode_func = ForcedOscillator;
      break;
    }
  }

  auto time_range = func::Range<double>::Fixed(0.0, 100.0, 0.001);

  auto initial_tensor = tensor::Tensor<double>::Vector(2);
  initial_tensor(0) = 0.0;
  initial_tensor(1) = 1.0;

  tensor::Tensor<double> result = ODESolver<double>::Builder()
                                      .Method(Method::kRK4)
                                      .SystemFunction(ode_func)
                                      .InitialConditions(initial_tensor)
                                      .CoordinatesRange(time_range)
                                      .Build()
                                      .Solve();

  ode::Print(result, time_range.Start(), time_range.Step());
}
