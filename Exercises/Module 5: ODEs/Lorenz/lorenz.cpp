#include "../odes.hpp"

tensor::Tensor<double> LorenzSystem(double t, tensor::Tensor<double> const &y) {
  auto dydt = tensor::Tensor<double>::Vector(3);

  dydt(0) = 10 * (y(1) - y(0));
  dydt(1) = -y(0) * y(2) + 28 * y(0) - y(1);
  dydt(2) = y(0) * y(1) - (8.0 / 3.0) * y(2);

  return dydt;
}
void PrintUsage(char const *name) {
  std::cerr << "Usage: " << name
            << " <ode_method (1: Euler, 2: RK2, 3: RK4) \n";
}

int main(int argc, char const **argv) {
  if (argc < 2) {
    PrintUsage(argv[0]);
    return 1;
  }

  // 1.0 1.0 1.0
  auto initial_tensor = tensor::Tensor<double>::Ones(3);

  auto time_range = func::Range<double>::Fixed(0.0, 10.0, 0.01);

  ode::Function<double> ode_func = LorenzSystem;

  tensor::Tensor<double> result;

  switch (std::stoi(argv[1])) {
    case 1:
      result = ode::Euler(initial_tensor, time_range, ode_func);
      break;
    case 2:
      result = ode::Midpoint(initial_tensor, time_range, ode_func);
      break;
    case 3:
      // TODO: Implement RK4
      break;
  }

  ode::Print(result, time_range.Start(), time_range.Step());

  return 0;
}
