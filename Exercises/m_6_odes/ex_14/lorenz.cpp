#include "../ode_resolver.hpp"
#include "../../m_3_interpolation/range.hpp"
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
  auto initial_tensor = tensor::Tensor<double>::One(3);

  auto time_range = func::Range<double>::Fixed(0.0, 10.0, 0.01);

  // Note method is 1 to 3 corresponding to Euler, Midpoint, RK4
  auto ode_method = static_cast<ode::Method>(std::stoi(argv[1]));

  auto solver = ode::ODESolver<double>::Builder()
                    .InitialConditions(initial_tensor)
                    .CoordinatesRange(time_range)
                    .SystemFunction(LorenzSystem)
                    .Method(ode_method)
                    .Build();

  tensor::Tensor<double> result = solver->Solve();

  ode::Print(result, time_range.Start(), time_range.Step());

  return 0;
}
