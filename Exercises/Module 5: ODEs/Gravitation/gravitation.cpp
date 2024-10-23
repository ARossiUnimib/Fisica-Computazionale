#include <cmath>

#include "../odes.hpp"

static tensor::Tensor<double> masses = tensor::Tensor<double>::One(3);

// System of ODEs representing the 3-body problem
tensor::Tensor<double> NewtonGrav(double t, const tensor::Tensor<double> &y) {
  // Create a tensor to hold the derivatives
  auto dydt = tensor::Tensor<double>::Vector(3 * 3 * 2);

  // Gravitational constant (set to 1 for simplicity)
  const double G = 1.0;

  // Helper lambda to compute the gravitational force between two bodies
  auto grav_accell = [&](int i, int j) -> tensor::Tensor<double> {
    auto r_ij = tensor::Tensor<double>::Vector(3);
    for (int k = 0; k < 3; k++) r_ij(k) = y(3 * i + k) - y(3 * j + k);

    double dist = sqrt(r_ij.Norm());
    double dist_cubed = std::pow(dist, 3);

    auto accell = tensor::Tensor<double>::Vector(3);

    for (int k = 0; k < 3; ++k)
      accell(k) = -G * masses(j) * r_ij(k) / dist_cubed;

    return accell;
  };

  // Calculate derivatives for each body
  for (int i = 0; i < 3; ++i) {
    // Position derivatives (velocity)
    for (int k = 0; k < 3; ++k) dydt(3 * i + k) = y(3 * (i + 3) + k);

    // Velocity derivatives (acceleration due to gravity)
    auto acc = tensor::Tensor<double>::Vector(3);

    auto acc1 = grav_accell(i, (i + 1) % 3);
    auto acc2 = grav_accell(i, (i + 2) % 3);

    acc = acc1 + acc2;

    for (int k = 0; k < 3; ++k) dydt(3 * (i + 3) + k) = acc(k);
  }

  return dydt;
}

int main(int argc, char const *argv[]) {
  auto value = argc != 2 ? 0 : std::stoi(argv[1]);

  if (!value) {
    std::cerr << "No es provided" << std::endl;
    return 1;
  }

  auto initial_conditions = tensor::Tensor<double>::Vector(3 * 3 * 2);

  // Initial positions

  // r1 = 1,0,0
  initial_conditions(0) = 1.0;

  // r2 = -1, 0, 0
  initial_conditions(3) = -1.0;

  // r3 = 0,0,0

  // Initial velocities
  if (value == 1) {
    masses(0) = masses(1) = masses(2) = 0.3;

    // Initial velocities

    // v1 = 0, 0.4, 0.0
    initial_conditions(10) = 0.15;
    initial_conditions(11) = -0.15;

    // v2 = 0, 0.0, 0.0
    initial_conditions(13) = -0.15;
    initial_conditions(14) = 0.15;

    // v3 = 0, 0, 0
  } else {
    // Masses of the bodies
    masses(0) = 1.6;
    masses(1) = masses(2) = 0.4;

    // v1 = 0, 0.4, 0.0
    initial_conditions(10) = 0.4;

    // v2 = 0, -0.8, 0.7
    initial_conditions(13) = -0.8;
    initial_conditions(14) = 0.7;

    // v3 = 0, -0.8, 0.7
    initial_conditions(16) = -0.8;
    initial_conditions(17) = 0.7;
  }

  auto time_range = func::Range<double>::Fixed(0.0, 100.0, 0.01);

  ode::Function<double> ode_func = NewtonGrav;

  tensor::Tensor<double> results =
      ode::Euler(initial_conditions, time_range, ode_func);

  ode::Print(results, time_range.Start(), time_range.Step());

  // TODO: calculate total energy of the system

  return 0;
}
