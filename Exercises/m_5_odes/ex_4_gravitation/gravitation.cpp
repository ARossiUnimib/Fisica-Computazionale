#include <cmath>

#include "../ode_resolver.hpp"
#include "../odes.hpp"

static tensor::Tensor<double> masses = tensor::Tensor<double>::One(3);

tensor::Tensor<double> NewtonGrav(double t, tensor::Tensor<double> const& y);

void CalculateTotalEnergy(tensor::Tensor<double> const& result,
                          func::Range<double>& time_range);

int main(int argc, char const* argv[]) {
  auto _conds = argc < 3 ? 0 : std::stoi(argv[1]);
  auto _energy = argc != 3 ? 0 : std::stoi(argv[2]);

  if (!_conds) {
    // print usage <initial_conds (1 or 2)> <energy 1 (default: 0)>
    std::cout << "Usage: gravitation <initial_conds (1 or 2)> <energy 0 or 1"
                 "(default: 0)>\n";
    return EXIT_FAILURE;
  }

  // Initial positions
  auto initial_conditions = tensor::Tensor<double>::Vector(3 * 3 * 2);

  // r1 = 1,0,0
  initial_conditions(0) = 1.0;
  // r2 = -1, 0, 0
  initial_conditions(3) = -1.0;
  // r3 = 0,0,0

  if (_conds == 1) {
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
    masses(0) = 1.6;
    masses(1) = masses(2) = 0.4;
    // Initial velocities
    // v1 = 0, 0.4, 0.0
    initial_conditions(10) = 0.4;
    // v2 = 0, -0.8, 0.7
    initial_conditions(13) = -0.8;
    initial_conditions(14) = 0.7;
    // v3 = 0, -0.8, 0.7
    initial_conditions(16) = -0.8;
    initial_conditions(17) = -0.7;
  }

  auto solver_builder = ode::ODESolver<double>::Builder()
                            .InitialConditions(initial_conditions)
                            .SystemFunction(NewtonGrav)
                            .Method(ode::Method::kRK4);

  // TODO: _h = argv[3] ? std::stod(argv[3]) : 0.01;
  double h = 0.001;

  auto time_range = func::Range<double>::Fixed(0.0, 10.0, h);
  tensor::Tensor<double> results =
      solver_builder.CoordinatesRange(time_range).Build()->Solve();

  if (!_energy) {
    ode::Print(results, time_range.Start(), time_range.Step());
  } else {
    CalculateTotalEnergy(results, time_range);
  }

  return 0;
}

// System of ODEs representing the 3-body problem
tensor::Tensor<double> NewtonGrav(double t, const tensor::Tensor<double>& y) {
  // Create a tensor to hold the derivatives
  auto dydt = tensor::Tensor<double>::Vector(3 * 3 * 2);

  // Gravitational constant (set to 1 for simplicity)
  const double kG = 1.0;

  // Helper lambda to compute the gravitational acceleration between two bodies
  auto grav_accell = [&](int i, int j) -> tensor::Tensor<double> {
    auto r_ij = tensor::Tensor<double>::Vector(3);
    for (int k = 0; k < 3; k++) {
      r_ij(k) = y(3 * i + k) - y(3 * j + k);
    }

    double dist = sqrt(r_ij.Norm());
    double dist_cubed = std::pow(dist, 3);

    auto accell = tensor::Tensor<double>::Vector(3);

    for (int k = 0; k < 3; ++k) {
      accell(k) = -kG * masses(j) * r_ij(k) / dist_cubed;
    }

    return accell;
  };

  // Calculate derivatives for each body
  for (int i = 0; i < 3; ++i) {
    // Position derivatives (velocity)
    for (int k = 0; k < 3; ++k) {
      dydt(3 * i + k) = y(3 * (i + 3) + k);
    }

    // Velocity derivatives (acceleration due to gravity)
    auto acc = tensor::Tensor<double>::Vector(3);

    // The position in the array of the other two bodies is (i + 1) % 3
    // and (i + 2) % 3

    auto acc1 = grav_accell(i, (i + 1) % 3);
    auto acc2 = grav_accell(i, (i + 2) % 3);

    acc = acc1 + acc2;

    for (int k = 0; k < 3; ++k) {
      dydt(3 * (i + 3) + k) = acc(k);
    }
  }

  return dydt;
}

void CalculateTotalEnergy(tensor::Tensor<double> const& result,
                          func::Range<double>& time_range) {
  for (int i = 0; i < result.Rows(); ++i) {
    auto r1 = tensor::Tensor<double>::Vector(3);
    auto r2 = tensor::Tensor<double>::Vector(3);
    auto r3 = tensor::Tensor<double>::Vector(3);
    auto v1 = tensor::Tensor<double>::Vector(3);
    auto v2 = tensor::Tensor<double>::Vector(3);
    auto v3 = tensor::Tensor<double>::Vector(3);

    // transofrm r and v in two matrixes
    // r = [r1_x, r2_x, r3_x]
    //     [r1_y, r2_y, r3_y]
    //     [r1_z, r2_z, r3_z]
    // TODO:
    auto r = tensor::Tensor<double>::SMatrix(3);

    for (int j = 0; j < 3; ++j) {
      r1(j) = result(i, 3 * 0 + j);
      r2(j) = result(i, 3 * 1 + j);
      r3(j) = result(i, 3 * 2 + j);
      v1(j) = result(i, 3 * (0 + 3) + j);
      v2(j) = result(i, 3 * (1 + 3) + j);
      v3(j) = result(i, 3 * (2 + 3) + j);
    }

    // Kinetic energy
    auto ke1 = 0.5 * masses(0) * v1.Norm();
    auto ke2 = 0.5 * masses(1) * v2.Norm();
    auto ke3 = 0.5 * masses(2) * v3.Norm();

    // Potential energy
    auto pe1 = 0.5 * masses(0) * masses(1) / (r1 - r2).Norm();
    pe1 += 0.5 * masses(0) * masses(2) / (r1 - r3).Norm();

    auto pe2 = 0.5 * masses(1) * masses(0) / (r2 - r1).Norm();
    pe2 += 0.5 * masses(1) * masses(2) / (r2 - r3).Norm();

    auto pe3 = 0.5 * masses(2) * masses(0) / (r3 - r1).Norm();
    pe3 += 0.5 * masses(2) * masses(1) / (r3 - r2).Norm();

    auto kinetic_energy = ke1 + ke2 + ke3;
    auto total_energy = ke1 + ke2 + ke3 - (pe1 + pe2 + pe3);
    std ::cout << time_range.Nodes()[i] << " " << kinetic_energy << " "
               << -(pe1 + pe2 + pe3) << std::endl;
  }
}
