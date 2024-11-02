#include <cmath>

#include "../ode_resolver.hpp"
#include "../odes.hpp"

static tensor::Tensor<double> masses = tensor::Tensor<double>::One(3);

tensor::Tensor<double> NewtonGrav(double t, tensor::Tensor<double> const& y);

void CalculateTotalEnergy(tensor::Tensor<double> const& result,
                          func::Range<double>& time_range);

int main(int argc, char const* argv[]) {
  int _conds = argc == 1 ? 0 : std::stoi(argv[1]);
  double h = argc != 3 ? 0.01 : std::stod(argv[2]);
  int _energy = argc != 4 ? 0 : std::stoi(argv[3]);

  if (!_conds) {
    std::cout << "Usage: gravitation <initial_conds (1 or 2)>  <step (default: "
                 "0.01)> <energy 0 or 1 (default: 0)\n";
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

  auto time_range = func::Range<double>::Fixed(0.0, 100.0, h);

  tensor::Tensor<double> results = ode::ODESolver<double>::Builder()
                                       .InitialConditions(initial_conditions)
                                       .SystemFunction(NewtonGrav)
                                       .Method(ode::Method::kRK4)
                                       .CoordinatesRange(time_range)
                                       .Build()
                                       ->Solve();

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

  // Helper lambda to compute the gravitational acceleration between two bodies
  auto grav_accell = [&](int i, int j) -> tensor::Tensor<double> {
    auto r_ij = tensor::Tensor<double>::Vector(3);
    for (int k = 0; k < 3; k++) {
      r_ij(k) = y(3 * i + k) - y(3 * j + k);
    }

    double dist = r_ij.Norm();
    double dist_cubed = std::pow(dist, 3);

    auto accell = tensor::Tensor<double>::Vector(3);

    for (int k = 0; k < 3; ++k) {
      accell(k) = -masses(j) * r_ij(k) / dist_cubed;
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
    // NOTE: copying results in new tensors is not necessary but the code is
    // more understandable this way
    auto r = tensor::Tensor<double>::SMatrix(3);
    auto v = tensor::Tensor<double>::SMatrix(3);

    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        r(j, k) = result(i, 3 * j + k);
        v(j, k) = result(i, 3 * (j + 3) + k);
      }
    }

    double ke = 0.0;
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 3; ++k) {
        ke += 0.5 * masses(j) * v(j, k) * v(j, k);
      }
    }

    double pe = 0.0;
    for (int j = 0; j < 3; ++j) {
      for (int k = j + 1; k < 3; ++k) {
        pe -= masses(j) * masses(k) / sqrt((r.Row(j) - r.Row(k)).NormSquared());
      }
    }

    auto total_energy = ke + pe;

    std ::cout << time_range.Nodes()[i] << " " << ke << " " << pe << " "
               << total_energy << std::endl;
  }
}
