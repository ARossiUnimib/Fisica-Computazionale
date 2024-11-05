#include <cmath>

#include "../ode_resolver.hpp"
#include "../odes.hpp"

// masses passed to NewtonGrav function without being function arguments
static tensor::Tensor<double> masses = tensor::Tensor<double>::One(3);

// The second order differential equation for the gravitational force is split
// in a system of first order differential equations, to be solved by the
// ODESolver
tensor::Tensor<double> NewtonGrav(double t, tensor::Tensor<double> const& y);

// Helper function for the second point of the exercise
void CalculateTotalEnergy(tensor::Tensor<double> const& result,
                          func::Range<double>& time_range);

int main(int argc, char const* argv[]) {
  // Interpreting input arguments
  int _conds = argc == 1 ? 0 : std::stoi(argv[1]);
  double h = argc != 3 ? 0.01 : std::stod(argv[2]);
  int _energy = argc != 4 ? 0 : std::stoi(argv[3]);

  if (!_conds) {
    std::cout << "Usage: gravitation <initial_conds (1 or 2)>  <step (default: "
                 "0.01)> <energy 0 or 1 (default: 0)\n";
    return EXIT_FAILURE;
  }

  // NOTE: the positions and the velocities are stored in this way:
  // r1_x r1_y ... r3_z v1_x v1_y v1_z ...
  // Position 1 of axis x will be y(3 * 0 + 0)
  // Velocity 1 of axis x will be y(3 * (0 + 3) + 0)
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

tensor::Tensor<double> NewtonGrav(double t, const tensor::Tensor<double>& y) {
  auto dydt = tensor::Tensor<double>::Vector(3 * 3 * 2);

  // Helper lambda to compute the gravitational acceleration between two bodies
  // i and j indicates body 0, 1 or 2
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

  for (int i = 0; i < 3; ++i) {
    for (int k = 0; k < 3; ++k) {
      // Update positions
      dydt(3 * i + k) = y(3 * (i + 3) + k);
    }

    // Velocity derivatives (acceleration due to gravity)
    auto acc = tensor::Tensor<double>::Vector(3);

    // The position in the array of the other two bodies is (i + 1) % 3
    // and (i + 2) % 3
    auto acc1 = grav_accell(i, (i + 1) % 3);
    auto acc2 = grav_accell(i, (i + 2) % 3);

    for (int k = 0; k < 3; ++k) {
      // Update velocities
      dydt(3 * (i + 3) + k) = (acc1 + acc2)(k);
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
        pe -= masses(j) * masses(k) / (r.Row(j) - r.Row(k)).Norm();
      }
    }

    double total_energy = ke + pe;

    std ::cout << time_range.Nodes()[i] << " " << ke << " " << pe << " "
               << total_energy << std::endl;
  }
}
