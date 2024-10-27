#include <cmath>

#include "../ode_resolver.hpp"
#include "../odes.hpp"

static tensor::Tensor<double> masses = tensor::Tensor<double>::One(3);

tensor::Tensor<double> NewtonGrav(double t, const tensor::Tensor<double>& y);

void CalculateTotalEnergy(std::vector<tensor::Tensor<double>>& result_list,
                          func::Range<double>& time_range);

int main(int argc, char const* argv[]) {
  auto _conds = argc != 2 ? 0 : std::stoi(argv[1]);

  if (!_conds) {
    std::cerr << "No es provided" << std::endl;
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
    initial_conditions(17) = 0.7;
  }

  auto solver_builder = ode::ODESolver<double>::Builder()
                            .InitialConditions(initial_conditions)
                            .SystemFunction(NewtonGrav)
                            .Method(ode::Method::kRK4);

  auto time_range = func::Range<double>::Fixed(0.0, 100.0, 0.1);

  {
    tensor::Tensor<double> results =
        solver_builder.CoordinatesRange(time_range).Build()->Solve();

    // ode::Print(results, time_range.Start(), time_range.Step());
  }
#if true

  // NOTE: smaller steps cause the total energy to fluctuate (example 0.01),
  // meanwhile higher steps cause an increasing trend in the energy (example
  // 0.1)
  // RK methods are not symplectic, so the total energy is not conserved
  // when the objects are near, the total energy becomes more negative requiring
  // a very high step precision, when at large distances the energy becomes more
  // stable
  // A solution could be do implement a variable step size method to avoid this
  // behavior

  auto h_range = func::Range<double>::FixedNum(0.01, 0.1, 15);

  auto result_list = std::vector<tensor::Tensor<double>>{};

  for (const auto& h : h_range) {
    time_range = func::Range<double>::Fixed(0.0, 100.0, h);

    result_list.push_back(
        solver_builder.CoordinatesRange(time_range).Build()->Solve());
  }

  CalculateTotalEnergy(result_list, time_range);

  return 0;
}
#endif

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

void CalculateTotalEnergy(std::vector<tensor::Tensor<double>>& result_list,
                          func::Range<double>& time_range) {
  for (int i = 0; i < result_list[0].Rows(); i++) {
    std::cout << time_range.Nodes()[i];
    for (const auto& results : result_list) {
      double total_energy = 0;
      for (int j = 0; j < 3; j++) {
        // kinetic energy
        double kinetic_energy = 0.5 * masses(j) *
                                (std::pow(results(i, 3 * j + 3), 2) +
                                 std::pow(results(i, 3 * j + 4), 2) +
                                 std::pow(results(i, 3 * j + 5), 2));
        // potential energy
        double potential_energy = 0;
        for (int k = 0; k < 3; k++) {
          if (j != k) {
            double r = 0;
            for (int l = 0; l < 3; l++) {
              r += std::pow(results(i, 3 * j + l) - results(i, 3 * k + l), 2);
            }
            r = std::sqrt(r);
            potential_energy += -masses(j) * masses(k) / r;
          }
        }
        total_energy += kinetic_energy + potential_energy;
      }

      std::cout << "\t" << total_energy;
    }

    std::cout << "\n";
  }

  std::cout << std::endl;
}
