#include <cmath>

#include "../ode_resolver.hpp"
#include "../odes.hpp"

#define PT_1 1
#define PT_2 2

// ODE of the ideal oscillator
tensor::Tensor<double> Oscillator(double t, tensor::Tensor<double> const &y);

// Exact solution of the ode
double ExactSolution(double t) { return sin(t); }

// Utlity function for printing error when using the executable
void PrintError(char const *name);

int main(int argc, char const *argv[]) {
  int _point = argc < 2 ? -1 : std::stoi(argv[1]);
  int _method = argc < 3 ? -1 : std::stoi(argv[2]);

  if (argc < 3 || _point > 2 || _point <= 0 || _method > 3 || _method <= 0) {
    PrintError(argv[0]);
    return 1;
  }

  // 1.0 1.0 1.0
  auto initial_tensor = tensor::Tensor<double>::Vector(2);
  initial_tensor(0) = 0.0;
  initial_tensor(1) = 1.0;

  auto method = static_cast<ode::Method>(_method);

  auto solver_builder =
      ode::ODESolver<double>::Builder()
          .InitialConditions(initial_tensor)
          .SystemFunction(Oscillator)
          // Note that _method corresponds to the Method values
          .Method(method);

  if (_point == PT_1) {
    auto time_range = func::Range<double>::Fixed(0.0, 100 * 3.14, 0.1);

    auto solver =
        solver_builder.Method(method).CoordinatesRange(time_range).Build();

    tensor::Tensor<double> solution = solver->Solve();
    ode::Print(solution, time_range.Start(), time_range.Step());
  }

  if (_point == PT_2) {
    // Set range of h to plot
    auto h_range = func::Range<double>::FixedNum(0.001, 1, 1000);

    for (const auto h : h_range) {
      // Set range of the solution given step h
      auto time_range = func::Range<double>::Fixed(0.0, 10, h);

      auto solver = solver_builder.SystemFunction(Oscillator)
                        .CoordinatesRange(time_range)
                        .Build();

      // Solve the current system with step h
      tensor::Tensor<double> h_solution = solver->Solve();

      // Get last position of h
      size_t max_h_pos = time_range.Nodes().size() - 1;

      // Save difference between last numerical solved and exact solution
      auto value = h_solution(max_h_pos, 0) -
                   ExactSolution(time_range.Nodes()[max_h_pos]);
      std::cout << h << " " << std::abs(value) << std::endl;
    }
  }

  return 0;
}

void PrintError(char const *name) {
  std::cerr
      << "Usage: " << name
      << " <es_point (1, 2)> <ode_method (1: Euler, 2: Midpoint, 3: RK4)>\n";
}

tensor::Tensor<double> Oscillator(double t, tensor::Tensor<double> const &y) {
  auto dydt = tensor::Tensor<double>::Vector(2);

  // ideal armonic oscillator
  dydt(0) = y(1);
  dydt(1) = -y(0);
  return dydt;
}
