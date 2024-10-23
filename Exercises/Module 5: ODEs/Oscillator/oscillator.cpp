#include <cmath>

#include "../odes.hpp"

#define PT_1 1
#define PT_2 2

#define EULER 1
#define MIDPOINT 2
#define RK4 3

tensor::Tensor<double> Oscillator(double t, tensor::Tensor<double> const &y); 

void PrintError(char const *name); 

double ExactSolution(double t) { return sin(t); }


using MethodFunc = tensor::Tensor<double> (*)(tensor::Tensor<double> const &,
                                              func::Range<double> const &,
                                              ode::Function<double>);

int main(int argc, char const *argv[]) {
  int point = argc < 2 ? -1 : std::stoi(argv[1]);
  int method = argc < 3 ? -1 : std::stoi(argv[2]);

  if (argc < 3 || point > 2 || point <= 0 || method > 3 || method <= 0) {
    PrintError(argv[0]);
    return 1;
  }

  // 1.0 1.0 1.0
  auto initial_tensor = tensor::Tensor<double>::Vector(2);
  initial_tensor(0) = 0.0;
  initial_tensor(1) = 1.0;

  auto time_range = func::Range<double>::Fixed(0.0, 100.0, 0.001);

  ode::Function<double> ode_func = Oscillator;
  MethodFunc method_func = nullptr;
  tensor::Tensor<double> results;

  switch (method) {
    case EULER: {
      method_func = reinterpret_cast<MethodFunc>(ode::Euler<double>);
      break;
    }
    case MIDPOINT: {
      method_func = reinterpret_cast<MethodFunc>(ode::Midpoint<double>);
      break;
    }
  }

  if (point == PT_1) {
    results = method_func(initial_tensor, time_range, ode_func);
    ode::Print(results, time_range.Start(), time_range.Step());
  }

  if (point == PT_2) {
    auto h_range = func::Range<double>::FixedNum(0.001, 1, 1000);

    func::Range<double> time_range;

    for (const auto h : h_range) {
      time_range = func::Range<double>::Fixed(0.0, 10, h);

      auto max_h = time_range.Nodes().size() - 1;

      auto mat = method_func(initial_tensor, time_range, ode_func);

      auto value = mat(max_h, 0) - ExactSolution(time_range.Nodes()[max_h]);

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

  // real armonic oscillator
  dydt(0) = y(1);
  dydt(1) = -y(0);
  return dydt;
}
