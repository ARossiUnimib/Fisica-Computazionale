#include "../spline.hpp"

#include <cstdlib>

#include "../range.hpp"

double RungeFunction(double x) { return 1 / (1 + 25 * x * x); }

int main(int argc, char const *argv[]) {
  using namespace func;

  int order = argv[1] ? atoi(argv[1]) : 1;
  int step = argv[2] ? atoi(argv[2]) : 1;

  if (order < 1 || order > 2 || step < 1 || step > 2) {
    std::cout << "Usage: " << argv[0] << " <order (1,2)> <step (1, 2)>"
              << std::endl;
    return EXIT_FAILURE;
  }

  double step_value = step == 1 ? 2.0 / 11 : 2.0 / 12;

  auto sample = Range<double>::Fixed(-1, 1, step_value);
  FunctionData<double> f_data(sample, RungeFunction);
  Spline<double> linear_spline(f_data, order);

  auto output_range = Range<double>::Fixed(-1, 1, 2.0 / 100.0);

  for (auto x : output_range) {
    std::cout << x << " " << linear_spline(x) << " " << RungeFunction(x)
              << "\n";
  }
  std::cout << std::endl;

  return EXIT_SUCCESS;
}
