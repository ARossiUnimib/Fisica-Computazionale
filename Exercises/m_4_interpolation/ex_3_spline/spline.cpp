#include "../spline.hpp"

#include <cstdlib>

#include "../range.hpp"

double RungeFunction(double x) { return 1 / (1 + 25 * x * x); }

int main() {
  auto x_range = func::Range<double>::Fixed(-1, 1, 2.0 / 11.0);
  std::vector<double> runge_values;

  for (auto x : x_range) {
    runge_values.push_back(RungeFunction(x));
  }

  func::FunctionData f_data(x_range.Nodes(), runge_values);

  Spline<double> linear_spline(f_data, 1);
  Spline<double> quad_spline(f_data, 2);

  auto output_range = func::Range<double>::Fixed(-1, 1, 2.0 / 100.0);

    int i = 0;
  for (auto x : output_range) {
    std::cout << x << " " << quad_spline(x) << " " << RungeFunction(x) << "\n";
        i++;
  }
  std::cout << std::endl;

  return EXIT_SUCCESS;
}
