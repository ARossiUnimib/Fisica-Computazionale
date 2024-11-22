#include <cmath>

#include "time_evolution.hpp"

Complex Phi0(ComplexTensor const& x, ComplexTensor const& y, double sigma_2,
             ComplexTensor const& p) {
  return std::exp(-(x - y).NormSquared() / sigma_2 + 1i * (p.Dot(x))(0));
}

int main() {
  // N1=N2=N
  int N = 64;
  double L = 1.0;

  double a = L / N;

  auto a_vec = ComplexTensor::One(2);
  a_vec(0) = a_vec(1) = a;

  double m = 8 / L;
  double tau = 0.002;
  double dt = tau / m;

  auto f_i = [&L, &m](ComplexTensor const& x) {
    // FIXME: why do we need negative sign here?
    auto p0 = ComplexTensor::FromData({-(2 * M_PI) * 6 / L, 0});
    auto y = ComplexTensor::Vector(2);
    y(0) = L / 4;
    y(1) = L / 2;
    double sigma = 0.5 / m;

    return Phi0(x, y, sigma * sigma, p0);
  };

  PotentialFunc V_func = [&m, &L, &a](ComplexTensor const& x) {
    if (x(0).real() <= L / 2 + a && x(0).real() >= L / 2) {
      return 32.0 * m;
    };
    return 0.0;
  };

  // trasform f_i into a matrix filled with values that for i,j are represented
  // as f(x,y)
  ComplexTensor f_i_mat = ComplexTensor::SMatrix(N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      f_i_mat(i, j) = f_i(ComplexTensor::FromData({i * a, j * a}));
    }
  }

  f_i_mat = f_i_mat / f_i_mat.Norm();

  /*
int n_frames = 100;
std::vector<ComplexTensor> times(n_frames);
times.push_back(f_i_mat);

ComplexTensor curr_func = ComplexTensor(f_i_mat);
for (double t = 0; t < 0.2; t += tau) {
  auto new_value = EvolveState(curr_func, V_func, a_vec, m, tau);

  curr_func = new_value;
}

for (int i = 0; i < f_i_mat.Rows(); i++) {
  for (int j = 0; j < f_i_mat.Cols(); j++) {
    std::cout << i * a << " " << j * a << " " << std::norm(f_i_mat(i, j))
              << " " << std::norm(curr_func(i, j)) << "\n";
  }
}
  */
  int n = 2000;
  std::vector<ComplexTensor> frames;

  // Initial state.
  ComplexTensor curr_func = ComplexTensor(f_i_mat);

  // Evolve the state for `n` frames.
  for (int frame = 0; frame < n; ++frame) {
    frames.push_back(curr_func);
    curr_func = EvolveState(curr_func, V_func, a_vec, m, tau);
  }

  // Print all frames side by side for each grid point.
  for (int i = 0; i < f_i_mat.Rows(); ++i) {
    for (int j = 0; j < f_i_mat.Cols(); ++j) {
      // Print position for the current grid point.
      std::cout << i * a << " " << j * a;

      // Print norm of the initial and evolved states for all frames.
      for (const auto& frame : frames) {
        std::cout << " " << std::norm(frame(i, j));
      }

      // End the line for this grid point.
      std::cout << "\n";
    }
  }
}
