#include <cmath>

#include "time_evolution.hpp"

Complex Phi0(ComplexTensor const& x, ComplexTensor const& y, double sigma_2,
             ComplexTensor const& p) {
    // time normalization factor
return std::exp(-(x - y).NormSquared() / sigma_2 + 1i * (p.Dot(x))(0)) * 1.0 / (2 * M_PI * sigma_2);
}

int main() {
  // N1=N2=N
  int N = 64;
  double L = 1;

  double a = L / N;

  auto a_vec = ComplexTensor::One(2);
  double m = 8 / L;
  double tau = 0.01;
  double dt = tau / m;

  auto f_i = [&L, &m](ComplexTensor const& x) {
    auto p0 = ComplexTensor::FromData({(2 * M_PI) * 6 / L, 0});
    auto y = ComplexTensor::Vector(2);
    y(0) = L / 4;
    y(1) = L / 2;
    double sigma = 0.5 / m;

    return Phi0(x, y, sigma * sigma, p0);
  };

  // TODO: maybe use a function instead
  ComplexTensor V_mat = ComplexTensor::SMatrix(N);

  // trasform f_i into a matrix filled with values that for i,j are represented
  // as f(x,y)
  ComplexTensor f_i_mat = ComplexTensor::SMatrix(N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      f_i_mat(i, j) = f_i(ComplexTensor::FromData({i * a, j * a}));
    }
  }

  ComplexTensor curr_func = ComplexTensor(f_i_mat);
  for (double t = 0; t < tau * 400; t += tau) {
    auto new_value = EvolveState(curr_func, V_mat, a_vec, m, tau);
    curr_func = new_value;
  }

  for (int i = 0; i < f_i_mat.Rows(); i++) {
    for (int j = 0; j < f_i_mat.Cols(); j++) {
      std::cout << i * a << " " << j * a << " "
                << std::abs(f_i_mat(i, j)) * std::abs(f_i_mat(i, j)) << " "
                << std::abs(curr_func(i, j)) * std::abs(curr_func(i, j)) << "\n";
    }
  }
}
