#pragma once

#include <complex>
#include <functional>

#include "../../m_2_matrices/tensor.hpp"

using Complex = std::complex<double>;
using ComplexTensor = tensor::Tensor<Complex>;
using PotentialFunc = std::function<Complex(ComplexTensor const&)>;
using Coords = std::vector<std::pair<int, int>>;
using namespace std::complex_literals;

ComplexTensor GenerateKappa(Complex m, ComplexTensor const& a) {
  auto k = ComplexTensor::Vector(a.Rows());

  for (int i = 0; i < k.Rows(); i++) {
    k(i) = 1.0 / (2.0 * m * m * a(i) * a(i));
  }

  return k;
}



// NOTE: Assuming same sample rate on all axis
ComplexTensor GenerateHamiltonian(ComplexTensor const& phi,
                                  PotentialFunc const& V_func,
                                  ComplexTensor const& a, Complex m) {
  ComplexTensor k = GenerateKappa(m, a);

  auto H = ComplexTensor::SMatrix(phi.Rows());
  const int kN = H.Rows();

  for (int i = 0; i < kN; i++) {
    for (int j = 0; j < kN; j++) {
      Complex sum = 0;
      for (int mu = 0; mu < 2; mu++) {
        if (!(i - 1 < 0 || i + 1 == kN || j - 1 < 0 || j + 1 == kN)) {
          sum -= k(mu) * (phi(i + 1, j) + phi(i - 1, j) + phi(i, j + 1) +
                          phi(i, j - 1));
        }
      }

      Complex sum_k = 0;
      for (int i = 0; i < k.Rows(); i++) {
        sum_k += k(i);
      }

      auto x = ComplexTensor::FromData({(Complex)i * a(0), (Complex)j * a(1)});
      Complex lambda = 2.0 * sum_k + V_func(x) / m;
      H(i, j) = sum + lambda * phi(i, j);
    }
  }

  return H;
}

ComplexTensor EvolveState(ComplexTensor const& phi, PotentialFunc const& V_func,
                          ComplexTensor const& a, Complex m, Complex t) {
  auto H_func = [&V_func, &m, &a](ComplexTensor const& phi) {
    return GenerateHamiltonian(phi, V_func, a, m);
  };

  ComplexTensor phi_1 = 1i * H_func(phi);

  ComplexTensor phi_2 = phi_1 + 1i * H_func(phi_1) * t / 2.0;
  ComplexTensor phi_3 = phi_1 + phi_2 * t / 2.0;
  ComplexTensor phi_4 = phi_1 + phi_3 * t;

  Complex c_2 = 2.0;
  ComplexTensor new_phi =
      phi + t / 6.0 * (phi_1 + c_2 * phi_2 + c_2 * phi_3 + phi_4);

  return new_phi / new_phi.Norm();
}
