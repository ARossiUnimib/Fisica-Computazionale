#pragma once

#include <complex>
#include <functional>

#include "../../m_2_matrices/tensor.hpp"

using Complex = std::complex<double>;
using ComplexTensor = tensor::Tensor<Complex>;
using Coords = std::vector<std::pair<int, int>>;
using namespace std::complex_literals;

ComplexTensor GenerateKappa(Complex m, ComplexTensor const& a) {
  auto k = ComplexTensor::Vector(a.Rows());

  for (int i = 0; i < k.Rows(); i++) {
    k(i) = 1.0 / (2.0 * m * m * a(i) * a(i));
  }

  return k;
}

// NOTE: We do not pass t since for the exercises V is time-independant
// we could pass V at the current timestep and it would be the same
ComplexTensor GenerateLambda(Complex m, ComplexTensor const& k,
                             ComplexTensor const& V_func) {
  const int kN = V_func.Rows();
  Complex sum = 0;

  for (int i = 0; i < k.Rows(); i++) {
    sum += k(i);
  }

  ComplexTensor result = ComplexTensor::SMatrix(kN);
  for (int i = 0; i < kN; i++) {
    for (int j = 0; j < kN; j++) {
      result(i, j) = 2.0 * sum + V_func(i, j) / m;
    }
  }

  return result;
}

// NOTE: Assuming same sample rate on all axis
ComplexTensor GenerateHamiltonian(ComplexTensor const& phi,
                                  ComplexTensor const& V_func,
                                  ComplexTensor const& a, Complex m) {
  ComplexTensor k = GenerateKappa(m, a);
  ComplexTensor lambda = GenerateLambda(m, k, V_func);

  auto H = ComplexTensor::SMatrix(phi.Rows());
  const int kN = H.Rows();

  for (int i = 0; i < kN; i++) {
    for (int j = 0; j < kN; j++) {
      Complex sum = 0;
      for (int mu = 0; mu < k.Rows(); mu++) {
        // TODO: recheck correctness
        if (i - 1 >= 0 && i + 1 < kN) {
          sum -= k(mu) * (phi(i - 1, j) + phi(i + 1, j));
        }
        if (j - 1 >= 0 && j + 1 < kN) {
          sum -= k(mu) * (phi(i, j - 1) + phi(i, j + 1));
        }
      }
      H(i, j) = sum + lambda(i, j) * phi(i, j);
    }
  }

  return H;
}

ComplexTensor EvolveState(ComplexTensor const& phi, ComplexTensor const& V_func,
                          ComplexTensor const& a, Complex m, Complex t) {
  auto H_func = [&V_func, &m, &a](ComplexTensor const& phi) {
    return GenerateHamiltonian(phi, V_func, a, m);
  };

  ComplexTensor phi_1 = H_func(phi);

  ComplexTensor phi_2 = phi_1 + H_func(phi_1) * t / 2.0;
  ComplexTensor phi_3 = phi_1 + phi_2 * t / 2.0;
  ComplexTensor phi_4 = phi_1 + phi_3 * t;

  Complex c_2 = 2.0;
  return phi + t / 6.0 * (phi_1 + c_2 * phi_2 + c_2 * phi_3 + phi_4);
}
