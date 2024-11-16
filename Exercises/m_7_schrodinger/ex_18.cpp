#include <cmath>
#include <complex>

#include "../m_2_matrices/tensor.hpp"

using Complex = std::complex<double>;
using ComplexTensor = tensor::Tensor<Complex>;
using namespace std::complex_literals;

Complex Phi0(ComplexTensor x, ComplexTensor y, double sigma_2,
             ComplexTensor p) {

  return std::exp(-(x - y).NormSquared() / sigma_2 + 1i * (p.Dot(x))(0));
}

int main() {
    // N1=N2=N
    int N = 64;
    double L = 1;
    double a = 1;

    auto y = ComplexTensor::Vector(2);
    y(0) = L/4;
    y(1) = L/2;
    double asigma = 4;
    double p0L = (2*M_PI)*6;
    double mL = 8;
    double t = 0.001;
}
