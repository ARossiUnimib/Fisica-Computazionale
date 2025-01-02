#pragma once

#include "../m_2_matrices/tensor.hpp"
#include "../m_2_matrices/tensor_utils.hpp"
#include "../utils.hpp"
#include "function_data.hpp"

template <typename T>
class Spline {
 public:
  Spline(func::FunctionData<T> func_data, int order)
      : func_data_(func_data),
        func_values_(std::vector<T>()),
        order_(order),
        matrix_(tensor::Tensor<T>::SMatrix((order+1) * (func_data_.Size() - 1))) {
    FillMatrix();
    GenerateCoefficients();
  }

  T operator()(T x) const {
    // get index for the specified x in the interval
    int i = 0;

    // value should be between func_data.X(i) and func_data_(X(i+1))
    while (x >= func_data_.X(i + 1) && i < func_data_.Size()) {
      i++;
    }

    assert(i < func_data_.Size());


    // TODO: fix last element of the spline


    // NOTE: coeffs are stored as a_1 b_1 c_1 a_2 b_2 ...
    // linear spline
    switch (order_) {
      case 1: {
        T a = coeffs_(2 * i);
        T b = coeffs_(2 * i + 1);
        return a + b * x;
      }
      // quadratic spline
      case 2: {
        T a = coeffs_(3 * i);
        T b = coeffs_(3 * i + 1);
        T c = coeffs_(3 * i + 2);
        return a + b * x + c * x * x;
      }
      default:
        LOG_ERROR("Order not implemented");
    }
  }

 private:
  void GenerateCoefficients() {
    auto func_tensor = tensor::Tensor<T>::FromData(func_values_);
    coeffs_ = tensor::GaussianElimination(matrix_, func_tensor);
  }

  // Calculate the matrix for the coefficients computation
  void FillMatrix() {
    int n = func_data_.Size();
    switch (order_) {
      // LINEAR
      case 1: {
        // block number
        int i = 0;
        // coordinates number (x1, x2, f1, f2...)
        int j = 0;
        while (i < 2 * (n - 1)) {
          matrix_(i, i) = matrix_(i + 1, i) = 1;

          // x1 x2 x2 x3 x3 ...
          matrix_(i, i + 1) = func_data_.X(j);
          matrix_(i + 1, i + 1) = func_data_.X(j + 1);

          // f1 f2 f2 f3 f3 ....

          func_values_.push_back(func_data_.F(j));
          func_values_.push_back(func_data_.F(j + 1));

          i += 2;
          j++;
        }

        break;
      }

      // QUADRATIC
      case 2: {
        int i = 0;
        int j = 0;

        while (i < 3 * (n - 1)) {
          matrix_(i, i) = matrix_(i + 1, i) = matrix_(i + 2, i + 1) = 1;

          // x_j and x_{j+1}
          auto x_1 = func_data_.X(j);
          auto x_2 = func_data_.X(j + 1);

          matrix_(i, i + 1) = x_1;
          matrix_(i, i + 2) = x_1 * x_1;
          matrix_(i + 1, i + 1) = x_2;
          matrix_(i + 1, i + 2) = x_2 * x_2;

          matrix_(i + 2, i + 2) = 2 * x_2;

          func_values_.push_back(func_data_.F(j));
          func_values_.push_back(func_data_.F(j + 1));
          func_values_.push_back(0);

          if (i < 3 * (n - 2)) {
            matrix_(i + 2, i + 4) = -1;
            matrix_(i + 2, i + 5) = -2 * x_2;
          }

          i += 3;
          j++;
        }

        break;
      }
      default:
        LOG_ERROR("Order not implemented");
    }
  }

  // Func values rearranged to fit the vector
  // for the matrix
  std::vector<T> func_values_;
  // Function values
  func::FunctionData<T> func_data_;
  // Order of spline interpolation
  int order_;
  // Matrix to calculate coefficients
  tensor::Tensor<T> matrix_;
  tensor::Tensor<T> coeffs_;
};
