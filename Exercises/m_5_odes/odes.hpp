#pragma once

#include <functional>
#include <iostream>

#include "../m_3_matrices/tensor.hpp"
#include "../m_4_interpolation/range.hpp"

namespace ode {

/**
 * @brief function which describes the behaviour of the system of every first
 * order derivatives
 *
 * Given a vector of derivatives, Function describes the ode system completely
 * passing the initial value through a RK method
 *
 * @tparam T precision
 */
template <typename T>
using Function = std::function<tensor::Tensor<T>(T, tensor::Tensor<T> const &)>;

/**
 * @brief Prints the solution of the ODE system in the order given by
 * ode::Function
 */
template <typename T>
void Print(tensor::Tensor<T> const &results, T t0, T h);

/* ------------------- IMPLEMENTATION --------------------------------------- */

template <typename T>
void Print(tensor::Tensor<T> const &results, T t0, T h) {
  std::cout << "t";
  for (int i = 0; i < results.Cols(); ++i) {
    std::cout << "\ty" << (i + 1);
  }
  std::cout << std::endl;

  T t = t0;
  for (int row = 0; row < results.Rows(); ++row) {
    std::cout << t;
    for (int col = 0; col < results.Cols(); ++col) {
      std::cout << "\t" << results(row, col);
    }
    std::cout << std::endl;

    t += h;
  }
}

}  // namespace ode
