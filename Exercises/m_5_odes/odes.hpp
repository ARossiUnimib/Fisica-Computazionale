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

/**
 * @brief Euler method for solving ODE
 *
 * @tparam T precision
 * @param y0 initial conditions
 * @param time_range range of time values
 * @param system_func function which describes the behaviour of the system
 * @return tensor::Tensor<T> result of the ode system
 */
template <typename T>
tensor::Tensor<T> Euler(tensor::Tensor<T> const &y0, func::Range<T> &time_range,
                        ode::Function<T> const &system_func);

/**
 * @brief Midpoint method for solving ODE
 *
 * @tparam T precision
 * @param y0 initial conditions
 * @param time_range range of time values
 * @param system_func function which describes the behaviour of the system
 * @return tensor::Tensor<T> result of the ode system
 */
template <typename T>
tensor::Tensor<T> Midpoint(const tensor::Tensor<T> &y0,
                           func::Range<T> &time_range, Function<T> system_func);

/* ------------------------- IMPLEMENTATION --------------------------------- */

template <typename T>
tensor::Tensor<T> Euler(tensor::Tensor<T> const &y0, func::Range<T> &time_range,
                        ode::Function<T> const &system_func) {
  using namespace tensor;

  Tensor<T> y = y0;
  Tensor<T> result = Tensor<T>::Matrix(time_range.Nodes().size(), y0.Rows());

  int step_index = 0;
  for (const auto &t : time_range) {
    // Store the current state in the result tensor
    for (int i = 0; i < y.Rows(); ++i) {
      result(step_index, i) = y(i);
    }

    Tensor<T> dydt = system_func(t, y);

    for (int i = 0; i < y.Rows(); ++i) {
      y(i) += time_range.Step() * dydt(i);
    }

    step_index++;
  }

  return result;
}

template <typename T>
tensor::Tensor<T> Midpoint(const tensor::Tensor<T> &y0,
                           func::Range<T> &time_range,
                           Function<T> system_func) {
  using namespace tensor;

  Tensor<T> y = y0;
  Tensor<T> result = Tensor<T>::Matrix(time_range.Nodes().size(), y0.Rows());

  int step_index = 0;
  for (const auto &t : time_range) {
    // Store the current state in the result tensor
    for (int i = 0; i < y.Rows(); ++i) {
      result(step_index, i) = y(i);
    }

    Tensor<T> k_1 = system_func(t, y);
    Tensor<T> k_2 = system_func(t + time_range.Step() / 2,
                                y + k_1 * (time_range.Step() / 2));

    for (int i = 0; i < y.Rows(); ++i) {
      y(i) += time_range.Step() * k_2(i);
    }

    step_index++;
  }

  return result;
}

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
