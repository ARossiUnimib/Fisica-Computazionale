#pragma once

#include <cassert>
#include <memory>

#include "../m_3_matrices/tensor.hpp"
#include "odes.hpp"

enum class Method { kNone = 0, kEuler = 1, kMidpoint = 2, kRK4 = 3 };

template <typename T>
class ODESolverBuilder;

template <typename T>
class ODESolver {
 private:
  template <typename U>
  friend class ODESolverBuilder;

 public:
  static ODESolverBuilder<T> Builder();

  /**
   * @brief Solves the system
   *
   * Returns the solution
   */
  tensor::Tensor<T> Solve();

 private:
  // Steps of different type

  tensor::Tensor<T> EulerStep(T coord_step, tensor::Tensor<T> prev_condition);
  tensor::Tensor<T> MidpointStep(T coord_step,
                                 tensor::Tensor<T> prev_condition);
  tensor::Tensor<T> RK4Step(T coord_step, tensor::Tensor<T> prev_condition);

 private:
  ODESolver<T>(Method method, tensor::Tensor<T> initial_conds,
               func::Range<T> coords_range, ode::Function<T> system_func)
      : method_(method),
        initial_conds_(initial_conds),
        coords_range_(coords_range),
        system_func_(system_func) {}

  ODESolver<T>() = default;

 private:
  Method method_;
  tensor::Tensor<T> initial_conds_;
  func::Range<T> coords_range_;
  ode::Function<T> system_func_;
};

template <typename T>
class ODESolverBuilder {
 public:
  ODESolverBuilder<T>() : temp_() {}

  ODESolverBuilder<T> Method(Method method) {
    temp_.method_ = method;
    return *this;
  }

  ODESolverBuilder<T> InitialConditions(
      tensor::Tensor<T> const &initial_conds) {
    temp_.initial_conds_ = initial_conds;
    return *this;
  }

  ODESolverBuilder<T> CoordinatesRange(func::Range<T> const &coords_range) {
    temp_.coords_range_ = coords_range;
    return *this;
  }

  ODESolverBuilder<T> SystemFunction(ode::Function<T> const &system_func) {
    temp_.system_func_ = system_func;
    return *this;
  }

  ODESolver<T> Build() {
    // Sanity check for the builder
    assert(temp_.method_ != Method::kNone);
    assert(temp_.initial_conds_.Rows() > 0);
    assert(temp_.coords_range_.Nodes().size() > 0);
    assert(temp_.system_func_ != nullptr);

    return temp_;
  }

  // For large objects a pointer should be prioritized
  std::unique_ptr<ODESolver<T>> BuildPtr() {
    return std::make_unique<ODESolver<T>>(temp_);
  }

 private:
  ODESolver<T> temp_;
};

/* -------------- IMPLEMENTATION -------------------------------------------- */

template <typename T>
ODESolverBuilder<T> ODESolver<T>::Builder() {
  return ODESolverBuilder<T>();
}

template <typename T>
tensor::Tensor<T> ODESolver<T>::Solve() {
  using namespace tensor;

  Tensor<T> y = initial_conds_;
  Tensor<T> result =
      Tensor<T>::Matrix(coords_range_.Nodes().size(), initial_conds_.Rows());

  int step_index = 0;
  for (const auto &t : coords_range_) {
    // Store the current state in the result tensor
    for (int i = 0; i < y.Rows(); ++i) {
      result(step_index, i) = y(i);
    }

    // FIX: don't check every loop
    switch (method_) {
      case Method::kEuler: {
        y = y + EulerStep(t, y) * coords_range_.Step();
        break;
      }

      case Method::kMidpoint: {
        y = y + MidpointStep(t, y) * coords_range_.Step();
        break;
      }

      case Method::kRK4: {
        y = y + RK4Step(t, y) * coords_range_.Step();
        break;
      }

      default:
        assert(0);
    }

    step_index++;
  }

  return result;
}

template <typename T>
tensor::Tensor<T> ODESolver<T>::EulerStep(T coord_step,
                                          tensor::Tensor<T> prev_condition) {
  return system_func_(coord_step, prev_condition);
}

template <typename T>
tensor::Tensor<T> ODESolver<T>::MidpointStep(T coord_step,
                                             tensor::Tensor<T> prev_condition) {
  tensor::Tensor<T> k_1 = system_func_(coord_step, prev_condition);
  tensor::Tensor<T> k_2 =
      system_func_(coord_step + coords_range_.Step() / 2,
                   prev_condition + k_1 * (coords_range_.Step() / 2));

  return k_2;
}

template <typename T>
tensor::Tensor<T> ODESolver<T>::RK4Step(T coord_step,
                                        tensor::Tensor<T> prev_condition) {
  tensor::Tensor<T> k_1 = system_func_(coord_step, prev_condition);
  tensor::Tensor<T> k_2 =
      system_func_(coord_step + coords_range_.Step() / 2,
                   prev_condition + k_1 * (coords_range_.Step() / 2));

  tensor::Tensor<T> k_3 =
      system_func_(coord_step + coords_range_.Step() / 2,
                   prev_condition + k_2 * (coords_range_.Step() / 2));

  tensor::Tensor<T> k_4 =
      system_func_(coord_step + coords_range_.Step(),
                   prev_condition + k_3 * coords_range_.Step());

  return (k_1 + k_2 * 2.0 + k_3 * 2.0 + k_4) * (1.0 / 6.0);
}
