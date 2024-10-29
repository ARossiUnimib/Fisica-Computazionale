#pragma once

#include <cassert>
#include <functional>
#include <memory>

#include "../m_3_matrices/tensor.hpp"
#include "../m_4_interpolation/range.hpp"
#include "odes.hpp"

namespace ode {

enum class Method { kNone = 0, kEuler = 1, kMidpoint = 2, kRK4 = 3 };

template <typename T>
class ODESolverBuilder;

template <typename T>
class ODESolver {
 private:
  // The builder can modify the private members of the class directly
  template <typename U>
  friend class ODESolverBuilder;

 public:
  /**
   * @brief Creates a builder for the ODESolver
   *
   * Instead of using a constructor a builder pattern is used to allow for a
   * more modular way of creating the object
   */
  static ODESolverBuilder<T> Builder();

  /**
   * @brief Solves the ode system
   *
   * Returns the solution
   */
  tensor::Tensor<T> Solve();

 private:
  /**
   * @brief Returns the function that implements the integrator method
   *
   * @param method The method to be used
   */
  ode::Function<T> GetMethodFunction(Method method);

  // Steps of different type

  tensor::Tensor<T> EulerStep(T coord_step,
                              tensor::Tensor<T> const &prev_condition);

  tensor::Tensor<T> MidpointStep(T coord_step,
                                 tensor::Tensor<T> const &prev_condition);

  tensor::Tensor<T> RK4Step(T coord_step,
                            tensor::Tensor<T> const &prev_condition);

 private:
  ODESolver(Method method, tensor::Tensor<T> initial_conds,
            func::Range<T> coords_range, ode::Function<T> system_func)
      : method_(method),
        initial_conds_(initial_conds),
        coords_range_(coords_range),
        system_func_(system_func) {}

  // ODESolver<T> &operator=(ODESolver<T> const &) = delete;
  // ODESolver<T> &operator=(ODESolver<T> &&) = delete;

  ODESolver()
      : method_(Method::kNone),
        initial_conds_(tensor::Tensor<T>::Vector(0)),
        coords_range_(func::Range<T>::Fixed(0, 0)),
        system_func_(ode::Function<T>{}) {}

 private:
  Method method_;
  tensor::Tensor<T> initial_conds_;
  func::Range<T> coords_range_;
  ode::Function<T> system_func_;
};

template <typename T>
class ODESolverBuilder {
  // HACK: This is a workaround to allow the empty constructor to be private
  struct empty_constr_ODESolver : public ODESolver<T> {};

 public:
  ODESolverBuilder()
      : temp_(std::move(std::make_shared<empty_constr_ODESolver>())) {}

  ODESolverBuilder<T> Method(Method method) {
    temp_->method_ = method;
    return *this;
  }

  ODESolverBuilder<T> InitialConditions(
      tensor::Tensor<T> const &initial_conds) {
    temp_->initial_conds_ = initial_conds;
    return *this;
  }

  ODESolverBuilder<T> CoordinatesRange(func::Range<T> const &coords_range) {
    temp_->coords_range_ = coords_range;
    return *this;
  }

  ODESolverBuilder<T> SystemFunction(ode::Function<T> const &system_func) {
    temp_->system_func_ = system_func;
    return *this;
  }

  std::unique_ptr<ODESolver<T>> Build() {
    // Sanity check for the builder
    assert(temp_->method_ != Method::kNone);
    assert(temp_->initial_conds_.Rows() > 0);

    assert(temp_->coords_range_.Nodes().size() > 0 &&

           // NOTE: for now we only support fixed step size, probably forever
           temp_->coords_range_.NodesType() ==
               func::Range<T>::NodeType::kFixed);

    assert(temp_->system_func_ != nullptr);

    return std::make_unique<ODESolver<T>>(std::move(*temp_));
  }

 private:
  // Since the builder can have multiple parts being created at different times,
  // we rely on shared pointers to keep track of the object state
  std::shared_ptr<ODESolver<T>> temp_;
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

  auto method_func = GetMethodFunction(method_);

  for (const auto &t : coords_range_) {
    // Store the current state in the result tensor
    for (int i = 0; i < y.Rows(); ++i) {
      result(step_index, i) = y(i);
    }

    // update the system through the chosen method
    y = y + method_func(t, y) * coords_range_.Step();

    step_index++;
  }

  return result;
}

// Binds class pointer "this", and the two arguments of the method to be usable
// as a member class function pointer
#define BIND_METHOD_FUNC(method)                                \
  std::bind(&ODESolver<T>::method, this, std::placeholders::_1, \
            std::placeholders::_2)

template <typename T>
ode::Function<T> ODESolver<T>::GetMethodFunction(Method method) {
  switch (method) {
    case Method::kEuler:
      return BIND_METHOD_FUNC(EulerStep);
    case Method::kMidpoint:
      return BIND_METHOD_FUNC(MidpointStep);
    case Method::kRK4:
      return BIND_METHOD_FUNC(RK4Step);
    default:
      assert(0);
  }
}

#undef BIND_METHOD_FUNC

template <typename T>
tensor::Tensor<T> ODESolver<T>::EulerStep(
    T coord_step, tensor::Tensor<T> const &prev_condition) {
  return system_func_(coord_step, prev_condition);
}

template <typename T>
tensor::Tensor<T> ODESolver<T>::MidpointStep(
    T coord_step, tensor::Tensor<T> const &prev_condition) {
  tensor::Tensor<T> k_1 = system_func_(coord_step, prev_condition);
  tensor::Tensor<T> k_2 =
      system_func_(coord_step + coords_range_.Step() / 2,
                   prev_condition + k_1 * (coords_range_.Step() / 2));

  return k_2;
}

template <typename T>
tensor::Tensor<T> ODESolver<T>::RK4Step(
    T coord_step, tensor::Tensor<T> const &prev_condition) {
  tensor::Tensor<T> k_1 = system_func_(coord_step, prev_condition);

  T h = coords_range_.Step();

  tensor::Tensor<T> k_2 =
      system_func_(coord_step + h / 2, prev_condition + k_1 * (h / 2));

  tensor::Tensor<T> k_3 =
      system_func_(coord_step + h / 2, prev_condition + k_2 * (h / 2));

  tensor::Tensor<T> k_4 =
      system_func_(coord_step + h, prev_condition + k_3 * h);

  return (k_1 + k_2 * 2.0 + k_3 * 2.0 + k_4) * (1.0 / 6.0);
}

}  // namespace ode
