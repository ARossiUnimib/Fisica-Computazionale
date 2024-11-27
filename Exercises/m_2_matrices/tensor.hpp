#pragma once

#include "tensor_obj.hpp"
#include "tensor_utils.hpp"

namespace tensor {

// Action via a scalar
template <typename T>
Tensor<T> operator*(T const &d, Tensor<T> const &a) {
  return a * d;
}

}  // namespace tensor
