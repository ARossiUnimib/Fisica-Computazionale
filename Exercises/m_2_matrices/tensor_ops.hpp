#pragma once

#include "tensor_obj.hpp"

// Define commutative operations for Tensor class

namespace tensor {

// Action via a scalar
template <typename T>
Tensor<T> operator*(T const &d, Tensor<T> const &a) {
  return a * d;
}

}
