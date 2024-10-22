#pragma once

#include "../Module 3: Matrices/tensor_utils.hpp"
#include "function_data.hpp"

namespace Interpolation
{

/**
 * @brief Compute the Vandermonde matrix
 *
 * @tparam T
 * @param values
 * @return Tensor<T>
 */
template <typename T> Tensor<T> VandermondeMatrix(const std::vector<T> &values);

/**
 * @brief Retrieve polynomial coefficients from direct interpolation using
 * Vandermonde matrix
 *
 * @tparam T
 * @param f
 * @return std::vector<T> coefficients of size f.Size()
 */
template <typename T>
std::vector<T> DirectCoefficients(const FunctionData<T> &f);

/**
 * @brief Direct interpolation using Vandermonde matrix
 *
 * @tparam T
 * @param f
 * @param range Values to be extrapolated from polynomial
 * @return FunctionData<T> polynomial data
 */
template <typename T>
FunctionData<T> DirectPolynomial(const FunctionData<T> &f,
                                 const Utils::Range<T> &range);

/**
 * @brief Compute the Newton coefficients of the function
 * @note use macro _SLOW_ALGORITHM to use the algorithm given during the course
 *
 * @tparam T
 * @param f
 * @return std::vector<T>
 */
template <typename T>
std::vector<T> NewtonCoefficients(const FunctionData<T> &f);

/**
 * @brief Compute the Newton polynomial of the function
 *
 * @tparam T
 * @param f
 * @param range
 * @return FunctionData<T>
 */
template <typename T>
FunctionData<T> NewtonPolynomial(const FunctionData<T> &f,
                                 const Utils::Range<T> &range);

} // namespace Interpolation

#include "interpolation.inl"
