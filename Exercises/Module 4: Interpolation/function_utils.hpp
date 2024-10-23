#pragma once

#include <iostream>

#include "function_data.hpp"
#include "range.hpp"

namespace func {

/**
 * @brief Create a polynomial function
 *
 * @tparam T
 * @param coefficients
 * @param x values
 * @return FunctionData<T>
 */
template <typename T>
FunctionData<T> Polynomial(std::vector<T> const &coefficients,
                           func::Range<T> const &x);

/**
 * @brief Print the given function
 *
 * @param f the function
 */
template <typename T>
void Print(func::FunctionData<T> const &f);

/* ---------------------------- IMPLEMENTATION ------------------------------ */

template <typename T>
FunctionData<T> Polynomial(std::vector<T> const &coefficients,
                           func::Range<T> const &x) {
    FunctionData<double> p;
    for (auto &&x : x) {
        double sum = 0;
        for (int j = 0; j < coefficients.size(); j++) {
            sum += coefficients[j] * pow(x, j);
        }
        p.Add(x, sum);
    }
    return p;
}

template <typename T>
void Print(func::FunctionData<T> const &f) {
    for (int i = 0; i < f.Size(); i++) {
        std::cout << f.X(i) << "\t" << f.F(i) << "\n";
    }

    std::cout << std::endl;
}

}  // namespace func
