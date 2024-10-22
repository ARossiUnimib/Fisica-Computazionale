#pragma once

#include "../Module 3: Matrices/tensor.hpp"
#include "../utils.hpp"
#include <functional>
#include <iostream>

namespace ODEUtils
{

template <typename T>
using ODEFunction = std::function<Tensor<T>(T, const Tensor<T> &)>;

template <typename T>
Tensor<T> Euler(const Tensor<T> &y0, const Utils::Range<T> &time_range,
                ODEFunction<T> system_func)
{
    Tensor<T> y = y0;
    Tensor<T> result = Tensor<T>::Matrix(time_range.nodes.size(), y0.Rows());

    int step_index = 0;
    for (const auto &t : time_range)
    {
        // Store the current state in the result tensor
        for (int i = 0; i < y.Rows(); ++i)
        {
            result(step_index, i) = y(i);
        }

        Tensor<T> dydt = system_func(t, y);

        for (int i = 0; i < y.Rows(); ++i)
        {
            y(i) += time_range.Step * dydt(i);
        }

        step_index++;
    }

    return result;
}

template <typename T>
Tensor<T> Midpoint(const Tensor<T> &y0, const Utils::Range<T> &time_range,
                   ODEFunction<T> system_func)
{
    Tensor<T> y = y0;
    Tensor<T> result = Tensor<T>::Matrix(time_range.nodes.size(), y0.Rows());

    int step_index = 0;
    for (const auto &t : time_range)
    {
        // Store the current state in the result tensor
        for (int i = 0; i < y.Rows(); ++i)
        {
            result(step_index, i) = y(i);
        }

        Tensor<T> k_1 = system_func(t, y);
        Tensor<T> k_2 = system_func(t + time_range.Step / 2,
                                    y + k_1 * (time_range.Step / 2));

        for (int i = 0; i < y.Rows(); ++i)
        {
            y(i) += time_range.Step * k_2(i);
        }

        step_index++;
    }

    return result;
}

template <typename T> void Print(const Tensor<T> &results, T t0, T h)
{
    std::cout << "t";
    for (int i = 0; i < results.Cols(); ++i)
    {
        std::cout << "\ty" << (i + 1);
    }
    std::cout << std::endl;

    T t = t0;
    for (int row = 0; row < results.Rows(); ++row)
    {
        std::cout << t;
        for (int col = 0; col < results.Cols(); ++col)
        {
            std::cout << "\t" << results(row, col);
        }
        std::cout << std::endl;

        t += h;
    }
}

} // namespace ODEUtils
