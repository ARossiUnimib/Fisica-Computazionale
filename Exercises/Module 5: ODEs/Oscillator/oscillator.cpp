#include "../../Module 3: Matrices/tensor_utils.hpp"
#include "../../utils.hpp"
#include "../odes.hpp"
#include <cmath>

Tensor<double> Oscillator(double t, const Tensor<double> &y)
{
    auto dydt = TensorBuilder<double>::Vector(2).Build();
    // real armonic oscillator
    dydt(0) = y(1);
    dydt(1) = -y(0);
    return dydt;
}

double ExactSolution(double t)
{
    return sin(t);
}

#define PT_1 1
#define PT_2 2

#define EULER 1
#define MIDPOINT 2
#define RK4 3

using MethodFunc =
    std::function<Tensor<double>(const Tensor<double>, const Utils::Range<double> &, ODEUtils::ODEFunction<double>)>;

void PrintError(const char *name)
{
    std::cerr << "Usage: " << name << " <es_point (1, 2)> <ode_method (1: Euler, 2: Midpoint, 3: RK4)>\n";
}

int main(int argc, char const *argv[])
{
    int point = argc < 2 ? -1 : std::stoi(argv[1]);
    int method = argc < 3 ? -1 : std::stoi(argv[2]);

    if (argc < 3 || point > 2 || point <= 0 || method > 3 || method <= 0)
    {
        PrintError(argv[0]);
        return 1;
    }

    // 1.0 1.0 1.0
    auto initial_tensor = TensorBuilder<double>::Vector(2).Build();
    initial_tensor(0) = 0.0;
    initial_tensor(1) = 1.0;

    auto time_range = Utils::Range<double>::Fixed(0.0, 100.0, 0.001);

    using namespace ODEUtils;

    ODEFunction<double> ode_func = Oscillator;
    MethodFunc method_func;

    Tensor<double> results;

    switch (method)
    {
    case EULER: {
        method_func = Euler<double>;
        break;
    }
    case MIDPOINT: {
        method_func = Midpoint<double>;
        break;
    }
    }

    if (point == PT_1)
    {
        results = method_func(initial_tensor, time_range, ode_func);
        Print(results, time_range.Start, time_range.Step);
    }

    if (point == PT_2)
    {
        auto h_range = Utils::Range<double>::FixedNum(0.001, 1, 1000);

        Utils::Range<double> time_range;

        for (const auto h : h_range)
        {
            time_range = Utils::Range<double>::Fixed(0.0, 10, h);

            auto max_h = time_range.nodes.size() - 1;

            auto mat = method_func(initial_tensor, time_range, ode_func);

            auto value = mat(max_h, 0) - ExactSolution(time_range.nodes[max_h]);

            std::cout << h << " " << std::abs(value) << std::endl;
        }
    }

    return 0;
}
