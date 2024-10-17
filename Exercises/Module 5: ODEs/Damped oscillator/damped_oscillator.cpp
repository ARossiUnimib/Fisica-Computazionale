#include <cmath>
#include <string>

#include "../../Module 3: Matrices/tensor_utils.hpp"
#include "../odes.hpp"

double g_Gamma = -1;
double g_A = -1;

Tensor<double> RealOscillator(double t, const Tensor<double> &y)
{
    auto dydt = TensorBuilder<double>::Vector(2).Build();

    // real armonic oscillator
    dydt(0) = y(1);
    dydt(1) = -sin(y(0));

    return dydt;
}

Tensor<double> DampedOscillator(double t, const Tensor<double> &y)
{
    auto dydt = TensorBuilder<double>::Vector(2).Build();

    // real armonic oscillator
    dydt(0) = y(1);
    dydt(1) = -sin(y(0)) - g_Gamma * y(1);

    return dydt;
}

Tensor<double> ForcedOscillator(double t, const Tensor<double> &y)
{
    auto dydt = TensorBuilder<double>::Vector(2).Build();

    // real armonic oscillator
    dydt(0) = y(1);
    dydt(1) = -sin(y(0)) - g_Gamma * y(1) + g_A * sin((2.0 / 3.0) * t);

    return dydt;
}

void PrintUsage(const char *name)
{
    std::cerr << "Usage: " << name << " <es_point (1, 2, 3)> <gamma (es 2, 3)> <A (es 3)> \n";
}

int main(int argc, const char *argv[])
{
    int point = argc < 2 ? -1 : std::stoi(argv[1]);
    double gamma = argc < 3 ? -1 : std::stod(argv[2]);
    double A = argc < 4 ? -1 : std::stod(argv[3]);

    if (argc < 2 || point > 3 || point <= 0)
    {
        PrintUsage(argv[0]);
        return 1;
    }

    auto initial_tensor = TensorBuilder<double>::Vector(2).Build();
    initial_tensor(0) = 0.0;
    initial_tensor(1) = 1.0;

    auto time_range = Utils::Range<double>::Fixed(0.0, 100.0, 0.001);

    using namespace ODEUtils;

    ODEFunction<double> ode_func;

    switch (point)
    {
    case 1:
        ode_func = RealOscillator;
        break;

    case 2:

        if (gamma >= 2 || gamma <= 0)
        {
            PrintUsage(argv[0]);
            return 1;
        }

        g_Gamma = gamma;

        ode_func = DampedOscillator;

        break;

    case 3:

        if (gamma >= 2 || gamma <= 0 || A >= 2 || A <= 0)
        {
            PrintUsage(argv[0]);
            return 1;
        }

        g_Gamma = gamma;
        g_A = A;

        ode_func = ForcedOscillator;

        break;
    }

    Tensor<double> result = Midpoint(initial_tensor, time_range, ode_func);

    Print(result, time_range.Start, time_range.Step);
}
