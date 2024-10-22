#include "../../utils.hpp"
#include "../odes.hpp"

Tensor<double> LorenzSystem(double t, const Tensor<double> &y)
{
    auto dydt = Tensor<double>::Vector(3);

    dydt(0) = 10 * (y(1) - y(0));
    dydt(1) = -y(0) * y(2) + 28 * y(0) - y(1);
    dydt(2) = y(0) * y(1) - (8.0 / 3.0) * y(2);

    return dydt;
}
void PrintUsage(const char *name)
{
    std::cerr << "Usage: " << name
              << " <ode_method (1: Euler, 2: RK2, 3: RK4) \n";
}

int main(int argc, char const *argv[])
{
    if (argc < 2)
    {
        PrintUsage(argv[0]);
        return 1;
    }

    // 1.0 1.0 1.0
    auto initial_tensor = Tensor<double>::Ones(3);

    auto time_range = Utils::Range<double>::Fixed(0.0, 10.0, 0.01);

    {
        using namespace ODEUtils;

        ODEFunction<double> ode_func = LorenzSystem;

        Tensor<double> result;

        switch (std::stoi(argv[1]))
        {
        case 1:
            result = Euler(initial_tensor, time_range, ode_func);
            break;
        case 2:
            result = Midpoint(initial_tensor, time_range, ode_func);
            break;
        case 3:
            // TODO: Implement RK4
            break;
        }

        Print(result, time_range.Start, time_range.Step);
    }

    return 0;
}
