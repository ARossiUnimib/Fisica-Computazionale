#include "../../utils.hpp"
#include "../function_data.hpp"
#include "../interpolation.hpp"

double runge_function(double x)
{
    return 1.0 / (1.0 + 25.0 * x * x);
}

// NOTE: direct poly becomes ill conditioned in chebeyshev points and it
// explodes faster than newton polys in fixed points due to worse handling of
// round off errors and stability

#ifdef CHEBYSHEV
#define FUNCTION Chebyshev
#else
#define FUNCTION FixedNum
#endif

int main()
{

    // Data sample for interpolation
    auto interp_sample = Utils::Range<double>::FUNCTION(-1.0, 1.0, 50);

    std::vector<double> runge_values;

    for (auto x : interp_sample)
    {

#ifdef PRINT_RUNGE
        std::cout << x << "\t" << runge_function(x) << std::endl;
        continue;
    }
#else

        runge_values.push_back(runge_function(x));
    }

    // Print runge sampled data

    FunctionData<double> runge_data = {interp_sample.nodes, runge_values};

    // Testing sample to visualize generated polynomial
    auto testing_sample = Utils::Range<double>::FUNCTION(-1.0, 1.0, 150);

    auto newton_interpolation = Interpolation::NewtonPolynomial(runge_data, testing_sample);
    auto direct_interpolation = Interpolation::DirectPolynomial(runge_data, testing_sample);

    // Print
    for (int i = 0; i < newton_interpolation.Size(); i++)
    {
        std::cout << newton_interpolation.X()[i] << "\t" << newton_interpolation.F()[i] << "\t"
                  << direct_interpolation.F()[i] << std::endl;
    }

#endif

    return 0;
}
