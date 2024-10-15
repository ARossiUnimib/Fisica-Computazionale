#include <cassert>
#include <cmath>
#include <vector>

#include "../../Module 3: Matrices/tensor_utils.hpp"
#include "../../utils.hpp"
#include "../function_data.hpp"
#include "../interpolation.hpp"

template <typename T> void DataPrint(const FunctionData<T> &new_poly, const FunctionData<T> &inv_poly)
{
    for (int i = 0; i < new_poly.Size(); i++)
    {
        std::cout << new_poly.X(i) << "\t" << new_poly.F(i) << "\t" << inv_poly.F(i) << std::endl;
    }
}

#define START 0
#define END 30
#define STEPS 0.25

void CalculateInterpolation(FunctionData<double>&& f, const Utils::Range<double> &slice_range)
{
    auto _x = Utils::slice(f.X(), slice_range);
    auto _f = Utils::slice(f.F(), slice_range);

    f = FunctionData<double>(_x, _f);

    auto range = Utils::Range<double>::Fixed(START, END, STEPS);

    FunctionData<double> newton_poly = Interpolation::NewtonPolynomial(f, range);

    FunctionData<double> direct_poly = Interpolation::DirectPolynomial(f, range);

    DataPrint(newton_poly, direct_poly);
}

int main()
{
    auto x = std::vector<double>{0, 10, 15, 20, 22.5, 30};
    auto f = std::vector<double>{0, 227.04, 362.78, 517.35, 602.97, 901.67};

    FunctionData<double> function(x, f);

    Utils::Range<double> slice_range;

#ifdef PART_1
    slice_range = {2, 3};
#elif PART_2
    slice_range = {1, 3};
#elif PART_3
    slice_range = {0, 4};
#else
    #error "Please define PART_[1-3]"
#endif

    CalculateInterpolation(std::move(function), slice_range);

    return 0;
}
