#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../../utils.hpp"
#include "../function_data.hpp"
#include "../interpolation.hpp"

using Func = func::FunctionData<double>;

void DataPrint(Func const &new_poly, Func const &inv_poly) {
    for (int i = 0; i < new_poly.Size(); i++) {
        std::cout << new_poly.X(i) << "\t" << new_poly.F(i) << "\t"
                  << inv_poly.F(i) << std::endl;
    }
}

void CalculateInterpolation(Func &&f, int start, int end) {
    // Select the range of data for the selected exercise
    auto _x = utils::slice(f.X(), start, end);
    auto _f = utils::slice(f.F(), start, end);

    f = {_x, _f};

    auto range = func::Range<double>::Fixed(0, 30, 0.25);

    Func newton_poly = interp::NewtonPolynomial(f, range);
    Func direct_poly = interp::DirectPolynomial(f, range);

    DataPrint(newton_poly, direct_poly);
}

int main(int argc, char const *argv[]) {
    auto es_num = argc != 2 ? -1 : std::stoi(argv[1]);

    if (es_num == -1) {
        std::cerr << "Usage: " << argv[0] << " <es_num>" << std::endl;
        return 1;
    }

    auto x = std::vector<double>{0, 10, 15, 20, 22.5, 30};
    auto f = std::vector<double>{0, 227.04, 362.78, 517.35, 602.97, 901.67};

    Func function = {x, f};

    std::array<double, 2> range;

    switch (es_num) {
        case 1:
            range = {2, 3};
            break;
        case 2:
            range = {1, 3};
            break;
        case 3:
            range = {0, 4};
            break;
        default:
            std::cerr << "Invalid exercise number" << std::endl;
            return 1;
    }

    CalculateInterpolation(std::move(function), range[0], range[1]);

    return 0;
}
