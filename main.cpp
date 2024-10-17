#include <cstdint>
#include <iostream>
#include <cmath>

#include "grid.h"

double func(double const x) {
    return 2 * pow(x, 2) + 0.1 * x + 3;
}
double func2D(double const x, double const y) {
    return -5 * pow(x, 2) * y + 2 * x * y + 10;
}
// Gaussian quadrature for 1 param functions
double quadGauss1D(double (*func)(double), double const start, double const end) {
    double result = 0.0;
    constexpr uint32_t steps = 1.e3;
    constexpr double x[2] = {-0.57735,0.57735}; // Points
    constexpr double w[2] = {1, 1}; // Weights

    const double dx = (end - start) / steps; // Size of 1 step

    for (int i = 0; i < steps; i++) {
        const double x1 = i * dx + start; // Current start
        const double x2 = (i + 1) * dx + start; // Current end

        for(int j = 0; j < 2; j++) {
            result += w[j] * func((x2 - x1) / 2 * x[j] + ((x1 + x2) / 2));
        }
    }
    return result * dx / 2;
}
// Gaussian quadrature for 2 param functions
double quadGauss2D(double (*func)(double, double), double const xStart, double const xEnd, double const yStart, double const yEnd) {
    double result = 0.0;
    constexpr uint32_t steps = 1.e3;
    constexpr double x[2] = {-0.57735,0.57735}; // Points
    constexpr double w[2] = {1, 1}; // Weights

    const double dx = (xEnd - xStart) / steps; // Size of 1 step
    const double dy = (yEnd - yStart) / steps;

    for (int ix = 0; ix < steps; ix++) {
        const double x1 = ix * dx + xStart; // Current start
        const double x2 = (ix + 1) * dx + xStart; // Current end
        for(int iy = 0; iy < steps; iy++) {
            const double y1 = (iy * dy + yStart);
            const double y2 = (iy + 1) * dy + yStart;

            for(int j = 0; j < 2; j++) {
                result += w[j] * func((x2 - x1) / 2 * x[j] + ((x1 + x2) / 2), (y1 - y2) / 2 * x[j] + ((y1 + y2) / 2));
            }
        }

    }
    return result * dx * dy / 2;
}

int main() {
    std::cout << "Grid creation:" << std::endl;
    GlobalData gd;
    gd.getData("../gridFiles/data.txt");
    gd.createGrid();
    gd.printData();
    gd.printGridNodes();
    gd.printGridElems();

    std::cout << "Quadratures:" << std::endl;
    std::cout << quadGauss1D(&func,-1.0,1.0) << std::endl;
    std::cout << quadGauss2D(&func2D,-1.0,1.0,-1.0, 1.0) << std::endl;

    return 0;
}
