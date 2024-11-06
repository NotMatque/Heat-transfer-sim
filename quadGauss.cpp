#include "quadGauss.h"

double quadGauss1D(double (*func)(double), double const start, double const end) {
    double result = 0.0;
    constexpr unsigned int steps = 1.e3;
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

double quadGauss2(double (*func)(double, double), double const xStart, double const xEnd, double const yStart, double const yEnd) {
    double result = 0.0;
    //constexpr uint32_t steps = 1.e3;
    constexpr double x[2] = {-0.57735,0.57735}; // Points
    constexpr double w[2] = {1, 1}; // Weights

    const double dx = (xEnd - xStart) / 2; // Size of 1 step
    const double dy = (yEnd - yStart) / 2;

    for (int ix = 0; ix < 2; ix++) {
        const double x1 = ix * dx + xStart; // Current start
        const double x2 = (ix + 1) * dx + xStart; // Current end
        for(int iy = 0; iy < 2; iy++) {
            const double y1 = (iy * dy + yStart);
            const double y2 = (iy + 1) * dy + yStart;

            for(int j = 0; j < 2; j++) {
                result += w[j] * func((x2 - x1) / 2 * x[j] + ((x1 + x2) / 2), (y1 - y2) / 2 * x[j] + ((y1 + y2) / 2));
            }
        }

    }
    return result * dx * dy / 2;
}

double quadGauss3(double (*func)(double, double), double xStart, double xEnd, double yStart, double yEnd) {
    double result = 0.0;
    //constexpr uint32_t steps = 1.e3;
    constexpr double x[3] = {-0.77459, 0, 0.77459}; // Points
    constexpr double w[3] = {0.55556, 0.88889, 0.55556}; // Weights

    const double dx = (xEnd - xStart) / 3; // Size of 1 step
    const double dy = (yEnd - yStart) / 3;

    for (int ix = 0; ix < 3; ix++) {
        const double x1 = ix * dx + xStart; // Current start
        const double x2 = (ix + 1) * dx + xStart; // Current end
        for(int iy = 0; iy < 3; iy++) {
            const double y1 = (iy * dy + yStart);
            const double y2 = (iy + 1) * dy + yStart;

            for(int j = 0; j < 3; j++) {
                result += w[j] * func((x2 - x1) / 2 * x[j] + ((x1 + x2) / 2), (y1 - y2) / 2 * x[j] + ((y1 + y2) / 2));
            }
        }
    }
    return result * dx * dy / 2;
}

