#pragma once

#ifndef GAUSS_H
#define GAUSS_H

#include <cstdint>

// Gaussian quadrature for 1 param functions
double quadGauss1D(double (*func)(double), double start, double end);
// Gaussian quadrature for 2 param functions using 2 integration points
double quadGauss2(double (*func)(double, double), double xStart, double xEnd, double yStart, double yEnd);
// Gaussian quadrature for 2 param functions using 3 integration points
double quadGauss3(double (*func)(double, double), double xStart, double xEnd, double yStart, double yEnd);

#endif //GAUSS_H
