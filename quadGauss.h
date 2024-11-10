#pragma once

#ifndef GAUSS_H
#define GAUSS_H

#include <cstdint>

// Kwadratura Gaussa dla funkcji o 1 parametrze
double quadGauss1D(double (*func)(double), double start, double end);
// Kwadratura Gaussa o 2 punktach całkowania dla funkcji o 2 parametrach
double quadGauss2(double (*func)(double, double), double xStart, double xEnd, double yStart, double yEnd);
// Kwadratura Gaussa o 3 punktach całkowania dla funkcji o 2 parametrach
double quadGauss3(double (*func)(double, double), double xStart, double xEnd, double yStart, double yEnd);

#endif //GAUSS_H
