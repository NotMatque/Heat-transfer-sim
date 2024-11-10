//
// Created by Mateusz on 04.11.2024.
//

#pragma once
#include <iostream>
#include <cmath>
#include <cstdint>


// Square matrix nxn
class SquareMatrix {
    double** matrix;
    unsigned int size;
public:
    SquareMatrix(unsigned int n = 2);
    ~SquareMatrix();
    unsigned int getSize(); // Zwraca wielkość macierzy
    double det() const; // Zwraca wyznacznik macierzy
    SquareMatrix inverse() const; // Zwraca macierz odwrotną
    friend std::ostream& operator<<(std::ostream& os, const SquareMatrix& n);
    SquareMatrix operator+(const SquareMatrix& n) const;
    SquareMatrix operator-(const SquareMatrix& n) const;
    SquareMatrix operator*(const double& n) const;
    double* operator[](unsigned int i) const;
};