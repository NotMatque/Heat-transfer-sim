#pragma once
#include <iostream>
#include <memory>
#include <cmath>
#include <cstdint>


// Square matrix nxn
class SquareMatrix {
    std::shared_ptr<double[]> matrix;
    unsigned int size;
public:
    SquareMatrix(unsigned int n = 2);
    ~SquareMatrix() = default;
    unsigned int getSize(); // Returns matrix size
    double det() const; // Returns determinant
    SquareMatrix inverse() const; // Returns inverse matrix
    friend std::ostream& operator<<(std::ostream& os, const SquareMatrix& n);
    SquareMatrix operator+(const SquareMatrix& n) const;
    SquareMatrix operator-(const SquareMatrix& n) const;
    SquareMatrix operator*(const double& n) const;
    double& operator()(unsigned int row, unsigned int col) const;
};