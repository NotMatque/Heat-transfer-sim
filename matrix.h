//
// Created by Mateusz on 04.11.2024.
//

#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <cmath>
#include <cstdint>


// Square matrix nxn
struct Matrix {
    double** matrix;
    double** invMatrix;
    unsigned int size;

    Matrix(unsigned int n = 2);
    ~Matrix();
    double det() const;
    void inverse();
    friend std::ostream& operator<<(std::ostream& os, const Matrix& n);
};



#endif //MATRIX_H
