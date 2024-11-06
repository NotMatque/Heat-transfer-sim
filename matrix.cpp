//
// Created by Mateusz on 04.11.2024.
//

#include "matrix.h"

Matrix::Matrix(unsigned int const n) {
    this->size = n;
    matrix = new double*[n];
    invMatrix = new double*[n];
    for (unsigned int i = 0; i < n; i++) {
        matrix[i] = new double[n] {0};
        invMatrix[i] = new double[n] {0};
    }
}
Matrix::~Matrix() {
    for(unsigned int i = 0; i < size; i++) {
        delete[] matrix[i];
        delete[] invMatrix[i];
    }
    delete[] matrix;
    delete[] invMatrix;
}
double Matrix::det() const {
    if (size == 1) {
        return matrix[0][0];
    }
    if (size == 2) {
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    }
    // KOD NIE TESTOWANY!
    double det = 0;
    for (int k = 0; k < size; ++k) {
        Matrix subMatrix(size - 1);
        for (int i = 1; i < size; ++i) {
            int subCol = 0;
            for (int j = 0; j < size; ++j) {
                if (j != k)
                    subMatrix.matrix[i - 1][subCol++] = matrix[i][j];
            }
        }
        det += matrix[0][k] * subMatrix.det() * pow(-1, k);
    }
    return det;
}
void Matrix::inverse() {
    if(det() != 0 && size == 2) {
        double temp = 1. / (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);
        invMatrix[0][0] = matrix[1][1] * temp;
        invMatrix[0][1] = -matrix[0][1] * temp;
        invMatrix[1][0] = -matrix[1][0] * temp;
        invMatrix[1][1] = matrix[0][0] * temp;
    }
}
std::ostream & operator<<(std::ostream &os, const Matrix &n) {
    os << "Matrix:" << std::endl;
    for(uint32_t i = 0; i < n.size; i++) {
        for(uint32_t j = 0; j < n.size; j++) {
            os << n.matrix[i][j] << " ";
        }
        os << std::endl;
    }
    return os;
}
