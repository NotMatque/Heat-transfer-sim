#include "matrix.h"


SquareMatrix::SquareMatrix(unsigned int const rows): size(rows){
        matrix = std::make_shared<double[]>(rows * rows);
}

unsigned int SquareMatrix::getSize() {
    return size;
}

double SquareMatrix::det() const {
    if (size == 1) {
        return matrix[0];
    }
    if (size == 2) {
        return matrix[0] * matrix[1 * size + 1] - matrix[1] * matrix[1 * size];
    }
    // KOD NIE TESTOWANY!
    double det = 0;
    for (int k = 0; k < size; ++k) {
        SquareMatrix subMatrix(size - 1);
        for (int i = 1; i < size; ++i) {
            int subCol = 0;
            for (int j = 0; j < size; ++j) {
                if (j != k) {
                    subMatrix.matrix[(i - 1) * subMatrix.getSize() + subCol] = matrix[i * size + j];
                    subCol++;
                }

            }
        }
        det += matrix[k] * subMatrix.det() * pow(-1, k);
    }
    return det;
}
SquareMatrix SquareMatrix::inverse() const {
    if(det() != 0 && size == 2) {
        SquareMatrix result(size);
        double temp = 1. / (matrix[0] * matrix[1 * size + 1] - matrix[1] * matrix[1 * size]);
        result(0, 0) = matrix[1 * size + 1] * temp;
        result(0, 1) = -matrix[1] * temp;
        result(1, 0) = -matrix[1 * size] * temp;
        result(1, 1) = matrix[0] * temp;
        return result;
    }
    else {
        std::cerr << "ERROR!\nCan't inverse the matrix!" << std::endl;
        exit(1);
    }
}

SquareMatrix SquareMatrix::operator+(const SquareMatrix &n) const {
    if(size != n.size) {
        std::cerr << "ERROR!\nCan't add matrices of different sizes!\n" << std::endl;
        exit(1);
    }
    SquareMatrix result(size);
    for(unsigned int i = 0; i < size; i++)
        for(unsigned int j = 0; j < size; j++)
            result.matrix[i * size + j] = matrix[i * size + j] + n.matrix[i * size + j];
    return result;
}

SquareMatrix SquareMatrix::operator-(const SquareMatrix &n) const {
    if(size != n.size) {
        std::cerr << "ERROR!\nCan't subtract matrices of different sizes!\n" << std::endl;
        exit(1);
    }
    SquareMatrix result(size);
    for(unsigned int i = 0; i < size; i++)
        for(unsigned int j = 0; j < size; j++)
            result.matrix[i * size + j] = matrix[i * size + j] - n.matrix[i * size + j];
    return result;
}

SquareMatrix SquareMatrix::operator*(const double& n) const {
    SquareMatrix result(size);
    for(unsigned int i = 0; i < size; i++)
        for(unsigned int j = 0; j < size; j++)
            result.matrix[i * size + j] = matrix[i * size + j] * n;
    return result;
}

double & SquareMatrix::operator()(const unsigned int row, const unsigned int col) const {
    return matrix[row * size + col];
}

std::ostream & operator<<(std::ostream &os, const SquareMatrix &n) {
    for(uint32_t i = 0; i < n.size; i++) {
        for(uint32_t j = 0; j < n.size; j++) {
            os << n.matrix[i * n.size + j] << " ";
        }
        os << std::endl;
    }
    return os;
}
