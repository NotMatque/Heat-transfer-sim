#include "matrix.h"


SquareMatrix::SquareMatrix(unsigned int const rows): size(rows), matrix(std::make_shared<std::shared_ptr<double[]>[]>(rows)) {
    for (unsigned int i = 0; i < rows; i++)
        matrix[i] = std::make_shared<double[]>(rows);
}

unsigned int SquareMatrix::getSize() {
    return size;
}

double SquareMatrix::det() const {
    if (size == 1) {
        return matrix[0][0];
    }
    if (size == 2) {
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    }
    // KOD NIE TESTOWANY!
    double det = 0;
    for (int k = 0; k < size; ++k) {
        SquareMatrix subMatrix(size - 1);
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
SquareMatrix SquareMatrix::inverse() const {
    if(det() != 0 && size == 2) {
        SquareMatrix result(size);
        double temp = 1. / (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]);
        result(0, 0) = matrix[1][1] * temp;
        result(0, 1) = -matrix[0][1] * temp;
        result(1, 0) = -matrix[1][0] * temp;
        result(1, 1) = matrix[0][0] * temp;
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
            result.matrix[i][j] = matrix[i][j] + n.matrix[i][j];
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
            result.matrix[i][j] = matrix[i][j] - n.matrix[i][j];
    return result;
}

SquareMatrix SquareMatrix::operator*(const double& n) const {
    SquareMatrix result(size);
    for(unsigned int i = 0; i < size; i++)
        for(unsigned int j = 0; j < size; j++)
            result.matrix[i][j] = matrix[i][j] * n;
    return result;
}

double & SquareMatrix::operator()(unsigned int row, unsigned int col) const {
    return matrix[row][col];
}

std::ostream & operator<<(std::ostream &os, const SquareMatrix &n) {
    for(uint32_t i = 0; i < n.size; i++) {
        for(uint32_t j = 0; j < n.size; j++) {
            os << n.matrix[i][j] << " ";
        }
        os << std::endl;
    }
    return os;
}
