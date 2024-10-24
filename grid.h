#pragma once

#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <cmath>

#define N_NODES_PER_ELEMENT 4
#define N_INTEGRATION_POINTS 4

inline double dN1_dKsi(double const etha) {
    return -0.25 * (1 - etha);
}
inline double dN2_dKsi(double const etha) {
    return 0.25 * (1 - etha);
}
inline double dN3_dKsi(double const etha) {
    return 0.25 * (1 + etha);
}
inline double dN4_dKsi(double const etha) {
    return -0.25 * (1 + etha);
}

inline double dN1_dEtha(double const ksi) {
    return -0.25 * (1 - ksi);
}
inline double dN2_dEtha(double const ksi) {
    return -0.25 * (1 + ksi);
}
inline double dN3_dEtha(double const ksi) {
    return 0.25 * (1 + ksi);
}
inline double dN4_dEtha(double const ksi) {
    return 0.25 * (1 - ksi);
}


struct Matrix;
struct Node;
struct Element;
struct Grid;
class GlobalData;

// Square matrix nxn
struct Matrix {
    double** matrix;
    double** invMatrix;
    uint32_t size;

    explicit Matrix(uint32_t);
    ~Matrix();
    double det() const;
    void inverse();
    friend std::ostream& operator<<(std::ostream& os, const Matrix& n);
};

// Węzeł siatki o określonych koordynatach x, y
struct Node {
    uint32_t id;
    double x;
    double y;

    Node();
    Node(uint32_t id, double x, double y);
    friend std::ostream& operator<<(std::ostream& os, const Node& n);
};

// Element składający się z (N_NODES_PER_ELEMENT =) 4 węzłów
// Węzły:
//   4---3
//  /   /
// 1---2
struct Element {
    uint32_t id;
    Node *nodes[N_NODES_PER_ELEMENT]; // Tablica ze wskaźnikami na węzły
    Matrix* jac; // Macierze Jakobiego dla 4 punktów całkowania

    Element();
    ~Element();
    void calculateJacobeans() const;
    friend std::ostream& operator<<(std::ostream& os, const Element& e);
};

// Siatka MES składająca się z elementów
// Elementy:
//     /---/---/
//    / 2 / 4 /
//   /---/---/
//  / 1 / 3 /
// /---/---/
struct Grid {
    uint32_t nNodes; // Liczba węzłów
    uint32_t nElems; // Liczba elementów
    uint32_t height; // Wysokość siatki
    uint32_t width; // Szerokość siatki
    Node *nodes; // Tablica węzłów
    Element *elems; // Tablica elementów

    Grid(uint32_t _nNodes, uint32_t _nElems, uint32_t _height, uint32_t _width);
    ~Grid();
    void generateGrid() const;
};

class GlobalData {
    double simTime;
    double simStepTime;
    double conductivity;
    double alpha;
    double tot;
    double initTemp;
    double density;
    double specificHeat;
    uint32_t gridHeight;
    uint32_t gridWidth;
    uint32_t nNodes;
    uint32_t nElems;
    Grid *grid;

    void checkData(std::fstream*, std::string, std::string);
public:
    GlobalData();
    void getOnlyData(const std::string &path);
    void printData() const;
    void createGrid();
    void printGridNodes() const;
    void printGridElems() const;
};
#endif //GRID_H
