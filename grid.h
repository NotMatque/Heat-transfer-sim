#pragma once

#include <iostream>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "matrix.h"

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

inline const double dN_dKsi_4[4][4] { //ksi, etha = +/- 1 / sqrt(3)
    {-0.394338, 0.394338, 0.105662, -0.105662},
    {-0.394338, 0.394338, 0.105662, -0.105662},
    {-0.105662, 0.105662, 0.394338, -0.394338},
    {-0.105662, 0.105662, 0.394338, -0.394338}
};
inline const double dN_dEtha_4[4][4] { //ksi, etha = +/- 1 / sqrt(3)
    {-0.394338, -0.105662, 0.105662, 0.394338},
    {-0.105662, -0.394338, 0.394338, 0.105662},
    {-0.394338, -0.105662, 0.105662, 0.394338},
    {-0.105662, -0.394338, 0.394338, 0.105662}
};
inline const double dN_dKsi_9[9][4] {
    {-0.443649, 0.443649, 0.056351, -0.056351},
    {-0.25, 0.25, 0.25, -0.25},
    {-0.056351, 0.056351, 0.443649, -0.443649},
    {-0.443649, 0.443649, 0.056351, -0.056351},
    {-0.25, 0.25, 0.25, -0.25},
    {-0.056351, 0.056351, 0.443649, -0.443649},
    {-0.443649, 0.443649, 0.056351, -0.056351},
    {-0.25, 0.25, 0.25, -0.25},
    {-0.056351, 0.056351, 0.443649, -0.443649}
};
inline const double dN_dEtha_9[9][4] {
    {-0.443649, -0.056351, 0.056351, 0.443649},
    {-0.443649, -0.056351, 0.056351, 0.443649},
    {-0.443649, -0.056351, 0.056351, 0.443649},
    {-0.25, -0.25, 0.25, 0.25},
    {-0.25, -0.25, 0.25, 0.25},
    {-0.25, -0.25, 0.25, 0.25},
    {-0.056351, -0.443649, -0.443649, 0.056351},
    {-0.056351, -0.443649, -0.443649, 0.056351},
    {-0.056351, -0.443649, -0.443649, 0.056351}
};

struct Node;
struct Element;
struct Grid;
class GlobalData;

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
    Node **nodes; // Tablica ze wskaźnikami na węzły
    Matrix* jac; // Macierze Jakobiego dla 2 lub 3 punktów całkowania
    unsigned int nIntergPoints;

    Element();
    ~Element();
    void setNIntergPoints(unsigned int nIntegrPoints);
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
    void getAllData(const std::string &path);
    void printData() const;
    void createGrid();
    void printGridNodes() const;
    void printGridElems() const;
};
