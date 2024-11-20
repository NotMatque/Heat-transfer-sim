#pragma once

#include <iostream>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <cmath>


#include "matrix.h"

#define N_NODES_PER_ELEMENT 4
#define N_INTEGRATION_POINTS 4

inline double dN1_dKsi(double const etha) { return -0.25 * (1 - etha);}
inline double dN2_dKsi(double const etha) { return 0.25 * (1 - etha);}
inline double dN3_dKsi(double const etha) { return 0.25 * (1 + etha);}
inline double dN4_dKsi(double const etha) { return -0.25 * (1 + etha);}
inline double (*dN_dKsi[4])(double) {dN1_dKsi, dN2_dKsi, dN3_dKsi, dN4_dKsi};

inline double dN1_dEta(double const ksi) { return -0.25 * (1 - ksi);}
inline double dN2_dEta(double const ksi) { return -0.25 * (1 + ksi);}
inline double dN3_dEta(double const ksi) { return 0.25 * (1 + ksi);}
inline double dN4_dEta(double const ksi) { return 0.25 * (1 - ksi);}
inline double (*dN_dEta[4])(double) {dN1_dEta, dN2_dEta, dN3_dEta, dN4_dEta};


// Węzeł siatki o określonych koordynatach x, y
struct Point {
    uint32_t id;
    double x;
    double y;

    Point();
    Point(uint32_t id, double x, double y);
    friend std::ostream& operator<<(std::ostream& os, const Point& n);
};

// Element składający się z (N_NODES_PER_ELEMENT =) 4 węzłów
// ID lokalne węzłów:
//   4---3
//  /   /
// 1---2
class Element {
    unsigned int nIntegrPoints; // Ilość punktów całkowania
    Point* integrPoints; // Tablica z punktami całkowania
    double* integrPointWeight;
    SquareMatrix* jMatrix; // Macierze Jakobiego dla 2 lub 3 punktów całkowania
public:
    uint32_t id;
    Point* nodes[4]; // Tablica ze wskaźnikami na węzły // TODO: shared_ptr
    SquareMatrix hMatrix;

    Element();
    ~Element();
    void setIntergPoints(unsigned int nIntegrPoints);
    void calculateJacobians() const;
    void calculateH(double) const;
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
    Point *nodes; // Tablica węzłów // TODO: shared_ptr
    Element *elems; // Tablica elementów

    Grid(uint32_t _nNodes, uint32_t _nElems, uint32_t _height, uint32_t _width);
    ~Grid();
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

    static void checkDataTag(std::fstream*, std::string, std::string);
public:
    GlobalData();
    void getOnlyData(const std::string &path);
    void getAllData(const std::string &path);
    void printData() const;
    void printGridNodes() const;
    void printGridElems() const;
};
