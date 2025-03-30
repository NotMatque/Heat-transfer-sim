#pragma once

#include <array>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <vector>
#include <sstream>

#include "matrix.h"
#include "substanceData.h"

#define N_NODES_PER_ELEMENT 4
#define FLOATING_POINT_PRECISION 12

inline double N1(double const ksi, double const eta) { return 0.25 * (1 - ksi) * (1 - eta); }
inline double N2(double const ksi, double const eta) { return 0.25 * (1 + ksi) * (1 - eta); }
inline double N3(double const ksi, double const eta) { return 0.25 * (1 + ksi) * (1 + eta); }
inline double N4(double const ksi, double const eta) { return 0.25 * (1 - ksi) * (1 + eta); }
inline double (*nFunc[4])(double, double){N1, N2, N3, N4};

inline double dN1_dKsi(double const eta) { return -0.25 * (1 - eta); }
inline double dN2_dKsi(double const eta) { return 0.25 * (1 - eta); }
inline double dN3_dKsi(double const eta) { return 0.25 * (1 + eta); }
inline double dN4_dKsi(double const eta) { return -0.25 * (1 + eta); }
inline double (*dN_dKsi[4])(double){dN1_dKsi, dN2_dKsi, dN3_dKsi, dN4_dKsi};

inline double dN1_dEta(double const ksi) { return -0.25 * (1 - ksi); }
inline double dN2_dEta(double const ksi) { return -0.25 * (1 + ksi); }
inline double dN3_dEta(double const ksi) { return 0.25 * (1 + ksi); }
inline double dN4_dEta(double const ksi) { return 0.25 * (1 - ksi); }
inline double (*dN_dEta[4])(double){dN1_dEta, dN2_dEta, dN3_dEta, dN4_dEta};

// Simple point
struct Point {
    double x, y;
    Point(): x(0), y(0) {}
    Point(double const x, double const y) : x(x), y(y) {}
    friend std::ostream &operator<<(std::ostream &, const Point &);
};
// Grid node on x and y coords
struct Node : Point {
    uint32_t id;
    bool isOnEdge;

    Node();
    Node(uint32_t, double, double);
    friend std::ostream &operator<<(std::ostream &, const Node &);
};

// Element composed of (N_NODES_PER_ELEMENT =) 4 nodes
// local node ID:
//   4---3
//  /   /
// 1---2
struct Element {
    std::shared_ptr<SubstanceData> substance; // Substance the elements is made of
    unsigned int nIntegrPoints; // Number of integration points
    Point *integrPoints; // Array of integration points
    double *integrPointWeights; // Array of integration point weights
    Point *sideIntegrPoints; // Array of side integration points
    double *sideIntegrPointWeights; // Array of side integration point weights
    SquareMatrix *jMatrix; // Array of Jacobi matrices for all integration points
public:
    uint32_t id; // Id of this element
    std::array<std::shared_ptr<Node>,4> nodes {};
    SquareMatrix hMatrix; //
    SquareMatrix hbcMatrix;
    SquareMatrix cMatrix; // Heat capacity matrix
    std::array<double, 4> pVector{};

    Element();
    ~Element();
    void setIntergPoints(unsigned int nIntegrPoints);
    void calculateJacobians() const;
    void calculateH() const;
    void calculateHbc() const;
    void calculateP(double);
    void calculateC() const;
    friend std::ostream &operator<<(std::ostream &os, const Element &e);
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
    std::vector<std::shared_ptr<Node>> nodes; // Tablica węzłów
    Element *elems; // Tablica elementów
    SquareMatrix hMatrixGlobal; // Macierz H globalna
    SquareMatrix cMatrixGlobal;
    double *pVector, *tVector;

    Grid(uint32_t _nNodes, uint32_t _nElems, double initTemp);
    ~Grid();
    void calculateHMatrixGlobal() const;
    void calculatePVectorGlobal(double) const;
    void calculateCMatrixGlobal() const;
    void calculateTVectorTransient(double);
    void calculateTVectorStaticState() const;

    void clearAllCalculations();
};

class GlobalData {
    double simTime;
    double simStepTime;
    double ambientTemp;
    double initTemp;
    unsigned int nSubstances;
    std::vector<std::shared_ptr<SubstanceData>> substances;
    uint32_t nNodes;
    uint32_t nElems;
    Grid *grid;

    static void checkDataTag(std::fstream *, std::string, std::string);
    void saveToFile(unsigned int) const;
    void checkIfDataIsLoaded() const;
public:
    GlobalData();

    void getAllDataFromDir(const std::string &directory);
    void getSimParamsFromFile(const std::string &path);
    void getSubstancesFromFile(const std::string &path);
    void getNodesFromFile(const std::string &path);
    void getElementsFromFile(const std::string &path);

    void runSimulationTransient() const;
    void runSimulationStaticState() const;

    void printData() const;
    void printGridNodes() const;
    void printGridElems() const;
};
