#pragma once

#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <cstdint>
#include <fstream>
#include <iomanip>

#define N_NODES_PER_ELEMENT 4

// Węzeł siatki o określonych koordynatach x, y
struct Node {
    uint32_t id;
    double x;
    double y;
    Node();
    Node(uint32_t id, double x, double y);
    friend std::ostream& operator<<(std::ostream& os, const Node& n);
};


// Element składający się z N_NODES_PER_ELEMENT węzłów
struct Element {
    uint32_t id;
    Node *nodes[N_NODES_PER_ELEMENT]; // Tablica ze wskaźnikami na węzły
    ~Element();
    friend std::ostream& operator<<(std::ostream& os, const Element& e);
};


struct Grid {
    uint32_t nNodes; // Liczba węzłów
    uint32_t nElems; // Liczba elementów
    uint32_t height; // Wysokość siatki
    uint32_t width; // Szerokość siatki
    Node *nodes; // Tablica węzłów
    Element *elems; // Tablica elementów

    Grid(uint32_t _nNodes, uint32_t _nElems, uint32_t _height, uint32_t _width);
    ~Grid();
    void createGrid() const;
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
public:
    GlobalData();
    void getData(const std::string &path);
    void printData() const;
    void createGrid();
    void printGridNodes() const;
    void printGridElems() const;
};
#endif //GRID_H
