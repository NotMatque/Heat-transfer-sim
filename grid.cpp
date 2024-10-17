#include "grid.h"

Node::Node() {
    id = 0;
    x = 0.0;
    y = 0.0;
}
Node::Node(const uint32_t id, const double x, const double y) {
    this->id = id;
    this->x = x;
    this->y = y;
}
std::ostream& operator<<(std::ostream& os, const Node& n){
    os << n.id << ": (" << n.x << ", " << n.y << ")";
    return os;
}

Element::~Element() {
    delete[] nodes;
}
std::ostream& operator<<(std::ostream& os, const Element& e) {
    os << e.id << "\n";
    for(uint32_t i = 0; i < N_NODES_PER_ELEMENT; i++) {
        os << *e.nodes[i] << "\n";
    }
    return os;
}

Grid::Grid(uint32_t const _nNodes, uint32_t const _nElems, uint32_t const _height, uint32_t const _width) {
    nNodes = _nNodes;
    nElems = _nElems;
    height = _height;
    width = _width;
    if(height*width != nNodes) {
        std::cerr << "ERROR: Grid size is not square!" << std::endl;
        exit(1);
    }
    nodes = new Node[nNodes];
    elems = new Element[nElems];
}
Grid::~Grid() {
    delete[] nodes;
    delete[] elems;
}
void Grid::createGrid() const { // BARDZO BRZYDKIE ROZWIĄZANIE!!!!!!!!!!!!!!!!!!!
    // Creating nodes
    uint32_t h = 0, w = 0;
    for(uint32_t i = 0; i < nNodes; i++) {
        nodes[i].id = i;
        nodes[i].x = w * 0.1; // Zmienić na step
        nodes[i].y = h++ * 0.1; // Zmienić na step
        if(h >= height) {
            h = 0;
            w++;
        }
    }
    // Creating elements
    for(uint32_t i = 0; i < nElems; i++) {
        elems[i].id = i;
        elems[i].nodes[0] = &nodes[i];
        elems[i].nodes[1] = &nodes[i + 1];
        elems[i].nodes[2] = &nodes[i + height];
        elems[i].nodes[3] = &nodes[i + 1 + height];
    }
}

GlobalData::GlobalData() {
    simTime = 0.0;
    simStepTime = 0.0;
    conductivity = 0.0;
    alpha = 0.0;
    tot = 0.0;
    initTemp = 0.0;
    density = 0.0;
    specificHeat = 0.0;
    gridHeight = 0;
    gridWidth = 0;
    nNodes = 0;
    nElems = 0;
    grid = nullptr;
}
// Wczytuje dane z pliku
void GlobalData::getData(const std::string &path) {
    std::fstream file;
    file.open(path);
    if (!file.is_open()) {
        std::cerr << "ERROR: Could not open file " << path << std::endl;
        file.close();
        exit(1);
    }

    std::string ignoreMe;
    file >> ignoreMe >> simTime; //SimulationTime
    file >> ignoreMe >> simStepTime; //SimulationStepTime
    file >> ignoreMe >> conductivity; //Conductivity
    file >> ignoreMe >> alpha; //Alfa
    file >> ignoreMe >> tot; //Tot
    file >> ignoreMe >> initTemp; //InitialTemp
    file >> ignoreMe >> density; //Density
    file >> ignoreMe >> specificHeat; //SpecificHeat
    file >> ignoreMe >> nNodes;
    file >> ignoreMe >> nElems;
    file >> ignoreMe >> gridHeight;
    file >> ignoreMe >> gridWidth;

    file.close();
}
// Wypisuje do konsoli wczytane dane
void GlobalData::printData() const {
    std::cout << "\nData\nSimulation time: " << simTime << " [s]\n";
    std::cout << "Simulation step time: " << simStepTime << " [s]\n";
    std::cout << "Conductivity: " << conductivity << " [W/(mK)]\n";
    std::cout << "Alpha: " << alpha << "\n";
    std::cout << "Tot: " << tot << "\n";
    std::cout << "Initial temperature: " << initTemp << " [K]\n";
    std::cout << "Density: " << density << " [kg/m^3]\n";
    std::cout << "Specific heat: " << specificHeat << "\n";
    std::cout << "Number of nodes: " << nNodes << "\n";
    std::cout << "Number of elements: " << nElems << "\n";
    std::cout << "Grid height: " << gridHeight << "\n";
    std::cout << "Grid width: " << gridWidth << "\n";
}
// Tworzy siatkę MES
void GlobalData::createGrid() {
    grid = new Grid(nNodes, nElems, gridHeight, gridWidth);
    grid->createGrid();
}
// Wypisuje koordynaty wszystkich węzłów siatki MES
void GlobalData::printGridNodes() const{
    std::cout << "Nodes:\n";
    for(uint32_t i = 0; i < nNodes; i++)
        std::cout << grid->nodes[i] << "\n";
}
// Wypisuje wszystkie dane elementów siatki MES
void GlobalData::printGridElems() const{
    std::cout << "Elements:\n";
    for(uint32_t i = 0; i < nElems; i++)
        std::cout << grid->elems[i] << "\n";
}
