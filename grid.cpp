#include "grid.h"


Matrix::Matrix(uint32_t const n) {
    this->size = n;
    matrix = new double*[n];
    invMatrix = new double*[n];
    for (uint32_t i = 0; i < n; i++) {
        matrix[i] = new double[n];
        invMatrix[i] = new double[n];
    }
}
Matrix::~Matrix() {
    for(uint32_t i = 0; i < size; i++) {
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

Element::Element() {
    id = 0;
    jac = new Matrix[4] { Matrix(2), Matrix(2), Matrix(2), Matrix(2)};
    for(int i = 0; i < N_NODES_PER_ELEMENT; i++) {
        nodes[i] = nullptr;
    }
}
Element::~Element() {
    delete[] jac;
    delete[] nodes;
}
void Element::calculateJacobeans() const {
    // Liczenie macierzy jakobiego dla pc1
    double ksi = -1./sqrt(3.);
    double etha = -1./sqrt(3.);
    jac[0].matrix[0][0] = dN1_dKsi(etha) * nodes[0]->x + dN2_dKsi(etha) * nodes[1]->x + dN3_dKsi(etha) * nodes[2]->x + dN4_dKsi(etha) * nodes[3]->x;
    jac[0].matrix[0][1] = dN1_dKsi(etha) * nodes[0]->y + dN2_dKsi(etha) * nodes[1]->y + dN3_dKsi(etha) * nodes[2]->y + dN4_dKsi(etha) * nodes[3]->y;
    jac[0].matrix[1][0] = dN1_dEtha(ksi) * nodes[0]->x + dN2_dEtha(ksi) * nodes[1]->x + dN3_dEtha(ksi) * nodes[2]->x + dN4_dEtha(ksi) * nodes[3]->x;
    jac[0].matrix[1][1] = dN1_dEtha(ksi) * nodes[0]->y + dN2_dEtha(ksi) * nodes[1]->y + dN3_dEtha(ksi) * nodes[2]->y + dN4_dEtha(ksi) * nodes[3]->y;
    // Liczenie macierzy jakobiego dla pc2
    ksi = -ksi;
    jac[1].matrix[0][0] = dN1_dKsi(etha) * nodes[0]->x + dN2_dKsi(etha) * nodes[1]->x + dN3_dKsi(etha) * nodes[2]->x + dN4_dKsi(etha) * nodes[3]->x;
    jac[1].matrix[0][1] = dN1_dKsi(etha) * nodes[0]->y + dN2_dKsi(etha) * nodes[1]->y + dN3_dKsi(etha) * nodes[2]->y + dN4_dKsi(etha) * nodes[3]->y;
    jac[1].matrix[1][0] = dN1_dEtha(ksi) * nodes[0]->x + dN2_dEtha(ksi) * nodes[1]->x + dN3_dEtha(ksi) * nodes[2]->x + dN4_dEtha(ksi) * nodes[3]->x;
    jac[1].matrix[1][1] = dN1_dEtha(ksi) * nodes[0]->y + dN2_dEtha(ksi) * nodes[1]->y + dN3_dEtha(ksi) * nodes[2]->y + dN4_dEtha(ksi) * nodes[3]->y;
    // Liczenie macierzy jakobiego dla pc3
    etha = -etha;
    jac[2].matrix[0][0] = dN1_dKsi(etha) * nodes[0]->x + dN2_dKsi(etha) * nodes[1]->x + dN3_dKsi(etha) * nodes[2]->x + dN4_dKsi(etha) * nodes[3]->x;
    jac[2].matrix[0][1] = dN1_dKsi(etha) * nodes[0]->y + dN2_dKsi(etha) * nodes[1]->y + dN3_dKsi(etha) * nodes[2]->y + dN4_dKsi(etha) * nodes[3]->y;
    jac[2].matrix[1][0] = dN1_dEtha(ksi) * nodes[0]->x + dN2_dEtha(ksi) * nodes[1]->x + dN3_dEtha(ksi) * nodes[2]->x + dN4_dEtha(ksi) * nodes[3]->x;
    jac[2].matrix[1][1] = dN1_dEtha(ksi) * nodes[0]->y + dN2_dEtha(ksi) * nodes[1]->y + dN3_dEtha(ksi) * nodes[2]->y + dN4_dEtha(ksi) * nodes[3]->y;
    // Liczenie macierzu jakobiego dla pc4
    ksi = -ksi;
    jac[3].matrix[0][0] = dN1_dKsi(etha) * nodes[0]->x + dN2_dKsi(etha) * nodes[1]->x + dN3_dKsi(etha) * nodes[2]->x + dN4_dKsi(etha) * nodes[3]->x;
    jac[3].matrix[0][1] = dN1_dKsi(etha) * nodes[0]->y + dN2_dKsi(etha) * nodes[1]->y + dN3_dKsi(etha) * nodes[2]->y + dN4_dKsi(etha) * nodes[3]->y;
    jac[3].matrix[1][0] = dN1_dEtha(ksi) * nodes[0]->x + dN2_dEtha(ksi) * nodes[1]->x + dN3_dEtha(ksi) * nodes[2]->x + dN4_dEtha(ksi) * nodes[3]->x;
    jac[3].matrix[1][1] = dN1_dEtha(ksi) * nodes[0]->y + dN2_dEtha(ksi) * nodes[1]->y + dN3_dEtha(ksi) * nodes[2]->y + dN4_dEtha(ksi) * nodes[3]->y;
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
void Grid::generateGrid() const { // BARDZO BRZYDKIE ROZWIĄZANIE!!!!!!!!!!!!!!!!!!!
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
        elems[i].calculateJacobeans();
    }
}

void GlobalData::checkData(std::fstream* file, std::string const curr, const std::string expected) {
    if(curr != expected) {
        std::cerr << "ERROR:\nExpected: " << expected << "\nFound: " << curr << std::endl;
        file->close();
        exit(1);
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
// Wczytuje dane z pliku dane (bez węzłów)
void GlobalData::getOnlyData(const std::string &path) {
    std::fstream file;
    file.open(path);
    if (!file.is_open()) {

    }

    std::string ignoreMe;
    file >> ignoreMe;
    checkData(&file, ignoreMe, "SimulationTime");
    file >> simTime; //SimulationTime
    file >> ignoreMe;
    checkData(&file, ignoreMe, "SimulationStepTime");
    file >> simStepTime; //SimulationStepTime
    file >> ignoreMe;
    checkData(&file, ignoreMe, "Conductivity");
    file >> conductivity; //Conductivity
    file >> ignoreMe;
    checkData(&file, ignoreMe, "Alfa");
    file  >> alpha; //Alfa
    file >> ignoreMe;
    checkData(&file, ignoreMe, "Tot");
    file  >> tot; //Tot
    file >> ignoreMe;
    checkData(&file, ignoreMe, "InitialTemp");
    file  >> initTemp; //InitialTemp
    file >> ignoreMe;
    checkData(&file, ignoreMe, "Density");
    file  >> density; //Density
    file >> ignoreMe;
    checkData(&file, ignoreMe, "SpecificHeat");
    file  >> specificHeat; //SpecificHeat
    file >> ignoreMe;
    checkData(&file, ignoreMe, "NodesNumber");
    file  >> nNodes;
    file >> ignoreMe;
    checkData(&file, ignoreMe, "ElementsNumber");
    file  >> nElems;
    file >> ignoreMe;
    checkData(&file, ignoreMe, "GridHeight");
    file  >> gridHeight;
    file >> ignoreMe;
    checkData(&file, ignoreMe, "GridWidth");
    file  >> gridWidth;

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
    grid->generateGrid();
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
