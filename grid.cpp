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

Element::Element() {
    id = 0;
    nIntegrPoints = 0;
    hMatrix = SquareMatrix(4);
}
Element::~Element() {
    delete[] jac;
    delete[] nodes;
}
void Element::setNIntergPoints(unsigned int nIntegrPoints) {
    if(nIntegrPoints != 2 && nIntegrPoints != 3) {
        std::cerr << "ERROR!\nWrong number of integration Points\n" << std::endl;
        exit(1);
    }
    this->nIntegrPoints = nIntegrPoints;
    jac = new SquareMatrix[nIntegrPoints * nIntegrPoints];
    for(int i = 0; i < nIntegrPoints * nIntegrPoints; i++)
        jac[i] = SquareMatrix(2);
}
void Element::calculateJacobeans() const {
    if(nIntegrPoints == 2) {
        // Liczenie macierzy jakobiego dla pc1
        for(unsigned int i = 0; i < nIntegrPoints * nIntegrPoints; i++) {
            jac[i][0][0] = dN_dKsi_4[i][0] * nodes[0]->x + dN_dKsi_4[i][1] * nodes[1]->x + dN_dKsi_4[i][2] * nodes[2]->x + dN_dKsi_4[i][3] * nodes[3]->x;
            jac[i][0][1] = dN_dKsi_4[i][0] * nodes[0]->y + dN_dKsi_4[i][1] * nodes[1]->y + dN_dKsi_4[i][2] * nodes[2]->y + dN_dKsi_4[i][3] * nodes[3]->y;
            jac[i][1][0] = dN_dEtha_4[i][0] * nodes[0]->x + dN_dEtha_4[i][1] * nodes[1]->x + dN_dEtha_4[i][2] * nodes[2]->x + dN_dEtha_4[i][3] * nodes[3]->x;
            jac[i][1][1] = dN_dEtha_4[i][0] * nodes[0]->y + dN_dEtha_4[i][1] * nodes[1]->y + dN_dEtha_4[i][2] * nodes[2]->y + dN_dEtha_4[i][3] * nodes[3]->y;
        }
    }
    if(nIntegrPoints == 3) {
        for(unsigned int i = 0; i < nIntegrPoints * nIntegrPoints; i++) {
            double ksi = 0, etha = 0;
            switch(i) {
                case 0:
                    ksi = -0.77459;
                    etha = -0.77459;
                    break;
                case 1:
                    ksi = -0.77459;
                    etha = 0.;
                    break;
                case 2:
                    ksi = -0.77459;
                    etha = 0.77459;
                    break;
                case 3:
                    ksi = 0.;
                    etha = -0.77459;
                    break;
                case 4:
                    ksi = 0.;
                    etha = 0.;
                    break;
                case 5:
                    ksi = 0.;
                    etha = 0.77459;
                    break;
                case 6:
                    ksi = 0.77459;
                    etha = -0.77459;
                    break;
                case 7:
                    ksi = 0.77459;
                    etha = 0.;
                    break;
                case 8:
                    ksi = 0.77459;
                    etha = 0.77459;
                    break;
            }
            jac[i][0][0] = dN1_dKsi(etha) * nodes[0]->x + dN2_dKsi(etha) * nodes[1]->x + dN3_dKsi(etha) * nodes[2]->x + dN4_dKsi(etha) * nodes[3]->x;
            jac[i][0][1] = dN1_dKsi(etha) * nodes[0]->y + dN2_dKsi(etha) * nodes[1]->y + dN3_dKsi(etha) * nodes[2]->y + dN4_dKsi(etha) * nodes[3]->y;
            jac[i][1][0] = dN1_dEtha(ksi) * nodes[0]->x + dN2_dEtha(ksi) * nodes[1]->x + dN3_dEtha(ksi) * nodes[2]->x + dN4_dEtha(ksi) * nodes[3]->x;
            jac[i][1][1] = dN1_dEtha(ksi) * nodes[0]->y + dN2_dEtha(ksi) * nodes[1]->y + dN3_dEtha(ksi) * nodes[2]->y + dN4_dEtha(ksi) * nodes[3]->y;
        }
    }
}
void Element::calculateH(double conductivity) const {
    calculateJacobeans();
    for(unsigned int i = 0; i < nIntegrPoints * nIntegrPoints; i++) {
        double dN_dx[4]; // Dla PCi
        double dN_dy[4];

        // Przemnażamy macierz Jakobiego przez odwrotność Jakobianu
        const double Jacobian = jac[i].det();
        SquareMatrix invJacobianMatrix(jac[i].getSize());
        invJacobianMatrix = jac[i].inverse();

        // Wyliczamy dNi_dx i dNi_dy dla PCi
        for(int j = 0; j < 4; j++) {
            dN_dx[j] = invJacobianMatrix[0][0] * dN_dKsi_4[i][j] + invJacobianMatrix[0][1] * dN_dEtha_4[i][j];
            dN_dy[j] = invJacobianMatrix[1][0] * dN_dKsi_4[i][j] + invJacobianMatrix[1][1] * dN_dEtha_4[i][j];
        }

        // Wyznaczamy H_PCi
        SquareMatrix hMatrixIntegrPointI(4);
        for(int j = 0; j < 4; j++)
            for(int k = 0; k < 4; k++) {
                hMatrixIntegrPointI[j][k] = (dN_dx[j] * dN_dx[k] + dN_dy[j] * dN_dy[k]) * conductivity * Jacobian;
            }
        std::cout << hMatrixIntegrPointI << std::endl;

        for(int j = 0; j < 4; j++)
            for(int k = 0; k < 4; k++)
                hMatrix[j][k] += hMatrixIntegrPointI[j][k] * 1 * 1; // DO ZMIANY! WAGI Z GAUSSA
    }
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
        std::cerr << "ERROR:\nReading from file\nExpected: " << expected << "\nFound: " << curr << std::endl;
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
        std::cerr << "ERROR: Could not open file " << path << std::endl;
        file.close();
        exit(1);
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
// Wczytuje wszystkie dane z pliku
void GlobalData::getAllData(const std::string &path) {
    getOnlyData(path);

    std::fstream file;
    std::string ignoreMe;
    file.open(path);
    if (!file.is_open()) {
        std::cerr << "ERROR: Could not open file " << path << std::endl;
        file.close();
        exit(1);
    }

    for(int i = 0; i < 24; i++)
        file >> ignoreMe;

    // Wczytywanie węzłów
    grid = new Grid(nNodes, nElems, gridHeight, gridWidth);
    file >> ignoreMe;
    checkData(&file, ignoreMe, "*Node");
    for(uint32_t i = 0; i < nNodes; i++) {
        file >> ignoreMe;
        std::string tempX, tempY;
        file >> tempX;
        tempX.pop_back();
        grid->nodes[i].x = stod(tempX);
        file >> tempY;
        grid->nodes[i].y = stod(tempY);
        grid->nodes[i].id = i;
    }

    file >> ignoreMe;
    checkData(&file, ignoreMe, "*Element,");
    file >> ignoreMe;
    checkData(&file, ignoreMe, "type=DC2D4");
    for(uint32_t i = 0; i < nElems; i++) {
        file >> ignoreMe; // Nr elementu
        file >> ignoreMe; // Node nr 1
        ignoreMe.pop_back();
        grid->elems[i].nodes[0] = &grid->nodes[stoi(ignoreMe) - 1];
        file >> ignoreMe; // Node nr 2
        ignoreMe.pop_back();
        grid->elems[i].nodes[1] = &grid->nodes[stoi(ignoreMe) - 1];
        file >> ignoreMe; // Node nr 3
        ignoreMe.pop_back();
        grid->elems[i].nodes[2] = &grid->nodes[stoi(ignoreMe) - 1];
        file >> ignoreMe; // Node nr 4
        grid->elems[i].nodes[3] = &grid->nodes[stoi(ignoreMe) - 1];
        grid->elems[i].id = i;
    }

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
