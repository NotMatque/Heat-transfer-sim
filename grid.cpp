#include "grid.h"

Point::Point() {
    id = 0;
    x = 0.0;
    y = 0.0;
}
Point::Point(const uint32_t id, const double x, const double y) {
    this->id = id;
    this->x = x;
    this->y = y;
}
std::ostream& operator<<(std::ostream& os, const Point& n){
    os << n.id << ": (" << n.x << ", " << n.y << ")";
    return os;
}

Element::Element() {
    id = 0;
    for(int i = 0; i < 4; i++)
        nodes[i] = nullptr;
    hMatrix = SquareMatrix(4);
    nIntegrPoints = 0;
    integrPoints = nullptr;
    integrPointWeight = nullptr;
    jac = nullptr;
}
Element::~Element() {
    delete[] integrPoints;
    delete[] integrPointWeight;
    delete[] jac;
}
void Element::setIntergPoints(unsigned int nIntegrPoints) {
    if(nIntegrPoints != 2 and nIntegrPoints != 3) {
        std::cerr << "ERROR!\nWrong number of integration Points\n" << std::endl;
        exit(1);
    }

    // Tworzenie tablicy punktów całkowania
    this->nIntegrPoints = nIntegrPoints;
    integrPoints = new Point[nIntegrPoints * nIntegrPoints];
    integrPointWeight = new double[nIntegrPoints * nIntegrPoints];

    // Tworzenie tablicy punktów całkowania
    switch (this->nIntegrPoints) {
        case 2: {
            // integrPoints[i]:
            //   3---2
            //  /   /
            // 0---1
            double sq1_3 = sqrt(1./3.);
            integrPoints[0].x = -sq1_3; integrPoints[0].y = -sq1_3;  integrPointWeight[0] = 1;
            integrPoints[1].x = sq1_3;  integrPoints[1].y = -sq1_3;  integrPointWeight[1] = 1;
            integrPoints[2].x = sq1_3;  integrPoints[2].y = sq1_3;    integrPointWeight[2] = 1;
            integrPoints[3].x = -sq1_3; integrPoints[3].y = sq1_3;  integrPointWeight[3] = 1;
            break;
        }
        case 3: {
            // integrPoints[i]:
            //     6---7---8
            //    /   /   /
            //   3---4---5
            //  /   /   /
            // 0---1---2
            double sq3_5 = sqrt(3./5.);
            integrPoints[0].x = -sq3_5; integrPoints[0].y = -sq3_5;  integrPointWeight[0] = 25./81.; // 5/9 * 5/9
            integrPoints[1].x = 0;      integrPoints[1].y = -sq3_5;  integrPointWeight[1] = 40./81; // 5/9 * 8/9
            integrPoints[2].x = sq3_5;  integrPoints[2].y = -sq3_5;  integrPointWeight[2] = 25./81.;
            integrPoints[3].x = -sq3_5; integrPoints[3].y = 0;       integrPointWeight[3] = 40./81;
            integrPoints[4].x = 0;      integrPoints[4].y = 0;       integrPointWeight[4] = 64./81; // 8/9 * 8/9
            integrPoints[5].x = sq3_5;  integrPoints[5].y = 0;       integrPointWeight[5] = 40./81;
            integrPoints[6].x = -sq3_5; integrPoints[6].y = sq3_5;   integrPointWeight[6] = 25./81.;
            integrPoints[7].x = 0;      integrPoints[7].y = sq3_5;   integrPointWeight[7] = 40./81;
            integrPoints[8].x = sq3_5;  integrPoints[8].y = sq3_5;   integrPointWeight[8] = 25./81.;
            break;
        }
    }

    // Tworzenie macierzy Jakobiego dla każdego PC
    jac = new SquareMatrix[nIntegrPoints * nIntegrPoints];
    for(int i = 0; i < nIntegrPoints * nIntegrPoints; i++)
        jac[i] = SquareMatrix(2);
}
void Element::calculateJacobians() const {
    for(unsigned int i = 0; i < nIntegrPoints * nIntegrPoints; i++) {
        jac[i][0][0] = dN_dKsi[0](integrPoints[i].y) * nodes[0]->x + dN_dKsi[1](integrPoints[i].y) * nodes[1]->x + dN_dKsi[2](integrPoints[i].y) * nodes[2]->x + dN_dKsi[3](integrPoints[i].y) * nodes[3]->x;
        jac[i][0][1] = dN_dKsi[0](integrPoints[i].y) * nodes[0]->y + dN_dKsi[1](integrPoints[i].y) * nodes[1]->y + dN_dKsi[2](integrPoints[i].y) * nodes[2]->y + dN_dKsi[3](integrPoints[i].y) * nodes[3]->y;
        jac[i][1][0] = dN_dEta[0](integrPoints[i].x) * nodes[0]->x + dN_dEta[1](integrPoints[i].x) * nodes[1]->x + dN_dEta[2](integrPoints[i].x) * nodes[2]->x + dN_dEta[3](integrPoints[i].x) * nodes[3]->x;
        jac[i][1][1] = dN_dEta[0](integrPoints[i].x) * nodes[0]->y + dN_dEta[1](integrPoints[i].x) * nodes[1]->y + dN_dEta[2](integrPoints[i].x) * nodes[2]->y + dN_dEta[3](integrPoints[i].x) * nodes[3]->y;
    }
}
void Element::calculateH(double conductivity) const {
    calculateJacobians();
    for(unsigned int i = 0; i < nIntegrPoints * nIntegrPoints; i++) {
        double dN_dx[4]; // Dla PCi
        double dN_dy[4];

        // Przemnażamy macierz Jakobiego przez odwrotność Jakobianu
        const double Jacobian = jac[i].det();
        SquareMatrix invJacobianMatrix(jac[i].getSize());
        invJacobianMatrix = jac[i].inverse();

        // Wyliczamy dNi_dx i dNi_dy dla PCi
        for(int j = 0; j < 4; j++) {
            dN_dx[j] = invJacobianMatrix[0][0] * dN_dKsi[j](integrPoints[i].y) + invJacobianMatrix[0][1] * dN_dEta[j](integrPoints[i].x);
            dN_dy[j] = invJacobianMatrix[1][0] * dN_dKsi[j](integrPoints[i].y) + invJacobianMatrix[1][1] * dN_dEta[j](integrPoints[i].x);
            std::cout << "PC" << i << ": " << dN_dx[j] << " " << dN_dy[j] << std::endl;
        }

        // Wyznaczamy H_PCi
        SquareMatrix hMatrixIntegrPointI(4);
        for(int j = 0; j < 4; j++)
            for(int k = 0; k < 4; k++) {
                hMatrixIntegrPointI[j][k] = (dN_dx[j] * dN_dx[k] + dN_dy[j] * dN_dy[k]) * conductivity * Jacobian;
            }
        //std::cout << hMatrixIntegrPointI << std::endl;

        for(int j = 0; j < 4; j++)
            for(int k = 0; k < 4; k++)
                hMatrix[j][k] += hMatrixIntegrPointI[j][k] * integrPointWeight[i]; // DO ZMIANY! WAGI Z GAUSSA
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
    nodes = new Point[nNodes];
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
        elems[i].calculateJacobians();
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
