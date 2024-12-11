#include "grid.h"

Node::Node(): Point() {
    id = 0;
    x = 0.0;
    y = 0.0;
    isOnEdge = false;
}

Node::Node(const uint32_t id, const double x, const double y): Point(x, y) {
    this->id = id;
    this->isOnEdge = false;
}

std::ostream & operator<<(std::ostream & os, const Point & p) {
    os << "(" << p.x << ", " << p.y << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const Node& n){
    os << n.id << ": (" << n.x << ", " << n.y << ") " << (n.isOnEdge ? "on edge" : "");
    return os;
}

Element::Element() { // TODO: Przy tworzeniu się elementu powinna być od razu ustalana ilość PC
    id = 0;
    for(int i = 0; i < 4; i++)
        nodes[i] = nullptr;
    hMatrix = SquareMatrix(4);
    hbcMatrix = SquareMatrix(4);
    nIntegrPoints = 0;
    integrPoints = nullptr;
    integrPointWeights = nullptr;
    sideIntegrPoints = nullptr;
    sideIntegrPointWeights = nullptr;
    jMatrix = nullptr;
}
Element::~Element() {
    delete[] integrPoints;
    delete[] integrPointWeights;
    delete[] jMatrix;
}
void Element::setIntergPoints(unsigned int nIntegrPoints) {
    if(nIntegrPoints != 2 and nIntegrPoints != 3 and nIntegrPoints != 4) {
        std::cerr << "ERROR!\nWrong number of integration Points\n" << std::endl;
        exit(1);
    }

    // Tworzenie tablicy punktów całkowania
    this->nIntegrPoints = nIntegrPoints;
    integrPoints = new Point[nIntegrPoints * nIntegrPoints];
    integrPointWeights = new double[nIntegrPoints * nIntegrPoints];
    sideIntegrPoints = new Point[4 * nIntegrPoints];
    sideIntegrPointWeights = new double[4 * nIntegrPoints];

    // Tworzenie tablicy punktów całkowania
    switch (this->nIntegrPoints) {
        case 2: {
            // integrPoints[i]:
            //   2---3
            //  /   /
            // 0---1
            double sq1_3 = sqrt(1./3.);
            double coords[2] = {-sq1_3, sq1_3};
            for(int i = 0; i < pow(nIntegrPoints, 2); i++) {
                integrPoints[i].x = coords[i % 2];
                integrPoints[i].y = coords[i / 2];
                integrPointWeights[i] = 1;
            }
            sideIntegrPoints[0].x = coords[0]; sideIntegrPoints[0].y = -1; sideIntegrPointWeights[0] = 1;
            sideIntegrPoints[1].x = coords[1]; sideIntegrPoints[1].y = -1; sideIntegrPointWeights[1] = 1;

            sideIntegrPoints[2].x = 1; sideIntegrPoints[2].y = coords[0]; sideIntegrPointWeights[2] = 1;
            sideIntegrPoints[3].x = 1; sideIntegrPoints[3].y = coords[1]; sideIntegrPointWeights[3] = 1;

            sideIntegrPoints[4].x = coords[1]; sideIntegrPoints[4].y = 1; sideIntegrPointWeights[4] = 1;
            sideIntegrPoints[5].x = coords[0]; sideIntegrPoints[5].y = 1; sideIntegrPointWeights[5] = 1;

            sideIntegrPoints[6].x = -1; sideIntegrPoints[6].y = coords[1]; sideIntegrPointWeights[6] = 1;
            sideIntegrPoints[7].x = -1; sideIntegrPoints[7].y = coords[0]; sideIntegrPointWeights[7] = 1;
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
            double coords[3] = {-sq3_5, 0., sq3_5};
            double weights[3] = {5./9., 8./9., 5./9.};
            unsigned int index = 0;
            for(int i = 0; i < nIntegrPoints; i++)
                for(int j = 0; j < nIntegrPoints; j++) {
                    integrPoints[index].x = coords[j];
                    integrPoints[index].y = coords[i];
                    integrPointWeights[index] = weights[j] * weights[i];
                    index++;
                }
            break;
        }
        case 4: {
            // integrPoints[i]:
            //       12--13--14--15
            //      /   /   /   /
            //     8---9---10--11
            //    /   /   /   /
            //   4---5---6---7
            //  /   /   /   /
            // 0---1---2---3
            double p1 = sqrt((3./7. - (2./7. * sqrt(6./5.))));
            double p2 = sqrt((3./7. + (2./7. * sqrt(6./5.))));
            double coords[4] = {-p2, -p1, p1, p2};
            double w1 = (18. + sqrt(30)) / 36.;
            double w2 = (18. - sqrt(30)) / 36.;
            double weights[4] = {w2, w1, w1, w2};
            unsigned int index = 0;
            for(int i = 0; i < nIntegrPoints; i++)
                for(int j = 0; j < nIntegrPoints; j++) {
                    integrPoints[index].x = coords[j];
                    integrPoints[index].y = coords[i];
                    integrPointWeights[index] = weights[j] * weights[i];
                    index++;
                }
            break;
        }
    }

    // Tworzenie macierzy Jakobiego dla każdego PC
    jMatrix = new SquareMatrix[nIntegrPoints * nIntegrPoints] (2);
}
void Element::calculateJacobians() const {
    for(unsigned int i = 0; i < nIntegrPoints * nIntegrPoints; i++) {
        jMatrix[i](0, 0) = dN_dKsi[0](integrPoints[i].y) * nodes[0]->x + dN_dKsi[1](integrPoints[i].y) * nodes[1]->x + dN_dKsi[2](integrPoints[i].y) * nodes[2]->x + dN_dKsi[3](integrPoints[i].y) * nodes[3]->x;
        jMatrix[i](0, 1) = dN_dKsi[0](integrPoints[i].y) * nodes[0]->y + dN_dKsi[1](integrPoints[i].y) * nodes[1]->y + dN_dKsi[2](integrPoints[i].y) * nodes[2]->y + dN_dKsi[3](integrPoints[i].y) * nodes[3]->y;
        jMatrix[i](1, 0) = dN_dEta[0](integrPoints[i].x) * nodes[0]->x + dN_dEta[1](integrPoints[i].x) * nodes[1]->x + dN_dEta[2](integrPoints[i].x) * nodes[2]->x + dN_dEta[3](integrPoints[i].x) * nodes[3]->x;
        jMatrix[i](1, 1) = dN_dEta[0](integrPoints[i].x) * nodes[0]->y + dN_dEta[1](integrPoints[i].x) * nodes[1]->y + dN_dEta[2](integrPoints[i].x) * nodes[2]->y + dN_dEta[3](integrPoints[i].x) * nodes[3]->y;
    }
}
void Element::calculateH(double const conductivity) const {
    calculateJacobians();
    for(unsigned int i = 0; i < nIntegrPoints * nIntegrPoints; i++) {
        double dN_dx[4]; // Dla PCi
        double dN_dy[4];

        // Przemnażamy macierz Jakobiego przez odwrotność Jakobianu
        const double Jacobian = jMatrix[i].det();
        SquareMatrix invJacobianMatrix(jMatrix[i].getSize());
        invJacobianMatrix = jMatrix[i].inverse();

        // Wyliczamy dNi_dx i dNi_dy dla PCi
        for(int j = 0; j < 4; j++) {
            dN_dx[j] = invJacobianMatrix(0, 0) * dN_dKsi[j](integrPoints[i].y) + invJacobianMatrix(0, 1) * dN_dEta[j](integrPoints[i].x);
            dN_dy[j] = invJacobianMatrix(1, 0) * dN_dKsi[j](integrPoints[i].y) + invJacobianMatrix(1, 1) * dN_dEta[j](integrPoints[i].x);
            //std::cout << "PC" << i << ": " << dN_dx[j] << " " << dN_dy[j] << std::endl;
        }

        // Wyznaczamy H_PCi
        SquareMatrix hMatrixIntegrPointI(4);
        for(int j = 0; j < 4; j++)
            for(int k = 0; k < 4; k++) {
                hMatrixIntegrPointI(j, k) = (dN_dx[j] * dN_dx[k] + dN_dy[j] * dN_dy[k]) * conductivity * Jacobian;
            }
        //std::cout << hMatrixIntegrPointI << std::endl;

        for(int j = 0; j < 4; j++)
            for(int k = 0; k < 4; k++)
                hMatrix(j, k) += hMatrixIntegrPointI(j, k) * integrPointWeights[i];
    }
}
void Element::calculateHbc(double const alpha) const {
    for(int side = 0; side < 4; side++) {
        double detJ = sqrt(pow(nodes[side]->x - nodes[(side + 1) % 4]->x,2)
            + pow(nodes[side]->y - nodes[(side + 1) % 4]->y,2)) / 2.;
        SquareMatrix sideHbcMatrix(4);
        for(int pci = 0; pci < nIntegrPoints; pci++) {
            double Nvec[4];
            for(int i = 0; i < 4; i++) {
                Nvec[i] = Nfunc[i](sideIntegrPoints[nIntegrPoints * side + pci].x,
                    sideIntegrPoints[nIntegrPoints * side + pci].y);
                //std::cout << Nvec[i] << std::endl;
            }


            for(int i = 0; i < 4; i++)
                for(int j = 0; j < 4; j++) {
                    sideHbcMatrix(i, j) += Nvec[i] * Nvec[j]
                    * sideIntegrPointWeights[nIntegrPoints * side + pci] * alpha * detJ;
                }
        }
        //std::cout << sideHbcMatrix << std::endl;
        for(int i = 0; i < 4; i++)
            for(int j = 0; j < 4; j++)
                hbcMatrix(i, j) += sideHbcMatrix(i, j);
    }
    //std::cout << hbcMatrix << std::endl;
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
    hMatrix = SquareMatrix(nNodes);
}
Grid::~Grid() {
    delete[] nodes;
    delete[] elems;
}
void Grid::calculateHMatrixGlobal(double conductivity, double alpha) const {
    for(int i = 0; i < nElems; i++) {
        elems[i].setIntergPoints(2); // TODO: Przenieść do konstruktora
        elems[i].calculateH(conductivity);
        elems[i].calculateHbc(alpha);

        for(int j = 0; j < 4; j++)
            for(int k = 0; k < 4; k++)
                hMatrix(elems[i].nodes[j]->id, elems[i].nodes[k]->id) += elems[i].hMatrix(j, k) + elems[i].hbcMatrix(j, k);
    }
    std::cout << hMatrix << std::endl;
}

void GlobalData::checkDataTag(std::fstream* file, std::string const curr, const std::string expected) {
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

    std::string tempString;
    // === SimulationTime ===
    file >> tempString;
    checkDataTag(&file, tempString, "SimulationTime");
    file >> tempString;
    simTime = stod(tempString);

    // === SimulationStepTime ===
    file >> tempString;
    checkDataTag(&file, tempString, "SimulationStepTime");
    file >> tempString;
    simStepTime = stod(tempString);

    // === Conductivity ===
    file >> tempString;
    checkDataTag(&file, tempString, "Conductivity");
    file >> tempString;
    conductivity = stod(tempString);

    // === Alpha ===
    file >> tempString;
    checkDataTag(&file, tempString, "Alfa");
    file  >> tempString;
    alpha = stod(tempString);

    // === Tot ===
    file >> tempString;
    checkDataTag(&file, tempString, "Tot");
    file  >> tempString;
    tot = stod(tempString);

    // === Initial temperature ===
    file >> tempString;
    checkDataTag(&file, tempString, "InitialTemp");
    file  >> tempString;
    initTemp = stod(tempString);

    // === Density ===
    file >> tempString;
    checkDataTag(&file, tempString, "Density");
    file  >> tempString;
    density = stod(tempString);

    // === Specific Heat ===
    file >> tempString;
    checkDataTag(&file, tempString, "SpecificHeat");
    file  >> tempString;
    specificHeat = stod(tempString);

    // === Number of Nodes ===
    file >> tempString;
    checkDataTag(&file, tempString, "NodesNumber");
    file  >> tempString;
    nNodes = stoi(tempString);

    // === Number of Elements ===
    file >> tempString;
    checkDataTag(&file, tempString, "ElementsNumber");
    file  >> tempString;
    nElems = stoi(tempString);

    // === Height of the grid ===
    file >> tempString;
    checkDataTag(&file, tempString, "GridHeight");
    file  >> tempString;
    gridHeight = stoi(tempString);

    // === Width of the grid
    file >> tempString;
    checkDataTag(&file, tempString, "GridWidth");
    file  >> tempString;
    gridWidth = stoi(tempString);

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
    checkDataTag(&file, ignoreMe, "*Node");
    for(uint32_t i = 0; i < nNodes; i++) {
        file >> ignoreMe; // node id
        std::string tempX, tempY;
        file >> tempX; // node x
        tempX.pop_back();
        grid->nodes[i].x = stod(tempX);
        file >> tempY; // node y
        grid->nodes[i].y = stod(tempY);
        grid->nodes[i].id = i;
    }

    // Wczytywanie elementów
    file >> ignoreMe;
    checkDataTag(&file, ignoreMe, "*Element,");
    file >> ignoreMe;
    checkDataTag(&file, ignoreMe, "type=DC2D4");
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

    // Ustawianie węzłów na granicy
    file >> ignoreMe;
    checkDataTag(&file, ignoreMe, "*BC");
    for(int i = 0; i < 11; i++) {
        file >> ignoreMe;
        ignoreMe.pop_back();
        grid->nodes[stoi(ignoreMe) - 1].isOnEdge = true;
    }
    file >> ignoreMe;
    grid->nodes[stoi(ignoreMe) - 1].isOnEdge = true;

    file.close();
}
// Wypisuje do konsoli wczytane dane
void GlobalData::printData() const {
    std::cout << "\nGlobal Data:\nSimulation time: " << simTime << " [s]\n";
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

void GlobalData::gridCalculateHMatrixGlobal() const {
    grid->calculateHMatrixGlobal(this->conductivity, this->alpha);
}
