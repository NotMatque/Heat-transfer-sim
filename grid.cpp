
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

Element::Element() {
    substance = nullptr;
    id = 0;
    nodes = new Node*[4];
    for (int i = 0; i < 4; i++)
        nodes[i] = nullptr;
    hMatrix = SquareMatrix(4);
    hbcMatrix = SquareMatrix(4);
    cMatrix = SquareMatrix(4);
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
void Element::setIntergPoints(unsigned int const nIntegrPoints) {
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
            for (int i = 0; i < 8; i++) {
                if (i < 2) { // Dół (y = -1)
                    sideIntegrPoints[i].x = coords[i];
                    sideIntegrPoints[i].y = -1;
                } else if (i < 4) { // Prawa strona (x = 1)
                    sideIntegrPoints[i].x = 1;
                    sideIntegrPoints[i].y = coords[i - 2];
                } else if (i < 6) { // Góra (y = 1)
                    sideIntegrPoints[i].x = coords[5 - i];
                    sideIntegrPoints[i].y = 1;
                } else { // Lewa strona (x = -1)
                    sideIntegrPoints[i].x = -1;
                    sideIntegrPoints[i].y = coords[7 - i];
                }
                sideIntegrPointWeights[i] = 1;
            }
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

            for (int i = 0; i < 4 * nIntegrPoints; i++) {
                if (i < 3) { // Dół (y = -1)
                    sideIntegrPoints[i].x = coords[i];
                    sideIntegrPoints[i].y = -1;
                } else if (i < 6) { // Prawa strona (x = 1)
                    sideIntegrPoints[i].x = 1;
                    sideIntegrPoints[i].y = coords[i - 3];
                } else if (i < 9) { // Góra (y = 1)
                    sideIntegrPoints[i].x = coords[8 - i];
                    sideIntegrPoints[i].y = 1;
                } else { // Lewa strona (x = -1)
                    sideIntegrPoints[i].x = -1;
                    sideIntegrPoints[i].y = coords[11 - i];
                }
                sideIntegrPointWeights[i] = weights[i % 3];
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

            for (int i = 0; i < 4 * nIntegrPoints; ++i) {
                if (i < 4) { // Dół (y = -1)
                    sideIntegrPoints[i].x = coords[i];
                    sideIntegrPoints[i].y = -1;
                } else if (i < 8) { // Prawa strona (x = 1)
                    sideIntegrPoints[i].x = 1;
                    sideIntegrPoints[i].y = coords[i - 4];
                } else if (i < 12) { // Góra (y = 1)
                    sideIntegrPoints[i].x = coords[11 - i];
                    sideIntegrPoints[i].y = 1;
                } else { // Lewa strona (x = -1)
                    sideIntegrPoints[i].x = -1;
                    sideIntegrPoints[i].y = coords[15 - i];
                }
                sideIntegrPointWeights[i] = weights[i % 4];
            }
            break;
        }
        default:
            break;
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
void Element::calculateH() const {
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
        }

        // Wyznaczamy H_PCi
        SquareMatrix hMatrixIntegrPointI(4);
        for(int j = 0; j < 4; j++)
            for(int k = 0; k < 4; k++) {
                hMatrixIntegrPointI(j, k) = (dN_dx[j] * dN_dx[k] + dN_dy[j] * dN_dy[k]) * substance->conductivity * Jacobian;
            }

        for(int j = 0; j < 4; j++)
            for(int k = 0; k < 4; k++)
                hMatrix(j, k) += hMatrixIntegrPointI(j, k) * integrPointWeights[i];
    }
}
void Element::calculateHbc() const {
    for(int side = 0; side < 4; side++) {
        if(!nodes[side]->isOnEdge or !nodes[(side + 1) % 4]->isOnEdge)
            continue;

        double detJ = sqrt(pow(nodes[side]->x - nodes[(side + 1) % 4]->x,2)
            + pow(nodes[side]->y - nodes[(side + 1) % 4]->y,2)) / 2.;

        SquareMatrix sideHbcMatrix(4);
        for(int pci = 0; pci < nIntegrPoints; pci++) {
            double Nvec[4];
            for(int i = 0; i < 4; i++) {
                Nvec[i] = nFunc[i](sideIntegrPoints[nIntegrPoints * side + pci].x,
                    sideIntegrPoints[nIntegrPoints * side + pci].y);
                //std::cout << Nvec[i] << std::endl;
            }

            for(int i = 0; i < 4; i++)
                for(int j = 0; j < 4; j++) {
                    sideHbcMatrix(i, j) += Nvec[i] * Nvec[j]
                    * sideIntegrPointWeights[nIntegrPoints * side + pci] * substance->alpha * detJ;
                }
        }
        for(int i = 0; i < 4; i++)
            for(int j = 0; j < 4; j++)
                hbcMatrix(i, j) += sideHbcMatrix(i, j);
    }
}
void Element::calculateP(double ambTemperature) {
    for(int side = 0; side < 4; side++) {
        if(!nodes[side]->isOnEdge or !nodes[(side + 1) % 4]->isOnEdge)
            continue;

        double detJ = sqrt(pow(nodes[side]->x - nodes[(side + 1) % 4]->x,2)
            + pow(nodes[side]->y - nodes[(side + 1) % 4]->y,2)) / 2.;
        std::array<double, 4> sidePVector {};
        for(int pci = 0; pci < nIntegrPoints; pci++) {
            for(int i = 0; i < 4; i++) {
                sidePVector[i] += nFunc[i](sideIntegrPoints[nIntegrPoints * side + pci].x, sideIntegrPoints[nIntegrPoints * side + pci].y)
                * sideIntegrPointWeights[nIntegrPoints * side + pci] * ambTemperature;
            }
        }
        for(int i = 0; i < 4; i++) {
            sidePVector[i] *= substance->alpha * detJ;
            pVector[i] += sidePVector[i];
        }
    }
}
void Element::calculateC() const{
    for(unsigned int i = 0; i < nIntegrPoints * nIntegrPoints; i++) {
        double detJ = jMatrix[i].det();

        SquareMatrix cMatrixIntegrPointI(4);
        for(int j = 0; j < 4; j++)
            for(int k = 0; k < 4; k++)
                cMatrixIntegrPointI(j,k) = substance->specificHeat * substance->density * detJ * nFunc[j](integrPoints[i].x, integrPoints[i].y) * nFunc[k](integrPoints[i].x, integrPoints[i].y);

        for(int j = 0; j < 4; j++)
            for(int k = 0; k < 4; k++)
                cMatrix(j, k) += cMatrixIntegrPointI(j, k) * integrPointWeights[i];
    }
}

std::ostream& operator<<(std::ostream& os, const Element& e) {
    os << e.id << "\n";
    for(uint32_t i = 0; i < N_NODES_PER_ELEMENT; i++) {
        os << *e.nodes[i] << "\n";
    }
    return os;
}

Grid::Grid(uint32_t const _nNodes, uint32_t const _nElems, double const initTemp) {
    nNodes = _nNodes;
    nElems = _nElems;

    nodes = new Node[nNodes];
    elems = new Element[nElems];
    for(int i = 0; i < nElems; i++)
        elems[i].setIntergPoints(3);
    hMatrix = SquareMatrix(nNodes);
    cMatrix = SquareMatrix(nNodes);
    pVector = new double[nNodes];
    tVector = new double[nNodes];
    for(int i = 0; i < nNodes; i++)
        tVector[i] = initTemp;
}
Grid::~Grid() {
    delete[] nodes;
    delete[] elems;
}
void Grid::calculateHMatrixGlobal() const {
    for(int i = 0; i < nElems; i++) {
        elems[i].calculateH();
        elems[i].calculateHbc();

        for(int j = 0; j < 4; j++)
            for(int k = 0; k < 4; k++)
                hMatrix(elems[i].nodes[j]->id, elems[i].nodes[k]->id) += elems[i].hMatrix(j, k) + elems[i].hbcMatrix(j, k);
    }
}
void Grid::calculatePVectorGlobal(double ambientTemp) const {
    // Creating global vector {P}
    for(int i = 0; i < nElems; i++) {
        elems[i].calculateP(ambientTemp);
        for(int j = 0; j < 4; j++)
            pVector[elems[i].nodes[j]->id] += elems[i].pVector[j];
    }
}
void Grid::calculateCMatrixGlobal() const {
    for(int i = 0; i < nElems; i++) {
        elems[i].calculateC();

        for(int j = 0; j < 4; j++)
            for(int k = 0; k < 4; k++)
                cMatrix(elems[i].nodes[j]->id, elems[i].nodes[k]->id) += elems[i].cMatrix(j, k);
    }
}
void Grid::calculateTVectorTransient(double stepTime) {
    // Preparing data
    SquareMatrix hTransientMatrix(nNodes);
    for(size_t i = 0; i < hTransientMatrix.getSize(); i++)
        for(size_t j = 0; j < hTransientMatrix.getSize(); j++)
            hTransientMatrix(i, j) = hMatrix(i, j) + cMatrix(i, j) / stepTime;

    for(int i = 0; i < nNodes; i++) {
        for(int j = 0; j < nNodes; j++)
            pVector[i] += cMatrix(i, j) / stepTime * tVector[j];
    }

    // Preparation for Crout (not implemented)
    auto *change = new unsigned int[nNodes];
    for(int i = 0; i < nNodes; i++)
        change[i] = i;

    // Gaussian elimination
    for(int i = 0; i < nNodes - 1; i++) {
        int iMax = 0; // TODO: Gauss-Crout or better

        for(int j = i + 1; j < nNodes; j++) {
            double ratio = hTransientMatrix(j,i) / hTransientMatrix(i,i);
            for(int k = i; k < nNodes; k++)
                hTransientMatrix(j,k) -= ratio * hTransientMatrix(i,k);
            pVector[j] -= ratio * pVector[i];
        }
    }

    // Calculating values of tVector
    for(int i = nNodes - 1; i >= 0; i--) {
        tVector[i] = pVector[i];
        for(int j = i + 1; j < nNodes; j++)
            tVector[i] -= hTransientMatrix(i,j) * tVector[j];
        tVector[i] /= hTransientMatrix(i,i);
    }
    delete[] change;
}
void Grid::calculateTVectorStaticState() const{
    // Gaussian elimination
    for(int i = 0; i < nNodes - 1; i++) {
        int iMax = 0; // TODO: Gauss-Crout or better

        for(int j = i + 1; j < nNodes; j++) {
            double ratio = hMatrix(j,i) / hMatrix(i,i);
            for(int k = i; k < nNodes; k++)
                hMatrix(j,k) -= ratio * hMatrix(i,k);
            pVector[j] -= ratio * pVector[i];
        }
    }

    // Calculating values of tVector
    for(int i = (int)nNodes - 1; i >= 0; i--) {
        tVector[i] = pVector[i];
        for(int j = i + 1; j < nNodes; j++)
            tVector[i] -= hMatrix(i,j) * tVector[j];
        tVector[i] /= hMatrix(i,i);
    }
}

void Grid::clearAllCalculations() {
    hMatrix = SquareMatrix(nNodes);
    cMatrix = SquareMatrix(nNodes);
    for(int i = 0; i < nNodes; i++)
        pVector[i] = 0;
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
    ambientTemp = 0.0;
    initTemp = 0.0;
    nSubstances = 0;
    substances = nullptr;
    nNodes = 0;
    nElems = 0;
    grid = nullptr;
}

void GlobalData::getAllDataFromDir(const std::string &directory) {
    getSimParamsFromFile(directory + "sim_params.txt");

    for(unsigned int i = 1; i <= nSubstances; i++)
        getSubstanceDataFromFile(i, directory + "substance_" + std::to_string(i) + ".txt");
    getNodesFromFile(directory + "nodes.csv");
    getElementsFromFile(directory + "elements.csv");
}
void GlobalData::getSimParamsFromFile(const std::string &path) {
    std::fstream simParamsFile;
    simParamsFile.open(path);

    if (!simParamsFile.is_open()) {
        std::cerr << "ERROR: Could not open file " << path << std::endl;
        exit(1);
    }

    std::string ignoreMe, readData;

    simParamsFile >> ignoreMe >> readData;
    checkDataTag(&simParamsFile, ignoreMe, "SimulationTime");
    this->simTime = stod(readData);

    simParamsFile >> ignoreMe >> readData;
    checkDataTag(&simParamsFile, ignoreMe, "SimulationStepTime");
    this->simStepTime = stod(readData);

    simParamsFile >> ignoreMe >> readData;
    checkDataTag(&simParamsFile, ignoreMe, "AmbientTemp");
    this->ambientTemp = stod(readData);

    simParamsFile >> ignoreMe >> readData;
    checkDataTag(&simParamsFile, ignoreMe, "InitialTemp");
    this->initTemp = stod(readData);

    simParamsFile >> ignoreMe >> readData;
    checkDataTag(&simParamsFile, ignoreMe, "SubstancesNumber");
    this->nSubstances = stoi(readData);
    substances = new SubstanceData[stoi(readData)];

    simParamsFile >> ignoreMe >> readData;
    checkDataTag(&simParamsFile, ignoreMe, "NodesNumber");
    this->nNodes = stoi(readData);

    simParamsFile >> ignoreMe >> readData;
    checkDataTag(&simParamsFile, ignoreMe, "ElementsNumber");
    this->nElems = stoi(readData);

    simParamsFile.close();
}
void GlobalData::getSubstanceDataFromFile(unsigned int substanceNumber, const std::string &path) {
    std:std::fstream substanceFile;
    substanceFile.open(path);

    if (!substanceFile.is_open()) {
        std::cerr << "ERROR: Could not open file " << path << std::endl;
        exit(1);
    }

    std::string ignoreMe, readData;

    substanceFile >> ignoreMe >> readData;
    checkDataTag(&substanceFile, ignoreMe, "Name");
    substances[substanceNumber - 1].name = readData;

    substanceFile >> ignoreMe >> readData;
    checkDataTag(&substanceFile, ignoreMe, "Conductivity");
    substances[substanceNumber - 1].conductivity = stod(readData);

    substanceFile >> ignoreMe >> readData;
    checkDataTag(&substanceFile, ignoreMe, "Alpha");
    substances[substanceNumber - 1].alpha = stod(readData);

    substanceFile >> ignoreMe >> readData;
    checkDataTag(&substanceFile, ignoreMe, "Density");
    substances[substanceNumber - 1].density = stod(readData);

    substanceFile >> ignoreMe >> readData;
    checkDataTag(&substanceFile, ignoreMe, "SpecificHeat");
    substances[substanceNumber - 1].specificHeat = stod(readData);

    substanceFile.close();
}
void GlobalData::getNodesFromFile(const std::string &path) {
    if(!grid)
        grid = new Grid(nNodes, nElems, initTemp);

    std::fstream nodesFile;
    nodesFile.open(path);

    if (!nodesFile.is_open()) {
        std::cerr << "ERROR: Could not open file " << path << std::endl;
        exit(1);
    }

    std::string ignoreMe, readX, readY, readIsOnEdge;

    nodesFile >> ignoreMe >> ignoreMe >> ignoreMe >> ignoreMe;
    for(unsigned int i = 0; i < nNodes; i++) {
        nodesFile >> ignoreMe >> readX >> readY >> readIsOnEdge;
        grid->nodes[i].id = i;
        readX.pop_back();
        grid->nodes[i].x = stod(readX);
        readY.pop_back();
        grid->nodes[i].y = stod(readY);
        grid->nodes[i].isOnEdge = stoi(readIsOnEdge);
    }

    nodesFile.close();
}
void GlobalData::getElementsFromFile(const std::string &path) {
    if(!grid)
        grid = new Grid(nNodes, nElems, initTemp);

    std::fstream elementsFile;
    elementsFile.open(path);

    if (!elementsFile.is_open()) {
        std::cerr << "ERROR: Could not open file " << path << std::endl;
        exit(1);
    }

    std::string ignoreMe, readNode1, readNode2, readNode3, readNode4, readSubstanceNum;

    elementsFile >> ignoreMe >> ignoreMe >> ignoreMe >> ignoreMe >> ignoreMe >> ignoreMe;
    for(unsigned int i = 0; i < nElems; i++) {
        elementsFile >> ignoreMe >> readNode1 >> readNode2 >> readNode3 >> readNode4 >> readSubstanceNum;
        grid->elems[i].id = i;
        readNode1.pop_back();
        grid->elems[i].nodes[0] = &grid->nodes[stoi(readNode1) - 1];
        readNode2.pop_back();
        grid->elems[i].nodes[1] = &grid->nodes[stoi(readNode2) - 1];
        readNode3.pop_back();
        grid->elems[i].nodes[2] = &grid->nodes[stoi(readNode3) - 1];
        readNode4.pop_back();
        grid->elems[i].nodes[3] = &grid->nodes[stoi(readNode4) - 1];
        grid->elems[i].substance = &substances[stoi(readSubstanceNum) - 1];
    }
    elementsFile.close();
}
void GlobalData::checkIfDataIsLoaded() const {
    if(grid == nullptr) {
        std::cerr << "ERROR: grid is null" << std::endl;
        exit(1);
    }
}

// Wypisuje do konsoli wczytane dane
void GlobalData::printData() const {
    std::cout << "\nGlobal Data:\nSimulation time: " << simTime << " [s]\n";
    std::cout << "Simulation step time: " << simStepTime << " [s]\n";
    std::cout << "Ambient temperature: " << ambientTemp << "[C]\n";
    std::cout << "Initial temperature: " << initTemp << " [C]\n";
    std::cout << "Substances: " << nSubstances << "\n";
    for(unsigned int i = 0; i < nSubstances; i++) {
        std::cout << substances[i] << "\n";
    }

    std::cout << "Number of nodes: " << nNodes << "\n";
    std::cout << "Number of elements: " << nElems << "\n";
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
void GlobalData::saveToFile(unsigned int const step) const {
    std::fstream file;
    std::string path = "../Results/step" + std::to_string(step) + ".txt";

    file.open(path, std::ios::out | std::ios::trunc);

    if (!file.is_open()) {
        std::cerr << "ERROR: Could not open file " << path << std::endl;
        exit(1);
    }

    for (uint32_t i = 0; i < nNodes; i++) {
        file << std::setprecision(FLOATING_POINT_PRECISION) << grid->tVector[i] << std::endl;
    }

    file.close();
}

void GlobalData::runSimulationTransient() const {
    checkIfDataIsLoaded();

    std::cout << "Running transient simulation...\n";
    const unsigned int steps = static_cast<int>(simTime / simStepTime);
    for(unsigned int i = 1; i <= steps; i++) {
        grid->clearAllCalculations();
        grid->calculateHMatrixGlobal();
        grid->calculatePVectorGlobal(this->ambientTemp);
        grid->calculateCMatrixGlobal();
        grid->calculateTVectorTransient(this->simStepTime);
        saveToFile(i);
    }
    std::cout << "Check Results in /Results directory\n";
}
void GlobalData::runSimulationStaticState() const {
    checkIfDataIsLoaded();

    std::cout << "Running static state simulation...\n";
    grid->calculateHMatrixGlobal();
    grid->calculatePVectorGlobal(this->ambientTemp);
    grid->calculateTVectorStaticState();
    saveToFile(0);
    std::cout << "Check Results in /Results/step0\n";
}