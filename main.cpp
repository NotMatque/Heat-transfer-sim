#include <cstdint>
#include <iostream>
#include <cmath>

#include "grid.h"
#include "quadGauss.h"
#include "matrix.h"

using namespace std;

void testFunc() {
    cout << "testFunc" << endl;
    Element elem = Element();
    elem.nodes[0] = new Node(1, 0., 0.);
    elem.nodes[0]->isOnEdge = true;
    elem.nodes[1] = new Node(1, 0.025, 0.);
    elem.nodes[1]->isOnEdge = true;
    elem.nodes[2] = new Node(1, 0.025, 0.025);
    elem.nodes[2]->isOnEdge = true;
    elem.nodes[3] = new Node(1, 0., 0.025);
    elem.nodes[3]->isOnEdge = true;

    elem.setIntergPoints(2);
    std::cout << "Integration points:" << endl;
    for(int i = 0; i < pow(elem.nIntegrPoints, 2); i++) {
        std::cout << elem.integrPoints[i] << " w:" << elem.integrPointWeights[i] << endl;
    }
    std::cout << "Side ingration points:" << endl;
    for(int i = 0; i < 4 * elem.nIntegrPoints; i++) {
        std::cout << elem.sideIntegrPoints[i] << " w:" << elem.sideIntegrPointWeights[i] << endl;
    }

    elem.calculateH(300.);
    std::cout << "H Matrix:" << endl << elem.hMatrix << endl;
    elem.calculateHbc(25.);
    std::cout << "H BC Matrix:" << endl << elem.hbcMatrix << endl;
}

void mainProgram() {
    GlobalData gData;
    gData.getAllData("../Data/Test1_4_4.txt");
    gData.printData();
    gData.printGridNodes();
    gData.printGridElems();

    gData.gridCalculateHMatrixGlobal();
}

int main() {
    //testFunc();
    mainProgram();
    return 0;
}
