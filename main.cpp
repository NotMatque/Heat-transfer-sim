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
    elem.nodes[0] = new Point(1, 0., 0.);
    elem.nodes[1] = new Point(1, 0.25, 0.);
    elem.nodes[2] = new Point(1, 0.25, 0.25);
    elem.nodes[3] = new Point(1, 0., 0.25);
}

int main() {
    GlobalData gData;
    gData.getAllData("../Data/Test1_4_4.txt");
    gData.printData();
    gData.printGridNodes();
    gData.printGridElems();

    gData.gridCalculateHMatrixGlobal();

    return 0;
}
