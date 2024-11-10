#include <cstdint>
#include <iostream>
#include <cmath>

#include "grid.h"
#include "quadGauss.h"
#include "matrix.h"

using namespace std;

void funcTest() {
    Element el;
    el.nodes[0] = new Point(1, 0.01, -0.01);
    el.nodes[1] = new Point(2, 0.025, 0);
    el.nodes[2] = new Point(3, 0.025, 0.025);
    el.nodes[3] = new Point(4, 0, 0.025);

    el.setIntergPoints(3);
    el.calculateH(30);
    std::cout << el.hMatrix << std::endl;
}

int main() {
    funcTest();
    return 0;
}
