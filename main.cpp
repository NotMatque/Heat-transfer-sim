#include <cstdint>
#include <iostream>
#include <cmath>

#include "grid.h"
#include "quadGauss.h"
#include "matrix.h"

using namespace std;

void funcTest() {
    Element el;
    el.nodes[0] = new Node(1, 0.01, -0.01);
    el.nodes[1] = new Node(2, 0.025, 0);
    el.nodes[2] = new Node(3, 0.025, 0.025);
    el.nodes[3] = new Node(4, 0, 0.025);

    el.setNIntergPoints(2);
    el.calculateH(30);
    std::cout << el.hMatrix << std::endl;
}

int main() {
    funcTest();
    return 0;
}
