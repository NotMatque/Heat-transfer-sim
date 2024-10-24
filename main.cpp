#include <cstdint>
#include <iostream>
#include <cmath>

#include "grid.h"
#include "gauss.h"



int main() {
    Node* nodes = new Node[4];
    nodes[0].x = 0.;
    nodes[0].y = 0.;

    nodes[1].x = 4.;
    nodes[1].y = 0.;

    nodes[2].x = 4.;
    nodes[2].y = 4.;

    nodes[3].x = 0.;
    nodes[3].y = 4.;

    Element elem;
    elem.nodes[0] = &nodes[0];
    elem.nodes[1] = &nodes[1];
    elem.nodes[2] = &nodes[2];
    elem.nodes[3] = &nodes[3];

    std::cout << elem << std::endl;
    std::cout << std::endl;

    elem.calculateJacobeans();
    std::cout << elem.jac[0] << std::endl;
    std::cout << elem.jac[1] << std::endl;
    std::cout << elem.jac[2] << std::endl;
    std::cout << elem.jac[3] << std::endl;

    return 0;
}
