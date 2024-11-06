#include <cstdint>
#include <iostream>
#include <cmath>

#include "grid.h"
#include "quadGauss.h"
#include "matrix.h"

int main() {
    Element el;
    el.nodes[0] = new Node(0, 0, 0);
    el.nodes[1] = new Node(0, 4, 0);
    el.nodes[2] = new Node(0, 4, 4);
    el.nodes[3] = new Node(0, 0, 4);

    el.setNIntergPoints(2);

    el.calculateJacobeans();

    for(int i =0; i < 4; i++)
        std::cout << el.jac[i] << std::endl;

    return 0;

}
