#include "grid.h"

using namespace std;

void mainProgram() {
    GlobalData gData;
    gData.getAllDataFromDir("../Data/");
    gData.printData();
    gData.printGridNodes();
    gData.printGridElems();

    gData.runSimulationTransient();
    //gData.runSimulationStaticState();
}

int main() {
    mainProgram();
    return 0;
}
