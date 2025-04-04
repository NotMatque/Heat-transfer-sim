## Heat transfer simulation 2D 🔥

**First version** of this program was made during college course "_Finite element method_"
at _AGH University of Krakow_.

Since then, Since then, I have worked on it independently,
driven by my fascination with finite element analysis,
and I have come to see it as a neat little passion project 🌻 of mine.

### Main features 👀

This program simulates 2D heat transfer on a given FEM grid over a specified time step.

Program outputs temperatures in each grid nodes that can be transformed into _.vlk_ format for easier visualization (may not be yet implemented 😜)

### How to use? 📰
In [Data](./Data) folder You can find 4 files that program takes as inputs for the simulation:

**1. File [sim_params](./Data/sim_params.txt)**\
In this file the program keeps all the data that control the simulation.

Definitions are described below:
```
SimulationTime <time of the simulation in seconds>
SimulationStepTime <step time in seconds>
AmbientTemp <constant ambient temperature in Celcius or Kalvin>
InitialTemp <initial temperature of the object  in Celcius or Kalvin>
SubstancesNumber <number of substances the object consists of>
NodesNumber <number of nodes in FEM grid>
ElementsNumber <number of square elements in FEM grid>
```
**2. File [substances](./Data/substances.csv)**\
In this file the program keeps all the data about substances simulated object is made of:

```
<substance name>, <conductivity>, <equivalent convection coefficient>, <density>, <specific heat>
```

**3. File [nodes](Data/nodes.csv)**\
In this file the program keeps all the data about grid nodes like:

```
<node id>, <x position>, <y position>, <information if the node is on edge of the grid>
```

**4. File [elements](Data/elements.csv)**\
In this file the program keeps all the data about square FEM elements:

```
<element number>, <no of upper right node>, <no of upper left node>, <no of lower left node>, <no of lower right node>, <id of substance the element is made of>
```

Now, with proper data the program outputs the calculations in each step in [Results](Results) directory 😉


