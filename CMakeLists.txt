cmake_minimum_required(VERSION 3.29)
project(Heat_transfer_sim)

set(CMAKE_CXX_STANDARD 20)

add_executable(Heat_transfer_sim Src/main.cpp
        Src/grid.cpp
        Src/grid.h
        Src/quadGauss.cpp
        Src/quadGauss.h
        Src/matrix.cpp
        Src/matrix.h
        Src/substanceData.h
)
