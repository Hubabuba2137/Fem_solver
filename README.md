# Fem_solver

Simple Finite Element Method solver, solving non-stationary heat transfer. Using triangular shape functions with Robin boundary condition. 
Loading configuration and mesh from .txt file.


## How to build 

Tou need a c/c++ compiler as well as CMAKE.

    mkdir build

    cd build

    cmake ..

    cmake --build .

    ./Debug/RAYIMGUI_temp.exe

## To do

    -reading .nas files

    -ability to freely define B.Cs(Robin and Neumann)

    -saving temperature at each step to text file for further visualisation

    -optimising gauss elimination