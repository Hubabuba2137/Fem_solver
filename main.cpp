#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>

#include "src/fem_solve.cpp"

#include "src/nas_to_mesh.cpp"

/*//first time:
mkdir build
cd build
cmake ..
cmake --build .
./Debug/FEM.exe

//building:
cmake --build .
./Debug/FEM.exe
*/

int main()
{
    //nas_to_mesh("mesh2.nas", "mesh2.txt");
    //std::string file_name = "Test1_4_4.txt";
    std::string file_name = "mesh2.txt";
    solve(file_name, 1);
    
    return 0;
}