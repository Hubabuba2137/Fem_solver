#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>

#include "src/fem_solve.cpp"

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
    std::string file_name = "Test1_4_4.txt";
    solve(file_name);
    
    return 0;
}