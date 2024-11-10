// Execute with: 
// g++ testefil_1.cpp src/utils.cpp -I include -o testefil_1.exe; ./testefil_1.exe

#include <iomanip>
#include <iostream>
#include "utils.hpp" // Header file

int main()
{
    Ising_lattice lattice = Ising_lattice(3, 3, 1);

}