// Execute with: 
// g++ testefil_1.cpp src/utils.cpp -I include -o testefil_1.exe; ./testefil_1.exe

#include <iomanip>
#include <iostream>
#include "utils.hpp" // Header file
using namespace std;

int main()
{
    Ising_lattice lattice = Ising_lattice(2, 3, 1);
    lattice.Lattice[0][0] = -1;
    lattice.Lattice[1][1] = -1;
    cout<< lattice.energy_tot() << endl;
}