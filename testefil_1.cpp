// Execute with: 
// g++ testefil_1.cpp src/utils.cpp -I include -o testefil_1.exe; ./testefil_1.exe

#include <iomanip>
#include <iostream>
#include <chrono> // Clock stuff
#include "utils.hpp" // Header file

using namespace std;

int main()
{
    unsigned int L = 4;
    double T = 3;
    double J = 1;
    unsigned int seed = chrono::system_clock::now().time_since_epoch().count();

    Ising_lattice lattice = Ising_lattice(L, T, J, seed);

    lattice.Lattice[0][0] = -1;
    lattice.Lattice[1][1] = -1;
    lattice.Lattice[2][2] = -1;
    lattice.Lattice[3][3] = -1;
    lattice.E = lattice.energy_tot();
    lattice.print();
    cout << lattice.energy_tot() << endl << endl;

    for (int i=0; i<=8; i++)
    {
        lattice.MCMC_cycle();    
    }
}