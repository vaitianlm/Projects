// Include guard
#ifndef __utils_hpp__  
#define __utils_hpp__

#include <vector>



class Ising_lattice
{
public:
    unsigned int L; // Length of lattice
    double T; // Temperature
    double J; // Coupling constant

    std::vector <std::vector<bool>> Lattice; // The lattice
    std::vector <double> Bf; // Possible Boltzmann factors


    // Contructor
    Ising_lattice(unsigned int length, double temp, double coup_const);
};

// Function for implementing periodic boundary conditions
int pbc(int i, unsigned int L, int add);







#endif  // end of include guard __utils_hpp__