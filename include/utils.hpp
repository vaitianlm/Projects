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
    double E; // Energy
    int M; // magnetisation

    std::vector <std::vector<int>> Lattice; // The lattice
    std::vector <double> Bf; // Possible Boltzmann factors


    // Contructor
    Ising_lattice(unsigned int length, double temp, double coup_const);

    // Function for implementing periodic boundary conditions
    int pbc(int i, int add);

    // Computes total magnetisation of lattice
    int magnetisation();

    // Computes total energy of lattice
    double energy_tot();

    // Computes difference in energy after spin flip
    double energy_diff(int i, int j);

    // Does one MCMC step on spin L_ij
    void MCMC_step(double r, int i, int j);

    // Does one MCMC cycle
    void MCMC_cycle();
};








#endif  // end of include guard __utils_hpp__