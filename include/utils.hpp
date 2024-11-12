// Include guard
#ifndef __utils_hpp__  
#define __utils_hpp__

#include <vector>
#include <random> 
#include <iostream>

class Ising_lattice
{
public:
    unsigned int L; // Length of lattice
    double T; // Temperature
    double J; // Coupling constant
    double E; // Energy
    int M; // magnetisation
    unsigned int Seed; // Seed for RNG

    std::mt19937 generator; // RNG
    std::uniform_int_distribution<int> uniform_L; // Uniform distrubution of integers between 0 and L-1
    std::uniform_real_distribution<double> uniform_01; // niform distrubution of real numbers between 0 and 1

    std::vector <std::vector<int>> Lattice; // The lattice
    std::vector <double> Bf; // Possible Boltzmann factors

    // Contructor
    Ising_lattice(unsigned int length, double temp, double coup_const, unsigned int seed);

    // Function for implementing periodic boundary conditions
    int pbc(int i, int add);

    // Prints lattice to terminal
    void print();

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