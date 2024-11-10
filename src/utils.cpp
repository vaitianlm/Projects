#include <cmath>
#include "utils.hpp"
#include <iostream>

// Contructor for Ising_lattice class
Ising_lattice::Ising_lattice(unsigned int length, double temp, double coup_const)
{
    // Assinging values
    L = length;
    T = temp;
    J = coup_const;

    // Initialising lattice
    Lattice.resize(length);
    for (int i=0; i<=length-1; i++)
    {
        Lattice[i].assign(length, true);
    }

    // Calculating Boltzmann factors
    for (int i=-8; i <= 8; i+=4)
    {
        Bf.push_back(exp(-i*coup_const/temp));
    }
}

// Function for implementing periodic boundary conditions
int pbc(int i, unsigned int L, int add)
{
    return (i + L + add) % L;
}