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
        Lattice[i].assign(length, 1);
    }

    M = magnetisation();
    E = energy_tot();

    // Calculating Boltzmann factors
    for (int i=-8; i <= 8; i+=4)
    {
        Bf.push_back(exp(-i*coup_const/temp));
    }
}

// Function for implementing periodic boundary conditions
int Ising_lattice::pbc(int i, int add)
{
    return (i + L + add) % L;
}

// Computes magnetisation
int Ising_lattice::magnetisation()
{
    int mag = 0;
    for (int i = 0; i <= L-1; i++)
    {
        for (int j = 0; j <= L-1; j++)
        {
            mag += Lattice[i][j];
        }
    }

    return mag;
}

// Computes total energy of lattice
double Ising_lattice::energy_tot()
{
    double energy = 0;
    for (int i=0; i <= L-1; i++)
    {
        for (int j=0; j <= L-1; j++)
        {
            energy -= Lattice[i][j]*(Lattice[pbc(i, -1)][j] + Lattice[i][pbc(j, -1)]);
        }
    }
    return energy;
}

// Computes dE (tested)
double Ising_lattice::energy_diff(int i, int j)
{
    double dE = 2*Lattice[i][j] *(Lattice[i][pbc(j, 1)] + Lattice[i][pbc(j, -1)] 
                                    + Lattice[pbc(i, 1)][j] + Lattice[pbc(i, -1)][j]);
    return dE;
}

// Does one MCMC step on spin L_ij
void Ising_lattice::MCMC_step(double r, int i, int j)
{
    if (energy_diff(i, j) <= r)
    {
        Lattice[i][j] *= -1;
        E += energy_diff(i, j);
        M += 2*Lattice[i][j];
    }
}

