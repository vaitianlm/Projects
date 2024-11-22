#include <cmath>
#include "utils.hpp"
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace std;

// Contructor for Ising_lattice class
Ising_lattice::Ising_lattice(unsigned int length, double temp, double coup_const, unsigned int seed)
{
    // Assigning values
    L = length;
    T = temp;
    J = coup_const;

    // Initialising lattice
    Lattice.resize(length);
    for (int i=0; i<=length-1; i++)
    {
        Lattice[i].assign(length, 1);
    }

    // Finding magnetisation and energy of the lattice
    M = magnetisation();
    E = energy_tot();

    // Seeding RNG
    generator.seed(seed);

    // Setting boundaries on uniform distributions
    uniform_L = std::uniform_int_distribution<int>(0, L-1);
    uniform_01 = uniform_real_distribution<double>(0,1);
    
    // Calculating Boltzmann factors
    for (int i=-8; i <= 8; i+=4)
    {
        Bf.push_back(exp(-i*coup_const/temp));
    }
}

// Function for randomising spins
void Ising_lattice::randomise_spins()
{
    // Uniform distribution of 0's and 1's
    std::uniform_int_distribution<int> uniform_int_01(0, 1);

    for (int i=0; i<= L-1; i++)
    {
        for (int j=0; j<= L-1; j++)
        {
            Lattice[i][j] = 1 - 2*uniform_int_01(generator); // Either 1 or -1
        }
    }
    // Evaluating magnetisation and energy of lattice
    M = magnetisation();
    E = energy_tot();
}

// Function for implementing periodic boundary conditions
int Ising_lattice::pbc(int i, int add)
{
    return (i + L + add) % L;
}

// Prints lattice to terminal
void Ising_lattice::print()
{
    for (int i = 0; i <= L-1; i++) 
    { 
        for (int j = 0; j <= L-1; j++) 
        { 
            cout << Lattice[i][j] << " "; 
        } 
        cout << endl;
    }
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

// Computes difference in energy after flipping spin L_ij
double Ising_lattice::energy_diff(int i, int j)
{
    double dE = 2*Lattice[i][j] *(Lattice[i][pbc(j, 1)] + Lattice[i][pbc(j, -1)] 
                                    + Lattice[pbc(i, 1)][j] + Lattice[pbc(i, -1)][j]);
    return dE;
}

// Does one MCMC step on spin L_ij
void Ising_lattice::MCMC_step(double r, int i, int j)
{
    double dE = energy_diff(i, j);
    int ind = dE/4 + 2; // Maps energy to an indice of Ising_lattice.Bf
    double w = Bf[ind];
    if (r < w)
    {
        Lattice[i][j] *= -1;
        E += dE;
        M += 2*Lattice[i][j];
    }
}

// Does one MCMC cycle
void Ising_lattice::MCMC_cycle()
{
    for (int k=0; k<=pow(L, 2)-1; k++)
    {
        double r = uniform_01(generator);
        int i = uniform_L(generator);
        int j = uniform_L(generator);

        MCMC_step(r, i, j);
    }
}

// Writes samples to file
void Ising_lattice::write_samples(string initial_config, int cycles, vector<vector<double>>& samples)
{
    stringstream s1; 
    s1<<T; 
    string T_string = s1.str(); // Converting the float value to string 

    stringstream s2; 
    s2<<L; 
    string L_string = s2.str(); // Converting the float value to string 

    ofstream ofile;
    string filename = "T" + T_string + "_L" +  L_string + "_N" + to_string(cycles) + "_" + initial_config + "_samples.txt";
    ofile.open(filename);

    // Parameters for formatting output
    int width = 12;
    int prec  = 4;

    // Writing to file
    for (int i = 0; i <= cycles-1; i++)
    {
        // Writing a line with the current x and y values (nicely formatted) to file
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << samples[i][0]
            << std::setw(width) << std::setprecision(prec) << std::scientific << samples[i][1]
            << std::setw(width) << std::setprecision(prec) << std::scientific << samples[i][2]
            << std::setw(width) << std::setprecision(prec) << std::scientific << samples[i][3]
            << std::setw(width) << std::setprecision(prec) << std::scientific << samples[i][4]
            << std::setw(width) << std::setprecision(prec) << std::scientific << samples[i][5]
            << std::endl;
    }  
    ofile.close();
}

// Writes calculated system values to file
void write_values(unsigned int L, int kmax, int cycles, std::vector<std::vector<double>>& values)
{
    ofstream ofile;
    string filename = "L" +  to_string(L) + "_N" + to_string(cycles) + "_" + to_string(kmax+1) + "_values.txt";
    ofile.open(filename);

    // Parameters for formatting output
    int width = 12;
    int prec  = 4;
    
    // Writing to file
    for (int i = 0; i <= kmax; i++)
    {
        // Writing a line with the current x and y values (nicely formatted) to file
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << values[i][0]
            << std::setw(width) << std::setprecision(prec) << std::scientific << values[i][1]
            << std::setw(width) << std::setprecision(prec) << std::scientific << values[i][2]
            << std::setw(width) << std::setprecision(prec) << std::scientific << values[i][3]
            << std::setw(width) << std::setprecision(prec) << std::scientific << values[i][4]
            << std::setw(width) << std::setprecision(prec) << std::scientific << values[i][5]
            << std::setw(width) << std::setprecision(prec) << std::scientific << values[i][6]
            << std::endl;
    }
    ofile.close();
}