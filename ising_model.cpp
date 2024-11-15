// Build: 
// g++ ising_model.cpp src/utils.cpp -I include -o ising_model.exe
// Execute:
// ./ising_model.exe <temperature> <burn in cycles> <initial config [o/u]> <output filename> <MCMC cycles>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <chrono> // Clock stuff
#include "utils.hpp" // Header file
using namespace std;

int main(int argc, char* argv[])
{

    if (argc != 6)  // Expect 3 command-line arguments
    {
        std::string executable_name = argv[0];
        std::cerr << "Error: Wrong number of input arguments." << std::endl;
        std::cerr << "Usage: " << executable_name << " <temperature> <burn in cycles> <initial config [o/u]> <output filename> <MCMC cycles>" << std::endl;

        return 1;   
    }
    
    std::string initial_config = argv[3];
    if (initial_config != "o" and initial_config != "u")
    {
        std::string executable_name = argv[0];
        std::cerr << "Error: third argument must be either o or u" << std::endl;
        std::cerr << "Usage: " << executable_name << " <temperature> <burn in cycles> <initial config [o/u]> <output filename> <MCMC cycles>" << std::endl;

        return 1;
    }
     
    double T = atof(argv[1]);   
    int burn_in = atoi(argv[2]);
    std::string filename = argv[4]; 
    int cycles = atoi(argv[5]);

    unsigned int L = 20;    // Length of lattice
    double J = 1;           // Value of coupling constant
    // Taking seed for RNG from clock
    unsigned int seed = chrono::system_clock::now().time_since_epoch().count();

    // Initialising lattice
    Ising_lattice lattice = Ising_lattice(L, T, J, seed);
    if (initial_config == "u")
    {
        lattice.randomise_spins();
    }

    // Preparing MCMC cycles
    // Burn-in
    for (int i=1; i <= burn_in; i++)
    {
        lattice.MCMC_cycle();    
    }

    vector<vector<double>> data(cycles, vector<double> (3, 0)); // vector to store data from simulation

    // MCMC cycles
    double eps_sum = 0;
    for (int i=0; i<= cycles-1; i++)
    {
        lattice.MCMC_cycle();    

        data[i][0] = i+1;
        data[i][1] = lattice.E/pow(L, 2);

        eps_sum += lattice.E/pow(L, 2);
        data[i][2] = eps_sum/(i+1);
    }

    // Writing to file
    ofstream ofile;
    ofile.open(filename);

    // Parameters for formatting output
    int width = 12;
    int prec  = 4;

    // Writing to file
    for (int i = 0; i <= cycles-1; i++)
    {
        // Writing a line with the current x and y values (nicely formatted) to file
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << data[i][0]
            << std::setw(width) << std::setprecision(prec) << std::scientific << data[i][1]
            << std::setw(width) << std::setprecision(prec) << std::scientific << data[i][2]
            << std::endl;
    }  

    ofile.close();
}