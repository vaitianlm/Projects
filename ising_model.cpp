// Build: 
// g++ ising_model.cpp -fopenmp src/utils.cpp -I include -o ising_model.exe -fopenmp
// Execute:
// ./ising_model.exe <min temperature> <max temperature> <lattice length> <initial config [o/u]> <store samples [true/false]> <MCMC cycles>


#include <iomanip>
#include <iostream>
#include <fstream>
#include <chrono> // Clock stuff
#include "utils.hpp" // Header file
using namespace std;

int main(int argc, char* argv[])
{
    // Handling command line arguments
    if (argc != 9)  // Expect 7 command-line arguments
    {
        string executable_name = argv[0];
        cerr << "Error: Wrong number of input arguments." << endl;
        cerr << "Usage: " << executable_name << "<min temperature> <max temperature> <#temp steps> <lattice length> <initial config [o/u]> <store samples [true/false]> <MCMC cycles> <burn-in cycles>" << std::endl;

        return 1;   
    }
    
    string initial_config = argv[5];
    if (initial_config != "o" and initial_config != "u")
    {
        string executable_name = argv[0];
        cerr << "Error: fifth argument must be either o or u" << endl;
        cerr << "Usage: " << executable_name << "<min temperature> <max temperature> <#temp steps> <lattice length> <initial config [o/u]> <store samples [true/false]> <MCMC cycles>" << std::endl;

        return 1;
    }
    
    std::string store_samples = argv[6];
    if (store_samples != "true" and store_samples != "false")
    {
        std::string executable_name = argv[0];
        std::cerr << "Error: sixth argument must be either true or false" << std::endl;
        std::cerr << "Usage: " << executable_name << "<min temperature> <max temperature> <#temp steps> <lattice length> <initial config [o/u]> <store samples [true/false]> <MCMC cycles>" << std::endl;

        return 1;
    }

    double T_min = atof(argv[1]);   // Minimum temperature [J/k]
    double T_max = atof(argv[2]);   // Maximum temperature [J/k]
    int T_num = atof(argv[3])-1;      // Number of temperature steps
    double T_step = 1;
    if (T_max != T_min)
    {
        T_step = (T_max-T_min)/T_num;   // Temperature stepsize
    }
    unsigned int L = atoi(argv[4]); // Length of lattice
    int cycles = atoi(argv[7]);     // Number of MCMC cycles to do
    int burn_in = atoi(argv[8]);    // Specify burn-in time
    double J = 1;                   // Value of coupling constant

    // Taking seed for RNG from clock
    unsigned int seed = chrono::system_clock::now().time_since_epoch().count();
    // unsigned int seed = 12345;

    // Or setting seed manually
    // unsigned int seed = 1;
    // Vector to store calculated system values
    vector<vector<double>> sys_vals(T_num+1, vector<double>(7)); 

    #pragma omp parallel for
    for(int k = 0; k <= T_num; k ++)
    {
        double T = T_min + k*T_step;
        
        // Giving each run of the loop a unique seed
        seed += k;

        // Initialising lattice
        Ising_lattice lattice = Ising_lattice(L, T, J, seed);

        // Resetting seed to not interfere with the rest of the loop
        seed -= k;

        // Randomise spin directions if specified
        if (initial_config == "u")
        {
            lattice.randomise_spins();
        }

        // Burn-in
        for (int i=1; i <= burn_in; i++)
        {
            lattice.MCMC_cycle();    
        }

        // MCMC cycles
        double eps = 0;
        double eps_sum = 0;
        double eps_sq_sum = 0;

        double m = 0;
        double m_sum_abs = 0;
        double m_sq_sum = 0;

        double eps_mean;
        double m_abs_mean;
        double heat_cap;
        double mag_sus;

        if (store_samples == "true" and T_min == T_max)
        {
            // Vector to store samples 
            vector<vector<double>> samples(cycles, vector<double> (6, 0)); 

            // Doing the MCMC cycles
            for (int i=0; i<= cycles-1; i++)
            {
                lattice.MCMC_cycle();    

                // Calculating energy/spin and related quantities
                eps = lattice.E/pow(L, 2);
                eps_sum += eps;
                eps_sq_sum += pow(eps, 2);

                // Calculating magnetisation/spin and related quantities
                m = lattice.M/pow(L, 2);
                m_sum_abs += abs(m);
                m_sq_sum += pow(m, 2);

                double eps_mean = eps_sum/(i+1);
                double m_abs_mean = m_sum_abs/(i+1);
                double heat_cap = pow(L, 2) /pow(T, 2) *(eps_sq_sum/(i+1) - pow(eps_mean, 2));
                double mag_sus = pow(L, 2) /T *(m_sq_sum/(i+1) - pow(m_abs_mean, 2));

                // Storing samples
                samples[i][0] = i+1;
                samples[i][1] = eps;
                samples[i][2] = eps_mean;
                samples[i][3] = m_abs_mean;
                samples[i][4] = heat_cap;
                samples[i][5] = mag_sus;
            }
            sys_vals[k][0] = T;
            sys_vals[k][1] = eps_sum/cycles;
            sys_vals[k][2] = m_sum_abs/cycles;
            sys_vals[k][3] = heat_cap;
            sys_vals[k][4] = mag_sus;
            sys_vals[k][5] = eps_sq_sum/cycles;
            sys_vals[k][6] = m_sq_sum/cycles;

            // ==============================================================
            // Writing to file
            lattice.write_samples(initial_config, cycles, samples);
        }
        else // Same as above, but doesn't store samples each cycle.
        {
            // Displaying progress in terminal
            if (k <= T_num/4)
            {
                cout << k << "/" << floor((T_num+1)/4) <<endl;
            }
            for (int i=0; i<= cycles-1; i++)
            {
                lattice.MCMC_cycle();    

                eps = lattice.E/pow(L, 2);
                eps_sum += eps;
                eps_sq_sum += pow(eps, 2);

                m = lattice.M/pow(L, 2);
                m_sum_abs += abs(m);
                m_sq_sum += pow(m, 2);
            }

            // Calculating system values (Only dependent on T and L, not fluctuations)
            double eps_mean = eps_sum/cycles;
            double m_abs_mean = m_sum_abs/cycles;
            double heat_cap = pow(L, 2) /pow(T, 2) *(eps_sq_sum/cycles - pow(eps_mean, 2));
            double mag_sus = pow(L, 2) /T *(m_sq_sum/cycles - pow(m_abs_mean, 2));

            sys_vals[k][0] = T;
            sys_vals[k][1] = eps_sum/cycles;
            sys_vals[k][2] = m_sum_abs/cycles;
            sys_vals[k][3] = heat_cap;
            sys_vals[k][4] = mag_sus;
            sys_vals[k][5] = eps_sq_sum/cycles;
            sys_vals[k][6] = m_sq_sum/cycles;
        }
    }
    // Writing to file
    write_values(L, T_num, cycles, sys_vals);
}