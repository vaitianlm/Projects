// Build with:
// g++ main.cpp src/Quantum_box.cpp -I include -larmadillo -o main.exe
// Execute with:
// ./main.exe <params/input_filename> <files/output_filename.bin>

#include <cmath>
#include "Quantum_box.hpp"
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{
    if (argc != 3)  // Expect 2 command-line argument
    {
        // Get the name of the executable file
        string executable_name = argv[0];

        cerr << "Error: Wrong number of input arguments." << endl;
        cerr << "Usage: " << executable_name << " <params/input_filename> <files/output_filename.bin>" << std::endl;

        // Exit program with non-zero return code to indicate a problem
        return 1;   
    }

    // Reading parameters from file
    string param_filename = argv[1];
    arma::vec params;
    bool ok = params.load(param_filename, arma::arma_ascii);

    // Checking for correct file loading
    if (ok == false)
    {
        cerr << "Problem with loading file. Did you include folder/?" << endl;
        return 1;
    }

    // Defining parameters
    double h = params(0);
    double dt = params(1);
    double T = params(2);
    double slits = params(3);
    double xc = params(4);
    double sigmax = params(5);
    double px = params(6);
    double yc = params(7);
    double sigmay = params(8);
    double py = params(9);

    // Initialising instance of Quantum_box
    Quantum_box double_slit = Quantum_box(T, dt, h, int(slits), xc, yc, px, py, sigmax, sigmay);
    
    // Running simulation
    double_slit.run_simulation();

    string filename = argv[2];
    double_slit.S.save(filename);

    string filename2 = "files/deviation.bin";
    double_slit.norm_dev.save(filename2);
}