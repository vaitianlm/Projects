// Execute with:
// g++ main.cpp src/Quantum_box.cpp -I include -larmadillo -o main.exe; ./main.exe

#include <cmath>
#include "Quantum_box.hpp"
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char* argv[])
{
    if (argc != 2)  // Expect 1 command-line argument
    {
        // Get the name of the executable file
        string executable_name = argv[0];

        cerr << "Error: Wrong number of input arguments." << endl;
        cerr << "Usage: " << executable_name << " <param_filename> " << std::endl;

        // Exit program with non-zero return code to indicate a problem
        return 1;   
    }

    // Reading parameters from file
    string param_filename = argv[1];
    arma::vec params;
    bool ok = params.load(param_filename, arma::arma_ascii);

    if (ok == false)
    {
        cerr << "Problem with loading file. Did you include folder/?" << endl;
        return 1;
    }

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

    Quantum_box double_slit = Quantum_box(T, dt, h, int(slits), xc, yc, px, py, sigmax, sigmay);
    
    double_slit.run_simulation();

    string filename = "files/double_slit_test.bin";
    bool success = double_slit.S.save(filename);

    string filename2 = "files/deviation.bin";
    double_slit.norm_dev.save(filename2);

    if (success == true)
    {
        cout << "Results written to file " << filename << endl;
    }
    else
    {
        cout << "Could not write to file" << endl;
    }   
}