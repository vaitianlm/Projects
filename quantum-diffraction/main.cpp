// Build with:
// g++ main.cpp src/Quantum_box.cpp -I include -larmadillo -o main.exe

// Execute with
// ./main.exe params/<input_filename.txt> files/<output_filename.bin> <track deviation [true/false]>

#include "Quantum_box.hpp"
using namespace std;

int main(int argc, char* argv[])
{
    if (argc != 4)  // Expect 3 command-line arguments
    {
        // Get the name of the executable file
        string executable_name = argv[0];

        cerr << "Error: Wrong number of input arguments." << endl;
        cerr << "Usage: " << executable_name << " params/<input_filename> files/<output_filename.bin> <track deviation [true/false]> "  << std::endl;

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
        cerr << "Problem with loading file. Did you include filepath?" << endl;
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

    // Writing to file
    string filename = argv[2];
    ok = double_slit.S.save(filename);

    if (ok == false)
    {
        cerr << "Problem with saving file. Does the specified folder for the output exist?" << endl;
        return 1;
    }

    // Writing deviation to file if wanted
    string write_deviation = argv[3];
    if (write_deviation != "false")
    {
        string filename2 = "files/deviation.bin";
        ok = double_slit.norm_dev.save(filename2);

        if (ok == false)
        {
            string filename2 = "deviation.bin";
            double_slit.norm_dev.save(filename2);
        }
        if (write_deviation != "true")
        {
            cout << "Third command line argument was neither true nor false. The deviation file was written anyways." << endl;
        }
    }
}