// Execute with:
// g++ main.cpp src/Quantum_box.cpp -I include -larmadillo -o main.exe; ./main.exe

#include <cmath>
#include "Quantum_box.hpp"
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace std;

int main()
{
    double h = 0.005;
    double dt = 2.5e-5;
    double T = dt;
    // double T = 0.008;
    double v0 = 0;
    double xc = 0.25;
    double sigmax = 0.05;
    double px = 200;
    double yc = 0.5;
    double sigmay = 0.05;
    double py = 0.0;

    Quantum_box double_slit = Quantum_box(T, dt, h, v0, xc, yc, px, py, sigmax, sigmay);
    
    double_slit.run_simulation();

    string filename = "files/double_slit_test.bin";
    bool success = double_slit.S.save(filename);

    if (success == true)
    {
        cout << "Results written to file " << filename << endl;
    }
    else
    {
        cout << "Could not write to file" << endl;
    }
    
}