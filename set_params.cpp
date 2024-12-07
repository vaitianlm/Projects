// Execute with:
// g++ set_params.cpp -larmadillo -o set_params.exe; ./set_params.exe

#include <armadillo>

int main()
{
    // Simulation Parameters
    double h = 0.005;
    double dt = 2.5e-5;
    double T = 0.008;
    double slits = 0;
    double xc = 0.25;
    double sigmax = 0.05;
    double px = 200;
    double yc = 0.5;
    double sigmay = 0.05;
    double py = 0.0;

    // Filling vector with params
    arma::vec params({h, dt, T, slits, xc, sigmax, px, yc, sigmay, py});

    // Writing to file
    std::string filename = "params/deviation.txt";
    params.save(filename, arma::arma_ascii);
}