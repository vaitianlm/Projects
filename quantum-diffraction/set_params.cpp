// Execute with:
// g++ set_params.cpp -larmadillo -o set_params.exe; ./set_params.exe

#include <armadillo>

int main()
{   // Choose filename
    std::string filename = "params/detector_2_slits.txt";

    // Simulation Parameters
    double h = 0.005;           // Stepsize x and y (space)
    double dt = 2.5e-5;         // Stepsize t (time)
    double T = 0.002;           // Total time of simulation
    double slits = 2;           // Number of slits. Must be one of {0, 1, 2, 3}. Value 0 removes the barrier.
    double xc = 0.25;           // x-value of the centre of the initial wavefunction (Gaussian)
    double sigmax = 0.05;       // Std in x-direction of initial wf
    double px = 200;            // Momentum in x-direction
    double yc = 0.5;            // y-value of initial wavefuntion centre
    double sigmay = 0.2;        // Std in y-direction of initial wf
    double py = 0.0;            // Momentum in y-direction

    // Filling vector with params
    arma::vec params({h, dt, T, slits, xc, sigmax, px, yc, sigmay, py});

    // Writing to file
    params.save(filename, arma::arma_ascii);
}