// Include guard
#ifndef __Quantum_box_hpp__  
#define __Quantum_box_hpp__

#include <armadillo>

class Quantum_box
{
public:
    double T;           // Total time of simulation
    double dt;          // Stepsize time
    double h;           // Stepsize space
    arma::cx_double r;  // i*dt/2h^2
    int M;              // Points along the spacial axes
    arma::sp_cx_mat V;  // Potential matrix
    arma::cx_vec u;     // Wave function
    arma::sp_cx_mat A;  // CN scheme matrix A
    arma::sp_cx_mat B;  // CN scheme matrix B
    arma::cx_cube S;    // 3D-matrix to store wavefunction at every timestep
    arma::vec norm_dev; // Stores deviation of wavefunction norm from 1

    // Contructor
    Quantum_box(double T_in, double dt_in, double h_in, int slits, double xc_in, 
                double yc_in, double px_in, double py_in, double sigmax_in, double sigmay_in);

    // Translates two indices i,j to single indice k  
    int vec_ind(int i, int j);

    // Generates vector containing diagonal elements of matrix A or B
    // Input 1 for A and -1 for B
    arma::cx_vec ab_vec(double pm);

    // Constructs matrices A or B for CN scheme
    // Input -r, ab_vec(1) for A and r, ab_vec(-1) for B
    arma::sp_cx_mat CN_matAB(arma::cx_double r, arma::cx_vec d);

    // Creates the wavefunction u(x,y,t=0)
    arma::cx_vec u_init(double xc,double yc, double px, double py, double sigmax, double sigmay);

    // Creates potential matrix V containing potantial at every point x,y
    // slits=0 makes potential 0 everywhere
    void slits_init(int slits);

    // Evolves wavefunction one timestep using Crank-Nicolson scheme
    void CN_step();

    // Stores wavefuntion at every timestep
    void save_wf(int t_ind);

    // Calculates and stores deviation of wavefuntion norm from 1
    void deviation(int t_ind);

    // Does the simulation and saves the wavefunction at every iteration
    void run_simulation();
};

#endif  // end of include guard __Quantum_box_hpp__